#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>
#include "kseq.h"
#include "kseq.cc.h"

// Input under kseq. zlib's gzread fills its whole buffer before it
// returns, so reading a pipe through it, records written slowly (by a
// running sequencer, say) are only seen once tens of kilobytes have
// arrived or the pipe closes. Here plain input returns whatever one read()
// delivers, and gzip input, any number of concatenated members, is decoded
// with inflate as it arrives. The format is decided from the first two
// bytes, on the first read. As with gzread, anything after a complete gzip
// member that is not another member (zero padding, say) ends the input.
struct Input
{
    int fd = -1;
    bool close_fd = false;
    bool detected = false;
    bool gzip = false;
    bool in_member = false;   // part way through a gzip member
    bool check_next = false;  // a member ended: is another one next?
    bool finished = false;    // no further members
    bool failed = false;      // a read error, corrupt data or a member cut short
    z_stream zs;
    unsigned char buf[65536];
    size_t pos = 0, len = 0;  // raw bytes read but not yet consumed

    ~Input()
    {
        if (gzip) inflateEnd(&zs);
        if (close_fd && fd >= 0) ::close(fd);
    }
};

static ssize_t raw_read(int fd, void *p, size_t n)
{
    ssize_t r;
    do {
        r = ::read(fd, p, n);
    } while (r < 0 && errno == EINTR);
    return r;
}

// After a gzip member ended, make at least two unconsumed bytes available
// and say whether they start another member: 1 yes, 0 no (or the input
// ended), -1 on a read error.
static int next_is_member(Input *in)
{
    if (in->zs.avail_in < 2) {
        size_t have = in->zs.avail_in;
        memmove(in->buf, in->zs.next_in, have);
        in->len = have;
        while (in->len < 2) {
            ssize_t r = raw_read(in->fd, in->buf + in->len, sizeof(in->buf) - in->len);
            if (r < 0) return -1;
            if (r == 0) break;
            in->len += r;
        }
        in->zs.next_in = in->buf;
        in->zs.avail_in = in->len;
        in->pos = in->len;
    }
    return in->zs.avail_in >= 2 && in->zs.next_in[0] == 0x1f && in->zs.next_in[1] == 0x8b;
}

static int input_read_raw(Input *in, void *out, unsigned out_len)
{
    if (!in->detected) {
        in->detected = true;
        while (in->len < 2) {
            ssize_t r = raw_read(in->fd, in->buf + in->len, sizeof(in->buf) - in->len);
            if (r < 0) return -1;
            if (r == 0) break;
            in->len += r;
        }
        if (in->len >= 2 && in->buf[0] == 0x1f && in->buf[1] == 0x8b) {
            memset(&in->zs, 0, sizeof(in->zs));
            if (inflateInit2(&in->zs, 15 + 16) != Z_OK) return -1;
            in->gzip = true;
        }
    }

    if (!in->gzip) {
        if (in->pos < in->len) {
            size_t n = std::min<size_t>(out_len, in->len - in->pos);
            memcpy(out, in->buf + in->pos, n);
            in->pos += n;
            return (int) n;
        }
        ssize_t r = raw_read(in->fd, out, out_len);
        return r < 0 ? -1 : (int) r;
    }

    if (in->finished) return 0;
    in->zs.next_out = (Bytef *) out;
    in->zs.avail_out = out_len;
    while (in->zs.avail_out == out_len) {
        if (in->check_next) {
            // Checked only once the text of the previous member has been
            // returned, so a pipe that pauses between members does not hold
            // it back.
            in->check_next = false;
            int next = next_is_member(in);
            if (next < 0) return -1;
            if (next == 0) {
                in->finished = true;
                return 0;
            }
            inflateReset(&in->zs);
        }
        if (in->zs.avail_in == 0) {
            if (in->pos >= in->len) {
                ssize_t r = raw_read(in->fd, in->buf, sizeof(in->buf));
                if (r < 0) return -1;
                if (r == 0) {
                    // the input ended; inside a member means it was cut short
                    return in->in_member ? -1 : 0;
                }
                in->pos = 0;
                in->len = r;
            }
            in->zs.next_in = in->buf + in->pos;
            in->zs.avail_in = in->len - in->pos;
            in->pos = in->len;
        }
        int ret = inflate(&in->zs, Z_NO_FLUSH);
        if (ret == Z_STREAM_END) {
            in->in_member = false;
            in->check_next = true;
        }
        else if (ret == Z_OK || ret == Z_BUF_ERROR) {
            in->in_member = true;
        }
        else {
            return -1;
        }
    }
    return (int) (out_len - in->zs.avail_out);
}

// kseq's read function: up to out_len bytes of text, 0 at the end of the
// input, -1 on a read error, a corrupt gzip stream or a gzip member cut
// short. A failure is also recorded, so that FastReader can drop a record
// that kseq returns after the stream failed part way through it.
static int input_read(Input *in, void *out, unsigned out_len)
{
    int n = input_read_raw(in, out, out_len);
    if (n < 0) in->failed = true;
    return n;
}

// A 64 KB parse buffer (kseq defaults to 16 KB) means fewer read and
// inflate calls per megabyte.
KSTREAM_INIT(Input *, input_read, 65536)
__KSEQ_TYPE(Input *)
__KSEQ_BASIC(static, Input *)
__KSEQ_READ(static)

struct FastReader::Impl
{
    Input input;
    kseq_t *seq = nullptr;

    ~Impl()
    {
        if (seq != nullptr) kseq_destroy(seq);
    }
};


FastReader::FastReader(std::string filename)
    : m_filename(filename), m_impl(new Impl)
{
    if (filename == "-") {
        m_impl->input.fd = fileno(stdin);
    }
    else {
        m_impl->input.fd = ::open(filename.c_str(), O_RDONLY);
        if (m_impl->input.fd < 0) {
            throw std::runtime_error("Cannot open " + filename + ": " + std::strerror(errno));
        }
        m_impl->input.close_fd = true;
    }
    m_impl->seq = kseq_init(&m_impl->input);
}


FastReader::~FastReader() = default;


int FastReader::read(Kraken2SequenceRequest& rec) {
    kseq_t *m_seq = m_impl->seq;
    int rtn;
    if ((rtn = kseq_read(m_seq)) < 0) {
        rec.Clear();
        m_last_status = rtn;
        if (rtn == -2) {
            std::cerr << m_filename << ": truncated quality string at record "
                      << (m_seq->name.s ? m_seq->name.s : "?")
                      << "; stopping." << std::endl;
        }
        else if (rtn == -3) {
            std::cerr << m_filename << ": error reading stream; stopping." << std::endl;
        }
        return rtn;
    }
    else if (m_impl->input.failed) {
        // The stream failed part way through this record (a gzip file cut
        // short, corrupt data): kseq returns what it had read, which may be
        // incomplete, so the record is dropped.
        std::cerr << m_filename << ": error reading stream at record "
                  << (m_seq->name.s ? m_seq->name.s : "?") << "; stopping." << std::endl;
        rec.Clear();
        m_last_status = -3;
        return -3;
    }
    else {
        // Once a file has given FASTQ records, a record without qualities
        // means the file ended part way through a record: kseq returns
        // the partial record as FASTA. An empty record in the middle of
        // the file is allowed.
        if (m_seq->qual.l > 0) {
            m_fastq = true;
        }
        else if (m_fastq && (m_seq->seq.l > 0 || ks_eof(m_seq->f))) {
            std::cerr << m_filename << ": record "
                      << (m_seq->name.s ? m_seq->name.s : "?")
                      << " has no quality string; the file appears to be truncated; stopping."
                      << std::endl;
            rec.Clear();
            m_last_status = -2;
            return -2;
        }
        m_last_status = rtn;
        // Only the fields the server uses are sent: id, sequence, qualities
        // and format. The full header line and a text rendering of the
        // record (proto fields `header` and `str_representation`) are
        // left empty to halve the request size.
        rec.set_id(m_seq->name.s, m_seq->name.l);
        rec.set_seq(m_seq->seq.s, m_seq->seq.l);
        if (m_seq->qual.l == 0)
        {
            rec.set_format(Kraken2SequenceRequest::FORMAT_FASTA);
        }
        else
        {
            rec.set_format(Kraken2SequenceRequest::FORMAT_FASTQ);
            rec.set_quals(m_seq->qual.s, m_seq->qual.l);
        }
    }
    return rtn;
}


static uint64_t record_bytes(const Kraken2SequenceRequest &rec)
{
    return rec.id().size() + rec.seq().size() + rec.quals().size();
}


int FastReader::read(std::vector<Kraken2SequenceRequest> &seqs, int batch_size,
                     uint64_t max_bytes)
{
    int rtn = 0;
    uint64_t bytes = 0;
    for(int i=0; i<batch_size && bytes < max_bytes; ++i) {
        Kraken2SequenceRequest rec;
        if (read(rec) >= 0) {
            rtn++;
            bytes += record_bytes(rec);
            seqs.push_back(std::move(rec));
        }
        else {
            break;
        }
    }
    return rtn;
}


int FastReader::read_pair(Kraken2SequenceRequest &rec, FastReader &mate_reader,
                          bool &mismatched)
{
    Kraken2SequenceRequest mate;
    bool got1 = read(rec) >= 0;
    bool got2 = mate_reader.read(mate) >= 0;
    if (got1 && got2) {
        rec.mutable_mate()->Swap(&mate);
        return 1;
    }
    // Different lengths, not a read error in one of the files.
    if (got1 != got2 && !failed() && !mate_reader.failed()) {
        mismatched = true;
    }
    rec.Clear();
    return 0;
}


int FastReader::read_pairs(std::vector<Kraken2SequenceRequest> &seqs, int batch_size,
                           FastReader &mate_reader, bool &mismatched,
                           uint64_t max_bytes)
{
    int rtn = 0;
    uint64_t bytes = 0;
    for(int i=0; i<batch_size && bytes < max_bytes; ++i) {
        Kraken2SequenceRequest rec, mate;
        bool got1 = read(rec) >= 0;
        bool got2 = mate_reader.read(mate) >= 0;
        if (got1 && got2) {
            bytes += record_bytes(rec) + record_bytes(mate);
            rec.mutable_mate()->Swap(&mate);
            rtn++;
            seqs.push_back(std::move(rec));
        }
        else {
            // Different lengths, not a read error in one of the files.
            if (got1 != got2 && !failed() && !mate_reader.failed()) {
                mismatched = true;
            }
            break;
        }
    }
    return rtn;
}


int FastReader::read_all(std::vector<Kraken2SequenceRequest> &seqs)
{
    int rtn = 0;
    while (true) {
        Kraken2SequenceRequest rec;
        if (read(rec) >= 0) {
            rtn++;
            seqs.push_back(std::move(rec));
        }
        else {
            break;
        }
    }
    return rtn;
}
