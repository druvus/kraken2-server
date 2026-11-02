#include <cerrno>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

#include <zlib.h>
#include "kseq.h"
#include "kseq.cc.h"

KSEQ_INIT(gzFile, gzread)

struct FastReader::Impl
{
    gzFile fp = nullptr;
    kseq_t *seq = nullptr;

    ~Impl()
    {
        if (seq != nullptr) kseq_destroy(seq);
        if (fp != nullptr) gzclose(fp);
    }
};


FastReader::FastReader(std::string filename)
    : m_filename(filename), m_impl(new Impl)
{
    // zlib owns the file descriptor and closes it in gzclose.
    errno = 0;
    m_impl->fp = (filename == "-") ? gzdopen(fileno(stdin), "r")
                                   : gzopen(filename.c_str(), "r");
    if (m_impl->fp == nullptr) {
        throw std::runtime_error(
            "Cannot open " + filename + ": " +
            (errno != 0 ? std::strerror(errno) : "out of memory"));
    }
    m_impl->seq = kseq_init(m_impl->fp);
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
