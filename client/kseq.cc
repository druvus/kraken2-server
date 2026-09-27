#include <cerrno>
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
    FILE *instream = (filename == "-") ? stdin : fopen(filename.c_str(), "r");
    if (instream == nullptr) {
        throw std::runtime_error(
            "Cannot open " + filename + ": " + std::strerror(errno));
    }
    m_impl->fp = gzdopen(fileno(instream), "r");
    if (m_impl->fp == nullptr) {
        throw std::runtime_error("Cannot open " + filename + " for reading.");
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


int FastReader::read(std::vector<Kraken2SequenceRequest> &seqs, int batch_size)
{
    int rtn = 0;
    seqs.reserve(batch_size);
    for(int i=0; i<batch_size; ++i) {
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


int FastReader::read_pairs(std::vector<Kraken2SequenceRequest> &seqs, int batch_size,
                           FastReader &mate_reader, bool &mismatched)
{
    int rtn = 0;
    seqs.reserve(batch_size);
    for(int i=0; i<batch_size; ++i) {
        Kraken2SequenceRequest rec, mate;
        bool got1 = read(rec) >= 0;
        bool got2 = mate_reader.read(mate) >= 0;
        if (got1 && got2) {
            rec.mutable_mate()->Swap(&mate);
            rtn++;
            seqs.push_back(std::move(rec));
        }
        else {
            if (got1 != got2) {
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
