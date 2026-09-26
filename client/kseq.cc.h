#pragma once

#include <string>
#include <vector>

#include "Kraken2.grpc.pb.h"

#include <zlib.h>
#include "kseq.h"

using kraken2proto::Kraken2SequenceRequest;

KSEQ_INIT(gzFile, gzread)

class FastReader
{
public:
    // Open a FASTA/FASTQ file, plain or gzip compressed. "-" reads stdin.
    // Throws std::runtime_error if the file cannot be opened.
    FastReader(std::string filename);
    ~FastReader();
    FastReader(const FastReader &) = delete;
    FastReader &operator=(const FastReader &) = delete;

    const std::string &filename() const { return m_filename; }
    // kseq status of the last read: >= 0 record length, -1 end of file,
    // -2 truncated quality string, -3 stream error.
    int last_status() const { return m_last_status; }
    bool failed() const { return m_last_status < -1; }
    // Read (at most) one sequence
    int read(Kraken2SequenceRequest&);
    // read (up to) batch_size sequences
    int read(std::vector<Kraken2SequenceRequest> &seqs, int batch_size);
    // read (up to) batch_size read pairs, one mate from this reader and one
    // from mate_reader. Stops at the end of the shorter input; if the two
    // inputs end at different points, mismatched is set to true.
    int read_pairs(std::vector<Kraken2SequenceRequest> &seqs, int batch_size,
                   FastReader &mate_reader, bool &mismatched);
    // read all sequences
    int read_all(std::vector<Kraken2SequenceRequest> &seqs);
private:
    std::string m_filename;
    gzFile m_fp;
    kseq_t *m_seq;
    int m_last_status = 0;
};
