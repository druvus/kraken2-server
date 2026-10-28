#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "Kraken2.grpc.pb.h"

using kraken2proto::Kraken2SequenceRequest;

// FASTA/FASTQ reader over kseq. The kseq types and functions live in
// kseq.cc only, behind FastReader::Impl, so including this header does not
// instantiate them.
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
    // Read up to batch_size sequences, stopping early once the batch holds
    // at least max_bytes of ids, sequence and quality data, so that a
    // batch of long reads stays small. Returns the number read.
    int read(std::vector<Kraken2SequenceRequest> &seqs, int batch_size,
             uint64_t max_bytes = UINT64_MAX);
    // Read up to batch_size read pairs (or max_bytes, as for read), one
    // mate from this reader and one from mate_reader. Stops at the end of
    // the shorter input; if the two inputs end at different points,
    // mismatched is set to true.
    int read_pairs(std::vector<Kraken2SequenceRequest> &seqs, int batch_size,
                   FastReader &mate_reader, bool &mismatched,
                   uint64_t max_bytes = UINT64_MAX);
    // read all sequences
    int read_all(std::vector<Kraken2SequenceRequest> &seqs);
private:
    struct Impl;
    std::string m_filename;
    std::unique_ptr<Impl> m_impl;
    int m_last_status = 0;
};
