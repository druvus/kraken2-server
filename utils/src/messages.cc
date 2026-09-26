#include "messages.h"

bool SequenceRequestToSequence(
    const kraken2proto::Kraken2SequenceRequest &req, kraken2::Sequence &seq)
{
    switch (req.format())
    {
    case kraken2proto::Kraken2SequenceRequest_SequenceFormat_FORMAT_AUTO_DETECT:
        seq.format = kraken2::SequenceFormat::FORMAT_AUTO_DETECT;
        break;
    case kraken2proto::Kraken2SequenceRequest_SequenceFormat_FORMAT_FASTQ:
        seq.format = kraken2::SequenceFormat::FORMAT_FASTQ;
        break;
    case kraken2proto::Kraken2SequenceRequest_SequenceFormat_FORMAT_FASTA:
        seq.format = kraken2::SequenceFormat::FORMAT_FASTA;
        break;
    default:
        seq.format = kraken2::SequenceFormat::FORMAT_AUTO_DETECT;
        break;
    }

    // kraken2's Sequence stores the read name in `header` and the rest of
    // the header line in `comment`. Only the name is used by the classifier.
    seq.header = req.id();
    seq.comment.clear();
    seq.seq = req.seq();
    seq.quals = req.quals();

    return true;
}

bool SequenceRequestToPair(
    const kraken2proto::Kraken2SequenceRequest &req,
    kraken2::Sequence &seq, kraken2::Sequence &seq2)
{
    SequenceRequestToSequence(req, seq);
    if (!req.has_mate())
    {
        return false;
    }
    SequenceRequestToSequence(req.mate(), seq2);
    return true;
}
