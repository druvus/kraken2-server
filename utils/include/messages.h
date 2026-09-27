#pragma once

#include "Kraken2.grpc.pb.h"
#include "seqreader.h"

// Copy a sequence request into a kraken2 Sequence. The request id becomes
// the Sequence header (name up to the first whitespace), matching kraken2.
bool SequenceRequestToSequence(
    const kraken2proto::Kraken2SequenceRequest &req, kraken2::Sequence &seq);

// Copy a request into seq and, when the request carries a mate, into seq2.
// Returns true when the request is a read pair.
bool SequenceRequestToPair(
    const kraken2proto::Kraken2SequenceRequest &req,
    kraken2::Sequence &seq, kraken2::Sequence &seq2);
