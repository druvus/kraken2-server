#include <doctest/doctest.h>

#include "messages.h"

using kraken2proto::Kraken2SequenceRequest;

TEST_CASE("SequenceRequestToSequence copies the fields the classifier uses")
{
    Kraken2SequenceRequest req;
    req.set_id("read1");
    req.set_seq("ACGT");
    req.set_quals("IIII");
    req.set_format(Kraken2SequenceRequest::FORMAT_FASTQ);

    kraken2::Sequence seq;
    CHECK(SequenceRequestToSequence(req, seq));
    CHECK(seq.header == "read1");
    CHECK(seq.comment.empty());
    CHECK(seq.seq == "ACGT");
    CHECK(seq.quals == "IIII");
    CHECK(seq.format == kraken2::FORMAT_FASTQ);
}

TEST_CASE("SequenceRequestToSequence maps every format")
{
    Kraken2SequenceRequest req;
    kraken2::Sequence seq;
    req.set_format(Kraken2SequenceRequest::FORMAT_FASTA);
    SequenceRequestToSequence(req, seq);
    CHECK(seq.format == kraken2::FORMAT_FASTA);
    req.set_format(Kraken2SequenceRequest::FORMAT_AUTO_DETECT);
    SequenceRequestToSequence(req, seq);
    CHECK(seq.format == kraken2::FORMAT_AUTO_DETECT);
}

TEST_CASE("SequenceRequestToSequence overwrites stale state when reused")
{
    kraken2::Sequence seq;
    Kraken2SequenceRequest a, b;
    a.set_id("a"); a.set_seq("AAAA"); a.set_quals("IIII");
    b.set_id("b"); b.set_seq("CC");
    SequenceRequestToSequence(a, seq);
    SequenceRequestToSequence(b, seq);
    CHECK(seq.header == "b");
    CHECK(seq.seq == "CC");
    CHECK(seq.quals.empty());
}

TEST_CASE("SequenceRequestToPair reports and fills the mate")
{
    Kraken2SequenceRequest req;
    req.set_id("frag/1");
    req.set_seq("ACGT");
    kraken2::Sequence seq, seq2;

    SUBCASE("single read")
    {
        CHECK_FALSE(SequenceRequestToPair(req, seq, seq2));
        CHECK(seq.header == "frag/1");
    }
    SUBCASE("paired read")
    {
        req.mutable_mate()->set_id("frag/2");
        req.mutable_mate()->set_seq("TTGG");
        CHECK(SequenceRequestToPair(req, seq, seq2));
        CHECK(seq.seq == "ACGT");
        CHECK(seq2.header == "frag/2");
        CHECK(seq2.seq == "TTGG");
    }
}
