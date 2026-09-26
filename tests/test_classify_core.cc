#include <doctest/doctest.h>

#include "classify_core.h"
#include "test_taxonomy.h"

using namespace kraken2server;

TEST_CASE("TrimPairInfo strips /1 and /2 only")
{
    CHECK(TrimPairInfo("read/1") == "read");
    CHECK(TrimPairInfo("read/2") == "read");
    CHECK(TrimPairInfo("read/3") == "read/3");
    CHECK(TrimPairInfo("read") == "read");
    CHECK(TrimPairInfo("/1") == "/1");  // too short to trim
    CHECK(TrimPairInfo("") == "");
}

TEST_CASE("AddHitlistString formats runs, ambiguous spans and markers")
{
    TestTaxonomy tax;
    std::ostringstream oss;

    SUBCASE("single run")
    {
        std::vector<taxid_t> taxa = {TestTaxonomy::A1, TestTaxonomy::A1, TestTaxonomy::A1};
        AddHitlistString(oss, taxa, tax.get());
        CHECK(oss.str() == "11:3");
    }
    SUBCASE("runs with unclassified and ambiguous spans")
    {
        std::vector<taxid_t> taxa = {TestTaxonomy::A1, TestTaxonomy::A1, 0, 0, 0,
                                     AMBIGUOUS_SPAN_TAXON, AMBIGUOUS_SPAN_TAXON, TestTaxonomy::B};
        AddHitlistString(oss, taxa, tax.get());
        CHECK(oss.str() == "11:2 0:3 A:2 20:1");
    }
    SUBCASE("mate pair border")
    {
        std::vector<taxid_t> taxa = {TestTaxonomy::A1, MATE_PAIR_BORDER_TAXON, 0};
        AddHitlistString(oss, taxa, tax.get());
        CHECK(oss.str() == "11:1 |:| 0:1");
    }
    SUBCASE("reading frame border")
    {
        std::vector<taxid_t> taxa = {0, READING_FRAME_BORDER_TAXON, TestTaxonomy::A};
        AddHitlistString(oss, taxa, tax.get());
        CHECK(oss.str() == "0:1 -:- 10:1");
    }
    SUBCASE("trailing marker")
    {
        std::vector<taxid_t> taxa = {TestTaxonomy::A, MATE_PAIR_BORDER_TAXON};
        AddHitlistString(oss, taxa, tax.get());
        CHECK(oss.str() == "10:1 |:|");
    }
}

TEST_CASE("ResolveTree picks the best supported taxon")
{
    TestTaxonomy tax;

    SUBCASE("clear winner")
    {
        taxon_counts_t hits = {{TestTaxonomy::A1, 5}, {TestTaxonomy::B, 1}};
        CHECK(ResolveTree(hits, tax.get(), 10, 0.0) == TestTaxonomy::A1);
    }
    SUBCASE("hits on an ancestor add to the descendant's score")
    {
        // A1 gets 2 + 2 (from A) = 4, B gets 3
        taxon_counts_t hits = {{TestTaxonomy::A1, 2}, {TestTaxonomy::A, 2}, {TestTaxonomy::B, 3}};
        CHECK(ResolveTree(hits, tax.get(), 10, 0.0) == TestTaxonomy::A1);
    }
    SUBCASE("tie resolves to the lowest common ancestor")
    {
        taxon_counts_t hits = {{TestTaxonomy::A1, 3}, {TestTaxonomy::A2, 3}};
        CHECK(ResolveTree(hits, tax.get(), 10, 0.0) == TestTaxonomy::A);
    }
    SUBCASE("confidence threshold moves the call up the tree")
    {
        // 6 of 10 minimizers hit within clade A
        taxon_counts_t hits = {{TestTaxonomy::A1, 3}, {TestTaxonomy::A2, 3}};
        CHECK(ResolveTree(hits, tax.get(), 10, 0.5) == TestTaxonomy::A);
        // 0.7 requires 7 hits, which not even root has, so no call
        taxon_counts_t hits2 = {{TestTaxonomy::A1, 3}, {TestTaxonomy::A2, 3}};
        CHECK(ResolveTree(hits2, tax.get(), 10, 0.7) == 0);
    }
    SUBCASE("confidence threshold with a single leaf")
    {
        // A1 alone has 4 of 10; at 0.5 the clade of A is still 4, root 4, so unclassified
        taxon_counts_t hits = {{TestTaxonomy::A1, 4}};
        CHECK(ResolveTree(hits, tax.get(), 10, 0.4) == TestTaxonomy::A1);
        taxon_counts_t hits2 = {{TestTaxonomy::A1, 4}};
        CHECK(ResolveTree(hits2, tax.get(), 10, 0.5) == 0);
    }
    SUBCASE("no hits")
    {
        taxon_counts_t hits;
        CHECK(ResolveTree(hits, tax.get(), 10, 0.0) == 0);
    }
}

TEST_CASE("DoubleStatToString uses fixed precision")
{
    CHECK(DoubleStatToString(1.0, 2) == "1.00");
    CHECK(DoubleStatToString(12.3456, 2) == "12.35");
    CHECK(DoubleStatToString(0.0, 1) == "0.0");
}

TEST_CASE("Stats reports handle empty streams without NaN")
{
    ClassificationStats empty = {0, 0, 0};
    struct timeval t0 = {0, 0}, t1 = {0, 0};
    std::string s = ReportStats(t0, t1, empty);
    CHECK(s.find("nan") == std::string::npos);
    CHECK(s.find("0 sequences (0.00 Mbp)") == 0);
    CHECK(ReportTotalStats(empty).find("0 sequences classified (0.00%)") != std::string::npos);
}

TEST_CASE("Stats reports show counts and percentages")
{
    ClassificationStats stats = {4, 4000000, 3};
    struct timeval t0 = {10, 0}, t1 = {12, 500000};
    std::string s = ReportStats(t0, t1, stats);
    CHECK(s.find("4 sequences (4.00 Mbp) processed in 2.50s") == 0);
    CHECK(s.find("3 sequences classified (75.00%)") != std::string::npos);
    CHECK(s.find("1 sequences unclassified (25.00%)") != std::string::npos);
}
