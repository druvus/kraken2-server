#include <doctest/doctest.h>

#include <sstream>
#include <vector>

#include "report_server.h"
#include "test_taxonomy.h"

static std::vector<std::string> Lines(const std::string &s)
{
    std::vector<std::string> out;
    std::istringstream is(s);
    std::string line;
    while (std::getline(is, line))
        out.push_back(line);
    return out;
}

TEST_CASE("ReportKrakenStyle writes header, unclassified line and clade tree")
{
    TestTaxonomy tax;
    kraken2::taxon_counters_t calls;
    calls[TestTaxonomy::A1].incrementReadCount();
    calls[TestTaxonomy::A1].incrementReadCount();
    calls[TestTaxonomy::B].incrementReadCount();

    std::ostringstream ss;
    kraken2::ReportKrakenStyle(ss, false, false, tax.get(), calls, 4, 1);
    auto lines = Lines(ss.str());

    REQUIRE(lines.size() == 6);
    CHECK(lines[0] == "% of Seqs\tClades\tTaxonomies\tRank\tTaxonomy ID\tScientific Name");
    // kraken2 prints the unclassified count in both the clade and taxon columns
    CHECK(lines[1] == " 25.00\t1\t1\tU\t0\tunclassified");
    CHECK(lines[2] == " 75.00\t3\t0\tR\t1\troot");
    // A (2 reads in clade) is listed before B (1 read)
    CHECK(lines[3] == " 50.00\t2\t0\tD\t10\t  A");
    CHECK(lines[4] == " 50.00\t2\t2\tS\t11\t    A1");
    CHECK(lines[5] == " 25.00\t1\t1\tD\t20\t  B");
}

TEST_CASE("ReportKrakenStyle with zero counts and kmer columns")
{
    TestTaxonomy tax;
    kraken2::taxon_counters_t calls;
    calls[TestTaxonomy::A2].incrementReadCount();
    calls[TestTaxonomy::A2].add_kmer(1);
    calls[TestTaxonomy::A2].add_kmer(1);
    calls[TestTaxonomy::A2].add_kmer(2);

    std::ostringstream ss;
    kraken2::ReportKrakenStyle(ss, true, true, tax.get(), calls, 1, 0);
    auto lines = Lines(ss.str());

    CHECK(lines[0] == "% of Seqs\tClades\tTaxonomies\tKmers\tDistinct Kmers\tRank\tTaxonomy ID\tScientific Name");
    // zero-count unclassified line is present because report_zeros is set
    CHECK(lines[1].rfind("  0.00\t0\t0\t0\t0\tU\t0\tunclassified", 0) == 0);
    // every taxon appears, including those with no reads
    REQUIRE(lines.size() == 7);
    // 3 kmers, 2 distinct
    CHECK(lines[4] == "100.00\t1\t1\t3\t2\tS\t12\t    A2");
}
