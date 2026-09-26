#include <doctest/doctest.h>

#include <sstream>

#include "classify_core.h"
#include "merged_taxonomy.h"
#include "report_server.h"
#include "test_taxonomy.h"

using namespace kraken2server;

namespace
{
    // A second synthetic taxonomy sharing root and B with TestTaxonomy:
    //   root(1) -> B(20) -> B1(21); root -> C(30)
    std::map<uint64_t, TaxonRecord> SecondTaxa()
    {
        std::map<uint64_t, TaxonRecord> r;
        r[1] = {1, 0, "root", "no rank"};
        r[20] = {20, 1, "B", "superkingdom"};
        r[21] = {21, 20, "B1", "species"};
        r[30] = {30, 1, "C", "superkingdom"};
        return r;
    }

    struct TempTaxonomy
    {
        std::string path;
        std::unique_ptr<kraken2::Taxonomy> taxonomy;
        explicit TempTaxonomy(const std::map<uint64_t, TaxonRecord> &records)
        {
            char tmpl[] = "/tmp/k2s_merged_XXXXXX";
            int fd = mkstemp(tmpl);
            close(fd);
            path = tmpl;
            WriteKrakenTaxonomy(records, path);
            taxonomy.reset(new kraken2::Taxonomy(path, false));
            taxonomy->GenerateExternalToInternalIDMap();
        }
        ~TempTaxonomy() { unlink(path.c_str()); }
    };
}

TEST_CASE("CollectTaxa reads every node with external parent, name and rank")
{
    TestTaxonomy tax;
    auto records = CollectTaxa(tax.get());
    REQUIRE(records.size() == 5);
    CHECK(records[1].parent_external_id == 0);
    CHECK(records[10].parent_external_id == 1);
    CHECK(records[11].parent_external_id == 10);
    CHECK(records[11].name == "A1");
    CHECK(records[11].rank == "species");
    CHECK(records[20].name == "B");
}

TEST_CASE("WriteKrakenTaxonomy round-trips through kraken2::Taxonomy")
{
    TestTaxonomy tax;
    auto records = CollectTaxa(tax.get());
    TempTaxonomy copy(records);
    auto again = CollectTaxa(*copy.taxonomy);
    CHECK(again == records);
    // structure is usable by the kraken2 tree functions
    auto a1 = copy.taxonomy->GetInternalID(11), a2 = copy.taxonomy->GetInternalID(12);
    auto a = copy.taxonomy->GetInternalID(10);
    CHECK(copy.taxonomy->LowestCommonAncestor(a1, a2) == a);
    CHECK(copy.taxonomy->IsAAncestorOfB(a, a1));
}

TEST_CASE("MergeTaxa unions consistent taxonomies and rejects conflicts")
{
    TestTaxonomy tax;
    std::map<uint64_t, TaxonRecord> merged;
    MergeTaxa(merged, CollectTaxa(tax.get()), "db1");
    MergeTaxa(merged, SecondTaxa(), "db2");
    CHECK(merged.size() == 7);
    CHECK(merged[21].parent_external_id == 20);
    CHECK(merged[30].parent_external_id == 1);

    auto conflicting = SecondTaxa();
    conflicting[11] = {11, 20, "A1", "species"};  // A1 under B instead of A
    CHECK_THROWS_WITH_AS(MergeTaxa(merged, conflicting, "db3"),
                         doctest::Contains("taxid 11"), std::runtime_error);
}

TEST_CASE("BuildMergedTaxonomy and MapToMerged")
{
    TestTaxonomy tax1;
    TempTaxonomy tax2(SecondTaxa());
    auto merged = BuildMergedTaxonomy({&tax1.get(), tax2.taxonomy.get()}, {"db1", "db2"});
    REQUIRE(merged->node_count() == 8);  // 7 taxa plus node 0

    auto map1 = MapToMerged(tax1.get(), *merged);
    auto map2 = MapToMerged(*tax2.taxonomy, *merged);
    // the same external id maps to the same merged internal id from both sides
    CHECK(map1[TestTaxonomy::B] == map2[tax2.taxonomy->GetInternalID(20)]);
    CHECK(map1[TestTaxonomy::ROOT] == map2[tax2.taxonomy->GetInternalID(1)]);
    CHECK(map1[0] == 0);
    // and external ids survive
    CHECK(merged->nodes()[map1[TestTaxonomy::A1]].external_id == 11);
    CHECK(merged->nodes()[map2[tax2.taxonomy->GetInternalID(21)]].external_id == 21);
    // B1 (only in db2) and A1 (only in db1) meet at root
    auto lca = merged->LowestCommonAncestor(map1[TestTaxonomy::A1], map2[tax2.taxonomy->GetInternalID(21)]);
    CHECK(merged->nodes()[lca].external_id == 1);
    // B1 and B meet at B
    lca = merged->LowestCommonAncestor(map1[TestTaxonomy::B], map2[tax2.taxonomy->GetInternalID(21)]);
    CHECK(merged->nodes()[lca].external_id == 20);
}

TEST_CASE("MergeTaxonVectors takes the LCA per position and keeps markers")
{
    TestTaxonomy tax;
    const auto &t = tax.get();
    std::vector<taxid_t> db1 = {TestTaxonomy::A1, 0, TestTaxonomy::A1, AMBIGUOUS_SPAN_TAXON, MATE_PAIR_BORDER_TAXON, TestTaxonomy::B};
    std::vector<taxid_t> db2 = {TestTaxonomy::A2, TestTaxonomy::B, 0, AMBIGUOUS_SPAN_TAXON, MATE_PAIR_BORDER_TAXON, TestTaxonomy::B};
    std::vector<taxid_t> out;
    MergeTaxonVectors({db1, db2}, t, out);
    CHECK(out == std::vector<taxid_t>{TestTaxonomy::A,   // LCA(A1, A2)
                                      TestTaxonomy::B,   // LCA(0, B) = B
                                      TestTaxonomy::A1,  // LCA(A1, 0) = A1
                                      AMBIGUOUS_SPAN_TAXON, MATE_PAIR_BORDER_TAXON,
                                      TestTaxonomy::B});

    SUBCASE("a single database passes through unchanged")
    {
        MergeTaxonVectors({db1}, t, out);
        CHECK(out == db1);
    }
    SUBCASE("three databases fold")
    {
        std::vector<taxid_t> db3 = {TestTaxonomy::B, 0, 0, AMBIGUOUS_SPAN_TAXON, MATE_PAIR_BORDER_TAXON, 0};
        MergeTaxonVectors({db1, db2, db3}, t, out);
        CHECK(out[0] == TestTaxonomy::ROOT);  // A and B meet at root
        CHECK(out[5] == TestTaxonomy::B);
    }
}

TEST_CASE("Merged hit counts resolve to the call kraken2's merge would make")
{
    // Read hits A1 in db1 and A2 in db2 at the same positions: merged
    // positions are A, so the call is A. With confidence 0, as merge uses.
    TestTaxonomy tax;
    const auto &t = tax.get();
    std::vector<taxid_t> db1(10, TestTaxonomy::A1), db2(10, TestTaxonomy::A2), out;
    MergeTaxonVectors({db1, db2}, t, out);
    taxon_counts_t hits;
    for (auto x : out) if (x) hits[x]++;
    CHECK(ResolveTree(hits, t, out.size(), 0.0) == TestTaxonomy::A);
}
