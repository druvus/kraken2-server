#pragma once

// Build one Kraken taxonomy from several loaded taxo.k2d taxonomies, for
// multi-database classification. See docs/MULTI_DB.md.

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "kraken2_data.h"
#include "taxonomy.h"

namespace kraken2server
{
    struct TaxonRecord
    {
        uint64_t external_id = 0;
        uint64_t parent_external_id = 0;  // 0 for the root
        std::string name;
        std::string rank;

        bool operator==(const TaxonRecord &o) const
        {
            return external_id == o.external_id && parent_external_id == o.parent_external_id &&
                   name == o.name && rank == o.rank;
        }
    };

    // All nodes of a loaded taxonomy as records keyed by external id.
    std::map<uint64_t, TaxonRecord> CollectTaxa(const kraken2::Taxonomy &taxonomy);

    // Add records to a union. Throws std::runtime_error if an external id
    // already present has a different parent (inconsistent taxonomies).
    void MergeTaxa(std::map<uint64_t, TaxonRecord> &merged,
                   const std::map<uint64_t, TaxonRecord> &records,
                   const std::string &source_name);

    // Write records as a Kraken taxonomy file (taxo.k2d format): nodes laid
    // out parents first with each node's children in one contiguous block.
    // Returns the number of nodes written, excluding the unused node 0.
    size_t WriteKrakenTaxonomy(const std::map<uint64_t, TaxonRecord> &records,
                               const std::string &path);

    // Union of the given taxonomies written to a temporary file and loaded,
    // with the external to internal id map generated.
    std::unique_ptr<kraken2::Taxonomy> BuildMergedTaxonomy(
        const std::vector<const kraken2::Taxonomy *> &taxonomies,
        const std::vector<std::string> &names);

    // For each internal id of `source`, the internal id of the same external
    // id in `merged` (0 for node 0).
    std::vector<kraken2::taxid_t> MapToMerged(const kraken2::Taxonomy &source,
                                              const kraken2::Taxonomy &merged);
}
