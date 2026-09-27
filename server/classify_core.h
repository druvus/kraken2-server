#pragma once

// Pure helpers used by the classifier. They have no dependency on gRPC or
// the server state, so they can be unit tested directly.

#include <cstdint>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <vector>

#include "kraken2_data.h"
#include "taxonomy.h"

namespace kraken2server
{
    using kraken2::taxid_t;
    using kraken2::Taxonomy;
    using kraken2::taxon_counts_t;
    using kraken2::TAXID_MAX;

    // Marker values inserted into the per-read taxon list, as in kraken2.
    static const taxid_t AMBIGUOUS_SPAN_TAXON = TAXID_MAX - 2;
    static const taxid_t MATE_PAIR_BORDER_TAXON = TAXID_MAX;
    static const taxid_t READING_FRAME_BORDER_TAXON = TAXID_MAX - 1;

    struct ClassificationStats
    {
        uint64_t total_sequences;
        uint64_t total_bases;
        uint64_t total_classified;
    };

    // Remove a trailing /1 or /2 from a read name.
    std::string TrimPairInfo(const std::string &id);

    // Format the per-read minimizer hit list, e.g. "562:13 A:2 |:| 0:4".
    void AddHitlistString(std::ostringstream &oss, const std::vector<taxid_t> &taxa,
                          const Taxonomy &taxonomy);

    // Choose the taxon for a read from its minimizer hit counts, applying the
    // confidence threshold. Adapted from kraken2 classify.cc. hit_counts is
    // not const because the upstream algorithm reads through operator[].
    taxid_t ResolveTree(taxon_counts_t &hit_counts, const Taxonomy &taxonomy,
                        size_t total_minimizers, double confidence_threshold);

    // Combine the per-database taxon vectors of one read into one vector in
    // the merged taxonomy: position by position, the lowest common ancestor
    // over databases, where LCA(x, 0) = x. Marker values pass through. All
    // vectors must have the same length. Mirrors kraken2's merge program.
    void MergeTaxonVectors(const std::vector<std::vector<taxid_t>> &per_database,
                           const Taxonomy &merged, std::vector<taxid_t> &out);

    std::string DoubleStatToString(double d, int precision);

    // Human readable summary of one stream, kraken2 style.
    std::string ReportStats(struct timeval time1, struct timeval time2,
                            const ClassificationStats &stats);

    // Human readable summary of the server history.
    std::string ReportTotalStats(const ClassificationStats &stats);
}
