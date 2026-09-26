#include "classify_core.h"

#include <cmath>
#include <iomanip>

namespace kraken2server
{

std::string TrimPairInfo(const std::string &id)
{
    size_t sz = id.size();
    if (sz <= 2)
        return id;
    if (id[sz - 2] == '/' && (id[sz - 1] == '1' || id[sz - 1] == '2'))
        return id.substr(0, sz - 2);
    return id;
}

// Adapted from kraken2 classify.cc.
void AddHitlistString(std::ostringstream &oss, const std::vector<taxid_t> &taxa,
                      const Taxonomy &taxonomy)
{
    auto last_code = taxa[0];
    auto code_count = 1;

    for (size_t i = 1; i < taxa.size(); i++)
    {
        auto code = taxa[i];

        if (code == last_code)
        {
            code_count += 1;
        }
        else
        {
            if (last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON)
            {
                if (last_code == AMBIGUOUS_SPAN_TAXON)
                {
                    oss << "A:" << code_count << " ";
                }
                else
                {
                    auto ext_code = taxonomy.nodes()[last_code].external_id;
                    oss << ext_code << ":" << code_count << " ";
                }
            }
            else
            { // mate pair/reading frame marker
                oss << (last_code == MATE_PAIR_BORDER_TAXON ? "|:| " : "-:- ");
            }
            code_count = 1;
            last_code = code;
        }
    }
    if (last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON)
    {
        if (last_code == AMBIGUOUS_SPAN_TAXON)
        {
            oss << "A:" << code_count << " ";
        }
        else
        {
            auto ext_code = taxonomy.nodes()[last_code].external_id;
            oss << ext_code << ":" << code_count;
        }
    }
    else
    { // mate pair/reading frame marker
        oss << (last_code == MATE_PAIR_BORDER_TAXON ? "|:|" : "-:-");
    }
}

// Adapted from kraken2 classify.cc.
taxid_t ResolveTree(taxon_counts_t &hit_counts, const Taxonomy &taxonomy,
                    size_t total_minimizers, double confidence_threshold)
{
    taxid_t max_taxon = 0;
    uint32_t max_score = 0;
    uint32_t required_score = ceil(confidence_threshold * total_minimizers);

    // Sum each taxon's LTR path, find taxon with highest LTR score
    for (auto &kv_pair : hit_counts)
    {
        taxid_t taxon = kv_pair.first;
        uint32_t score = 0;

        for (auto &kv_pair2 : hit_counts)
        {
            taxid_t taxon2 = kv_pair2.first;

            if (taxonomy.IsAAncestorOfB(taxon2, taxon))
            {
                score += kv_pair2.second;
            }
        }

        if (score > max_score)
        {
            max_score = score;
            max_taxon = taxon;
        }
        else if (score == max_score)
        {
            max_taxon = taxonomy.LowestCommonAncestor(max_taxon, taxon);
        }
    }

    // Reset max. score to be only hits at the called taxon
    max_score = hit_counts[max_taxon];
    // We probably have a call w/o required support (unless LCA resolved tie)
    while (max_taxon && max_score < required_score)
    {
        max_score = 0;
        for (auto &kv_pair : hit_counts)
        {
            taxid_t taxon = kv_pair.first;
            // Add to score if taxon in max_taxon's clade
            if (taxonomy.IsAAncestorOfB(max_taxon, taxon))
            {
                max_score += kv_pair.second;
            }
        }
        // Score is now sum of hits at max_taxon and w/in max_taxon clade
        if (max_score >= required_score)
            // Kill loop and return, we've got enough support here
            return max_taxon;
        else
            // Run up tree until confidence threshold is met
            // Run off tree if required score isn't met
            max_taxon = taxonomy.nodes()[max_taxon].parent_id;
    }

    return max_taxon;
}

void MergeTaxonVectors(const std::vector<std::vector<taxid_t>> &per_database,
                       const Taxonomy &merged, std::vector<taxid_t> &out)
{
    out.clear();
    if (per_database.empty()) return;
    const size_t n = per_database[0].size();
    out.reserve(n);
    for (size_t i = 0; i < n; i++)
    {
        taxid_t t = per_database[0][i];
        if (t == AMBIGUOUS_SPAN_TAXON || t == MATE_PAIR_BORDER_TAXON ||
            t == READING_FRAME_BORDER_TAXON)
        {
            out.push_back(t);
            continue;
        }
        for (size_t d = 1; d < per_database.size(); d++)
        {
            t = merged.LowestCommonAncestor(t, per_database[d][i]);
        }
        out.push_back(t);
    }
}

std::string DoubleStatToString(double d, int precision)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << d;
    return stream.str();
}

std::string ReportStats(struct timeval time1, struct timeval time2,
                        const ClassificationStats &stats)
{
    time2.tv_usec -= time1.tv_usec;
    time2.tv_sec -= time1.tv_sec;
    if (time2.tv_usec < 0)
    {
        time2.tv_sec--;
        time2.tv_usec += 1000000;
    }
    double seconds = time2.tv_usec;
    seconds /= 1e6;
    seconds += time2.tv_sec;

    uint64_t total_unclassified = stats.total_sequences - stats.total_classified;
    // Guard against empty streams and zero elapsed time.
    double denom_seqs = stats.total_sequences > 0 ? stats.total_sequences : 1.0;
    double minutes = seconds > 0 ? seconds / 60 : 1.0 / 60;

    return std::to_string(stats.total_sequences) + " sequences (" + DoubleStatToString(stats.total_bases / 1.0e6, 2) + " Mbp) processed in " + DoubleStatToString(seconds, 2) + "s (" + DoubleStatToString(stats.total_sequences / 1.0e3 / minutes, 2) + " Kseq/m, " + DoubleStatToString(stats.total_bases / 1.0e6 / minutes, 2) + " Mbp/m).\n" +
           "\t" + std::to_string(stats.total_classified) + " sequences classified (" + DoubleStatToString(stats.total_classified * 100.0 / denom_seqs, 2) + "%)\n" +
           "\t" + std::to_string(total_unclassified) + " sequences unclassified (" + DoubleStatToString(total_unclassified * 100.0 / denom_seqs, 2) + "%)\n";
}

std::string ReportTotalStats(const ClassificationStats &stats)
{
    uint64_t total_unclassified = stats.total_sequences - stats.total_classified;
    double denom_seqs = stats.total_sequences > 0 ? stats.total_sequences : 1.0;

    return std::to_string(stats.total_sequences) + " sequences (" + DoubleStatToString(stats.total_bases / 1.0e6, 2) + " Mbp) processed.\n" +
           std::to_string(stats.total_classified) + " sequences classified (" + DoubleStatToString(stats.total_classified * 100.0 / denom_seqs, 2) + "%).\n" +
           std::to_string(total_unclassified) + " sequences unclassified (" + DoubleStatToString(total_unclassified * 100.0 / denom_seqs, 2) + "%).\n";
}

} // namespace kraken2server
