#include "merged_taxonomy.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>

#include <unistd.h>

namespace kraken2server
{

using kraken2::Taxonomy;
using kraken2::TaxonomyNode;
using kraken2::taxid_t;

std::map<uint64_t, TaxonRecord> CollectTaxa(const Taxonomy &taxonomy)
{
    std::map<uint64_t, TaxonRecord> records;
    const TaxonomyNode *nodes = taxonomy.nodes();
    for (size_t i = 1; i < taxonomy.node_count(); i++)
    {
        TaxonRecord r;
        r.external_id = nodes[i].external_id;
        r.parent_external_id = nodes[i].parent_id ? nodes[nodes[i].parent_id].external_id : 0;
        r.name = taxonomy.name_data() + nodes[i].name_offset;
        r.rank = taxonomy.rank_data() + nodes[i].rank_offset;
        records[r.external_id] = r;
    }
    return records;
}

void MergeTaxa(std::map<uint64_t, TaxonRecord> &merged,
               const std::map<uint64_t, TaxonRecord> &records,
               const std::string &source_name)
{
    for (const auto &kv : records)
    {
        auto it = merged.find(kv.first);
        if (it == merged.end())
        {
            merged[kv.first] = kv.second;
        }
        else if (it->second.parent_external_id != kv.second.parent_external_id)
        {
            throw std::runtime_error(
                "Inconsistent taxonomies: taxid " + std::to_string(kv.first) +
                " has parent " + std::to_string(it->second.parent_external_id) +
                " in an earlier database but parent " +
                std::to_string(kv.second.parent_external_id) + " in " + source_name);
        }
    }
}

size_t WriteKrakenTaxonomy(const std::map<uint64_t, TaxonRecord> &records,
                           const std::string &path)
{
    // children per external id; the root is the node whose parent is 0 or
    // itself
    std::map<uint64_t, std::vector<uint64_t>> children;
    std::vector<uint64_t> roots;
    for (const auto &kv : records)
    {
        uint64_t parent = kv.second.parent_external_id;
        if (parent == 0 || parent == kv.first)
            roots.push_back(kv.first);
        else
            children[parent].push_back(kv.first);
    }
    if (roots.size() != 1)
    {
        throw std::runtime_error("Merged taxonomy must have exactly one root, found " +
                                 std::to_string(roots.size()));
    }

    // Breadth-first order gives every node a contiguous block of children
    // placed after the parent, which is what the Kraken taxonomy needs.
    std::vector<uint64_t> order;         // internal id - 1 -> external id
    std::map<uint64_t, uint64_t> internal;  // external id -> internal id
    order.push_back(roots[0]);
    internal[roots[0]] = 1;
    for (size_t i = 0; i < order.size(); i++)
    {
        auto it = children.find(order[i]);
        if (it == children.end()) continue;
        for (uint64_t child : it->second)
        {
            internal[child] = order.size() + 1;
            order.push_back(child);
        }
    }
    if (order.size() != records.size())
    {
        throw std::runtime_error("Merged taxonomy has nodes not reachable from the root");
    }

    std::string names, ranks;
    std::map<std::string, uint64_t> rank_offsets;
    std::vector<TaxonomyNode> nodes(order.size() + 1);
    memset(&nodes[0], 0, sizeof(TaxonomyNode));
    for (size_t i = 0; i < order.size(); i++)
    {
        const TaxonRecord &r = records.at(order[i]);
        TaxonomyNode &n = nodes[i + 1];
        memset(&n, 0, sizeof(n));
        n.external_id = r.external_id;
        n.parent_id = (i == 0) ? 0 : internal[r.parent_external_id];
        auto ch = children.find(r.external_id);
        if (ch != children.end() && !ch->second.empty())
        {
            n.first_child = internal[ch->second.front()];
            n.child_count = ch->second.size();
        }
        n.name_offset = names.size();
        names.append(r.name);
        names.push_back('\0');
        auto ro = rank_offsets.find(r.rank);
        if (ro == rank_offsets.end())
        {
            rank_offsets[r.rank] = ranks.size();
            ranks.append(r.rank);
            ranks.push_back('\0');
            ro = rank_offsets.find(r.rank);
        }
        n.rank_offset = ro->second;
    }

    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot write merged taxonomy to " + path);
    const char magic[] = "K2TAXDAT";
    uint64_t node_count = nodes.size(), name_len = names.size(), rank_len = ranks.size();
    f.write(magic, strlen(magic));
    f.write((const char *) &node_count, sizeof(node_count));
    f.write((const char *) &name_len, sizeof(name_len));
    f.write((const char *) &rank_len, sizeof(rank_len));
    f.write((const char *) nodes.data(), sizeof(TaxonomyNode) * node_count);
    f.write(names.data(), name_len);
    f.write(ranks.data(), rank_len);
    if (!f) throw std::runtime_error("Error writing merged taxonomy to " + path);
    return order.size();
}

std::unique_ptr<Taxonomy> BuildMergedTaxonomy(
    const std::vector<const Taxonomy *> &taxonomies,
    const std::vector<std::string> &names)
{
    std::map<uint64_t, TaxonRecord> merged;
    for (size_t i = 0; i < taxonomies.size(); i++)
    {
        MergeTaxa(merged, CollectTaxa(*taxonomies[i]),
                  i < names.size() ? names[i] : std::to_string(i));
    }

    char path[] = "/tmp/k2server_merged_taxo_XXXXXX";
    int fd = mkstemp(path);
    if (fd < 0) throw std::runtime_error("Cannot create temporary file for merged taxonomy");
    close(fd);
    WriteKrakenTaxonomy(merged, path);
    std::unique_ptr<Taxonomy> result(new Taxonomy(std::string(path), false));
    unlink(path);
    result->GenerateExternalToInternalIDMap();
    return result;
}

std::vector<taxid_t> MapToMerged(const Taxonomy &source, const Taxonomy &merged)
{
    std::vector<taxid_t> map(source.node_count(), 0);
    for (size_t i = 1; i < source.node_count(); i++)
    {
        map[i] = merged.GetInternalID(source.nodes()[i].external_id);
    }
    return map;
}

} // namespace kraken2server
