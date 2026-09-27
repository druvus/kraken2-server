#pragma once

// A small taxonomy written to a temporary file in the kraken2 taxo.k2d
// format, for tests that need a kraken2::Taxonomy.
//
//   internal  external  name   parent
//   1         1         root   -
//   2         10        A      root
//   3         20        B      root
//   4         11        A1     A
//   5         12        A2     A
//
// Children of a node occupy a contiguous block of internal ids, as the
// kraken2 taxonomy builder guarantees.

#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include <unistd.h>

#include "taxonomy.h"

class TestTaxonomy
{
public:
    enum Node : uint64_t { ROOT = 1, A = 2, B = 3, A1 = 4, A2 = 5 };

    TestTaxonomy()
    {
        char tmpl[] = "/tmp/k2s_taxo_XXXXXX";
        int fd = mkstemp(tmpl);
        close(fd);
        path_ = tmpl;
        Write(path_);
        taxonomy_ = new kraken2::Taxonomy(path_, false);
    }

    ~TestTaxonomy()
    {
        delete taxonomy_;
        unlink(path_.c_str());
    }

    kraken2::Taxonomy &get() { return *taxonomy_; }

private:
    std::string path_;
    kraken2::Taxonomy *taxonomy_;

    static void Write(const std::string &path)
    {
        // name and rank data are NUL separated super strings
        std::string names, ranks;
        auto add = [](std::string &data, const char *s) {
            uint64_t off = data.size();
            data.append(s);
            data.push_back('\0');
            return off;
        };
        uint64_t n_root = add(names, "root"), n_a = add(names, "A"), n_b = add(names, "B"),
                 n_a1 = add(names, "A1"), n_a2 = add(names, "A2");
        uint64_t r_none = add(ranks, "no rank"), r_king = add(ranks, "superkingdom"),
                 r_species = add(ranks, "species");

        std::vector<kraken2::TaxonomyNode> nodes(6);
        // parent, first_child, child_count, name_offset, rank_offset, external_id, godparent
        nodes[0] = {0, 0, 0, 0, 0, 0, 0};
        nodes[ROOT] = {0, A, 2, n_root, r_none, 1, 0};
        nodes[A] = {ROOT, A1, 2, n_a, r_king, 10, 0};
        nodes[B] = {ROOT, 0, 0, n_b, r_king, 20, 0};
        nodes[A1] = {A, 0, 0, n_a1, r_species, 11, 0};
        nodes[A2] = {A, 0, 0, n_a2, r_species, 12, 0};

        std::ofstream f(path, std::ios::binary);
        const char magic[] = "K2TAXDAT";
        uint64_t node_count = nodes.size(), name_len = names.size(), rank_len = ranks.size();
        f.write(magic, strlen(magic));
        f.write((const char *)&node_count, sizeof(node_count));
        f.write((const char *)&name_len, sizeof(name_len));
        f.write((const char *)&rank_len, sizeof(rank_len));
        f.write((const char *)nodes.data(), sizeof(kraken2::TaxonomyNode) * node_count);
        f.write(names.data(), name_len);
        f.write(ranks.data(), rank_len);
    }
};
