#ifndef KRAKEN2_SERVER_REPORTS_H_
#define KRAKEN2_SERVER_REPORTS_H_

#include <sstream>

#include "kraken2_headers.h"
#include "taxonomy.h"
#include "kraken2_data.h"
#include "readcounts.h"

namespace kraken2
{
    // Write a kraken2 style report to a string stream. The tree walk and
    // line formatting are the upstream functions from reports.cc; this
    // wrapper exists because upstream ReportKrakenStyle writes to a file.
    // A header line naming the columns is written first.
    void ReportKrakenStyle(std::ostringstream &ss, bool report_zeros, bool report_kmer_data,
                           Taxonomy &taxonomy, taxon_counters_t &call_counters, uint64_t total_seqs,
                           uint64_t total_unclassified);
}
#endif
