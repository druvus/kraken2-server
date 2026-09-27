#include "report_server.h"
#include "reports.h"

namespace kraken2
{
  void ReportKrakenStyle(std::ostringstream &ss, bool report_zeros, bool report_kmer_data,
                         Taxonomy &taxonomy, taxon_counters_t &call_counters, uint64_t total_seqs,
                         uint64_t total_unclassified)
  {
    taxon_counters_t clade_counters = GetCladeCounters(taxonomy, call_counters);

    ss << "% of Seqs\tClades\tTaxonomies\t";
    if (report_kmer_data)
    {
      ss << "Kmers\tDistinct Kmers\t";
    }
    ss << "Rank\tTaxonomy ID\tScientific Name\n";

    // Special handling of the unclassified sequences
    if (total_unclassified != 0 || report_zeros)
    {
      READCOUNTER rc(total_unclassified, 0);
      PrintKrakenStyleReportLine(ss, report_kmer_data, total_seqs, rc,
                                 rc, "U", 0, "unclassified", 0);
    }
    // DFS through the taxonomy, printing nodes as encountered
    KrakenReportDFS(1, ss, report_zeros, report_kmer_data, taxonomy,
                    clade_counters, call_counters, total_seqs, 'R', -1, 0);
  }
}
