#include <condition_variable>
#include <map>
#include <fstream>
#include <getopt.h>
#include <memory>
#include <thread>
#include <sysexits.h>

#include "classify_server.h"
#include "messages.h"

using kraken2server::AddHitlistString;
using kraken2server::ResolveTree;
using kraken2server::TrimPairInfo;
using kraken2server::ReportStats;
using kraken2server::ReportTotalStats;

using namespace std::chrono_literals; // ns, us, ms, s, h, etc.

// Minimizer token stream, adapted from kraken2 classify.cc. Minimizers are
// first tokenised, then looked up in one batch, then replayed to build the
// per-read taxon list.
enum MinTokKind { TOK_LOOKUP, TOK_SKIP, TOK_REPEAT, TOK_AMBIG,
                  TOK_BORDER_MATE, TOK_BORDER_FRAME };
struct MinToken { uint8_t kind; uint32_t key_idx; };

Kraken2ServerClassifier::Kraken2ServerClassifier(Options &options)
        : opts(options), pool(opts.thread_pool < 0 ? 0 : (size_t) opts.thread_pool) {
    // The pool handles classification tasks; the index loads in the
    // background.
    std::cout << "Created classification thread pool with "
              << pool.thread_count() << " thread(s)." << std::endl;
    loader = std::thread([this]() { LoadIndex(); });
}

Kraken2ServerClassifier::~Kraken2ServerClassifier(){
    // The loader touches members, so it must finish before they are
    // destroyed. Loading cannot be interrupted, so a shutdown during load
    // waits for it to complete.
    if (loader.joinable()) {
        loader.join();
    }
}

std::string Kraken2ServerClassifier::GetSummary() {
    std::lock_guard<std::mutex> lock(stats_mtx);
    return summary;
}

static std::string Basename(const std::string &path) {
    std::string p = path;
    while (p.size() > 1 && p.back() == '/') p.pop_back();
    auto slash = p.find_last_of('/');
    return slash == std::string::npos ? p : p.substr(slash + 1);
}

Index Kraken2ServerClassifier::LoadOne(const std::string &path) {
    Index idx;
    idx.path = path;
    idx.name = Basename(path);
    std::string options_filename = path + "/opts.k2d";
    std::string taxonomy_filename = path + "/taxo.k2d";
    std::string index_filename = path + "/hash.k2d";

    idx.options = {0};
    ifstream idx_opt_fs(options_filename);
    struct stat sb;
    if (stat(options_filename.c_str(), &sb) < 0)
        throw std::runtime_error("Unable to get filesize of " + options_filename);
    if (!idx_opt_fs.read((char *)&idx.options, sb.st_size))
        throw std::runtime_error("Unable to read " + options_filename);

    std::cerr << "[" << idx.name << "] Loading taxonomy..." << std::endl;
    idx.taxonomy.reset(new Taxonomy(taxonomy_filename, opts.use_memory_mapping));

    std::cerr << "[" << idx.name << "] Loading hash table..." << std::endl;
    switch (GetKVStoreCellType(index_filename)) {
    case CompactHash32:
        idx.hash.reset(new CompactHashTable<CompactHashCell>(index_filename, opts.use_memory_mapping));
        break;
    case CompactHash40:
        idx.hash.reset(new CompactHashTable<CompactHashCell40>(index_filename, opts.use_memory_mapping));
        break;
    default:
        throw std::runtime_error("Unable to determine width of compact hash cell in " + index_filename);
    }
    return idx;
}

void Kraken2ServerClassifier::LoadIndex() {
    index_available = false;
    std::cerr << "Loading database information..." << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(opts.wait));

    try {
        indexes.clear();
        for (const std::string &path : opts.db_paths) {
            indexes.push_back(LoadOne(path));
        }
        if (indexes.empty())
            throw std::runtime_error("No database given.");

        // All databases must share the scanner parameters, since one
        // minimizer stream serves them all.
        idx_opts = indexes[0].options;
        min_hash_any = idx_opts.minimum_acceptable_hash_value;
        for (size_t i = 1; i < indexes.size(); i++) {
            const IndexOptions &o = indexes[i].options;
            if (o.k != idx_opts.k || o.l != idx_opts.l ||
                o.spaced_seed_mask != idx_opts.spaced_seed_mask ||
                o.toggle_mask != idx_opts.toggle_mask || o.dna_db != idx_opts.dna_db ||
                o.revcom_version != idx_opts.revcom_version) {
                throw std::runtime_error(
                    "Database " + indexes[i].name + " has different k-mer, minimizer or "
                    "sequence type settings from " + indexes[0].name +
                    "; databases classified together must be built alike.");
            }
            min_hash_any = std::min(min_hash_any, o.minimum_acceptable_hash_value);
        }

        if (opts.use_translated_search && idx_opts.dna_db) {
            std::cerr << "Warning: --translated-search requested but the database "
                      << "is nucleotide; using nucleotide search." << std::endl;
        }
        opts.use_translated_search = !idx_opts.dna_db;
        if (opts.use_translated_search) {
            // kraken2 >= 2.1.3 initialises the codon lookup tables explicitly
            // once at startup rather than lazily in TranslateToAllFrames.
            // Without this call every codon translates to 'K'.
            initLookUpTables();
            std::cerr << "Protein database: using translated search." << std::endl;
        }

        if (indexes.size() == 1) {
            taxonomy = indexes[0].taxonomy.get();
            indexes[0].to_merged.resize(taxonomy->node_count());
            for (size_t i = 0; i < indexes[0].to_merged.size(); i++)
                indexes[0].to_merged[i] = i;
        } else {
            std::cerr << "Building merged taxonomy from " << indexes.size()
                      << " databases..." << std::endl;
            std::vector<const Taxonomy *> taxonomies;
            std::vector<std::string> names;
            for (const Index &idx : indexes) {
                taxonomies.push_back(idx.taxonomy.get());
                names.push_back(idx.name);
            }
            merged_taxonomy = kraken2server::BuildMergedTaxonomy(taxonomies, names);
            taxonomy = merged_taxonomy.get();
            for (Index &idx : indexes)
                idx.to_merged = kraken2server::MapToMerged(*idx.taxonomy, *taxonomy);
            std::cerr << "Merged taxonomy has " << taxonomy->node_count() - 1
                      << " taxa." << std::endl;
            if (opts.confidence_threshold > 0 || opts.minimum_hit_groups != 2) {
                std::cerr << "Note: with several databases the merged call follows "
                          << "kraken2's merge program: --confidence and --hit-groups "
                          << "are not applied to it (see docs/MULTI_DB.md)." << std::endl;
            }
        }
    }
    catch (const std::exception &ex) {
        std::cerr << "Unable to load index"
                  << ": " << ex.what() << std::endl;
        index_broken = true;
        return;
    }
    std::cerr << "Successfully loaded " << indexes.size() << " index"
              << (indexes.size() == 1 ? "" : "es") << "." << std::endl;
    index_available = true;
}

// Drain classified batches onto the gRPC stream until the queue is closed
// and empty. Batches finish in any order on the pool; they are held here
// and written in input order, as kraken2 does, so the client sees results
// in the same order as its reads. Runs on its own thread so writes overlap
// with classification.
static void ResultsHandler(
        ServerStream *stream,
        taxon_counters_t &stream_taxon_counters,
        ClassificationStats &stream_stats,
        ThreadSafeQueue<BatchResults> &results_queue) {
    std::map<uint64_t, BatchResults> pending;
    uint64_t next_sequence = 0;
    while (std::optional<BatchResults> res = results_queue.pop_wait()) {
        pending.emplace(res->sequence, std::move(*res));
        for (auto it = pending.find(next_sequence); it != pending.end();
             it = pending.find(++next_sequence)) {
            BatchResults &batch = it->second;
            // The client is configured to receive messages up to INT_MAX, and
            // the client limits request batches to 128 MB, so a batch of
            // results always fits.
            Kraken2SequenceStreamResult result;
            result.mutable_classifications()->Swap(&batch.k2results);
            stream->Write(result, WriteOptions().set_buffer_hint());
            // update stats for the stream
            stream_stats.total_bases += batch.stats.total_bases;
            stream_stats.total_classified += batch.stats.total_classified;
            stream_stats.total_sequences += batch.stats.total_sequences;
            // update taxon_counters for the stream
            for (auto &kv_pair : batch.taxon_counters) {
                stream_taxon_counters[kv_pair.first] += std::move(kv_pair.second);
            }
            pending.erase(it);
        }
    }
    if (!pending.empty()) {
        std::cerr << "Warning: " << pending.size()
                  << " result batch(es) never reached the writer." << std::endl;
    }
}

void Kraken2ServerClassifier::ProcessSequenceStream(
        ServerContext *context, ServerStream *stream, std::string &results) {
    std::cerr << "Starting stream handler." << std::endl;
    stream->SendInitialMetadata();

    // Stats for the whole stream
    taxon_counters_t stream_taxon_counters;
    ClassificationStats stream_stats = {0, 0, 0};

    struct timeval tv1, tv2;
    gettimeofday(&tv1, nullptr);

    // Queue and thread that aggregate batch results and write them to the
    // output stream.
    ThreadSafeQueue<BatchResults> results_queue;
    std::thread results_thread(ResultsHandler,
        stream, std::ref(stream_taxon_counters), std::ref(stream_stats),
        std::ref(results_queue));

    // Backpressure: limit the number of request batches held by this stream
    // (queued in the pool or being classified) so a fast client cannot
    // make the server buffer its whole input.
    const size_t max_in_flight = std::max<size_t>(4, 2 * pool.thread_count());
    std::mutex in_flight_mtx;
    std::condition_variable in_flight_cv;
    size_t in_flight = 0;

    // Classify while reads are still being received on the input stream
    uint64_t sequence = 0;
    while (!context->IsCancelled()) {
        // Each batch is owned by a shared_ptr so the pool can copy the task
        // object without copying the message.
        auto req = std::make_shared<Kraken2SequenceRequestMulti>();
        if (!stream->Read(req.get())) {
            break;
        }
        {
            std::unique_lock<std::mutex> lock(in_flight_mtx);
            in_flight_cv.wait(lock, [&] { return in_flight < max_in_flight; });
            in_flight++;
        }
        const uint64_t this_sequence = sequence++;
        pool.push([this, req, this_sequence, &results_queue, &in_flight_mtx, &in_flight_cv, &in_flight]() {
            try {
                ProcessBatch(*req, this_sequence, results_queue);
            }
            catch (const std::exception &ex) {
                std::cerr << "Error classifying batch: " << ex.what() << std::endl;
            }
            std::lock_guard<std::mutex> lock(in_flight_mtx);
            in_flight--;
            in_flight_cv.notify_one();
        });
    }

    // Wait for the outstanding batches, then let the results thread drain
    // the queue and finish.
    {
        std::unique_lock<std::mutex> lock(in_flight_mtx);
        in_flight_cv.wait(lock, [&] { return in_flight == 0; });
    }
    results_queue.close();
    results_thread.join();

    gettimeofday(&tv2, nullptr);
    // generate the report, and update servers total history
    GenerateReport(
        results, summary, opts, *taxonomy, tv1, tv2, stream_stats, total_stats,
        stream_taxon_counters, total_taxon_counters, stats_mtx);

    std::cerr << "Finished stream handler." << std::endl;
}

void Kraken2ServerClassifier::ProcessBatch(
    const Kraken2SequenceRequestMulti &reqs, uint64_t sequence,
    ThreadSafeQueue<BatchResults> &result_q) {

    MinimizerScanner scanner(
        idx_opts.k, idx_opts.l, idx_opts.spaced_seed_mask,
        idx_opts.dna_db, idx_opts.toggle_mask,
        idx_opts.revcom_version);
    vector<taxid_t> taxa;
    taxon_counts_t hit_counts;
    vector<string> translated_frames(6);

    BatchResults results = BatchResults();
    results.sequence = sequence;
    results.k2results.mutable_classes()->Reserve(reqs.seqs_size());

    kraken2::Sequence seq, seq2;
    for (auto &req : reqs.seqs()) {
        bool paired = SequenceRequestToPair(req, seq, seq2);
        // A pair counts as one fragment, as in kraken2.
        results.stats.total_sequences++;
        results.stats.total_bases += seq.seq.size();
        if (paired)
            results.stats.total_bases += seq2.seq.size();
        if (opts.minimum_quality_score > 0) {
            bool ok = MaskLowQualityBases(seq, opts.minimum_quality_score);
            if (paired)
                ok = MaskLowQualityBases(seq2, opts.minimum_quality_score) && ok;
            if (!ok) {
                // Malformed record from the client. Report it as
                // unclassified rather than terminating the server.
                results.k2results.mutable_classes()->Add(
                    UnclassifiedResult(seq, paired ? &seq2 : nullptr));
                continue;
            }
        }

        Kraken2SequenceResult classification = ClassifySequence(
            seq, paired ? &seq2 : nullptr, results.stats, scanner,
            taxa, hit_counts, translated_frames, results.taxon_counters);

        results.k2results.mutable_classes()->Add(std::move(classification));
    }

    result_q.push(std::move(results));
}

Kraken2SequenceResult Kraken2ServerClassifier::UnclassifiedResult(
    const Sequence &dna, const Sequence *dna2)
{
    Kraken2SequenceResult result;
    std::string id = dna.header;
    result.set_id(dna2 != nullptr ? TrimPairInfo(id) : id);
    result.set_classified(false);
    result.set_size(dna.seq.size());
    if (dna2 != nullptr)
    {
        result.set_paired(true);
        result.set_size2(dna2->seq.size());
    }
    result.set_hitlist("0:0");
    return result;
}

////////////////////////////////
// The following methods are adapted from the Kraken2 source code
// (classify.cc). Quick mode and text output have been removed; paired
// reads are handled per fragment rather than through a global option.
////////////////////////////////

Kraken2SequenceResult Kraken2ServerClassifier::ClassifySequence(
    Sequence &dna, Sequence *dna2,
    ClassificationStats &stats, MinimizerScanner &scanner,
    vector<taxid_t> &taxa, taxon_counts_t &hit_counts,
    vector<string> &tx_frames, taxon_counters_t &curr_taxon_counts)
{
    Taxonomy &taxonomy = *this->taxonomy;
    const bool multi_db = indexes.size() > 1;
    uint64_t *minimizer_ptr;
    taxid_t call = 0;
    taxa.clear();
    hit_counts.clear();
    const bool paired = dna2 != nullptr;
    auto frame_ct = opts.use_translated_search ? 6 : 1;
    int64_t minimizer_hit_groups = 0;

    // Pool threads are long lived, so thread_local scratch space avoids
    // reallocating per read.
    static thread_local std::vector<uint64_t> lookup_keys;
    static thread_local std::vector<MinToken> tok_stream;
    lookup_keys.clear();
    tok_stream.clear();

    // Pass 1: tokenise minimizers from each mate and each frame.
    for (int mate_num = 0; mate_num < 2; mate_num++)
    {
        if (mate_num == 1 && !paired)
            break;
        Sequence &cur = mate_num == 0 ? dna : *dna2;

        if (opts.use_translated_search)
        {
            TranslateToAllFrames(cur.seq, tx_frames);
        }
        // index of frame is 0 - 5 w/ tx search (or 0 if no tx search)
        for (int frame_idx = 0; frame_idx < frame_ct; frame_idx++)
        {
            if (opts.use_translated_search)
            {
                scanner.LoadSequence(tx_frames[frame_idx]);
            }
            else
            {
                scanner.LoadSequence(cur.seq);
            }
            uint64_t last_minimizer = UINT64_MAX;
            while ((minimizer_ptr = scanner.NextMinimizer()) != nullptr)
            {
                if (scanner.is_ambiguous())
                {
                    tok_stream.push_back({TOK_AMBIG, 0});
                }
                else if (*minimizer_ptr != last_minimizer)
                {
                    last_minimizer = *minimizer_ptr;
                    bool skip_lookup = min_hash_any &&
                        MurmurHash3(*minimizer_ptr) < min_hash_any;
                    if (skip_lookup)
                    {
                        tok_stream.push_back({TOK_SKIP, 0});
                    }
                    else
                    {
                        tok_stream.push_back({TOK_LOOKUP, (uint32_t) lookup_keys.size()});
                        lookup_keys.push_back(*minimizer_ptr);
                    }
                }
                else
                {
                    tok_stream.push_back({TOK_REPEAT, 0});
                }
            }
            if (opts.use_translated_search && frame_idx != 5)
                tok_stream.push_back({TOK_BORDER_FRAME, 0});
        }
        if (paired && mate_num == 0)
            tok_stream.push_back({TOK_BORDER_MATE, 0});
    }

    // Pass 2: one batched hash lookup per database for all distinct
    // minimizers. Each database's values are mapped into the (merged)
    // taxonomy; a key below a database's own hash threshold is no hit there.
    static thread_local std::vector<std::vector<taxid_t>> db_taxa;
    static thread_local std::vector<hvalue_t> lookup_vals;
    db_taxa.resize(indexes.size());
    for (size_t d = 0; d < indexes.size(); d++)
    {
        const Index &idx = indexes[d];
        std::vector<taxid_t> &vals = db_taxa[d];
        vals.assign(lookup_keys.size(), 0);
        lookup_vals.resize(lookup_keys.size());
        if (!lookup_keys.empty())
            idx.hash->GetBatch(lookup_keys.data(), lookup_vals.data(), lookup_keys.size());
        const uint64_t min_hash = idx.options.minimum_acceptable_hash_value;
        for (size_t k = 0; k < lookup_keys.size(); k++)
        {
            if (lookup_vals[k] == 0) continue;
            if (multi_db && min_hash > min_hash_any && MurmurHash3(lookup_keys[k]) < min_hash)
                continue;
            vals[k] = idx.to_merged[lookup_vals[k]];
        }
    }
    // With several databases, combine per key: the LCA over databases in
    // the merged taxonomy, as kraken2's merge program does per position.
    std::vector<taxid_t> &merged_vals = db_taxa[0];
    for (size_t d = 1; d < indexes.size(); d++)
        for (size_t k = 0; k < lookup_keys.size(); k++)
            merged_vals[k] = taxonomy.LowestCommonAncestor(merged_vals[k], db_taxa[d][k]);

    // Pass 3: replay tokens to build the taxon list and hit counts.
    {
        taxid_t last_taxon = 0;
        for (size_t ti = 0; ti < tok_stream.size(); ti++)
        {
            const MinToken &tok = tok_stream[ti];
            taxid_t taxon = 0;
            switch (tok.kind)
            {
            case TOK_AMBIG:
                taxa.push_back(AMBIGUOUS_SPAN_TAXON);
                continue;
            case TOK_BORDER_FRAME:
                taxa.push_back(READING_FRAME_BORDER_TAXON);
                continue;
            case TOK_BORDER_MATE:
                taxa.push_back(MATE_PAIR_BORDER_TAXON);
                continue;
            case TOK_SKIP:
                taxon = 0;
                last_taxon = 0;
                break;
            case TOK_LOOKUP:
                taxon = merged_vals[tok.key_idx];
                last_taxon = taxon;
                if (taxon)
                {
                    // Increment only for a DB hit on a new minimizer.
                    minimizer_hit_groups++;
                    curr_taxon_counts[taxon].add_kmer(lookup_keys[tok.key_idx]);
                }
                break;
            default: // TOK_REPEAT
                taxon = last_taxon;
                if (taxon)
                {
                    curr_taxon_counts[taxon].increaseKmerCount(1);
                }
                break;
            }
            taxa.push_back(taxon);
            if (taxon)
            {
                hit_counts[taxon]++;
            }
        }
    }

    auto total_kmers = taxa.size();
    if (paired) // account for the mate pair marker
        total_kmers--;
    if (opts.use_translated_search) // account for reading frame markers
        total_kmers -= paired ? 4 : 2;
    if (multi_db)
    {
        // kraken2's merge program recomputes the call from the merged hit
        // lists with confidence 0 and no hit-group filter (docs/MULTI_DB.md).
        call = ResolveTree(hit_counts, taxonomy, total_kmers, 0.0);
    }
    else
    {
        call = ResolveTree(hit_counts, taxonomy, total_kmers, opts.confidence_threshold);
        // Void a call made by too few minimizer groups
        if (call && minimizer_hit_groups < opts.minimum_hit_groups)
            call = 0;
    }

    if (call)
    {
        stats.total_classified++;
        curr_taxon_counts[call].incrementReadCount();
    }

    Kraken2SequenceResult result;
    result.set_id(paired ? TrimPairInfo(dna.header) : dna.header);
    if (call)
    {
        result.set_classified(true);
        result.set_tax_id(taxonomy.nodes()[call].external_id);
        // The scientific name is not sent per read; the client does not
        // use it and it would add to every response.
    }
    else
        result.set_classified(false);
    result.set_size(dna.seq.size());
    if (paired)
    {
        result.set_paired(true);
        result.set_size2(dna2->seq.size());
    }
    if (taxa.empty())
        result.set_hitlist("0:0");
    else
    {
        std::ostringstream hitlist;
        AddHitlistString(hitlist, taxa, taxonomy);
        result.set_hitlist(hitlist.str());
    }

    return result;
}

bool Kraken2ServerClassifier::MaskLowQualityBases(Sequence &dna, int minimum_quality_score)
{
    if (dna.format != FORMAT_FASTQ)
        return true;
    if (dna.seq.size() != dna.quals.size())
    {
        std::cerr << dna.header << ": Sequence length (" << dna.seq.size()
                  << ") != Quality string length (" << dna.quals.size()
                  << "); reporting as unclassified." << std::endl;
        return false;
    }
    for (size_t i = 0; i < dna.seq.size(); i++)
    {
        if ((dna.quals[i] - '!') < minimum_quality_score)
            dna.seq[i] = 'x';
    }
    return true;
}

void Kraken2ServerClassifier::GenerateReport(
        std::string &results, std::string &summary, Options &opts, Taxonomy &taxonomy,
        timeval &tv1, timeval &tv2, ClassificationStats &stats,
        ClassificationStats &total_stats, taxon_counters_t &taxon_counters, taxon_counters_t &total_taxon_counters,
        std::mutex &stats_mtx)
{
    std::ostringstream ss;
    auto total_unclassified = stats.total_sequences - stats.total_classified;
    ReportKrakenStyle(ss,
                      opts.report_zero_counts,
                      opts.report_kmer_data,
                      taxonomy,
                      taxon_counters,
                      stats.total_sequences,
                      total_unclassified);
    results.assign(ss.str());

    std::cerr << ReportStats(tv1, tv2, stats) << std::endl;

    if (opts.stats)
    {
        std::lock_guard<std::mutex> lock(stats_mtx);

        total_stats.total_sequences += stats.total_sequences;
        total_stats.total_classified += stats.total_classified;
        total_stats.total_bases += stats.total_bases;
        for (auto &kv_pair : taxon_counters)
        {
            total_taxon_counters[kv_pair.first] += std::move(kv_pair.second);
        }
        ss.str(std::string());
        total_unclassified = total_stats.total_sequences - total_stats.total_classified;
        ReportKrakenStyle(ss,
                          opts.report_zero_counts,
                          opts.report_kmer_data,
                          taxonomy,
                          total_taxon_counters,
                          total_stats.total_sequences,
                          total_unclassified);

        ss << "\n"
           << ReportTotalStats(total_stats);
        summary.assign(ss.str());
    }
}
