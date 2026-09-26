#include <condition_variable>
#include <fstream>
#include <getopt.h>
#include <memory>
#include <thread>
#include <sysexits.h>

#include "classify_server.h"
#include "messages.h"

using namespace std::chrono_literals; // ns, us, ms, s, h, etc.

// Minimizer token stream, adapted from kraken2 classify.cc. Minimizers are
// first tokenised, then looked up in one batch, then replayed to build the
// per-read taxon list.
enum MinTokKind { TOK_LOOKUP, TOK_SKIP, TOK_REPEAT, TOK_AMBIG,
                  TOK_BORDER_MATE, TOK_BORDER_FRAME };
struct MinToken { uint8_t kind; uint32_t key_idx; };


Kraken2ServerClassifier::Kraken2ServerClassifier(Options &options)
        : opts(options) {
    // start a thread pool to handle classification tasks and
    // start loading the index in the background.
    pool.reset(opts.thread_pool);
    std::cout << "Created classification thread pool with "
              << pool.get_thread_count() << " thread(s)." << std::endl;
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


void Kraken2ServerClassifier::LoadIndex() {
    index_available = false;
    std::cerr << "Loading database information..." << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(opts.wait));

    try {
        idx_opts = {0};
        ifstream idx_opt_fs(opts.options_filename);
        struct stat sb;
        if (stat(opts.options_filename.c_str(), &sb) < 0)
            throw std::runtime_error("Unable to get filesize of index file.");
        auto opts_filesize = sb.st_size;
        if (!idx_opt_fs.read((char *)&idx_opts, opts_filesize))
            throw std::runtime_error("Unable to read " + opts.options_filename);
        if (opts.use_translated_search && idx_opts.dna_db) {
            std::cerr << "Warning: --translated-search requested but the database "
                      << "is nucleotide; using nucleotide search." << std::endl;
        }
        opts.use_translated_search = !idx_opts.dna_db;

        std::cerr << "Loading taxonomy..." << std::endl;
        taxonomy.reset(new Taxonomy(opts.taxonomy_filename, opts.use_memory_mapping));

        std::cerr << "Loading hash table..." << std::endl;
        switch (GetKVStoreCellType(opts.index_filename)) {
        case CompactHash32:
            hash.reset(new CompactHashTable<CompactHashCell>(
                opts.index_filename, opts.use_memory_mapping));
            break;
        case CompactHash40:
            hash.reset(new CompactHashTable<CompactHashCell40>(
                opts.index_filename, opts.use_memory_mapping));
            break;
        default:
            throw std::runtime_error("Unable to determine width of compact hash cell.");
        }
    }
    catch (const std::exception &ex) {
        std::cerr << "Unable to load index"
                  << ": " << ex.what() << std::endl;
        index_broken = true;
        return;
    }
    std::cerr << "Successfully loaded index." << std::endl;
    index_available = true;
}

// Drain classified batches onto the gRPC stream until the queue is closed
// and empty. Runs on its own thread so writes overlap with classification.
static void ResultsHandler(
        ServerStream *stream,
        taxon_counters_t &stream_taxon_counters,
        ClassificationStats &stream_stats,
        ThreadSafeQueue<BatchResults> &results_queue) {
    while (std::optional<BatchResults> res = results_queue.pop_wait()) {
        // The client is configured to receive messages up to INT_MAX, and
        // the client limits request batches to 128 MB, so a batch of
        // results always fits.
        Kraken2SequenceStreamResult result;
        result.mutable_classifications()->Swap(&res->k2results);
        stream->Write(result, WriteOptions().set_buffer_hint());
        // update stats for the stream
        stream_stats.total_bases += res->stats.total_bases;
        stream_stats.total_classified += res->stats.total_classified;
        stream_stats.total_sequences += res->stats.total_sequences;
        // update taxon_counters for the stream
        for (auto &kv_pair : res->taxon_counters) {
            stream_taxon_counters[kv_pair.first] += std::move(kv_pair.second);
        }
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
    const size_t max_in_flight = std::max<size_t>(4, 2 * pool.get_thread_count());
    std::mutex in_flight_mtx;
    std::condition_variable in_flight_cv;
    size_t in_flight = 0;

    // Classify while reads are still being received on the input stream
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
        pool.push_task([this, req, &results_queue, &in_flight_mtx, &in_flight_cv, &in_flight]() {
            try {
                ProcessBatch(*req, results_queue);
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
    const Kraken2SequenceRequestMulti &reqs,
    ThreadSafeQueue<BatchResults> &result_q) {

    MinimizerScanner scanner(
        idx_opts.k, idx_opts.l, idx_opts.spaced_seed_mask,
        idx_opts.dna_db, idx_opts.toggle_mask,
        idx_opts.revcom_version);
    vector<taxid_t> taxa;
    taxon_counts_t hit_counts;
    vector<string> translated_frames(6);

    BatchResults results = BatchResults();
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
            seq, paired ? &seq2 : nullptr,
            *hash, *taxonomy, idx_opts, opts, results.stats, scanner,
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

void Kraken2ServerClassifier::AddHitlistString(
    ostringstream &oss, vector<taxid_t> &taxa, Taxonomy &taxonomy)
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

Kraken2SequenceResult Kraken2ServerClassifier::ClassifySequence(
    Sequence &dna, Sequence *dna2,
    KeyValueStore &hash, Taxonomy &taxonomy, IndexOptions &idx_opts,
    Options &opts, ClassificationStats &stats, MinimizerScanner &scanner,
    vector<taxid_t> &taxa, taxon_counts_t &hit_counts,
    vector<string> &tx_frames, taxon_counters_t &curr_taxon_counts)
{
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
                    bool skip_lookup = idx_opts.minimum_acceptable_hash_value &&
                        MurmurHash3(*minimizer_ptr) < idx_opts.minimum_acceptable_hash_value;
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

    // Pass 2: one batched hash lookup for all distinct minimizers.
    static thread_local std::vector<hvalue_t> lookup_vals;
    lookup_vals.resize(lookup_keys.size());
    if (!lookup_keys.empty())
        hash.GetBatch(lookup_keys.data(), lookup_vals.data(), lookup_keys.size());

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
                taxon = lookup_vals[tok.key_idx];
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
    call = ResolveTree(hit_counts, taxonomy, total_kmers, opts);
    // Void a call made by too few minimizer groups
    if (call && minimizer_hit_groups < opts.minimum_hit_groups)
        call = 0;

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


taxid_t Kraken2ServerClassifier::ResolveTree(taxon_counts_t &hit_counts,
                                             Taxonomy &taxonomy, size_t total_minimizers, Options &opts)
{
    taxid_t max_taxon = 0;
    uint32_t max_score = 0;
    uint32_t required_score = ceil(opts.confidence_threshold * total_minimizers);

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

std::string Kraken2ServerClassifier::ReportStats(struct timeval time1, struct timeval time2,
                                                 ClassificationStats &stats)
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

std::string Kraken2ServerClassifier::ReportTotalStats(ClassificationStats &stats)
{
    uint64_t total_unclassified = stats.total_sequences - stats.total_classified;
    double denom_seqs = stats.total_sequences > 0 ? stats.total_sequences : 1.0;

    return std::to_string(stats.total_sequences) + " sequences (" + DoubleStatToString(stats.total_bases / 1.0e6, 2) + " Mbp) processed.\n" +
           std::to_string(stats.total_classified) + " sequences classified (" + DoubleStatToString(stats.total_classified * 100.0 / denom_seqs, 2) + "%).\n" +
           std::to_string(total_unclassified) + " sequences unclassified (" + DoubleStatToString(total_unclassified * 100.0 / denom_seqs, 2) + "%).\n";
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

std::string Kraken2ServerClassifier::TrimPairInfo(std::string &id)
{
    size_t sz = id.size();
    if (sz <= 2)
        return id;
    if (id[sz - 2] == '/' && (id[sz - 1] == '1' || id[sz - 1] == '2'))
        return id.substr(0, sz - 2);
    return id;
}

std::string Kraken2ServerClassifier::DoubleStatToString(double d, const int precision)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << d;
    return stream.str();
}