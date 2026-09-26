#pragma once

#include <atomic>
#include <iomanip>
#include <future>
#include <memory>

// kraken2
#include "kraken2_data.h"
#include "taxonomy.h"
#include "kv_store.h"
#include "compact_hash.h"
#include "mmscanner.h"
#include "seqreader.h"
#include "aa_translate.h"
#include "utilities.h"

// kraken2 server
#include "classify_core.h"
#include "worker_pool.h"
#include "report_server.h"
#include "thread_safe_queue.h"
#include "Kraken2.grpc.pb.h"

using namespace kraken2;

using grpc::ServerContext;
using grpc::ServerReaderWriter;
using grpc::WriteOptions;

using kraken2proto::Kraken2Service;
using kraken2proto::Kraken2SequenceRequest;
using kraken2proto::Kraken2SequenceRequestMulti;
using kraken2proto::Kraken2SequenceResult;
using kraken2proto::Kraken2SequenceResultMulti;
using kraken2proto::Kraken2SequenceStreamResult;

typedef ServerReaderWriter<Kraken2SequenceStreamResult, Kraken2SequenceRequestMulti> ServerStream;

using kraken2server::AMBIGUOUS_SPAN_TAXON;
using kraken2server::MATE_PAIR_BORDER_TAXON;
using kraken2server::READING_FRAME_BORDER_TAXON;
using kraken2server::ClassificationStats;


struct Options {
    string db_path;
    string host = "localhost";
    int port = 8080;
    int max_queue = 0;
    int thread_pool = 0;  // 0 selects the number of hardware threads

    string index_filename;
    string taxonomy_filename;
    string options_filename;
    string report_filename = "latest_run.txt";
    bool report_kmer_data = false;
    bool report_zero_counts = false;
    bool use_translated_search = false;
    bool stats = true;
    double confidence_threshold = 0.0;
    int minimum_quality_score = 0;
    int minimum_hit_groups = 2;
    bool use_memory_mapping = false;
    int wait = 0;
};




struct BatchResults {
   // Position of the request batch in the stream, used to write results
   // back in input order.
   uint64_t sequence = 0;
   Kraken2SequenceResultMulti k2results;
   taxon_counters_t taxon_counters;
   ClassificationStats stats = {0, 0, 0};
};


class Kraken2ServerClassifier {

public:
    // Written by the loader thread, read by gRPC handler threads.
    std::atomic<bool> index_available{false};
    std::atomic<bool> index_broken{false};

    /**
     * @brief Construct a new Kraken 2 Server Classifier. Starts the thread
     *        pool and begins loading the database asynchronously. The
     *        database is loaded only once and reused for all requests.
     */
    Kraken2ServerClassifier(Options &options);
    ~Kraken2ServerClassifier();

    /**
     * @brief Load kraken2 index options, taxonomy and hash table.
     */
    void LoadIndex();

    /**
     * @brief Classify sequences in a input queue and populate the classification queue.
     */
    void ProcessSequenceStream(
        ServerContext *context, ServerStream *stream, std::string &results);

    /**
     * @brief Classifies the sequences (or pairs) in a request batch and
     *        pushes the results and per-batch counters onto result_q.
     */
    void ProcessBatch(
        const Kraken2SequenceRequestMulti &reqs, uint64_t sequence,
        ThreadSafeQueue<BatchResults> &result_q);

    /**
     * @brief Return a copy of the summary of historical classifications.
     */
    std::string GetSummary();

private:
    // Database and Historical Stats
    Options opts;
    std::unique_ptr<Taxonomy> taxonomy;
    std::unique_ptr<KeyValueStore> hash;
    IndexOptions idx_opts;
    taxon_counters_t total_taxon_counters;
    ClassificationStats total_stats = {0, 0, 0};
    std::string summary;
    std::mutex stats_mtx;
    WorkerPool pool;
    std::thread loader;

    /**
     * @brief Classify one fragment. dna2 is the mate for paired reads and
     *        nullptr for single-end reads.
     */
    Kraken2SequenceResult ClassifySequence(
        Sequence &dna, Sequence *dna2,
        KeyValueStore &hash, Taxonomy &taxonomy, IndexOptions &idx_opts,
        Options &opts, ClassificationStats &stats, MinimizerScanner &scanner,
        vector<taxid_t> &taxa, taxon_counts_t &hit_counts,
        vector<string> &tx_frames, taxon_counters_t &curr_taxon_counts);

    /**
     * @brief Replace bases below the quality threshold with 'x'. Returns
     *        false if the sequence and quality strings differ in length.
     */
    bool MaskLowQualityBases(Sequence &dna, int minimum_quality_score);

    /**
     * @brief Build an unclassified result for a read that could not be
     *        processed.
     */
    Kraken2SequenceResult UnclassifiedResult(const Sequence &dna, const Sequence *dna2);

    void GenerateReport(
        std::string &results, std::string &summary, Options &opts, Taxonomy &taxonomy,
        timeval &tv1, timeval &tv2, ClassificationStats &stats, ClassificationStats &total_stats,
        taxon_counters_t &taxon_counters, taxon_counters_t &total_taxon_counters, std::mutex &stats_mtx);

};
