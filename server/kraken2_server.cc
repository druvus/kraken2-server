#include <atomic>
#include <csignal>
#include <memory>
#include <sysexits.h>

#include <grpc/grpc.h>
#include <grpc++/server.h>
#include <grpc++/server_builder.h>
#include <grpc++/server_context.h>
#include <grpc++/security/server_credentials.h>

#include "cli.h"
#include "messages.h"
#include "classify_server.h"

using grpc::ResourceQuota;
using grpc::Server;
using grpc::ServerBuilder;
using grpc::ServerContext;
using grpc::ServerReader;
using grpc::ServerReaderWriter;
using grpc::Status;
using grpc::StatusCode;

using kraken2proto::Kraken2ReadyRequest;
using kraken2proto::Kraken2ReadyResult;
using kraken2proto::Kraken2SummaryRequest;
using kraken2proto::Kraken2SummaryResults;
using kraken2proto::Kraken2ShutdownRequest;
using kraken2proto::Kraken2ShutdownResult;
using kraken2proto::Kraken2SequenceRequest;
using kraken2proto::Kraken2SequenceStreamResult;
using kraken2proto::Kraken2Service;



// Shared between the signal handler, the RemoteShutdown RPC and RunServer.
// The promise may only be fulfilled once, so a flag guards against a second
// signal or a shutdown request arriving after Ctrl-C.
static std::promise<void> exit_promise;
static std::atomic<bool> exit_flag{false};
static std::promise<void> *exit_requested = &exit_promise;

static void RequestExit() {
    if (!exit_flag.exchange(true)) {
        exit_requested->set_value();
    }
}


class ServiceImpl final : public Kraken2Service::Service {

public:
    ServiceImpl(Options opts, Kraken2ServerClassifier *classifier)
    : options(opts), classifier(classifier)
    {}

    /**
     * @brief Endpoint to request a summary of the classification history on the server.
     */
    Status GetSummary(
           ServerContext *context, const Kraken2SummaryRequest *req,
            Kraken2SummaryResults *results) override {
        if (!classifier->index_available) {
            return IndexStatus();
        }

        // Only return summary if the server is recording history.
        if (options.stats) {
            results->set_summary(classifier->GetSummary());
        }
        // Else indicate to the user it is not available.
        else {
            results->set_summary("Summary not available on this server.");
        }
        return Status::OK;
    }

    /** 
     * @brief Endpoint to request ask server if it is ready.
     */
    Status ServerReady(
            ServerContext *context, const Kraken2ReadyRequest *req,
            Kraken2ReadyResult *results) override {
        results->set_ready(classifier->index_available && !classifier->index_broken);
        return IndexStatus();
    }

    /**
     * @brief Endpoint to initiate a remote shutdown of the server.
     */
    Status RemoteShutdown(
            ServerContext *context, const Kraken2ShutdownRequest *req,
            kraken2proto::Kraken2ShutdownResult *result) override {
        std::cerr << "Received shutdown request." << std::endl;
        RequestExit();
        result->set_successful(true);
        std::cerr << "Shutdown request made." << std::endl;
        return Status::OK;
    }

    /**
     * @brief Endpoint to classify a stream of sequences and return
     *        a stream of classifications as response.
     */
    Status ClassifyStream(
            ServerContext *context, ServerStream *reader_writer) override {
        if (!classifier->index_available) {
            return IndexStatus();
        }

        std::string results;

        classifier->ProcessSequenceStream(context, reader_writer, std::ref(results));

        // If connection is open, send a final message containing the summary.
        if (!context->IsCancelled()) {
            Kraken2SequenceStreamResult summary;
            summary.set_summary(results.c_str());
            reader_writer->Write(summary);
        }

        return Status::OK;
    }

private:
    Options options;
    Kraken2ServerClassifier *classifier;

    grpc::Status IndexNotLoaded = grpc::Status(grpc::StatusCode::UNAVAILABLE, "Index not loaded yet, please wait.");
    grpc::Status IndexError = grpc::Status(grpc::StatusCode::FAILED_PRECONDITION, "There was an error loading the index, the server will remain unavailable without intervention.");
    grpc::Status IndexLoaded = grpc::Status(grpc::StatusCode::OK, "Index loaded.");

    grpc::Status IndexStatus(){
        if(!classifier->index_available) {
            return classifier->index_broken ? IndexError : IndexNotLoaded;
        }
        return IndexLoaded;
    } 

};

void RunServer(Options opts, Kraken2ServerClassifier *classifier) {
    std::string server_address = opts.host + ":" + std::to_string(opts.port);
    ServiceImpl service(opts, classifier);
    // Sets the max number of concurrent requests
    ResourceQuota rq;
    if (opts.max_queue > 0){
        // +1 to account for the server thread itself. If one thread were to
        // be specified, all requests are rejected.
        rq.SetMaxThreads(opts.max_queue + 1);
    }
    // According to gRPC devs, does not fully track gRPC memory usage and no
    // docs on unit of memory allocation - not recommended
    // rq.Resize(new_memory_allocation);
    ServerBuilder builder;
    builder.AddListeningPort(server_address, grpc::InsecureServerCredentials());
    // don't use port if already in use
    builder.AddChannelArgument(GRPC_ARG_ALLOW_REUSEPORT, 0);
    builder.SetResourceQuota(rq);
    builder.RegisterService(&service);
    // allow 128Mb messages
    builder.SetMaxSendMessageSize(128 * 1024 * 1024);
    builder.SetMaxMessageSize(128 * 1024 * 1024);
    builder.SetMaxReceiveMessageSize(128 * 1024 * 1024);
    // Shared pointer so graceful shutdown can be invoked.
    std::shared_ptr<Server> server(builder.BuildAndStart());
    if (server == nullptr) {
        std::cout << "Failed to start server on " << server_address
                  << ". See above for more details." << std::endl;     
    } else {
        std::cout << "Server listening on " << server_address
                  << ". Press Ctrl-C to end." << std::endl;
        // handle interrupts
        auto handler = [](int s) { RequestExit(); };
        std::signal(SIGINT, handler);
        std::signal(SIGTERM, handler);
        std::signal(SIGQUIT, handler);

        // block until exit request is set
        auto f = exit_requested->get_future();
        f.wait();
        server->Shutdown();
    }
}


void ParseCommandLine(int argc, char **argv, Options &opts) {
    using cli::Parser;
    Parser parser("kraken2_server");
    parser.add({"db", 'd', true, "[path]", "Path to Kraken 2 database",
        [&](const std::string &v) {
            opts.db_path = v;
            opts.taxonomy_filename = v + "/taxo.k2d";
            opts.options_filename = v + "/opts.k2d";
            opts.index_filename = v + "/hash.k2d";
        }, true});
    parser.add({"max-requests", 'r', true, "[int]", "Max number of client requests processed concurrently (0 for default)",
        [&](const std::string &v) { opts.max_queue = Parser::ParseInt(v, "--max-requests", 0, 100000); }});
    parser.add({"thread-pool", 'x', true, "[int]", "Classification threads shared by all clients (0, the default, uses all hardware threads)",
        [&](const std::string &v) { opts.thread_pool = Parser::ParseInt(v, "--thread-pool", 0, 100000); }});
    parser.add({"no-stats", 's', false, "", "Do not track statistics of all processed sequences on this server. Saves memory long-term.",
        [&](const std::string &) { opts.stats = false; }});
    parser.add({"host-ip", 'i', true, "[addr]", "Server IP address (default: localhost)",
        [&](const std::string &v) { opts.host = v; }});
    parser.add({"port", 'p', true, "[int]", "Port number on which to listen for requests (0 - 65535, default 8080)",
        [&](const std::string &v) { opts.port = Parser::ParseInt(v, "--port", 0, 65535); }});
    parser.add({"report-kmer", 'k', false, "", "Include distinct k-mers in reports",
        [&](const std::string &) { opts.report_kmer_data = true; }});
    parser.add({"report-zero", 'z', false, "", "Include zero count taxons in reports",
        [&](const std::string &) { opts.report_zero_counts = true; }});
    parser.add({"translated-search", 't', false, "", "Use translated search (set automatically for protein databases)",
        [&](const std::string &) { opts.use_translated_search = true; }});
    parser.add({"confidence", 'c', true, "[double]", "Confidence score threshold (default: 0.0) (0 - 1)",
        [&](const std::string &v) { opts.confidence_threshold = Parser::ParseDouble(v, "--confidence", 0.0, 1.0); }});
    parser.add({"min-quality", 'q', true, "[int]", "Minimum base quality used in classification (default: 0), FASTQ input only",
        [&](const std::string &v) { opts.minimum_quality_score = Parser::ParseInt(v, "--min-quality", 0, 1000); }});
    parser.add({"hit-groups", 'g', true, "[int]", "Minimum number of hit groups (overlapping k-mers sharing the same minimizer) needed to make a call (default: 2)",
        [&](const std::string &v) { opts.minimum_hit_groups = Parser::ParseInt(v, "--hit-groups", 0, 1000000); }});
    parser.add({"memory-mapping", 'o', false, "", "Avoids loading database into RAM",
        [&](const std::string &) { opts.use_memory_mapping = true; }});
    parser.add({"wait", 'w', true, "[int]", "Delay database loading by this many seconds (for testing)",
        [&](const std::string &v) { opts.wait = Parser::ParseInt(v, "--wait", 0, 100000); }});

    cli::Result result = parser.parse(argc, argv);
    if (result.status == cli::Status::Help) {
        std::cerr << parser.usage();
        exit(0);
    }
    if (result.status == cli::Status::Error) {
        std::cerr << result.message << std::endl << std::endl << parser.usage();
        exit(EX_USAGE);
    }
}


int main(int argc, char **argv) {
    Options opts;
    ParseCommandLine(argc, argv, opts);
    std::unique_ptr<Kraken2ServerClassifier> classifier(new Kraken2ServerClassifier(opts));
    RunServer(opts, classifier.get());

    return classifier->index_available ? EX_OK : EX_IOERR;
}
