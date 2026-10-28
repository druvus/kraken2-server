#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <mutex>
#include <fstream>
#include <future>
#include <random>
#include <thread>
#include <sysexits.h>

#include <grpc/grpc.h>
#include <grpc++/channel.h>
#include <grpc++/client_context.h>
#include <grpc++/create_channel.h>
#include <grpc++/security/credentials.h>

#include "cli.h"
#include "utils.h"
#include "thread_safe_queue.h"
#include "Kraken2.grpc.pb.h"

#include "kseq.cc.h"

using namespace std::chrono_literals; // ns, us, ms, s, h, etc.

using grpc::Channel;
using grpc::ClientContext;
using grpc::ClientReaderWriter;
using grpc::ClientWriter;
using grpc::Status;
using grpc::WriteOptions;

using kraken2proto::Kraken2ReadyRequest;
using kraken2proto::Kraken2ReadyResult;
using kraken2proto::Kraken2SequenceRequest;
using kraken2proto::Kraken2SequenceRequestMulti;
using kraken2proto::Kraken2SequenceResult;
using kraken2proto::Kraken2SequenceResultMulti;
using kraken2proto::Kraken2SequenceStreamResult;
using kraken2proto::Kraken2Service;
using kraken2proto::Kraken2SummaryRequest;
using kraken2proto::Kraken2SummaryResults;
using kraken2proto::Kraken2ShutdownRequest;
using kraken2proto::Kraken2ShutdownResult;

// Command line options
struct Options
{
    std::string sequence;
    std::string sequence2;  // mate file for paired-end reads, empty for single-end
    std::string report_file;
    std::string host = "localhost";
    int port = 8080;
    bool shutdown = false;
    // TLS: enabled by --tls or any --tls-* option. tls_ca overrides the system
    // roots; tls_cert and tls_key present a client certificate.
    bool tls = false;
    std::string tls_ca;
    std::string tls_cert;
    std::string tls_key;
    std::string tls_server_name;  // expected certificate name when it differs from the host
};

typedef std::shared_ptr<ClientReaderWriter<Kraken2SequenceRequestMulti, Kraken2SequenceStreamResult>> ClientStream;


// Buffering is bounded in bytes as well as reads so that long reads do not
// make the client hold its whole input in memory: a file batch holds at
// most FASTQ_BATCH_SIZE reads or FASTQ_BATCH_BYTES of data, MAX_BATCHES of
// them are buffered ahead of the writer, and at most MAX_IN_FLIGHT_BYTES
// (or MAX_IN_FLIGHT reads) are sent and not yet answered. Worst case
// the client holds about MAX_BATCHES * FASTQ_BATCH_BYTES + MAX_IN_FLIGHT_BYTES.
#define ST_BATCH_SIZE 2000                            // reads in a gRPC message
#define MAX_IN_FLIGHT 64000                           // reads sent but not yet answered
#define MAX_IN_FLIGHT_BYTES (256ULL * 1024 * 1024)    // bytes sent but not yet answered
#define FASTQ_BATCH_SIZE 4000                         // reads in a file batch
#define FASTQ_BATCH_BYTES (8ULL * 1024 * 1024)        // bytes in a file batch
#define MAX_BATCHES 16                                // file batches buffered ahead of the writer

// Reads and bytes sent and not yet answered. The writer adds a message
// after sending it and waits while the limits are exceeded; the reader
// retires reads as results arrive, in order, and releases a message's
// bytes once all of its reads are answered.
class InFlight {
public:
    void add(uint64_t reads, uint64_t bytes) {
        std::lock_guard<std::mutex> lock(mtx_);
        messages_.push_back({reads, bytes});
        reads_ += reads;
        bytes_ += bytes;
    }
    void retire_one() {
        std::lock_guard<std::mutex> lock(mtx_);
        reads_--;
        if (!messages_.empty() && ++answered_in_head_ == messages_.front().reads) {
            bytes_ -= messages_.front().bytes;
            messages_.pop_front();
            answered_in_head_ = 0;
        }
        changed_.notify_all();
    }
    // Block until a message of the given size may be sent, or the stream
    // has closed. Returns true when it may be sent.
    bool wait_for_room(uint64_t reads, uint64_t bytes, std::atomic<bool> &closed) {
        std::unique_lock<std::mutex> lock(mtx_);
        bool waited = false;
        while (!closed && reads_ > 0 &&
               (reads_ + reads > MAX_IN_FLIGHT || bytes_ + bytes > MAX_IN_FLIGHT_BYTES)) {
            if (!waited) {
                waited = true;
                std::cerr << "Waiting before sending more. In-flight: " << reads_
                          << " reads, " << bytes_ / (1024 * 1024) << " MB." << std::endl;
            }
            changed_.wait_for(lock, 100ms);
        }
        return !closed;
    }
    void wake() { changed_.notify_all(); }
    uint64_t reads() {
        std::lock_guard<std::mutex> lock(mtx_);
        return reads_;
    }
private:
    struct Message { uint64_t reads, bytes; };
    std::mutex mtx_;
    std::condition_variable changed_;
    std::deque<Message> messages_;
    uint64_t reads_ = 0, bytes_ = 0;
    uint64_t answered_in_head_ = 0;
};

class SequenceClient {

public:
    SequenceClient(std::shared_ptr<Channel> channel) : sequence_stub(Kraken2Service::NewStub(channel)) {}


    /**
     * @brief Send sequences from a kseq file as a stream and receive classifications individually as a stream.
     *
     * @param sequence_name  path to reads (R1 for paired data)
     * @param sequence2_name path to mates (R2), empty for single-end data
     * @return EX_IOERR if sequences could not be read
     * @return EX_DATAERR if the two paired-end inputs have different lengths
     * @return EX_UNAVAILABLE if sequences could not be sent to server
     * @return EX_PROTOCOL if the server returned fewer results than reads sent
     * @return else gRPC status code
     */
    int ClassifySequences(const std::string &sequence_name,
                          const std::string &sequence2_name,
                          const std::string &report_file) {
        std::cerr << "Classifying sequence stream." << std::endl;
        int state = WaitForServer();
        if (state != 0) {return state;}

        ClientContext context;
        ClientStream stream(sequence_stub->ClassifyStream(&context));
        InFlight in_flight;

        // Bounded queue of read batches between the file reader and the
        // gRPC writer. The reader blocks when MAX_BATCHES are buffered and
        // closes the queue when the input is exhausted.
        ThreadSafeQueue<std::vector<Kraken2SequenceRequest>> batches_queue(MAX_BATCHES);

        // reads data from file into queue
        std::future<int> fastq_batches = std::async(
            std::launch::async, &SequenceClient::FastBatcher, this,
            std::ref(sequence_name), std::ref(sequence2_name), std::ref(batches_queue));

        // take data from queue and send over gRPC
        std::future<int> stream_batches = std::async(
            std::launch::async, &SequenceClient::StreamWriter, this,
            std::ref(in_flight), std::ref(batches_queue), std::ref(stream));

        // reading back results on gRPC stream
        std::future<int> recv_reads = std::async(
            std::launch::async, &SequenceClient::StreamReader, this,
            std::ref(in_flight), std::ref(report_file), std::ref(stream));

        // wait for things to finish in order
        fastq_batches.wait();
        stream_batches.wait();
        recv_reads.wait();
        std::cerr << "Done waiting" << std::endl;

        const int sent = stream_batches.get();
        const int received = recv_reads.get();
        std::cerr << "Sent    : " << sent << std:: endl;
        std::cerr << "Received: " << received << std::endl;

        // Handle the stream response
        Status status = stream->Finish();
        if (!status.ok()) {
            std::cerr << "Client RPC stream failed: " << status.error_message() << std::endl;
            return status.error_code();
        }
        if (received != sent) {
            // The server completed the stream but dropped results, which
            // it reports on its own log.
            std::cerr << "Error: the server returned " << received << " results for "
                      << sent << " reads; the output is incomplete." << std::endl;
            return EX_PROTOCOL;
        }
        if (fastq_batches.get() < 0) {
            return EX_IOERR;
        }
        if (pairs_mismatched) {
            return EX_DATAERR;
        }
        return status.error_code();
    }

    /**
     * @brief Request a summary of the classification history on the server.
     *
     * @return gRPC status code of request
     */
    int GetSummary() {
        ClientContext context;
        Kraken2SummaryRequest req;
        Kraken2SummaryResults response;

        Status status = sequence_stub->GetSummary(&context, req, &response);
        if (!status.ok())
        {
            std::cerr << "Could not retrieve Kraken2 server summary." << std::endl;
        }
        std::cout << response.summary() << std::endl;
        return status.error_code();
    }

    /**
     * @brief Shutdown the server remotely
     *
     * @return gRPC status code of request
     */
    int ShutdownServer() {
        ClientContext context;
        Kraken2ShutdownRequest req;
        Kraken2ShutdownResult response;
        Status status = sequence_stub->RemoteShutdown(&context, req, &response);
        if (!status.ok()) {
            std::cerr << "Failed to send shutdown request: " << status.error_message() << std::endl;
            return status.error_code();
        }
        if (response.successful()) {
            std::cerr << "Shutdown request processed." << std::endl;
        }
        else{
            std::cerr << "Shutdown request not processed correctly." << std::endl;
        }
        return status.error_code();
    }

    // Send one message, waiting first until it fits within the in-flight
    // limits. Returns the number of reads sent, 0 once the stream has failed.
    uint64_t SendMessage(Kraken2SequenceRequestMulti &req, uint64_t msg_size,
                         InFlight &in_flight, ClientStream &writer, bool buffer_hint) {
        const uint64_t reads = req.seqs_size();
        if (!in_flight.wait_for_room(reads, msg_size, stream_closed)) {
            return 0;
        }
        bool ok = buffer_hint ? writer->Write(req, WriteOptions().set_buffer_hint())
                              : writer->Write(req);
        if (!ok) {
            StreamFailed();
            return 0;
        }
        in_flight.add(reads, msg_size);
        return reads;
    }

    int StreamWriter(
            InFlight &in_flight,
            ThreadSafeQueue<std::vector<Kraken2SequenceRequest>> &batches,
            ClientStream &writer) {
        int seqs_sent = 0;
        try {
            // pop_wait returns an empty optional once the reader has closed
            // the queue and every batch has been taken.
            while (std::optional<std::vector<Kraken2SequenceRequest>> item = batches.pop_wait()) {
                std::vector<Kraken2SequenceRequest> batch = std::move(*item);

                // rebatch to smaller batches for stream, moving each read
                // into the message rather than copying it. File batches are
                // already bounded in bytes; the message limit is a safety net
                // for a single read of exceptional size.
                const uint64_t MAX_SIZE = 128 * 1024 * 1024;
                for(size_t i = 0; i < batch.size() && !stream_closed; i += ST_BATCH_SIZE) {
                    auto last = std::min(batch.size(), i + ST_BATCH_SIZE);
                    size_t bsize = last - i;
                    Kraken2SequenceRequestMulti req;
                    req.mutable_seqs()->Reserve(bsize);
                    for (size_t k = i; k < last; ++k) {
                        req.mutable_seqs()->Add(std::move(batch[k]));
                    }
                    uint64_t msg_size = req.ByteSizeLong();
                    if (msg_size > MAX_SIZE) {
                        // send one by one
                        for (auto &read : *req.mutable_seqs()) {
                            Kraken2SequenceRequestMulti single;
                            single.mutable_seqs()->Add(std::move(read));
                            uint64_t single_size = single.ByteSizeLong();
                            if (single_size > MAX_SIZE) {
                                std::cerr << single.seqs(0).id() << ": read is larger than the "
                                          << MAX_SIZE / (1024 * 1024) << " MB message limit; skipping." << std::endl;
                                continue;
                            }
                            seqs_sent += SendMessage(single, single_size, in_flight, writer, true);
                            if (stream_closed) {
                                break;
                            }
                        }
                    }
                    else {
                        seqs_sent += SendMessage(req, msg_size, in_flight, writer, false);
                    }
                }
                if (stream_closed) {
                    break;
                }
            }
        }
        catch (const std::exception &ex) {
            std::cerr << "Failed to send sequences"
                      << ": " << ex.what() << std::endl;
            StreamFailed();
        }

        if (stream_closed) {
            // Release the file reader, which may be blocked on a full queue.
            batches.close();
        }
        writer->WritesDone();
        return seqs_sent;
    }

    int StreamReader(
            InFlight &in_flight, const std::string &report_file,
            ClientStream &reader) {
        Kraken2SequenceStreamResult result;
        int n_reads = 0;
        try {
            while (reader->Read(&result)) {
                if (result.has_classifications()) {
                    for (auto &res : result.classifications().classes()){
                        n_reads++;
                        PrintClassification(res);
                        in_flight.retire_one();
                    }
                }
                else if (result.has_summary()) {
                    PrintSummary(result.summary(), report_file);
                }
                else {
                    std::cerr << "Result had neither classifications or summary :/" << std::endl;
                }
            }
        }
        catch (const std::exception &ex) {
            std::cerr << "Failed to receive responses"
                      << ": " << ex.what() << std::endl;
        }
        // No more results will arrive, whether the stream ended normally
        // or the server went away.
        stream_closed = true;
        in_flight.wake();
        return n_reads;
    }
    
    /**
     * @brief Read batches of reads (or read pairs when sequence2_file is
     *        non-empty) from disk onto the batches queue.
     *
     * @return number of batches read, or -1 if the input could not be read.
     */
    int FastBatcher(
            const std::string &sequence_file,
            const std::string &sequence2_file,
            ThreadSafeQueue<std::vector<Kraken2SequenceRequest>> &batches_queue) {
        int n_batches = 0;
        const bool paired = !sequence2_file.empty();
        // Whatever happens, the writer must be released.
        struct Closer {
            ThreadSafeQueue<std::vector<Kraken2SequenceRequest>> &q;
            ~Closer() { q.close(); }
        } closer{batches_queue};
        try {
            FastReader reader(sequence_file);
            std::unique_ptr<FastReader> mate_reader;
            if (paired) {
                mate_reader.reset(new FastReader(sequence2_file));
                std::cerr << "Reading read pairs from files: " << sequence_file
                          << " and " << sequence2_file << std::endl;
            }
            else {
                std::cerr << "Reading sequences from file: " << sequence_file << std::endl;
            }
            bool mismatched = false;
            while (true) {
                std::vector<Kraken2SequenceRequest> seqs;
                int n_reads;
                if (paired) {
                    n_reads = reader.read_pairs(seqs, FASTQ_BATCH_SIZE, *mate_reader, mismatched,
                                                FASTQ_BATCH_BYTES);
                }
                else {
                    n_reads = reader.read(seqs, FASTQ_BATCH_SIZE, FASTQ_BATCH_BYTES);
                }
                if (n_reads > 0) {
                    n_batches++;
                    // blocks while MAX_BATCHES are already buffered; false
                    // once the writer has closed the queue after a failure
                    if (!batches_queue.push(std::move(seqs))) {
                        break;
                    }
                }
                // A batch ends early at the byte budget; the input is done
                // when a reader has hit the end of its file or an error.
                if (reader.last_status() < 0 || (paired && mate_reader->last_status() < 0)) {
                    break;
                }
            }
            if (mismatched) {
                std::cerr << "Warning: paired-end inputs have different numbers of reads; "
                          << "only complete pairs were classified: "
                          << sequence_file << ", " << sequence2_file << std::endl;
                pairs_mismatched = true;
            }
            if (reader.failed() || (paired && mate_reader->failed())) {
                // The reader has already described the problem.
                return -1;
            }
        }
        catch (const std::exception &ex) {
            std::cerr << "Failed to read sequences from file: " << sequence_file
                      << ": " << ex.what() << std::endl;
            return -1;
        }
        return n_batches;
    }


private:

    // The gRPC service stub for the service defined in Kraken2.proto
    std::unique_ptr<kraken2proto::Kraken2Service::Stub> sequence_stub;

    // Set by FastBatcher when R1 and R2 have different read counts.
    std::atomic<bool> pairs_mismatched{false};

    // Set when the stream can no longer carry results: a write failed or
    // the reader saw the end of the stream. The writer stops sending and
    // stops waiting for in-flight results, which would never arrive.
    std::atomic<bool> stream_closed{false};

    void StreamFailed() {
        if (!stream_closed.exchange(true)) {
            std::cerr << "The stream to the server failed; not sending more reads." << std::endl;
        }
    }

    int WaitForServer() {
        // The server answers UNAVAILABLE with a message starting "Index not
        // loaded" while its database loads; wait for that indefinitely.
        // Any other UNAVAILABLE comes from the transport (connection refused,
        // TLS mismatch) and is retried only a few times.
        const int max_connect_attempts = 5;
        int connect_attempts = 0;
        while (true) {
            ClientContext context;
            Kraken2ReadyRequest req;
            Kraken2ReadyResult response;
            Status status = sequence_stub->ServerReady(&context, req, &response);
            if (status.ok()) {
                std::cerr << "Server responded as ready." << std::endl;
                return EX_OK;
            }
            if (status.error_code() != grpc::StatusCode::UNAVAILABLE) {
                std::cerr << "Server is in error state: "
                          << status.error_message() << std::endl;
                return status.error_code();
            }
            if (status.error_message().rfind("Index not loaded", 0) == 0) {
                std::cerr << "Server is not ready: " << status.error_message() << std::endl;
                std::cerr << "Waiting 10s..." << std::endl;
                std::this_thread::sleep_for(10s);
                continue;
            }
            connect_attempts++;
            std::cerr << "Cannot reach server (" << connect_attempts << "/" << max_connect_attempts
                      << "): " << status.error_message() << std::endl;
            if (connect_attempts >= max_connect_attempts) {
                std::cerr << "Giving up. Check the address and port, and whether the "
                          << "server requires TLS (--tls, --tls-ca)." << std::endl;
                return EX_UNAVAILABLE;
            }
            std::this_thread::sleep_for(2s);
        }
    }

    void PrintSummary(const std::string &summary, const std::string &report_file) {
        if (report_file != "") {
            try {
                std::ofstream summary_file(report_file, std::ofstream::out);
                summary_file << summary;
                summary_file.close();
            }
            catch (const std::exception &ex) {
                std::cerr << "Failed to write report file"
                          << ": " << ex.what() << std::endl;
            }
        }
    }

    /**
     * @brief Print the classification.
     *
     * @param classification
     */
    void PrintClassification(const Kraken2SequenceResult &classification) {
        std::string classified = classification.classified() ? "C" : "U";
        std::cout
            << classified << '\t'
            << classification.id() << '\t'
            << classification.tax_id() << '\t'
            << classification.size();
        // kraken2 reports paired fragments as len1|len2
        if (classification.paired()) {
            std::cout << '|' << classification.size2();
        }
        std::cout
            << '\t'
            << classification.hitlist() << '\n';
    }
};

void ParseCommandLine(int argc, char **argv, Options &opts) {
    using cli::Parser;
    Parser parser("kraken2_client");
    parser.add({"sequence", 's', true, "[path]", "Path to sequence file (*.fast(a|q)(.gz)), or - for stdin. Omit to request the server summary.",
        [&](const std::string &v) { opts.sequence = v; }});
    parser.add({"sequence2", '2', true, "[path]", "Path to mate file for paired-end reads (same order as --sequence)",
        [&](const std::string &v) { opts.sequence2 = v; }});
    parser.add({"report", 'r', true, "[path]", "Path to output report file",
        [&](const std::string &v) { opts.report_file = v; }});
    parser.add({"host-ip", 'i', true, "[addr]", "Server IP address (default: localhost)",
        [&](const std::string &v) { opts.host = v; }});
    parser.add({"port", 'p', true, "[int]", "Server port (default: 8080)",
        [&](const std::string &v) { opts.port = Parser::ParseInt(v, "--port", 0, 65535); }});
    parser.add({"shutdown", 'k', false, "", "Shutdown server",
        [&](const std::string &) { opts.shutdown = true; }});
    parser.add({"tls", 't', false, "", "Connect with TLS using the system CA roots",
        [&](const std::string &) { opts.tls = true; }});
    parser.add({"tls-ca", 0, true, "[path]", "PEM CA bundle to verify the server (implies --tls)",
        [&](const std::string &v) { opts.tls_ca = v; opts.tls = true; }});
    parser.add({"tls-cert", 0, true, "[path]", "PEM client certificate for mutual TLS (implies --tls)",
        [&](const std::string &v) { opts.tls_cert = v; opts.tls = true; }});
    parser.add({"tls-key", 0, true, "[path]", "PEM client private key for --tls-cert",
        [&](const std::string &v) { opts.tls_key = v; opts.tls = true; }});
    parser.add({"tls-server-name", 0, true, "[name]", "Name to verify in the server certificate when it differs from --host-ip",
        [&](const std::string &v) { opts.tls_server_name = v; opts.tls = true; }});

    cli::Result result = parser.parse(argc, argv);
    if (result.status == cli::Status::Help) {
        std::cerr << parser.usage();
        exit(0);
    }
    if (result.status == cli::Status::Error) {
        std::cerr << result.message << std::endl << std::endl << parser.usage();
        exit(EX_USAGE);
    }
    if (!opts.sequence2.empty() && opts.sequence.empty()) {
        std::cerr << "--sequence2 requires --sequence." << std::endl;
        exit(EX_USAGE);
    }
    if (opts.tls_cert.empty() != opts.tls_key.empty()) {
        std::cerr << "--tls-cert and --tls-key must be given together." << std::endl;
        exit(EX_USAGE);
    }
}

int main(int argc, char **argv) {
    Options opts;
    ParseCommandLine(argc, argv, opts);

    int rtn_code = 0;
    std::string server_address = opts.host + ":" + std::to_string(opts.port);

    std::cerr << "Connecting to server: " << server_address << "." << std::endl;

    // when we send messages from the client we break up sequence
    // batches to stay below 128Mb message size as set up in kraken2_server.cc
    // We need to similarly ensure the the client can receive more than the
    // default 4MB message size. Just set it to the max
    grpc::ChannelArguments ch_args;
    ch_args.SetMaxReceiveMessageSize(INT_MAX);
    std::shared_ptr<grpc::ChannelCredentials> credentials;
    if (opts.tls) {
        grpc::SslCredentialsOptions ssl_opts;
        try {
            if (!opts.tls_ca.empty()) ssl_opts.pem_root_certs = read_file(opts.tls_ca);
            if (!opts.tls_cert.empty()) {
                ssl_opts.pem_cert_chain = read_file(opts.tls_cert);
                ssl_opts.pem_private_key = read_file(opts.tls_key);
            }
        }
        catch (const std::exception &ex) {
            std::cerr << "Cannot configure TLS: " << ex.what() << std::endl;
            return EX_NOINPUT;
        }
        if (!opts.tls_server_name.empty()) {
            ch_args.SetSslTargetNameOverride(opts.tls_server_name);
        }
        credentials = grpc::SslCredentials(ssl_opts);
    } else {
        credentials = grpc::InsecureChannelCredentials();
    }
    std::shared_ptr<grpc::Channel> ch =
        grpc::CreateCustomChannel(server_address, credentials, ch_args);
    SequenceClient client(ch);

    if (opts.shutdown) {
        rtn_code = client.ShutdownServer();
    }
    else if (opts.sequence.empty()) {
        rtn_code = client.GetSummary();
    }
    else {
        const std::string filename(opts.sequence);
        const std::string filename2(opts.sequence2);
        const std::string report_file(opts.report_file);
        rtn_code = client.ClassifySequences(filename, filename2, report_file);
    }

    std::cerr << "Return code: " << rtn_code << std::endl;
    return rtn_code;
}
