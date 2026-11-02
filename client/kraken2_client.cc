#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cerrno>
#include <cstdlib>
#include <cstring>
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
    // Ping the server after this many seconds without data on a stream and
    // give up when no answer comes within KEEPALIVE_TIMEOUT_MS; 0 disables.
    int keepalive = 60;
};

// How long a ping to the server may go unanswered before the client gives
// up on the connection.
#define KEEPALIVE_TIMEOUT_MS 20000

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
// Reads that arrive slowly (a pipe fed by a running sequencer) are sent
// once nothing else has arrived for this long, rather than when a file
// batch is full.
#define FLUSH_INTERVAL 1s

// Reads that do not yet fill a file batch. The file reader adds to it and
// hands it to the batch queue once full; the writer takes it when the
// queue has been empty for FLUSH_INTERVAL.
struct PendingReads {
    std::mutex mtx;
    std::vector<Kraken2SequenceRequest> reads;
    uint64_t bytes = 0;

    // Add a read; returns the full batch to queue, or an empty vector.
    std::vector<Kraken2SequenceRequest> add(Kraken2SequenceRequest &&rec, uint64_t rec_bytes) {
        std::lock_guard<std::mutex> lock(mtx);
        reads.push_back(std::move(rec));
        bytes += rec_bytes;
        if (reads.size() >= FASTQ_BATCH_SIZE || bytes >= FASTQ_BATCH_BYTES) {
            return take_locked();
        }
        return {};
    }
    std::vector<Kraken2SequenceRequest> take() {
        std::lock_guard<std::mutex> lock(mtx);
        return take_locked();
    }
    // The next batch for the writer after the queue stayed empty: a batch
    // queued meanwhile comes first (the file reader queues a full batch
    // before it adds the next read here, so a non-empty pending batch is
    // always later than anything in the queue), else the pending reads.
    std::vector<Kraken2SequenceRequest> take_after(
            ThreadSafeQueue<std::vector<Kraken2SequenceRequest>> &queue) {
        std::lock_guard<std::mutex> lock(mtx);
        if (std::optional<std::vector<Kraken2SequenceRequest>> queued = queue.pop()) {
            return std::move(*queued);
        }
        return take_locked();
    }
private:
    std::vector<Kraken2SequenceRequest> take_locked() {
        std::vector<Kraken2SequenceRequest> out;
        out.swap(reads);
        bytes = 0;
        return out;
    }
};

static uint64_t RecordBytes(const Kraken2SequenceRequest &rec) {
    uint64_t n = rec.id().size() + rec.seq().size() + rec.quals().size();
    if (rec.has_mate()) {
        n += rec.mate().seq().size() + rec.mate().quals().size();
    }
    return n;
}

// Reads and bytes sent and not yet answered. The writer adds a message
// just before sending it (so its first result cannot arrive before it is
// counted) and waits while the limits are exceeded; the reader retires
// reads as results arrive, in order, and releases a message's bytes once
// all of its reads are answered.
class InFlight {
public:
    void add(uint64_t reads, uint64_t bytes) {
        std::lock_guard<std::mutex> lock(mtx_);
        messages_.push_back({reads, bytes});
        reads_ += reads;
        bytes_ += bytes;
    }
    // Take back the message added last, when sending it failed.
    void undo_last() {
        std::lock_guard<std::mutex> lock(mtx_);
        if (messages_.empty()) return;
        reads_ -= messages_.back().reads;
        bytes_ -= messages_.back().bytes;
        messages_.pop_back();
        changed_.notify_all();
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
     * @return EX_NOINPUT if an input file cannot be opened
     * @return EX_IOERR if sequences could not be read, or the output could
     *         not be written
     * @return EX_CANTCREAT if the report could not be written
     * @return EX_DATAERR if the two paired-end inputs have different
     *         lengths, or a read exceeds the message limit
     * @return EX_UNAVAILABLE if the server cannot be reached
     * @return EX_PROTOCOL if the server returned fewer results than reads sent
     * @return else the gRPC status code of a failed stream
     */
    int ClassifySequences(const std::string &sequence_name,
                          const std::string &sequence2_name,
                          const std::string &report_file) {
        // Open the inputs first, so a missing file is reported as such
        // rather than after waiting for the server.
        std::unique_ptr<FastReader> reader, mate_reader;
        try {
            reader.reset(new FastReader(sequence_name));
            if (!sequence2_name.empty()) {
                mate_reader.reset(new FastReader(sequence2_name));
            }
        }
        catch (const std::exception &ex) {
            std::cerr << ex.what() << std::endl;
            return EX_NOINPUT;
        }

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
        PendingReads pending;

        // reads data from file into queue
        std::future<int> fastq_batches = std::async(
            std::launch::async, &SequenceClient::FastBatcher, this,
            reader.get(), mate_reader.get(), std::ref(batches_queue), std::ref(pending));

        // take data from queue and send over gRPC
        std::future<int> stream_batches = std::async(
            std::launch::async, &SequenceClient::StreamWriter, this,
            std::ref(in_flight), std::ref(batches_queue), std::ref(pending), std::ref(stream));

        // reading back results on gRPC stream
        std::future<int> recv_reads = std::async(
            std::launch::async, &SequenceClient::StreamReader, this,
            std::ref(in_flight), std::ref(report_file), std::ref(stream),
            std::ref(batches_queue));

        // The stream is over once the reader has seen its end; that also
        // releases the writer. The file reader is waited for last, and only
        // briefly, as it may be blocked on a pipe that has no more data.
        recv_reads.wait();
        stream_batches.wait();
        const int sent = stream_batches.get();
        const int received = recv_reads.get();
        std::cerr << "Sent    : " << sent << std:: endl;
        std::cerr << "Received: " << received << std::endl;

        int code = EX_OK;
        Status status = stream->Finish();
        if (!status.ok()) {
            std::cerr << "Client RPC stream failed: " << status.error_message() << std::endl;
            code = status.error_code();
        }
        else if (read_too_large) {
            code = EX_DATAERR;
        }
        else if (received != sent) {
            // The server completed the stream but dropped results, which
            // it reports on its own log.
            std::cerr << "Error: the server returned " << received << " results for "
                      << sent << " reads; the output is incomplete." << std::endl;
            code = EX_PROTOCOL;
        }

        if (fastq_batches.wait_for(2s) != std::future_status::ready) {
            // The input is still open, for example a pipe with no data, but
            // the stream is over. Exit without waiting for it; the threads
            // cannot be joined while the read blocks.
            std::cerr << "The stream ended while the input was still open; not reading further."
                      << std::endl;
            std::cout.flush();
            std::cerr.flush();
            std::_Exit(code != EX_OK ? code : EX_SOFTWARE);
        }
        if (fastq_batches.get() < 0 && code == EX_OK) {
            code = EX_IOERR;
        }
        if (pairs_mismatched && code == EX_OK) {
            code = EX_DATAERR;
        }
        if (report_failed && code == EX_OK) {
            code = EX_CANTCREAT;
        }
        std::cout.flush();
        if (!std::cout) {
            std::cerr << "Error: could not write the classifications to standard output."
                      << std::endl;
            if (code == EX_OK) {
                code = EX_IOERR;
            }
        }
        return code;
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
        in_flight.add(reads, msg_size);
        bool ok = buffer_hint ? writer->Write(req, WriteOptions().set_buffer_hint())
                              : writer->Write(req);
        if (!ok) {
            in_flight.undo_last();
            StreamFailed();
            return 0;
        }
        return reads;
    }

    int StreamWriter(
            InFlight &in_flight,
            ThreadSafeQueue<std::vector<Kraken2SequenceRequest>> &batches,
            PendingReads &pending,
            ClientStream &writer) {
        int seqs_sent = 0;
        try {
            while (true) {
                std::vector<Kraken2SequenceRequest> batch;
                if (std::optional<std::vector<Kraken2SequenceRequest>> item =
                        batches.pop_wait_for(FLUSH_INTERVAL)) {
                    batch = std::move(*item);
                }
                else if (batches.drained()) {
                    // the file reader has queued everything and finished
                    break;
                }
                else {
                    // nothing arrived for FLUSH_INTERVAL: send what the
                    // file reader holds so far
                    batch = pending.take_after(batches);
                    if (batch.empty()) {
                        continue;
                    }
                }

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
                                // The server cannot receive it. Stop here
                                // rather than drop it, so the output keeps
                                // one line per input read up to this one.
                                std::cerr << single.seqs(0).id() << ": read is larger than the "
                                          << MAX_SIZE / (1024 * 1024) << " MiB message limit and "
                                          << "cannot be classified; stopping." << std::endl;
                                read_too_large = true;
                                stream_closed = true;
                                break;
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
            ClientStream &reader,
            ThreadSafeQueue<std::vector<Kraken2SequenceRequest>> &batches) {
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
                    // one flush per result message, so results reach a
                    // file or pipe as they arrive
                    std::cout.flush();
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
        // or the server went away. Release the writer and the file reader
        // if they are still waiting on each other.
        stream_closed = true;
        in_flight.wake();
        batches.close();
        return n_reads;
    }
    
    /**
     * @brief Read batches of reads (or read pairs when mate_reader is not
     *        null) onto the batches queue.
     *
     * @return number of batches read, or -1 if the input could not be read.
     */
    int FastBatcher(
            FastReader *reader_ptr,
            FastReader *mate_reader,
            ThreadSafeQueue<std::vector<Kraken2SequenceRequest>> &batches_queue,
            PendingReads &pending) {
        int n_batches = 0;
        const bool paired = mate_reader != nullptr;
        FastReader &reader = *reader_ptr;
        const std::string &sequence_file = reader.filename();
        // Whatever happens, the writer must be released.
        struct Closer {
            ThreadSafeQueue<std::vector<Kraken2SequenceRequest>> &q;
            ~Closer() { q.close(); }
        } closer{batches_queue};
        try {
            if (paired) {
                std::cerr << "Reading read pairs from files: " << sequence_file
                          << " and " << mate_reader->filename() << std::endl;
            }
            else {
                std::cerr << "Reading sequences from file: " << sequence_file << std::endl;
            }
            bool mismatched = false;
            bool stopped = false;
            // One read (or pair) at a time into the pending batch, which is
            // queued when full; the writer also takes it when reads arrive
            // slowly.
            while (true) {
                Kraken2SequenceRequest rec;
                bool got = paired ? reader.read_pair(rec, *mate_reader, mismatched) > 0
                                  : reader.read(rec) >= 0;
                if (!got) {
                    break;
                }
                uint64_t rec_bytes = RecordBytes(rec);
                std::vector<Kraken2SequenceRequest> full = pending.add(std::move(rec), rec_bytes);
                if (!full.empty()) {
                    n_batches++;
                    // blocks while MAX_BATCHES are already buffered; false
                    // once the queue has been closed after a failure
                    if (!batches_queue.push(std::move(full))) {
                        stopped = true;
                        break;
                    }
                }
            }
            if (!stopped) {
                std::vector<Kraken2SequenceRequest> rest = pending.take();
                if (!rest.empty()) {
                    n_batches++;
                    batches_queue.push(std::move(rest));
                }
            }
            if (mismatched) {
                std::cerr << "Warning: paired-end inputs have different numbers of reads; "
                          << "only complete pairs were classified: "
                          << sequence_file << ", " << mate_reader->filename() << std::endl;
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

    // Set by StreamWriter when a read exceeds the message limit.
    std::atomic<bool> read_too_large{false};

    // Set by PrintSummary when the report file cannot be written.
    std::atomic<bool> report_failed{false};

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
        if (report_file.empty()) {
            return;
        }
        std::ofstream summary_file(report_file, std::ofstream::out);
        if (summary_file) {
            summary_file << summary;
            summary_file.close();
        }
        if (!summary_file) {
            std::cerr << "Error: could not write the report to " << report_file
                      << ": " << std::strerror(errno) << std::endl;
            report_failed = true;
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
    parser.add({"keepalive", 0, true, "[seconds]", "Ping the server after this long without data and give up if it does not answer within 20 s; 0 disables (default: 60, minimum 10)",
        [&](const std::string &v) {
            opts.keepalive = Parser::ParseInt(v, "--keepalive", 0, 86400);
            if (opts.keepalive != 0 && opts.keepalive < 10) {
                throw std::invalid_argument("--keepalive must be 0 or at least 10 seconds, got '" + v + "'");
            }
        }});
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
    if (opts.keepalive > 0) {
        // Detect a server that stops answering (frozen, or a network that
        // drops without a reset) during a stream. The server accepts pings
        // every 10 s or more without data.
        ch_args.SetInt(GRPC_ARG_KEEPALIVE_TIME_MS, opts.keepalive * 1000);
        ch_args.SetInt(GRPC_ARG_KEEPALIVE_TIMEOUT_MS, KEEPALIVE_TIMEOUT_MS);
        // Current gRPC times a keepalive ping out with the general HTTP/2
        // ping timeout (60 s by default, no public macro) rather than the
        // keepalive timeout.
        ch_args.SetInt("grpc.http2.ping_timeout_ms", KEEPALIVE_TIMEOUT_MS);
        ch_args.SetInt(GRPC_ARG_KEEPALIVE_PERMIT_WITHOUT_CALLS, 0);
        ch_args.SetInt(GRPC_ARG_HTTP2_MAX_PINGS_WITHOUT_DATA, 0);
    }
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
