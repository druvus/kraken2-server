#include <atomic>
#include <cassert>
#include <chrono>
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

#include <zlib.h>
#include "kseq.h"
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
};

typedef std::shared_ptr<ClientReaderWriter<Kraken2SequenceRequestMulti, Kraken2SequenceStreamResult>> ClientStream;


#define ST_BATCH_SIZE 2000     // reads in a gRPC batch
#define MAX_IN_FLIGHT 64000    // total reads in gRPC system
#define FASTQ_BATCH_SIZE 4000  // reads to read at once from fastq
#define MAX_BATCHES 64         // number of stream batches to buffer from fastq

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
     * @return EX_UNAVAILABLE if sequences could nto be sent to server
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
        std::atomic<uint64_t> seqs_in_flight = 0;

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
            std::ref(seqs_in_flight), std::ref(batches_queue), std::ref(stream));

        // reading back results on gRPC stream
        std::future<int> recv_reads = std::async(
            std::launch::async, &SequenceClient::StreamReader, this,
            std::ref(seqs_in_flight), std::ref(report_file), std::ref(stream));

        // wait for things to finish in order
        fastq_batches.wait();
        stream_batches.wait();
        recv_reads.wait();
        std::cerr << "Done waiting" << std::endl;

        std::cerr << "Sent    : " << stream_batches.get() << std:: endl;
        std::cerr << "Received: " << recv_reads.get() << std::endl;
        assert(seqs_in_flight==0);

        // Handle the stream response
        Status status = stream->Finish();
        if (!status.ok()) {
            std::cerr << "Client RPC stream failed: " << status.error_message() << std::endl;
            return status.error_code();
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
            std::cerr << "Failed to send shutdown request." << std::endl;
        }
        if (response.successful()) {
            std::cerr << "Shutdown request processed." << std::endl;
        }
        else{
            std::cerr << "Shutdown request not processed correctly." << std::endl;
        }
        return status.error_code();
    }

    int StreamWriter(
            std::atomic<uint64_t> &seqs_in_flight,
            ThreadSafeQueue<std::vector<Kraken2SequenceRequest>> &batches,
            ClientStream &writer) {
        int seqs_sent = 0;
        try {
            // pop_wait returns an empty optional once the reader has closed
            // the queue and every batch has been taken.
            while (std::optional<std::vector<Kraken2SequenceRequest>> item = batches.pop_wait()) {
                std::vector<Kraken2SequenceRequest> batch = std::move(*item);
                bool show_msg = true;
                while (seqs_in_flight + batch.size() >= MAX_IN_FLIGHT) {
                    if (show_msg) {
                        show_msg = false;
                        std::cerr << "Waiting before sending more. In-flight: " << seqs_in_flight << "." << std::endl;
                    }
                    std::this_thread::sleep_for(10ms);
                }

                // rebatch to smaller batches for stream, moving each read
                // into the message rather than copying it
                const uint64_t MAX_SIZE = 128 * 1024 * 1024;
                for(size_t i = 0; i < batch.size(); i += ST_BATCH_SIZE) {
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
                            if (single.ByteSizeLong() > MAX_SIZE) {
                                std::cerr << "Read is too large! Skipping." << std::endl;
                                continue;
                            }
                            writer->Write(single, WriteOptions().set_buffer_hint());
                            seqs_in_flight.fetch_add(1);
                            seqs_sent++;
                        }
                    }
                    else {
                        writer->Write(req);
                        seqs_in_flight.fetch_add(bsize);
                        seqs_sent += bsize;
                    }
                }
            }
        }
        catch (const std::exception &ex) {
            std::cerr << "Failed to send sequences"
                      << ": " << ex.what() << std::endl;
            writer->WritesDone();
            return seqs_sent;
        }

        writer->WritesDone();
        return seqs_sent;
    }

    int StreamReader(
            std::atomic<uint64_t> &seqs_in_flight, const std::string &report_file,
            ClientStream &reader) {
        Kraken2SequenceStreamResult result;
        int n_reads = 0;
        try {
            while (reader->Read(&result)) {
                if (result.has_classifications()) {
                    for (auto &res : result.classifications().classes()){
                        n_reads++;
                        PrintClassification(res);
                        seqs_in_flight--;
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
            return n_reads;
        }
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
                    n_reads = reader.read_pairs(seqs, FASTQ_BATCH_SIZE, *mate_reader, mismatched);
                }
                else {
                    n_reads = reader.read(seqs, FASTQ_BATCH_SIZE);
                }
                if (n_reads > 0) {
                    n_batches++;
                    // blocks while MAX_BATCHES are already buffered
                    batches_queue.push(std::move(seqs));
                }
                if (n_reads < FASTQ_BATCH_SIZE) { break; }
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

    int WaitForServer() {
        // wait for server
        while (true) {
            ClientContext context;
            Kraken2ReadyRequest req;
            Kraken2ReadyResult response;
            Status status;
            try {
                status = sequence_stub->ServerReady(&context, req, &response);
                if (status.ok())
                {
                    std::cerr << "Server responded as ready." << std::endl;
                    break;
                }
            }
            // complete failure
            catch (const std::exception &ex) {
                std::cerr << "Server status check failed: "
                          << ex.what() << std::endl;
                return EX_UNAVAILABLE;
            }
            // not ready condition
            if (status.error_code() == grpc::StatusCode::UNAVAILABLE) {
                // server may come back
                std::cerr << "Server is not ready: " << status.error_message() << std::endl;
                std::cerr << "Waiting 10s..." << std::endl;
                std::this_thread::sleep_for(10s);
            }
            // unknown error
            else {
                std::cerr << "Server is in error state: "
                          << status.error_message() << std::endl;
                return status.error_code();
            }
        }
        return EX_OK;
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
    std::shared_ptr<grpc::Channel> ch =
        grpc::CreateCustomChannel(
            server_address,
            grpc::InsecureChannelCredentials(), ch_args);
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
