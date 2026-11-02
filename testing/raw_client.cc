// Minimal test client that sends hand-built sequence records to a
// kraken2_server without going through a FASTA/FASTQ parser. It exists to
// exercise server code paths that the regular client cannot reach, such as
// records whose quality string length differs from the sequence length.
//
// Input on stdin, one record per line, tab separated:
//   id  seq  quals  [mate_seq  mate_quals]
// An empty quals field sends the record as FASTA. All records are sent as
// one batch. Results are printed as: C/U  id  taxid  size[|size2]  hitlist
//
// Usage: raw_client <host:port> [--hold SECONDS] < records.tsv
//
// With --hold, the stream is kept open after the batch is sent, and the
// client fails (exit 3) unless the first result arrives within SECONDS.
// A server that holds results back until the client sends more would make
// a real client wait forever.

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <thread>
#include <string>
#include <vector>

#include <grpc++/channel.h>
#include <grpc++/client_context.h>
#include <grpc++/create_channel.h>
#include <grpc++/security/credentials.h>

#include "Kraken2.grpc.pb.h"

using kraken2proto::Kraken2SequenceRequest;
using kraken2proto::Kraken2SequenceRequestMulti;
using kraken2proto::Kraken2SequenceStreamResult;
using kraken2proto::Kraken2Service;

static void Fill(Kraken2SequenceRequest &rec, const std::string &id,
                 const std::string &seq, const std::string &quals) {
    rec.set_id(id);
    rec.set_seq(seq);
    if (quals.empty()) {
        rec.set_format(Kraken2SequenceRequest::FORMAT_FASTA);
    } else {
        rec.set_format(Kraken2SequenceRequest::FORMAT_FASTQ);
        rec.set_quals(quals);
    }
}

int main(int argc, char **argv) {
    int hold_seconds = -1;
    if (argc == 4 && std::string(argv[2]) == "--hold") {
        hold_seconds = std::atoi(argv[3]);
    } else if (argc != 2) {
        std::cerr << "Usage: raw_client <host:port> [--hold SECONDS] < records.tsv" << std::endl;
        return 2;
    }

    Kraken2SequenceRequestMulti batch;
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        std::vector<std::string> f;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, '\t')) f.push_back(field);
        while (f.size() < 5) f.push_back("");
        Kraken2SequenceRequest *rec = batch.add_seqs();
        Fill(*rec, f[0], f[1], f[2]);
        if (!f[3].empty()) {
            Fill(*rec->mutable_mate(), f[0], f[3], f[4]);
        }
    }

    auto channel = grpc::CreateChannel(argv[1], grpc::InsecureChannelCredentials());
    auto stub = Kraken2Service::NewStub(channel);
    grpc::ClientContext context;
    auto stream = stub->ClassifyStream(&context);
    stream->Write(batch);

    std::atomic<bool> got_result{false};
    std::thread reader([&]() {
        Kraken2SequenceStreamResult result;
        while (stream->Read(&result)) {
            if (!result.has_classifications()) continue;
            got_result = true;
            for (const auto &r : result.classifications().classes()) {
                std::cout << (r.classified() ? "C" : "U") << '\t' << r.id() << '\t'
                          << r.tax_id() << '\t' << r.size();
                if (r.paired()) std::cout << '|' << r.size2();
                std::cout << '\t' << r.hitlist() << '\n';
            }
        }
    });

    if (hold_seconds >= 0) {
        auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(hold_seconds);
        while (!got_result && std::chrono::steady_clock::now() < deadline) {
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
        if (!got_result) {
            std::cerr << "No result within " << hold_seconds
                      << " s while the stream was open." << std::endl;
            context.TryCancel();
            reader.join();
            return 3;
        }
    }
    stream->WritesDone();
    reader.join();
    grpc::Status status = stream->Finish();
    if (!status.ok()) {
        std::cerr << "RPC failed: " << status.error_message() << std::endl;
        return 1;
    }
    return 0;
}
