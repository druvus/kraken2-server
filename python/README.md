# kraken2-client (Python)

A small Python client for `kraken2_server`. It speaks the same gRPC protocol
as the C++ client and produces the same per-read output.

## Install

From a checkout of the repository (the gRPC stubs are generated from
`../protos/Kraken2.proto` during installation):

```
pip install --no-build-isolation ./python
```

Requirements: `grpcio`, `protobuf`; `grpcio-tools` at build time. Install
`grpcio-tools` into the same environment first (for example from
conda-forge) and pass `--no-build-isolation`, so the stubs are generated
with a `grpcio-tools` that matches the installed `grpcio`. Stubs generated
by a newer `grpcio-tools` refuse to load with an older `grpcio`.

## Command line

`k2client` takes the same options as `kraken2_client`:

```
k2client --port 8080 --sequence reads.fq.gz --report report.txt > classifications.txt
k2client --port 8080 --sequence reads_1.fq.gz --sequence2 reads_2.fq.gz
k2client --port 8080 --tls-ca server.crt --sequence reads.fq.gz
k2client --port 8080                # print the server summary
k2client --port 8080 --shutdown     # needs --allow-remote-shutdown on the server
```

## Library

```python
from kraken2_client import Client, TLS, read_fastx, read_pairs

with Client("localhost", 8080) as client:
    client.wait_ready()
    with open("report.txt", "w") as report:
        for hit in client.classify(read_fastx("reads.fq.gz"), report=report):
            print(hit.classified, hit.read_id, hit.tax_id, hit.length, hit.hitlist)

    # paired-end
    for hit in client.classify(read_pairs("reads_1.fq.gz", "reads_2.fq.gz")):
        ...

# TLS
client = Client("server.example.org", 8080, TLS(ca="ca.pem"))
```

`classify()` accepts any iterable of `(id, seq, quals)` or
`(id, seq, quals, mate_seq, mate_quals)` tuples, so reads can come from a
running sequencer as easily as from a file. Results are yielded in input
order while the stream is still being sent.

## Tests

```
KRAKEN2_SERVER=../build/server/kraken2_server \
KRAKEN2_CLIENT=../build/client/kraken2_client \
KRAKEN2_DB=../testing/virus-zymo-kraken2 \
pytest python/tests
```

The tests start a server on a free port and compare the Python client's
output with the C++ client's.
