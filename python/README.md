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

`k2client` takes the same options as `kraken2_client`, and `--wait N` to
give up after N seconds if the server's database is still loading
(default: wait as long as it takes). It exits with the same statuses for
input and output errors (65, 66, 73, 74; see the main README) and 7 when
a remote shutdown is refused. Any other failed request to the server
exits 69.

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

### Records and results

`classify()` accepts any iterable of `(id, seq, quals)` or
`(id, seq, quals, mate_seq, mate_quals)` tuples (`quals` empty for FASTA),
so reads can come from a running sequencer as easily as from a file.
Results are yielded in input order, as `Classification` objects, while
the stream is still being sent; `to_kraken_line()` formats one the way
`kraken2` does.

- **Slow sources.** An iterable that may block (standard input, a pipe,
  or any iterable other than `read_fastx` or `read_pairs` on regular
  files) is read on a separate thread, and records that arrive slowly are
  sent once nothing more has arrived for a second, so each result comes
  back within about a second of its record.
- **Memory.** Records are sent in messages of at most 2000 records or
  8 MB, and at most 256 MB (or 64,000 records) are awaiting results, so
  memory use is bounded whatever the read length. A record over the
  server's 128 MiB message limit raises `ValueError`.
- **Report.** If `report` is a writable text stream, the kraken2 style
  report for the stream is written to it at the end.
- **Flushing.** `on_message`, if given, is called after the results of
  each message from the server have been yielded; `k2client` uses it to
  flush its output, so results reach a file or pipe as they arrive.
- **Errors from the iterable**, such as the `ValueError` from
  `read_pairs` on inputs of different lengths, are raised by `classify()`
  after the results for the records already sent have been yielded and
  the report for them has been written.
- **Errors from the server** are `grpc.RpcError`: `INTERNAL` means the
  server could not classify a batch and the results are incomplete,
  `UNAVAILABLE` that the server went away or stopped answering, and
  `FAILED_PRECONDITION` that its database failed to load.

### Reading files

`read_fastx()` and `read_pairs()` read FASTA or FASTQ, plain or gzipped,
from files, named pipes (including process substitution such as
`<(zcat reads.fq.gz)`), or standard input (`-`, also gzipped). They open
their files at once, so a missing file raises `OSError` from the call
itself. A FASTQ file that ends part way through a record, a truncated or
corrupt gzip file, or a read error raises `InputError` (a `ValueError`)
when it is reached. Bytes that are not valid UTF-8 in read ids are kept:
they are decoded with `surrogateescape` and sent to the server unchanged.
Reading from `-` leaves standard input open.

### Connection

`Client(..., keepalive=60)` pings the server after that many seconds
without data on a stream and ends the stream with `UNAVAILABLE` when no
answer comes within 20 seconds (minimum 10, 0 disables); `k2client`
takes `--keepalive`. `wait_ready()` waits for the server's database to
load and raises `ServerNotReady` after its `timeout`.

## Tests

From the repository root, after building the C++ server and client in
`build/` and installing this package:

```
./testing/parity_test.sh        # downloads the test database and reads on first use
pytest python/tests
```

The tests start servers on free ports and compare the Python client's
output with the C++ client's. They are skipped when the server, the
database or the reads are missing. `KRAKEN2_SERVER`, `KRAKEN2_CLIENT`,
`KRAKEN2_DB` and `KRAKEN2_READS` override the default paths
(`build/server/kraken2_server`, `build/client/kraken2_client`,
`testing/virus-zymo-kraken2`, `testing/virus-zymo-kraken.reads.fastq.gz`).
Python 3.9 or later is needed.
