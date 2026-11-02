# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixed
- A stream hung forever with very long reads (about 30 Mb and more, such as
  assembled contigs): the server held results back with a buffering hint
  while the client waited for them before sending more. Results are now
  written as soon as they are next in order.
- A client that stops reading results could make the server hold
  classified batches without limit. A request batch now counts against
  the stream's limit until its results are written.
- A batch that failed on the server stopped the output ordering; it is
  now replaced by a marker so the stream ends with `INTERNAL` promptly.
  Exceptions of any type in classification or result writing fail the
  stream instead of terminating the server, and statistics count only
  results that reached the client.
- SIGTERM and Ctrl-C did not stop the server while a client held a stream
  open, for example one fed by a running sequencer. Open streams now get
  10 seconds and are then cancelled, and a second signal exits at once.
- A port already in use or unreadable TLS files were reported at once,
  but the server exited only after the database had finished loading. The
  database is now loaded after the server is listening, and these errors
  exit immediately with 69 and 78.
- Signals are taken on a dedicated thread with `sigwait` instead of a
  signal handler that was not async-signal-safe.
- `kraken2_client` exited 0 when the report could not be written, when
  standard output could not be written (a full disk), when a plain FASTQ
  file ended part way through a record (the partial record was classified
  as a shorter read), and when it skipped a read over the 128 MiB message
  limit. These now exit 73, 74, 74 and 65; the client stops at an
  oversized read so the output keeps one line per input read.
- `kraken2_client` waited on its input after the server went away, so a
  client reading a quiet pipe never exited. It now exits with the stream's
  status two seconds after the stream ends.
- `kraken2_client` reports a missing input file (exit 66) before
  connecting to the server, closes its input files, counts a message as
  in flight before sending it, and no longer reports a read error in one
  of two paired files as mismatched lengths.
- `k2client` ended in a traceback with exit 1 for a FASTQ file cut off
  part way through a record or a truncated gzip file, a missing input
  file, an unwritable report, and standard output that could not be
  written. It now exits 74, 66, 73 and 74 like `kraken2_client`. `read_fastx` raises `InputError` for the
  truncated files and opens its file at once, so a missing file is
  reported by the call; reading from `-` no longer closes standard input.
- The report passed to `classify()` is written for the records that were
  classified also when the record source fails part way.

### Changed
- The distinct k-mer estimate (a HyperLogLog sketch per taxon) is computed
  only with `--report-kmer`, which is the only report that shows it. The
  server uses about 10% less CPU on short paired reads.
- Per-thread scratch buffers larger than 64 MB are released after the
  read that grew them. After four 90 Mb contigs the server holds 0.8 GB
  instead of about 4.5 GB, and its peak drops from 6.4 to 5.3 GB.
- The parity test compares `--report-kmer` reports with
  `kraken2 --report-minimizer-data`.
- The Python client's in-flight tracker took about 11 µs per result with
  a full window of 64,000 records, because it removed from the front of a
  list; it now uses a deque (0.3 µs).

## [v0.2.1] - 2026-10-28
### Fixed
- The C++ client hung forever when the server went away mid-stream with
  more than 64,000 reads still to send, waiting for results that would
  never arrive. It now stops sending and exits with the gRPC status.
- The Python client failed with `RESOURCE_EXHAUSTED` and no output on long
  reads averaging more than about 32 kb, because it cut request messages
  by record count only and exceeded the server's 128 MB message limit.
  Messages are now also cut at 8 MB of sequence.
- `k2client` on paired-end inputs with different numbers of reads exited
  69 with "Exception iterating requests!" and no results, because the
  error was raised inside the gRPC request generator. It now classifies
  the complete pairs, prints the message and exits 65 like `kraken2_client`.
- A batch that could not be classified on the server was dropped with a
  log line, every later batch was held back, and the client exited 0 with
  fewer lines than reads. The server now fails the stream with `INTERNAL`
  (client exit 13), and the C++ client exits 76 should a completed stream
  ever carry fewer results than reads sent.
### Changed
- Both clients bound their buffering in bytes as well as reads: file
  batches of at most 8 MB, 16 batches read ahead, 256 MB in flight. The
  C++ client's peak memory on a 1.5 GB nanopore input drops from 1.1 GB to
  under 0.2 GB, and throughput against a 4-thread server improves because
  smaller batches pipeline better. Output is unchanged.
- `testing/make_longread_testdata.py` builds Illumina-style pairs and long
  reads from the test reads; the parity script and the Python tests cover
  them, together with the lost-server and failed-batch paths.
- The build notes and CI install `libsystemd` on Linux, which current
  conda-forge `libgrpc` builds link against without depending on.
- The GitHub workflow is a single Linux job (build, unit tests, parity
  suite, Python tests) on pull requests and pushes to `master`, cancelling
  superseded runs. macOS builds are covered by the GitLab conda packaging.
- The GitLab environments install `libsystemd` as well.

## [v0.2.0] - 2026-09-27
### Added
- Classification against several databases at once (`--db` repeated),
  matching `k2 classify --db a,b` in Kraken 2.17: a merged taxonomy is built
  from the databases' taxonomies and each read is called from the LCA of its
  hits across databases. Output is identical to `k2` on the test databases.
  `--strict-merge` applies `--confidence` and `--hit-groups` to the merged
  call, which `k2` does not. Design in `docs/MULTI_DB.md`.
- TLS. The server takes `--tls-cert` and `--tls-key`, and `--tls-ca` to
  require client certificates. The client takes `--tls`, `--tls-ca`,
  `--tls-cert`, `--tls-key` and `--tls-server-name`. Without these the
  server prints a warning and runs unencrypted as before.
- A Python client package (`python/`, `pip install --no-build-isolation
  ./python`) with a `k2client` command taking the same options as the C++
  client, and a library API that streams any iterable of records and yields
  results in input order. Its tests compare its output with the C++ client.
- A GitHub Actions workflow building and testing on Linux and macOS.
- `--allow-remote-shutdown` on the server. The `RemoteShutdown` RPC is
  refused with `PERMISSION_DENIED` unless it is given.
- Paired-end classification. The client accepts a mate file with `--sequence2`
  and the server classifies each pair as one fragment, matching `kraken2 --paired`
  (pooled minimizers, `|:|` marker in the hit list, `len1|len2` length column,
  `/1` and `/2` stripped from the read id). Pairing is decided per read, so
  single-end and paired-end clients can share a server.
### Changed
- Kraken2 submodule updated from 2.1.2 to 2.17.2 (commit 01fb1d9). The server
  now reads databases with 40-bit hash cells and uses batched hash lookups.
  Per-read output and reports are identical to `kraken2` 2.17.2.
- The taxonomy and hash table are loaded in the background after the server
  starts listening, as the client readiness check already assumed.
- The classification thread pool defaults to the number of hardware threads
  instead of one thread.
- The client no longer sends the full header line and a text copy of each
  record, and the server no longer sends the scientific name per read. Both
  fields were unused and roughly doubled the request size.
- Result batches are handed to the writer thread through a blocking queue
  instead of two busy-waiting loops, and read batches are moved rather than
  copied between the reader, the queue and the gRPC messages on both sides.
- The server limits the number of request batches buffered per stream, so a
  fast client cannot make the server hold its whole input in memory.
- The conda recipe and CI use the conda-forge `libgrpc` and `libprotobuf`
  packages instead of building gRPC from source.
- CI compiles and smoke tests the server and client on every push, runs the
  unit tests, and runs the kraken2 parity test with a cached test database.
- Results are returned in the order the reads were sent, as with `kraken2`,
  whatever the number of classification threads.
- The vendored BS::thread_pool header is replaced by a small worker pool
  built on the project's own queue.
- Warnings are errors for the project's own sources in CI (`-DENABLE_WERROR=ON`).
- Server and client share a table-driven option parser; each option is
  declared once, usage text is generated, and error messages name the
  option and accepted range.
- Unit tests (doctest, `-DBUILD_TESTS=ON`) for the pure classification
  helpers, report writer, request conversion and queue. The helpers moved
  from the classifier class into `server/classify_core.cc`.
- The report writer uses the upstream kraken2 functions instead of a local
  copy.
- Usage errors on the command line exit with status 64 (`EX_USAGE`) instead
  of 0.
### Removed
- The unused `thread_pool_light.hpp`. The unused proto fields `header`,
  `str_representation` and `name` are marked deprecated but kept for
  compatibility.
### Fixed
- Translated search against protein databases classified nothing after the
  Kraken2 update, because kraken2 2.1.3 and later require an explicit call to
  initialise the codon tables. The server now calls it when the database is
  a protein database, and the parity test covers translated search with a
  synthetic protein database built by `testing/make_protein_db.sh`.
- A FASTQ record whose quality string length differs from its sequence, with
  `--min-quality` set, is reported as unclassified instead of terminating the
  server. Covered by `testing/parity_test.sh` using the `raw_client` test tool. The client now reports truncated or unreadable records and exits
  with a non-zero status rather than treating them as end of input.
- The classification summary was read without locking while another stream
  could be rewriting it.
- Ctrl-C during database load waits for the loader thread instead of
  destroying the classifier underneath it. A second shutdown signal no longer
  throws.
- Statistics for an empty stream printed NaN percentages.
- The client `-i` / `-I` short options for `--host-ip` were not accepted.
- The client closed an uninitialised file handle when a reader was destroyed.
- The client exits with a non-zero status when an input file cannot be opened
  or when paired-end inputs contain different numbers of reads.

## [v0.1.8]
### Fixed
- Receive size of messages in client raised to accomodate larger requests.

## [v0.1.7]
### Fixed
- Several CLI options did not correctly require a value.

## [v0.1.6]
### Fixed
- Memory mapping CLI option was setting translated search mode not enabling memory mapping!

## [v0.1.5]
### Changed
- Update license file.

## [v0.1.4]
### Fixed
- Explicitly send metadata from server to avoid stall, remove hacky sleep.

## [v0.1.3]
### Added
- Client can read input from stdin.
### Fixed
- Work around a non-understood error.

## [v0.1.2]
### Fixed
- Batching of reads for gRPC stream was incorrect, leading to client errors.

## [v0.1.1]
### Added
- Perform classifications in shared thread pool across clients. Allows speed up of
  processing data from a single client.
### Removed
- The "batch" processing mode of the server and client. This was a stream->unary gRPC,
  that did not achieve any useful functionality beyond the stream->stream gRPC.

## [v0.1.0]
### Fixed
- Memory use of client when processing large input files.

## [v0.0.10]
### Fixed
- Loading of database when user does not provide trailing `/`.

## [v0.0.9]
### Fixed
- Segmentation fault in client when fasta/q comment sections are empty.

## [v0.0.8]
### Added
- --wait option to server fo testing.
### Changed
- client waits if server is active but not ready.

## [v0.0.7]
### Changed
- Server loads kraken2 database asynchronously and will inform clients it is unavailable until loaded.

## [v0.0.6]
### Fixed
- Stop multiple servers running from running on the same port.

## [v0.0.5]
### Added
- Remote shutdown RPC so client can stop server.
- Add IP address option to server.

## [v0.0.4]
### Changed
- kraken2 report data now output to file. 


## [v0.0.3]
First useful release.

### Added
- kraken2_server/client programs.
- conda packaging.
