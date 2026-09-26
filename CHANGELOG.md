# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.2.0]
### Added
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
- CI compiles and smoke tests the server and client on every push and runs
  the unit tests.
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
