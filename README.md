# Kraken2 Server

Kraken2 is a taxonomic sequence classification system. This project wraps
its classifier in a gRPC server so that a database is loaded once into a
persistent process and reads can be classified as they arrive, from the
same machine or from clients elsewhere. Two uses follow:

* Access to classification from low-resource machines, for example a
  sequencing laptop, by sending reads to a more powerful server.
* One-time loading of a database, a slow step compared to classification,
  with independent classification requests as data becomes available.

The server compiles the Kraken 2.17.2 sources (a git submodule) and
produces per-read output and reports identical to the `kraken2` program on
the test databases used in `testing/parity_test.sh`, for nucleotide and
protein databases, single-end and paired-end reads, and several databases
at once.

## Installation

The server and client are available through the nanoporetech conda channel:

```
mamba create -n kraken2 -c conda-forge -c nanoporetech kraken2-server
```

The Python client is installed from a checkout, see "Python client".

## Usage

Start a server:

```
kraken2_server --db <db_path>
```

where `<db_path>` is a directory containing a standard kraken2 database
(`hash.k2d`, `opts.k2d`, `taxo.k2d`). The server starts listening at once
and loads the database in the background; clients wait until it is ready.
Classification uses all hardware threads unless `--thread-pool` says
otherwise.

Classify reads with the client:

```
kraken2_client --port 8080 --sequence <reads.fq.gz> --report <report.txt> > classifications.txt
```

`<reads.fq.gz>` can be FASTQ or FASTA, plain or gzip compressed; `-` reads
standard input. Paired-end reads are classified as one fragment per pair,
as with `kraken2 --paired`, by giving the mate file as a second input:

```
kraken2_client --port 8080 --sequence <reads_1.fq.gz> --sequence2 <reads_2.fq.gz>
```

The two files must list mates in the same order. If they contain different
numbers of reads, only complete pairs are classified and the client exits
with status 65. Pairing is decided per read, so one server serves
single-end and paired-end clients at the same time.

Without `--sequence` the client prints the server's cumulative report over
every read it has classified since it started (unless the server runs with
`--no-stats`).

### Output

Per-read lines on standard output have the columns of `kraken2` output:

```
C  read_id  562  1683        562:13 0:4 A:2 562:9
C  frag1    1639 761|509     0:727 |:| 0:33 1639:16
```

classified flag, read id (with `/1` and `/2` stripped for pairs), taxonomy
id, sequence length (`len1|len2` for pairs), and the minimizer hit list with
`A` for ambiguous spans, `|:|` between mates and `-:-` between reading
frames. Results arrive in the order the reads were sent, whatever the
number of server threads.

The report written with `--report` is a kraken2 report with one header line
added before the standard columns. `--report-kmer` on the server adds the
k-mer and distinct k-mer columns.

### Server options

| option | meaning |
|---|---|
| `-d, --db PATH` | database directory; repeat for several databases |
| `-x, --thread-pool N` | classification threads shared by all clients (default: all hardware threads) |
| `-r, --max-requests N` | gRPC threads, bounds concurrent client requests (0: gRPC default) |
| `-i, --host-ip ADDR`, `-p, --port N` | listen address and port (default localhost:8080) |
| `-c, --confidence X` | confidence threshold 0 to 1 (default 0) |
| `-g, --hit-groups N` | minimum distinct minimizer hits for a call (default 2) |
| `-q, --min-quality N` | mask bases below this FASTQ quality (default 0, off) |
| `-k, --report-kmer`, `-z, --report-zero` | report k-mer columns, include zero-count taxa |
| `-o, --memory-mapping` | map the database instead of reading it into RAM |
| `-s, --no-stats` | do not keep the server-wide summary |
| `--strict-merge` | with several databases, apply `--confidence` and `--hit-groups` to the merged call |
| `--tls-cert`, `--tls-key`, `--tls-ca` | TLS, see "Security" |
| `--allow-remote-shutdown` | permit `kraken2_client --shutdown` |
| `-w, --wait N` | delay database loading N seconds (testing) |

Short options are accepted in either case (`-d` or `-D`). `--translated-search`
is accepted for compatibility; protein databases are detected from `opts.k2d`.

### Client options

| option | meaning |
|---|---|
| `-s, --sequence PATH` | reads, FASTA or FASTQ, plain or gzipped, `-` for stdin; omit to print the server summary |
| `-2, --sequence2 PATH` | mate file for paired-end reads |
| `-r, --report PATH` | write the kraken2 style report here |
| `-i, --host-ip ADDR`, `-p, --port N` | server address (default localhost:8080) |
| `-k, --shutdown` | ask the server to stop (needs `--allow-remote-shutdown` on the server) |
| `-t, --tls`, `--tls-ca`, `--tls-cert`, `--tls-key`, `--tls-server-name` | TLS, see "Security" |

Exit status is 0 on success, 64 for a usage error, 65 for mismatched
paired-end inputs, 69 when the server cannot be reached, 74 when an input
file cannot be read or parsed, or the gRPC status code of a failed request.
A client started while the server is still loading its database waits.

### Memory and threads

The hash table is read into RAM by default. `--memory-mapping` maps the
files instead, so the operating system page cache is shared with other
processes using the same database (for example the `kraken2` program) and
a large database starts serving before it is fully paged in. Each stream
buffers at most twice the thread count in request batches, so a fast
client cannot make the server hold its whole input in memory.

## Several databases

`--db` may be repeated. Reads are then looked up in every database and the
results combined the way `k2 classify --db a,b` does in Kraken 2.17: a
merged taxonomy is built from the databases' taxonomies, each minimizer's
taxa across databases are reduced to their lowest common ancestor, and the
read is called from the merged hits.

```
kraken2_server --db archaea --db viral --db univec
```

The databases must be built with the same k-mer and minimizer settings and
must all be nucleotide or all protein. Their taxonomies must agree on the
parent of every shared taxid. Memory use is the sum of the databases. With
several databases the merged call is made with confidence 0 and without the
hit-group filter, as kraken2's merge program does; `--confidence` and
`--hit-groups` therefore have no effect on the merged call and the output is
identical to `k2` on the test databases. Add `--strict-merge` to apply both
thresholds to the merged call instead; the result is then filtered the way
a single-database run would be, and can differ from `k2`. Unlike `k2`, the
server needs only the `.k2d` files, not `nodes.dmp` or `seqid2taxid.map`.
See `docs/MULTI_DB.md`.

## Python client

`python/` holds a Python package with the same protocol and output as the
C++ client, for use from scripts, notebooks and pipelines:

```
pip install --no-build-isolation ./python
k2client --port 8080 --sequence reads.fq.gz
```

```python
from kraken2_client import Client, read_fastx
with Client("localhost", 8080) as client:
    for hit in client.classify(read_fastx("reads.fq.gz")):
        print(hit.read_id, hit.tax_id)
```

`classify()` takes any iterable of records, so reads can be streamed from a
running sequencer. See `python/README.md`.

## Security

By default the server listens without TLS and refuses remote shutdown
requests. It is intended for a workstation or a trusted network. For anything
else, enable TLS:

```
kraken2_server --db <db_path> --tls-cert server.crt --tls-key server.key
kraken2_client --port 8080 --tls-ca server.crt --sequence reads.fq.gz
```

`--tls-ca` on the client names the CA bundle (or the self-signed server
certificate) used to verify the server; `--tls` alone uses the system roots.
`--tls-server-name` overrides the name checked in the certificate when it
differs from `--host-ip`. Giving `--tls-ca` to the server as well turns on
mutual TLS: clients must present a certificate signed by that CA with
`--tls-cert` and `--tls-key`.

`kraken2_client --shutdown` only works when the server was started with
`--allow-remote-shutdown`. Otherwise stop the server with Ctrl-C or SIGTERM.

## Building from source

The project needs `cmake` >= 3.13 and a C++17 compiler.

The Kraken2 sources are a git submodule (`kraken2/`, pinned to upstream
commit `01fb1d9`, Kraken 2.17.2) compiled into the server. Clone with
submodules, or fetch them afterwards:

```
git clone --recurse-submodules https://github.com/epi2me-labs/kraken2-server.git
cd kraken2-server
# or, in an existing clone:
git submodule update --init
```

gRPC and protobuf come from the conda-forge `libgrpc` and `libprotobuf`
packages, which is also what the conda recipe in `conda/` and the CI use:

```
mamba create -n k2build -c conda-forge cmake libgrpc libprotobuf zlib llvm-openmp cxx-compiler make
mamba activate k2build
mkdir build && cd build
cmake -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DCMAKE_BUILD_TYPE=Release ..
make -j 8
```

The executables are `build/server/kraken2_server` and
`build/client/kraken2_client`. Useful CMake options: `-DBUILD_TESTS=ON`
(unit tests, needs the `doctest` package), `-DBUILD_TEST_TOOLS=ON`
(`testing/raw_client`), `-DENABLE_WERROR=ON` (warnings are errors in the
project's own sources; on in CI). `conda build conda/` produces the conda
package.

Note that the older `grpc-cpp` conda package pins a protobuf without CMake
configuration files and does not work here.

### Fallback: building gRPC from source

If conda is not an option, gRPC and protobuf can be built from source (see
[gRPC Dependencies for C++](https://grpc.io/docs/languages/cpp/quickstart/))
and passed to CMake through `CMAKE_PREFIX_PATH`:

```
INSTALL_ROOT=$PWD  # or something else
export PROTO_DIR=$INSTALL_ROOT/proto-build
export PATH="$PROTO_DIR/bin:$PATH"
# this flag might be needed on some platforms
export LDFLAGS="-lrt"

mkdir -p $PROTO_DIR
git clone --recurse-submodules -b v1.83.0 --depth 1 --shallow-submodules https://github.com/grpc/grpc
mkdir -p grpc/cmake/build
pushd grpc/cmake/build
cmake -DgRPC_INSTALL=ON \
      -DgRPC_BUILD_TESTS=OFF \
      -DCMAKE_INSTALL_PREFIX=$PROTO_DIR \
      ../..
make -j
make install
popd

mkdir -p build && cd build
cmake -DCMAKE_PREFIX_PATH=${PROTO_DIR} -DCMAKE_BUILD_TYPE=Release ..
make -j 8
```

## Testing

Unit tests use [doctest](https://github.com/doctest/doctest). Configure with
`-DBUILD_TESTS=ON` and run `build/tests/unit_tests` or `ctest`. They cover
the classification helpers, the merged taxonomy, the report writer, request
conversion, the option parser, the queue and worker pool, and the
translation tables, using a small synthetic taxonomy and no database.

`testing/parity_test.sh` is the acceptance test. It compares server and
client output with the `kraken2` and `k2` programs from the same Kraken2
version on small databases:

* the virus-zymo nucleotide test database (downloaded on first use):
  single-end, paired-end, mismatched pairs, empty input, output order;
* a synthetic protein database built by `testing/make_protein_db.sh`:
  translated search, single-end and paired-end;
* two synthetic nucleotide databases sharing a species, built by
  `testing/make_multidb.sh`: multi-database classification against
  `k2 classify --db a,b`, and `--strict-merge`;
* malformed quality strings sent with `testing/raw_client`, TLS, mutual TLS,
  the refused remote shutdown, and clean exit on SIGTERM.

It needs `kraken2`, `kraken2-build`, `k2`, `seqkit`, `openssl` and `python3`
on the `PATH`; sections whose tools are missing are skipped with a notice.
`python/tests` checks the Python client against the C++ client with pytest.

Both CI configurations (`.github/workflows/build.yml` for the GitHub
mirror, `.gitlab-ci.yml` for the ONT GitLab) build with warnings as errors,
run the unit tests and the parity script on every push, and package on
tags.

## Relationship to upstream Kraken2

Most of the classifier is compiled unchanged from the submodule. Four
functions are adapted copies of upstream `classify.cc` and are marked as
such in comments: `ClassifySequence` and `MaskLowQualityBases` in
`server/classify_server.cc`, and `ResolveTree` and `AddHitlistString` in
`server/classify_core.cc`. When the submodule is moved to a new upstream
commit, diff those functions against upstream, rebuild, and run
`testing/parity_test.sh` with a `kraken2` of the same version; the script
is the definition of "still compatible".

## Benchmarks

Server throughput on 73,700 nanopore reads (198 Mbp, the virus-zymo test
reads repeated twenty times) with one client on the same machine, v0.2.0,
MacBook Pro M-series:

| server threads | server time / s | throughput / Gbp per min |
|---|---|---|
| 1 | 10.2 | 1.2 |
| 4 | 3.2 | 3.7 |
| 8 | 2.1 | 5.7 |

Earlier measurements with `testing/run_server.sh` on v0.1.x, kept for
reference:

**Single client** (MacBook Pro 14-inch 2021, M1 Max, 64Gb. macOS 13.2.1.
Clang 13.1.6. 1190.33 Mbp per client)

| clients | server threads | client time / s | server throughput / Mbp/m |
|---------|----------------|-----------------|---------------------------|
|       1 |              2 |            48.8 |                      1462 |
|       1 |              4 |            25.1 |                      2839 |
|       1 |              6 |            18.4 |                      3876 |
|       1 |              8 |            16.8 |                      4237 |

**Multiple clients** (Intel Xeon Gold 6230, Ubuntu 16.04.7, gcc 11.3.0.
991.94 Mbp per client)

| clients | server threads | client time / s | server throughput / Mbp/m |
|---------|----------------|-----------------|---------------------------|
|       1 |             64 |            19.8 |                      3011 |
|       2 |             64 |            21.6 |                      5491 |
|       4 |             64 |            24.1 |                      9802 |
|       8 |             64 |            28.9 |                     16335 |
|      16 |             64 |            62.2 |                     14958 |
