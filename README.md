# Kraken2 Server

Kraken2 is a taxonomic sequence classification system. This project builds
on the classification functionality to provide a server-client architecture
to allow two use cases:

* Access to the classification algorithms in low-resource settings by
  sending requests to remote, more powerful servers
* One-time loading of databases (a slow step compared to classification)
  into a persistent process, with subsequent independent classification
  requests as data becomes available

The software is currently a public beta release.

## Installation

The kraken2 server and client are available through our conda channel:

```
mamba create -n kraken2 -c conda-forge -c nanoporetech kraken2-server
```

## Usage

To start a server run:

```
kraken2_server --db <db_path>
```

where `<db_path>` is a directory containing a standard kraken2 database. The
server will wait for requests for clients and respond as necessary.

To classify reads run a client with:

```
kraken2_client --port 8080 --sequence <reads.fq.gz>
```

where `<reads.fq.gz>` can be FASTQ or FASTA either plain text or gzip compressed.
Use `--sequence -` to read from standard input.

Paired-end reads are classified as one fragment per pair, as with
`kraken2 --paired`, by giving the mate file as a second input:

```
kraken2_client --port 8080 --sequence <reads_1.fq.gz> --sequence2 <reads_2.fq.gz>
```

The two files must list mates in the same order. If they contain different
numbers of reads, only complete pairs are classified and the client exits
with a non-zero status. Pairing is decided per read, so a single server can
serve single-end and paired-end clients at the same time.

The per-read output has the same columns as the standard `kraken2` output
(classified flag, read id, taxonomy id, sequence length or `len1|len2` for
pairs, and the minimizer hit list), in the same order as the input reads.
The report file adds a header line before the standard kraken2 report
columns.


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
`--hit-groups` therefore have no effect. Output is identical to `k2` on the
test databases (`testing/parity_test.sh`). See `docs/MULTI_DB.md`.

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

See `python/README.md` for details.

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

The project can be built with `cmake` >3.13 and a C++17 compiler.

The Kraken2 sources are included as a git submodule (`kraken2/`, pinned to
upstream commit `01fb1d9`, Kraken 2.17.2) and are compiled into the server.
Clone with submodules, or fetch them afterwards:

```
git clone --recurse-submodules https://github.com/epi2me-labs/kraken2-server.git
cd kraken2-server
# or, in an existing clone:
git submodule update --init
```

The server-client architecture uses gRPC and protobuf to communicate. The
recommended way to obtain them is the conda-forge `libgrpc` and `libprotobuf`
packages, which is also what the conda recipe in `conda/` and the CI use:

```
mamba create -n k2build -c conda-forge cmake libgrpc libprotobuf zlib llvm-openmp cxx-compiler make
mamba activate k2build
mkdir build && cd build
cmake -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DCMAKE_BUILD_TYPE=Release ..
make -j 8
```

The server and client executables are written to:

```
build/server/kraken2_server
build/client/kraken2_client
```

Note that the older `grpc-cpp` conda package pins a protobuf without CMake
configuration files and does not work here. `conda build conda/` produces
the package that the conda channel distributes.

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

Unit tests use [doctest](https://github.com/doctest/doctest) (conda-forge
package `doctest`). Configure with `-DBUILD_TESTS=ON` and run
`build/tests/unit_tests` or `ctest`. They cover the pure classification
helpers (taxon resolution, hit list formatting, pair name trimming), the
report writer, the request conversion and the blocking queue, using a small
synthetic taxonomy and no database.

`testing/parity_test.sh` compares the server and client output with the
`kraken2` command line program (the same version as the submodule) on a small
database, for single-end, paired-end, mismatched and empty input. It needs
`kraken2` and `seqkit` on the `PATH` and built binaries in `build/`.
When `kraken2-build` is available the script also builds a tiny synthetic
protein database with `testing/make_protein_db.sh` and checks translated
search against `kraken2` the same way.
Configuring with `-DBUILD_TEST_TOOLS=ON` also builds `testing/raw_client`,
which sends hand-built records and lets the script check how the server
handles FASTQ records with a truncated quality string.

## Benchmarks

Benchmarking script from `testing/run_server.sh`

**Single client test**

*MacBook Pro 14-inch 2021, M1 Max, 64Gb. macOS 13.2.1. Clang 13.1.6. 1190.33 Mbp per client*

| clients | server threads | client time / s | server throughput / Mbp/m |
|---------|----------------|-----------------|---------------------------|
|       1 |              2 |            48.8 |                      1462 |
|       1 |              4 |            25.1 |                      2839 |
|       1 |              6 |            18.4 |                      3876 |
|       1 |              8 |            16.8 |                      4237 |

**Multi client test**

*Intel Xeon Gold 6230, Ubuntu 16.04.7, gcc 11.3.0. 991.94 Mbp per client*

| clients | server threads | client time / s | server throughput / Mbp/m |
|---------|----------------|-----------------|---------------------------|
|       1 |             64 |            19.8 |                      3011 |
|       2 |             64 |            21.6 |                      5491 |
|       4 |             64 |            24.1 |                      9802 |
|       8 |             64 |            28.9 |                     16335 |
|      16 |             64 |            62.2 |                     14958 |

