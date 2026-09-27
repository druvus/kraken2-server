#!/bin/bash
set -eo pipefail

# gRPC and protobuf come from the libgrpc and libprotobuf conda packages in
# the host environment, so CMake only needs to be pointed at $PREFIX.
echo "Building kraken2-server in: $PWD"
# A separate directory so a developer's build/ (copied with the source)
# cannot interfere.
rm -rf build-conda
mkdir -p build-conda
pushd build-conda
cmake ${CMAKE_ARGS} \
      -DCMAKE_PREFIX_PATH="${PREFIX}" \
      -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
      -DCMAKE_BUILD_TYPE=Release \
      ..
make -j"${CPU_COUNT:-4}"

mkdir -p "${PREFIX}/bin"
for bin in server/kraken2_server client/kraken2_client; do
    cp "$bin" "${PREFIX}/bin/$(basename "$bin")"
done
popd
