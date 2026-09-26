#!/bin/bash
# Build two tiny nucleotide kraken2 databases sharing one species, for
# testing multi-database classification.
#
# Needs kraken2-build and python3 on the PATH. Nothing is downloaded.
# The seqid2taxid.map and taxonomy/ directory are kept in each database so
# that `k2 classify --db A,B` can build its merged taxonomy from them.
#
# Usage: ./make_multidb.sh [dir]   (default: multidb; databases dir/dbA, dir/dbB)

set -euo pipefail

cd "$(dirname "$0")"
DIR=${1:-multidb}

if [ -f "$DIR/dbA/hash.k2d" ] && [ -f "$DIR/dbB/hash.k2d" ]; then
    echo "+++ Multi-database test databases in $DIR already built +++"
    exit 0
fi

echo "+++ Generating synthetic genomes and reads +++"
rm -rf "$DIR"
python3 make_multidb_testdata.py "$DIR"

for db in dbA dbB; do
    echo "+++ Building $db +++"
    mkdir -p "$DIR/$db"
    cp -r "$DIR/taxonomy" "$DIR/$db/taxonomy"
    kraken2-build --db "$DIR/$db" --add-to-library "$DIR/$db.fna" --no-masking > "$DIR/$db/build.log" 2>&1
    kraken2-build --db "$DIR/$db" --build --threads 2 >> "$DIR/$db/build.log" 2>&1
    ls "$DIR/$db"/*.k2d "$DIR/$db"/seqid2taxid.map > /dev/null
done
echo "+++ Done: $DIR/dbA $DIR/dbB +++"
