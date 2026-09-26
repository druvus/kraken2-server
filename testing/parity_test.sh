#!/bin/bash
# Compare kraken2_server/kraken2_client output with the kraken2 command line
# program on the same database and reads, for single-end and paired-end input.
#
# Requires kraken2 (same version as the kraken2/ submodule) and seqkit on the
# PATH, and built binaries in ../build (or set SERVER and CLIENT). The test
# database and reads are downloaded on first use, as in run_server.sh.
#
# Usage: ./parity_test.sh [port]

set -u

cd "$(dirname "$0")"
PORT=${1:-8123}
SERVER=${SERVER:-../build/server/kraken2_server}
CLIENT=${CLIENT:-../build/client/kraken2_client}
DB=virus-zymo-kraken2
READS=virus-zymo-kraken.reads.fastq.gz
WORK=$(mktemp -d)
FAILED=0

for tool in kraken2 seqkit "$SERVER" "$CLIENT"; do
    if ! command -v "$tool" > /dev/null; then
        echo "Missing: $tool"; exit 2
    fi
done

if [ ! -d $DB ]; then
    echo "+++ Downloading database +++"
    curl -sSL -o db.tar.gz https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/misc/virus-zymo-kraken.tar.gz
    tar -xzf db.tar.gz && rm db.tar.gz
fi
if [ ! -f $READS ]; then
    echo "+++ Downloading reads +++"
    curl -sSL -O https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/misc/$READS
fi

# Synthetic pairs: reads 1-1500 as mate 1, 1501-3000 as mate 2, renamed so
# mates share a name with /1 and /2 suffixes. A truncated mate 2 file
# exercises the mismatched-length path.
seqkit range -r 1:1500 $READS 2>/dev/null | seqkit replace -p '.*' -r 'frag{nr}/1' -o $WORK/r1.fq.gz 2>/dev/null
seqkit range -r 1501:3000 $READS 2>/dev/null | seqkit replace -p '.*' -r 'frag{nr}/2' -o $WORK/r2.fq.gz 2>/dev/null
seqkit range -r 1:1000 $WORK/r2.fq.gz -o $WORK/r2_short.fq.gz 2>/dev/null

check() {
    # check <name> <expected> <actual>
    if [ "$2" == "$3" ]; then
        echo "PASS  $1"
    else
        echo "FAIL  $1 (expected '$2', got '$3')"
        FAILED=1
    fi
}

# Number of lines that differ between two files, ignoring order of lines.
ndiff() {
    diff <(sort "$1") <(sort "$2") | grep -c '^[<>]'
}

echo "+++ Reference kraken2 +++"
kraken2 --db $DB --threads 4 --minimum-hit-groups 2 \
    --report $WORK/ref_single.report --output $WORK/ref_single.out $READS > $WORK/ref_single.log 2>&1
check "reference kraken2 single-end run" 0 $?
kraken2 --db $DB --threads 4 --minimum-hit-groups 2 --paired \
    --report $WORK/ref_paired.report --output $WORK/ref_paired.out $WORK/r1.fq.gz $WORK/r2.fq.gz > $WORK/ref_paired.log 2>&1
check "reference kraken2 paired-end run" 0 $?
for f in ref_single.out ref_single.report ref_paired.out ref_paired.report; do
    [ -s $WORK/$f ] || { echo "FAIL  reference output $f missing"; FAILED=1; }
done

echo "+++ Starting server on port $PORT +++"
"$SERVER" --db $DB --host-ip 127.0.0.1 --port $PORT > $WORK/server.log 2>&1 &
SERVER_PID=$!
trap 'kill $SERVER_PID 2> /dev/null' EXIT
sleep 3

echo "+++ Single-end +++"
"$CLIENT" --port $PORT --host-ip 127.0.0.1 --sequence $READS \
    --report $WORK/single.report > $WORK/single.out 2> $WORK/single.err
check "single-end client exit code" 0 $?
check "single-end per-read output identical" 0 "$(ndiff $WORK/single.out $WORK/ref_single.out)"
check "single-end output in input order" 0 "$(diff $WORK/single.out $WORK/ref_single.out | grep -c '^[<>]')"
check "single-end report identical" 0 "$(diff <(tail -n +2 $WORK/single.report) $WORK/ref_single.report | grep -c '^[<>]')"

echo "+++ Paired-end +++"
"$CLIENT" --port $PORT --host-ip 127.0.0.1 --sequence $WORK/r1.fq.gz --sequence2 $WORK/r2.fq.gz \
    --report $WORK/paired.report > $WORK/paired.out 2> $WORK/paired.err
check "paired-end client exit code" 0 $?
check "paired-end per-read output identical" 0 "$(ndiff $WORK/paired.out $WORK/ref_paired.out)"
check "paired-end output in input order" 0 "$(diff $WORK/paired.out $WORK/ref_paired.out | grep -c '^[<>]')"
check "paired-end report identical" 0 "$(diff <(tail -n +2 $WORK/paired.report) $WORK/ref_paired.report | grep -c '^[<>]')"
check "paired-end hit lists contain mate marker" 1500 "$(grep -c '|:|' $WORK/paired.out)"

echo "+++ Mismatched pair files +++"
"$CLIENT" --port $PORT --host-ip 127.0.0.1 --sequence $WORK/r1.fq.gz --sequence2 $WORK/r2_short.fq.gz \
    > $WORK/mism.out 2> $WORK/mism.err
check "mismatched pairs exit code" 65 $?
check "mismatched pairs classify complete pairs only" 1000 "$(wc -l < $WORK/mism.out | tr -d ' ')"

echo "+++ Empty input +++"
: | "$CLIENT" --port $PORT --host-ip 127.0.0.1 --sequence - > /dev/null 2>&1
check "empty input exit code" 0 $?

echo "+++ Shutdown +++"
"$CLIENT" --port $PORT --host-ip 127.0.0.1 --shutdown > /dev/null 2>&1
wait $SERVER_PID
check "server exit code" 0 $?
trap - EXIT

# Translated search against a tiny synthetic protein database, built on
# first use by make_protein_db.sh (needs kraken2-build and python3).
PROT_DB=tinyprot
if command -v kraken2-build > /dev/null && ./make_protein_db.sh $PROT_DB > $WORK/protdb.log 2>&1; then
    echo "+++ Translated search (protein database) +++"
    kraken2 --db $PROT_DB --minimum-hit-groups 2 \
        --report $WORK/ref_prot.report --output $WORK/ref_prot.out $PROT_DB/reads.fq > /dev/null 2>&1
    check "reference kraken2 protein run" 0 $?
    kraken2 --db $PROT_DB --minimum-hit-groups 2 --paired \
        --report $WORK/ref_prot_p.report --output $WORK/ref_prot_p.out $PROT_DB/reads_1.fq $PROT_DB/reads_2.fq > /dev/null 2>&1
    check "reference kraken2 protein paired run" 0 $?
    # the synthetic reads all derive from the library, so kraken2 should classify them
    check "reference classifies every synthetic protein read" 0 "$(awk '$1=="U"' $WORK/ref_prot.out | wc -l | tr -d ' ')"

    "$SERVER" --db $PROT_DB --host-ip 127.0.0.1 --port $PORT > $WORK/server_prot.log 2>&1 &
    SERVER_PID=$!
    trap 'kill $SERVER_PID 2> /dev/null' EXIT
    sleep 2
    "$CLIENT" --port $PORT --host-ip 127.0.0.1 --sequence $PROT_DB/reads.fq \
        --report $WORK/prot.report > $WORK/prot.out 2> $WORK/prot.err
    check "protein single-end client exit code" 0 $?
    check "protein single-end output identical" 0 "$(diff $WORK/prot.out $WORK/ref_prot.out | grep -c '^[<>]')"
    check "protein single-end report identical" 0 "$(diff <(tail -n +2 $WORK/prot.report) $WORK/ref_prot.report | grep -c '^[<>]')"
    "$CLIENT" --port $PORT --host-ip 127.0.0.1 --sequence $PROT_DB/reads_1.fq --sequence2 $PROT_DB/reads_2.fq \
        --report $WORK/prot_p.report > $WORK/prot_p.out 2> $WORK/prot_p.err
    check "protein paired-end client exit code" 0 $?
    check "protein paired-end output identical" 0 "$(diff $WORK/prot_p.out $WORK/ref_prot_p.out | grep -c '^[<>]')"
    check "protein paired-end report identical" 0 "$(diff <(tail -n +2 $WORK/prot_p.report) $WORK/ref_prot_p.report | grep -c '^[<>]')"
    check "server reports translated search" 1 "$(grep -c 'translated search' $WORK/server_prot.log)"
    "$CLIENT" --port $PORT --host-ip 127.0.0.1 --shutdown > /dev/null 2>&1
    wait $SERVER_PID
    check "protein server exit code" 0 $?
    trap - EXIT

    echo "+++ --translated-search with a nucleotide database +++"
    "$SERVER" --db $DB --translated-search --host-ip 127.0.0.1 --port $PORT > $WORK/server_warn.log 2>&1 &
    SERVER_PID=$!
    trap 'kill $SERVER_PID 2> /dev/null' EXIT
    sleep 3
    "$CLIENT" --port $PORT --host-ip 127.0.0.1 --shutdown > /dev/null 2>&1
    wait $SERVER_PID
    trap - EXIT
    check "warning when --translated-search given for nucleotide database" 1 "$(grep -c 'Warning: --translated-search' $WORK/server_warn.log)"
else
    echo "SKIP  translated search test (kraken2-build not available or protein database build failed)"
fi

# Server-side handling of records whose quality string length differs from
# the sequence length. kseq rejects such records in the regular client, so
# they are sent with the raw test client (built with -DBUILD_TEST_TOOLS=ON).
RAW_CLIENT=${RAW_CLIENT:-../build/testing/raw_client}
if [ -x "$RAW_CLIENT" ]; then
    echo "+++ Malformed quality strings with --min-quality +++"
    gzip -dc $READS | awk '/^@af6a6aee/{getline s; getline; getline q;
        printf "good\t%s\t%s\n", s, q;
        printf "badqual\t%s\t%s\n", s, substr(q, 1, length(q) - 10);
        printf "pair_badmate\t%s\t%s\t%s\t%s\n", s, q, s, substr(q, 1, 5);
        printf "after\t%s\t%s\n", s, q; exit }' > $WORK/records.tsv
    "$SERVER" --db $DB --host-ip 127.0.0.1 --port $PORT --min-quality 10 > $WORK/server_mq.log 2>&1 &
    SERVER_PID=$!
    trap 'kill $SERVER_PID 2> /dev/null' EXIT
    sleep 3
    "$RAW_CLIENT" 127.0.0.1:$PORT < $WORK/records.tsv > $WORK/raw.out 2> $WORK/raw.err
    check "raw client exit code" 0 $?
    check "well-formed read classified" C "$(awk '$2=="good"{print $1}' $WORK/raw.out)"
    check "truncated quality read unclassified" "U 0:0" "$(awk '$2=="badqual"{print $1, $5}' $WORK/raw.out)"
    check "pair with truncated mate unclassified" "U 0:0" "$(awk '$2=="pair_badmate"{print $1, $5}' $WORK/raw.out)"
    check "read after malformed records classified" C "$(awk '$2=="after"{print $1}' $WORK/raw.out)"
    check "malformed records logged" 2 "$(grep -c 'reporting as unclassified' $WORK/server_mq.log)"
    "$CLIENT" --port $PORT --host-ip 127.0.0.1 --shutdown > /dev/null 2>&1
    wait $SERVER_PID
    check "server survived malformed records" 0 $?
    trap - EXIT
else
    echo "SKIP  malformed quality test (raw_client not built; use -DBUILD_TEST_TOOLS=ON)"
fi

if [ $FAILED -eq 0 ]; then
    echo "All checks passed."
    rm -rf $WORK
else
    echo "Some checks failed. Outputs kept in $WORK"
fi
exit $FAILED
