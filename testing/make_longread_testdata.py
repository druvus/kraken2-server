#!/usr/bin/env python3
"""Build Illumina-style paired-end and long-read test sets from a FASTQ of
nanopore reads (the virus-zymo test reads), for parity_test.sh and the
Python client tests.

Usage: make_longread_testdata.py reads.fastq.gz outdir

Writes to outdir:
  ill_casava_{1,2}.fq.gz  150 bp pairs with CASAVA 1.8 names (same name in
                          both mates, "1:N:0:..." / "2:N:0:..." comments)
  ill_slash_{1,2}.fq      the same pairs with "/1" and "/2" name suffixes
  ill_short_2.fq.gz       the first 1000 mates of ill_casava_2, for the
                          mismatched-length test
  long_mixed.fq           302 long reads, 6 kb to 5 Mb, concatenated from
                          the source reads
  long_40k.fq             2100 reads of 40 kb or more, so that a batch of
                          2000 reads exceeds the 128 MB gRPC message limit

Some pairs have mates shorter than the k-mer size or runs of N, to cover
the empty hit list and ambiguous span paths. Every record is deterministic
for a given input.
"""
import gzip
import os
import random
import sys


def read_fastq(path):
    reads = []
    with gzip.open(path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                return reads
            seq = f.readline().strip()
            f.readline()
            qual = f.readline().strip()
            reads.append((seq, qual))


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def casava_names(n):
    base = "M00001:7:000000000-A1B2C:1:1101:%d:%d" % (10000 + n, 1000 + n)
    return base + " 1:N:0:ACGTACGT", base + " 2:N:0:ACGTACGT"


def slash_names(n):
    return "SRR000001.%d/1" % (n + 1), "SRR000001.%d/2" % (n + 1)


def write_pairs(reads, path1, path2, namer, limit=None):
    opener = gzip.open if path1.endswith(".gz") else open
    n = 0
    with opener(path1, "wt") as f1, opener(path2, "wt") as f2:
        for seq, qual in reads:
            if len(seq) < 500:
                continue
            if limit is not None and n >= limit:
                break
            r1s, r1q = seq[0:150], qual[0:150]
            r2s, r2q = revcomp(seq[250:400]), qual[250:400][::-1]
            if n % 61 == 5:    # ambiguous bases
                r1s = r1s[:40] + "N" * 10 + r1s[50:]
            if n % 50 == 7:    # mate 2 shorter than k
                r2s, r2q = r2s[:20], r2q[:20]
            if n % 97 == 3:    # both mates shorter than k
                r1s, r1q = r1s[:30], r1q[:30]
                r2s, r2q = r2s[:25], r2q[:25]
            name1, name2 = namer(n)
            f1.write("@%s\n%s\n+\n%s\n" % (name1, r1s, r1q))
            f2.write("@%s\n%s\n+\n%s\n" % (name2, r2s, r2q))
            n += 1
    return n


def write_long(reads, path, lengths, rng):
    """Each output read is a concatenation of source reads drawn at random
    until it is at least the requested length."""
    name = os.path.splitext(os.path.basename(path))[0]
    total = 0
    with open(path, "w") as f:
        for i, target in enumerate(lengths):
            parts, quals, length = [], [], 0
            while length < target:
                seq, qual = rng.choice(reads)
                parts.append(seq)
                quals.append(qual)
                length += len(seq)
            f.write("@%s_%d len=%d\n%s\n+\n%s\n" % (name, i, length, "".join(parts), "".join(quals)))
            total += length
    return len(lengths), total


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    reads = read_fastq(src)
    rng = random.Random(7)

    n = write_pairs(reads, os.path.join(out, "ill_casava_1.fq.gz"), os.path.join(out, "ill_casava_2.fq.gz"), casava_names)
    write_pairs(reads, os.path.join(out, "ill_slash_1.fq"), os.path.join(out, "ill_slash_2.fq"), slash_names)
    write_pairs(reads, os.path.join(out, "ill_short_1.fq.gz"), os.path.join(out, "ill_short_2.fq.gz"), casava_names, limit=1000)
    os.remove(os.path.join(out, "ill_short_1.fq.gz"))
    print("Illumina-style pairs:", n)

    lengths = [int(rng.lognormvariate(11.0, 0.9)) for _ in range(300)] + [2_000_000, 5_000_000]
    print("long_mixed: %d reads, %d bp" % write_long(reads, os.path.join(out, "long_mixed.fq"), lengths, rng))
    print("long_40k: %d reads, %d bp" % write_long(reads, os.path.join(out, "long_40k.fq"), [40_000] * 2100, rng))


if __name__ == "__main__":
    main()
