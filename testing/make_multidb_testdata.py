#!/usr/bin/env python3
"""Generate two small nucleotide reference sets that overlap in one species,
for testing multi-database classification.

Writes, into the given directory:
  taxonomy/nodes.dmp, names.dmp   root -> genus G (5) -> species 10, 20, 30
  dbA.fna                         genomes of species 10 and 20
  dbB.fna                         genomes of species 20 and 30
  reads.fq                        fragments from all three genomes
  reads_1.fq, reads_2.fq          the same fragments as mates

Genomes are random sequences from a fixed seed. Headers carry the taxid in
the kraken:taxid form that kraken2-build understands.
"""

import argparse
import os
import random

COMPLEMENT = str.maketrans("ACGT", "TGCA")
GENUS = 5
SPECIES = [(10, "Alpha synthetica"), (20, "Beta synthetica"), (30, "Gamma synthetica")]


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def write_taxonomy(outdir):
    taxdir = os.path.join(outdir, "taxonomy")
    os.makedirs(taxdir, exist_ok=True)
    with open(os.path.join(taxdir, "nodes.dmp"), "w") as nodes, \
         open(os.path.join(taxdir, "names.dmp"), "w") as names:
        nodes.write("1\t|\t1\t|\tno rank\t|\n")
        names.write("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        nodes.write("%d\t|\t1\t|\tgenus\t|\n" % GENUS)
        names.write("%d\t|\tSynthetica\t|\t\t|\tscientific name\t|\n" % GENUS)
        for taxid, name in SPECIES:
            nodes.write("%d\t|\t%d\t|\tspecies\t|\n" % (taxid, GENUS))
            names.write("%d\t|\t%s\t|\t\t|\tscientific name\t|\n" % (taxid, name))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("outdir")
    ap.add_argument("--genome-length", type=int, default=20000)
    ap.add_argument("--fragment-length", type=int, default=300)
    ap.add_argument("--fragments-per-genome", type=int, default=40)
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    write_taxonomy(args.outdir)

    genomes = {}
    for taxid, _ in SPECIES:
        genomes[taxid] = "".join(rng.choice("ACGT") for _ in range(args.genome_length))

    def write_db(filename, taxids):
        with open(os.path.join(args.outdir, filename), "w") as f:
            for taxid in taxids:
                f.write(">genome%d|kraken:taxid|%d| %s\n" % (taxid, taxid, dict(SPECIES)[taxid]))
                g = genomes[taxid]
                for i in range(0, len(g), 80):
                    f.write(g[i:i + 80] + "\n")

    write_db("dbA.fna", [10, 20])
    write_db("dbB.fna", [20, 30])

    qual = "I" * args.fragment_length
    half = args.fragment_length // 2
    n = 0
    with open(os.path.join(args.outdir, "reads.fq"), "w") as fq, \
         open(os.path.join(args.outdir, "reads_1.fq"), "w") as fq1, \
         open(os.path.join(args.outdir, "reads_2.fq"), "w") as fq2:
        for taxid, _ in SPECIES:
            g = genomes[taxid]
            for _ in range(args.fragments_per_genome):
                start = rng.randrange(0, len(g) - args.fragment_length)
                dna = g[start:start + args.fragment_length]
                if rng.random() < 0.5:
                    dna = revcomp(dna)
                n += 1
                name = "frag%d_tax%d" % (n, taxid)
                fq.write("@%s\n%s\n+\n%s\n" % (name, dna, qual))
                fq1.write("@%s/1\n%s\n+\n%s\n" % (name, dna[:half], qual[:half]))
                fq2.write("@%s/2\n%s\n+\n%s\n" % (name, revcomp(dna[half:]), qual[half:]))
    print("wrote 2 reference sets and %d reads to %s" % (n, args.outdir))


if __name__ == "__main__":
    main()
