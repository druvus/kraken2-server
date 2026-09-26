#!/usr/bin/env python3
"""Generate a tiny protein reference set and matching DNA reads for testing
translated search.

Writes, into the given directory:
  taxonomy/nodes.dmp, taxonomy/names.dmp   a root and three species
  library.faa                              one random protein per species,
                                           headers carry kraken:taxid|N
  reads.fq                                 DNA fragments back-translated from
                                           the proteins with random synonymous
                                           codons, half reverse complemented
  reads_1.fq, reads_2.fq                   the same fragments split into mates

Everything is derived from a fixed random seed, so the files are
reproducible. The sequences are synthetic; the test only needs the server
and the kraken2 program to agree on them.
"""

import argparse
import os
import random

CODONS = {
    "A": ["GCT", "GCC", "GCA", "GCG"], "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "N": ["AAT", "AAC"], "D": ["GAT", "GAC"], "C": ["TGT", "TGC"], "Q": ["CAA", "CAG"],
    "E": ["GAA", "GAG"], "G": ["GGT", "GGC", "GGA", "GGG"], "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"], "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "K": ["AAA", "AAG"], "M": ["ATG"], "F": ["TTT", "TTC"],
    "P": ["CCT", "CCC", "CCA", "CCG"], "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"], "W": ["TGG"], "Y": ["TAT", "TAC"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
}
AMINO_ACIDS = "".join(sorted(CODONS))
TAXA = [(10, "Alpha synthetica"), (20, "Beta synthetica"), (30, "Gamma synthetica")]
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def write_taxonomy(outdir):
    taxdir = os.path.join(outdir, "taxonomy")
    os.makedirs(taxdir, exist_ok=True)
    with open(os.path.join(taxdir, "nodes.dmp"), "w") as nodes, \
         open(os.path.join(taxdir, "names.dmp"), "w") as names:
        nodes.write("1\t|\t1\t|\tno rank\t|\n")
        names.write("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        for taxid, name in TAXA:
            nodes.write("%d\t|\t1\t|\tspecies\t|\n" % taxid)
            names.write("%d\t|\t%s\t|\t\t|\tscientific name\t|\n" % (taxid, name))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("outdir")
    ap.add_argument("--protein-length", type=int, default=600)
    ap.add_argument("--fragment-length", type=int, default=300)
    ap.add_argument("--fragments-per-protein", type=int, default=20)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    write_taxonomy(args.outdir)

    proteins = {}
    with open(os.path.join(args.outdir, "library.faa"), "w") as faa:
        for taxid, name in TAXA:
            protein = "M" + "".join(rng.choice(AMINO_ACIDS) for _ in range(args.protein_length - 1))
            proteins[taxid] = protein
            # kraken2 reads the taxid from a pipe-delimited token in the id
            faa.write(">prot%d|kraken:taxid|%d| %s\n%s\n" % (taxid, taxid, name, protein))

    aa_per_fragment = args.fragment_length // 3
    qual = "I" * args.fragment_length
    half = args.fragment_length // 2
    with open(os.path.join(args.outdir, "reads.fq"), "w") as fq, \
         open(os.path.join(args.outdir, "reads_1.fq"), "w") as fq1, \
         open(os.path.join(args.outdir, "reads_2.fq"), "w") as fq2:
        n = 0
        for taxid, protein in proteins.items():
            for _ in range(args.fragments_per_protein):
                start = rng.randrange(0, len(protein) - aa_per_fragment)
                dna = "".join(rng.choice(CODONS[a]) for a in protein[start:start + aa_per_fragment])
                if rng.random() < 0.5:
                    dna = revcomp(dna)
                n += 1
                name = "frag%d_tax%d" % (n, taxid)
                fq.write("@%s\n%s\n+\n%s\n" % (name, dna, qual))
                fq1.write("@%s/1\n%s\n+\n%s\n" % (name, dna[:half], qual[:half]))
                fq2.write("@%s/2\n%s\n+\n%s\n" % (name, revcomp(dna[half:]), qual[half:]))
    print("wrote %d proteins and %d reads to %s" % (len(proteins), n, args.outdir))


if __name__ == "__main__":
    main()
