#!/usr/bin/env python3
import sys
from os import path

from Bio import SeqIO
from Bio.Seq import Seq

# python remove_nonb.py ../LANL_alignment/HIV1_FLT_2018_genome_DNA.fasta ../results/alignments

# removes non subtype B sequences
def remove_nonb(orig_matrix):
    for seqid in orig_matrix:
        if seqid.startswith("B."):
            yield orig_matrix[seqid]

matrix_in = path.abspath(sys.argv[1])
out_dir = path.abspath(sys.argv[2])
out_prefix = path.join(out_dir, path.splitext(path.basename(matrix_in))[0]+"_subtypeB")

orig_matrix = SeqIO.index(matrix_in, "fasta")

with open(f"{out_prefix}.fa", "w") as f:
    SeqIO.write(remove_nonb(orig_matrix),f,"fasta-2line")