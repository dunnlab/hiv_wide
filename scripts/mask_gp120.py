#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from os import path
import random

from Bio import SeqIO
from Bio.Seq import Seq

random.seed(42)

# masks gp120 from alignment
# gp120 coordinates on hxb2: 6615, 6812
# to run: python mask_gp120.py ../results/alignments/HIV1_FLT_2018_genome_DNA_subtypeB.fa ../results/alignments

def get_hxb2_coords(orig_matrix, start, stop):
    """
    Get hxb2 coordinates for any given start and stop position
    """
    ungapped_hxb2_start = start  # 0-based
    ungapped_hxb2_stop = stop
    for seqid in orig_matrix:
        if "HXB2" in seqid:
            HXB2 = seqid
            break
    nogap = 0
    hxb2_start = 0
    hxb2_stop = 0
    for i, nt in enumerate(orig_matrix[HXB2]):
        if nogap == ungapped_hxb2_start:
            hxb2_start = i
        elif nogap == ungapped_hxb2_stop:
            hxb2_stop = i
        if hxb2_start > 0 and hxb2_stop > 0:
            break
        if nt != "-":
            nogap += 1
    hxb2_nt = Seq(
        "".join([x for x in orig_matrix[HXB2][hxb2_start:hxb2_stop] if x != "-"])
    )
    return (hxb2_start, hxb2_stop)

def mask_region(seqs_dict, mask_char, start, stop):
    for i in seqs_dict:
        seq = seqs_dict[i]
        yield seq[:start] + len(seq[start:stop]) * mask_char + seq[stop:]

matrix_in = path.abspath(sys.argv[1])
out_dir = path.abspath(sys.argv[2])
out_prefix = path.join(out_dir, path.splitext(path.basename(matrix_in))[0])
mask_character = "-"

orig_matrix = SeqIO.index(matrix_in, "fasta")
gp120_start, gp120_stop = get_hxb2_coords(orig_matrix, 6615, 6812)
shuff_ids_list = [x for x in orig_matrix]

with open(f"{out_prefix}_maskgp120.fa", "w") as masked_matrix:
    SeqIO.write(
        mask_region(
            orig_matrix,
            mask_character,
            gp120_start,
            gp120_stop,
            ),
            masked_matrix,
            "fasta-2line",
        )
