#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from os import path
import random

from Bio import SeqIO
from Bio.Seq import Seq

random.seed(42)


def get_HXB2_pol_coords(orig_matrix):
    """
    This is the ~1300nt region of pol that is commonly sequenced for clinical care (usually called "PRRT" for Protease + the beginning of Reverse Transcriptase):
    HXB2 positions 2253-3554
    """
    ungapped_pol_start = 2252  # 0-based
    ungapped_pol_stop = 3554
    for seqid in orig_matrix:
        if "HXB2" in seqid:
            HXB2 = seqid
            break
    nogap = 0
    pol_start = 0
    pol_stop = 0
    for i, nt in enumerate(orig_matrix[HXB2]):
        if nogap == ungapped_pol_start:
            pol_start = i
        elif nogap == ungapped_pol_stop:
            pol_stop = i
        if pol_start > 0 and pol_stop > 0:
            break
        if nt != "-":
            nogap += 1
    pol_nt = Seq(
        "".join([x for x in orig_matrix[HXB2][pol_start:pol_stop] if x != "-"])
    )
    pol = pol_nt.translate()
    if not pol.startswith("PQVTL") or not pol.endswith("QGQG"):
        print("Couldn't find HXB2's pol, check expected pol positions")
        print(pol_nt, pol)
        sys.exit(1)
    return (pol_start, pol_stop)


def mask_around(seq, mask_char, start, stop):
    left_pad = len(seq[:start]) * mask_char
    right_pad = len(seq[stop:]) * mask_char
    return left_pad + seq[start:stop] + right_pad


def generate_sequences(seqs_dict, id_list, n_to_mask, mask_char, start, stop):
    for i, seq_id in enumerate(id_list):
        if i <= n_to_mask:
            yield mask_around(seqs_dict[seq_id], mask_char, start, stop)
        else:
            yield seqs_dict[seq_id]

matrix_in = path.abspath(sys.argv[1])
out_dir = path.abspath(sys.argv[2])
out_prefix = path.join(out_dir, path.splitext(path.basename(matrix_in))[0])
mask_character = "-"

orig_matrix = SeqIO.index(matrix_in, "fasta")
pol_start, pol_stop = get_HXB2_pol_coords(orig_matrix)

shuff_ids_list = [x for x in orig_matrix]
random.shuffle(shuff_ids_list)
for pct_to_mask in range(1, 11):
    n_to_mask = round(len(shuff_ids_list) * pct_to_mask / 10)
    with open(f"{out_prefix}_mask{pct_to_mask*10:0>3}.fa", "w") as masked_matrix:
        SeqIO.write(
            generate_sequences(
                orig_matrix,
                shuff_ids_list,
                n_to_mask,
                mask_character,
                pol_start,
                pol_stop,
            ),
            masked_matrix,
            "fasta-2line",
        )
