#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from os import path
import random

from Bio import SeqIO

random.seed(42)


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

# could have gotten fancy but instead I just found the starts and stops of pol in the alignment by hand
# currently made up, waiting for muscle to finsh
pol_start = 1287
pol_stop = 8422
orig_matrix = SeqIO.index(matrix_in, "fasta-2line")

shuff_ids_list = [x for x in orig_matrix]
random.shuffle(shuff_ids_list)
for pct_to_mask in range(1, 11):
    n_to_mask = round(len(shuff_ids_list) * pct_to_mask / 10)
    with open(f"{out_prefix}_mask{pct_to_mask:0>3}.fa") as masked_matrix:
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
