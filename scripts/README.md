# Overview

Scripts to generate the masked alignments and trees are here.

# Script descriptions

 * `remove_nonb.py`: Filters alignment down to only subtypeB sequences
    * `python remove_nonb.py ../LANL_alignment/HIV1_FLT_2018_genome_DNA.fasta ../results/alignments/`
    * Arguments: path to full alignment, output directory
    * Filters out all sequences from alignment to leave only subtype B (`remove_nonb`)
    * Writes subtype B only alignment to output directory

 * `shuffle_and_mask.py`: Shuffles rows of alignment and outputs masked alignments
    * `python shuffle_and_mask.py ../results/alignments/HIV1_FLT_2018_genome_DNA_subtypeB.fa ../results/alignments/`
    * Arguments: path to alignment, output directory
    * Gets coordinates corresponding to pol/prrt region in HXB2 (`get_HXB2_pol_coords`)
    * Shuffles rows of alignment once
    * Masks alignments down to prrt region for levels 10% to 100% (`mask_around` and `generate_sequences`)
    * Writes masked alignments to output directory