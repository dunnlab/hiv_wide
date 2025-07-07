#!/bin/bash
#SBATCH -t 1:00:00 --mem=4G
#SBATCH -J mask_and_shuffle_intermediate
#SBATCH -e /gpfs/data/cbc/aguang/hiv_wide/logs/mask_and_shuffle-intermediate-%A-%a.err
#SBATCH -o /gpfs/data/cbc/aguang/hiv_wide/logs/mask_and_shuffle-intermediate-%A-%a.out

# To reproduce: change singularity_bindpath and workdir to appropriate directories (likely repository)

export SINGULARITY_BINDPATH="/gpfs/data/cbc/aguang/hiv_wide"

WORKDIR=/gpfs/data/cbc/aguang/hiv_wide
SINGULARITY_IMG=${WORKDIR}/metadata/rkantor_hiv.simg
ALIGNMENTS=${WORKDIR}/results/alignments

fa=HIV1_FLT_2018_genome_DNA_subtypeB.fa
python shuffle_and_mask_intermediate.py ${ALIGNMENTS}/${fa} ${ALIGNMENTS}