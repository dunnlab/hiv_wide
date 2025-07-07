#!/bin/bash
#SBATCH -t 1:00:00 --mem=4G
#SBATCH -J subset_mask_and_shuffle
#SBATCH -e /gpfs/data/cbc/aguang/hiv_wide/logs/subset_mask_and_shuffle-%A-%a.err
#SBATCH -o /gpfs/data/cbc/aguang/hiv_wide/logs/subset_mask_and_shuffle-%A-%a.out

# To reproduce: change singularity_bindpath and workdir to appropriate directories (likely repository)

export SINGULARITY_BINDPATH="/gpfs/data/cbc/aguang/hiv_wide"

WORKDIR=/oscar/data/cbc/aguang/hiv_wide
SINGULARITY_IMG=${WORKDIR}/metadata/rkantor_hiv.simg
ALIGNMENTS=${WORKDIR}/results/subset/alignments

mkdir -p ${ALIGNMENTS}

fa=brazil18858.padded.fasta
python ../shuffle_and_mask.py ${WORKDIR}/results/subset/${fa} ${ALIGNMENTS}
python shuffle_and_mask_intermediate.py ${WORKDIR}/results/subset/${fa} ${ALIGNMENTS}
