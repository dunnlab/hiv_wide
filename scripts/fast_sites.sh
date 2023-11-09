#!/bin/bash

#SBATCH -n 1
#SBATCH --mem 8G
#SBATCH -t 1:00:00

#module load anaconda/2022.05
#conda activate fisher

REPO=/gpfs/data/cbc/aguang/hiv_wide

fast_site_remover.py -s 500 -m $REPO/results/alignments/HIV1_FLT_2018_genome_DNA_subtypeB.fa -tr $REPO/results/trees/HIV1_FLT_2018_genome_DNA_subtypeB.fa.treefile 
