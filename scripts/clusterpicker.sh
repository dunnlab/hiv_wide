#!/bin/bash
#SBATCH -t 4-0 --mem=16G -c8
#SBATCH --array=0-10
#SBATCH -J clusterpicker
#SBATCH -e /gpfs/data/cbc/aguang/hiv_wide/logs/%J.err
#SBATCH -o /gpfs/data/cbc/aguang/hiv_wide/logs/%J.out

export SINGULARITY_BINDPATH="/gpfs/data/cbc/aguang/hiv_wide"

WORKDIR=/gpfs/data/cbc/aguang/hiv_wide
SINGULARITY_IMG=${WORKDIR}/metadadta/rkantor_hiv.simg
ALIGNMENTS=${WORKDIR}/results/alignments
TREES=${WORKDIR}/results/trees

masks=( 000 010 020 030 040 050 060 070 080 090 100 )
fa=HIV1_FLT_2018_genome_DNA_subtypeB_mask${masks[$SLURM_ARRAY_TASK_ID]}.fa
tree=HIV1_FLT_2018_genome_DNA_subtypeB_mask${masks[$SLURM_ARRAY_TASK_ID]}.treefile

java –jar /home/rstudio/ClusterPicker_1.2.jar $fa $tree 0.9 0.99 0.015 10


