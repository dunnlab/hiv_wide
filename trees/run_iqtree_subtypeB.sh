#!/bin/bash
#SBATCH -t 4-0 --mem=16G -c8
#SBATCH --array=0-9
#SBATCH -J iqtree
#SBATCH -e /gpfs/data/cbc/aguang/hiv_wide/logs/%J-%A-%a.err
#SBATCH -o /gpfs/data/cbc/aguang/hiv_wide/logs/%J-%A-%a.out

export SINGULARITY_BINDPATH="/gpfs/data/cbc/aguang/hiv_wide"

WORKDIR=/gpfs/data/cbc/aguang/hiv_wide
SINGULARITY_IMG=${WORKDIR}/metadata/rkantor_hiv.simg
ALIGNMENTS=${WORKDIR}/results/alignments

masks=( 010 020 030 040 050 060 070 080 090 100 )
fa=HIV1_FLT_2018_genome_DNA_subtypeB_mask${masks[$SLURM_ARRAY_TASK_ID]}.fa

singularity exec ${SINGULARITY_IMG} iqtree -nt 8 -mem 16G  -s ${ALIGNMENTS}/${fa} -m GTR+F+I+G4 -alrt 1000 -bb 1000 -wbt -wbtl
