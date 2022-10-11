#!/bin/bash
#SBATCH -J clusterpicker
#SBATCH --time=10:00:00
#SBATCH --mem=8G
#SBATCH --array=0-10
#SBATCH -e /gpfs/data/cbc/aguang/hiv_wide/logs/cluster-picker-%A-%a-%J.err
#SBATCH -o /gpfs/data/cbc/aguang/hiv_wide/logs/cluster-picker-%A-%a-%J.out

export SINGULARITY_BINDPATH="/gpfs/home/$USER,/gpfs/scratch/$USER,/gpfs/data/cbc"
module load java

BASE_DIR=/gpfs/data/cbc/aguang/hiv_wide
SINGULARITY_IMG=${BASE_DIR}/metadata/rkantor_hiv.simg
TREE_DIR=${BASE_DIR}/trees

newgrp cbcollab

masks=( "" _mask010 _mask020 _mask030 _mask040 _mask050 _mask060 _mask070 _mask080 _mask090 _mask100 )

alignment=HIV1_FLT_2018_genome_DNA${masks[$SLURM_ARRAY_TASK_ID]}.fa
tree=${alignment}.treefile

singularity exec $SINGULARITY_IMG java -jar ClusterPicker_command.jar $alignment $tree 0.9 0.99 0.015 10
