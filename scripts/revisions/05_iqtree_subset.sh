#!/bin/bash
#SBATCH -t 4-0 --mem=16G -c8
#SBATCH --array=0-999
#SBATCH -A cbc-condo
#SBATCH -J subset-iqtree
#SBATCH -e /gpfs/data/cbc/aguang/hiv_wide/logs/subset-iqtree/iqtree-%A-%a.err
#SBATCH -o /gpfs/data/cbc/aguang/hiv_wide/logs/subset-iqtree/iqtree-%A-%a.out

# To reproduce: change singularity_bindpath and workdir to appropriate directories (likely repository)

export SINGULARITY_BINDPATH="/gpfs/data/cbc/aguang/hiv_wide"

WORKDIR=/gpfs/data/cbc/aguang/hiv_wide
SINGULARITY_IMG=${WORKDIR}/metadata/rkantor_hiv_06132025.sif
ALIGNMENTS=${WORKDIR}/results/subset/alignments

cd $ALIGNMENTS
seeds=(*/)
seed=${seeds[$(( $SLURM_ARRAY_TASK_ID % 100 ))]%/} # values 0-99 for indexing
masks=( 010 020 030 040 050 060 070 080 090 100 )
mask=${masks[$(( $SLURM_ARRAY_TASK_ID % 10 ))]} # values 0-9 for indexing

fa=brazil18858.padded_${seed}_mask${mask}.fa

mkdir -p ${WORKDIR}/results/subset/trees/${seed}
singularity exec ${SINGULARITY_IMG} iqtree2 -nt 8 -mem 16G -seed ${seed} -s ${ALIGNMENTS}/${seed}/${fa} -m GTR+F+I+G4 -alrt 1000 -bb 1000 -wbt -wbtl -pre ${WORKDIR}/results/trees/${seed}/${fa}
