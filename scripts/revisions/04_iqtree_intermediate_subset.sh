#!/bin/bash
#SBATCH -t 4-0 --mem=16G -c8
#SBATCH --array=0-999
#SBATCH -J subset-iqtree-intermediate
#SBATCH -e /gpfs/data/cbc/aguang/hiv_wide/logs/subset-iqtree/iqtree-intermediate-%A-%a.err
#SBATCH -o /gpfs/data/cbc/aguang/hiv_wide/logs/subset-iqtree/iqtree-intermediate-%A-%a.out

# To reproduce: change singularity_bindpath and workdir to appropriate directories (likely repository)

export SINGULARITY_BINDPATH="/gpfs/data/cbc/aguang/hiv_wide"

WORKDIR=/gpfs/data/cbc/aguang/hiv_wide
#SINGULARITY_IMG=${WORKDIR}/metadata/rkantor_hiv.simg
SINGULARITY_IMG=${WORKDIR}/metadata/rkantor_hiv_06132025.sif
ALIGNMENTS=${WORKDIR}/results/subset/alignments

cd $ALIGNMENTS
seeds=(*/)
seed=${seeds[$(( $SLURM_ARRAY_TASK_ID % 100 ))]%/} # values 0-99 for indexing
masks=( 001 002 003 004 005 006 007 008 009 000 )
mask=${masks[$(( $SLURM_ARRAY_TASK_ID % 10 ))]} # values 0-9 for indexing

fa=brazil18858.padded_${seed}_mask${mask}.fa

# Check if the mask variable is exactly "000"
if [ "$mask" = "000" ]; then
  echo "Mask is 000. Copying over whole genome to ${ALIGNMENTS}/${seed}..."

  cp ${WORKDIR}/results/subset/brazil18858.padded.fa ${ALIGNMENTS}/${seed}/${fa}

fi

mkdir -p ${WORKDIR}/results/subset/trees/${seed}
singularity exec ${SINGULARITY_IMG} iqtree2 -nt 8 -mem 16G -seed ${seed} -s ${ALIGNMENTS}/${seed}/${fa} -m GTR+F+I+G4 -alrt 1000 -bb 1000 -wbt -wbtl -pre ${WORKDIR}/results/trees/${seed}/${fa}
