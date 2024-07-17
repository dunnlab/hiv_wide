#!/bin/bash
#SBATCH -N 1
#SBATCH -t 1-0 --mem=16G -c8
#SBATCH --array=0-999
#SBATCH -J iqtree
#SBATCH -e /gpfs/data/cbc/aguang/hiv_wide/logs/fixedparams/iqtree-%A-%a.err
#SBATCH -o /gpfs/data/cbc/aguang/hiv_wide/logs/fixedparams/iqtree-%A-%a.out

# To reproduce: change singularity_bindpath and workdir to appropriate directories (likely repository)

export SINGULARITY_BINDPATH="/gpfs/data/cbc/aguang/hiv_wide"

WORKDIR=/gpfs/data/cbc/aguang/hiv_wide
SINGULARITY_IMG=${WORKDIR}/metadata/rkantor_hiv.simg
ALIGNMENTS=${WORKDIR}/results/alignments

cd $ALIGNMENTS
seeds=(*/)
seed=${seeds[$(( $SLURM_ARRAY_TASK_ID / 10 ))]%/} # values 0-99 for indexing
echo $seed
masks=( 010 020 030 040 050 060 070 080 090 100 )
mask=${masks[$(( $SLURM_ARRAY_TASK_ID % 10 ))]} # values 0-9 for indexing
echo $mask

fa=HIV1_FLT_2018_genome_DNA_subtypeB_${seed}_mask${mask}.fa

mkdir -p ${WORKDIR}/results/trees/fixedparams/${seed}
singularity exec ${SINGULARITY_IMG} iqtree -nt 8 -mem 16G  -s ${ALIGNMENTS}/${seed}/${fa} -m "GTR{2.09420,4.76260,0.89580,0.98320,5.5322}+F{0.363,0.17600,0.23900,0.22200}+I+G4" -alrt 1000 -bb 1000 -wbt -wbtl -pre ${WORKDIR}/results/trees/fixedparams/${seed}/${fa} -a 0.51560 -i 0.10650

#singularity exec ${SINGULARITY_IMG} iqtree -nt 8 -mem 16G  -s ${ALIGNMENTS}/${seed}/${fa} -m GTR{2.09420,4.76260,0.89580,0.98320,5.5322}+I+G4 -alrt 1000 -bb 1000 -wbt -wbtl -pre ${WORKDIR}/results/trees/fixedparams/${seed}/${fa} -a 0.51560 -i 0.10650
