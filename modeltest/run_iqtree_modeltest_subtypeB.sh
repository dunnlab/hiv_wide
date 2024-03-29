#!/bin/bash
#SBATCH -t 2-0 --mem=16G -c16

#module load IQ-TREE/1.6.12-GCCcore-7.3.0-avx512
export SINGULARITY_BINDPATH="/gpfs/data/cbc/aguang/hiv_wide"

WORKDIR=/gpfs/data/cbc/aguang/hiv_wide
SINGULARITY_IMG=${WORKDIR}/metadata/rkantor_hiv.simg
ALIGNMENTS=${WORKDIR}/results/alignments
singularity exec ${SINGULARITY_IMG} iqtree -nt ${SLURM_CPUS_ON_NODE} -s ${ALIGNMENTS}/HIV1_FLT_2018_genome_DNA_subtypeB.fa -m TESTONLY
