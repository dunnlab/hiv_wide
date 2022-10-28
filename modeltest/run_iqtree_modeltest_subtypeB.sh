#!/bin/bash
#SBATCH -t 2-0 -p pi_dunn --mem=16G -c16

module load IQ-TREE/1.6.12-GCCcore-7.3.0-avx512
iqtree -nt ${SLURM_CPUS_ON_NODE} -s ../results/alignments/HIV1_FLT_2018_genome_DNA_subtypeB.fasta -m TESTONLY
