#!/bin/bash
rm -f jobs.txt
ALIGNMENT_PATH=../results/alignments
cd $ALIGNMENT_PATH
for fa in *.fa;
do
    echo "iqtree -nt 8 -mem 16G -o 'K.CD.97.97ZR_EQTB11.AJ249235' -s $fa -m GTR+F+I+G4 -alrt 1000 -bb 1000 -wbt -wbtl" >> jobs.txt
done

module load dSQ
dsq --job-file jobs.txt -C avx512 -t 4-0 -p pi_dunn,general -c8 --mem=16G --batch-file iqtree_hiv_wide.sh --output /dev/null --suppress-stats-file
sed -i '/^# DO NOT/i module load IQ-TREE\/1.6.12-GCCcore-10.2.0-avx512' iqtree_hiv_wide.sh
sbatch iqtree_hiv_wide.sh
