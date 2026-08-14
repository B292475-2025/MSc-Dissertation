#!/bin/sh
#$ -cwd
#$ -N Classifier
#$ -l h_rt=12:00:00
#$ -l h_vmem=64G
#$ -pe sharedmem 4

echo -e "This script trains the full length classifier used for taxonomic assignment"
echo "loading qiime2 environment"

source /exports/applications/apps/RL9/anaconda/2024.02/etc/profile.d/conda.sh
conda activate rachis-qiime2-2026.4

echo "Training Full-Length Classifier"
qiime feature-classifier fit-classifier-naive-bayes \
	--i-reference-reads silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
	--i-reference-taxonomy silva-138.2-ssu-nr99-tax-derep-uniq.qza \
	--o-classifier silva-138.2-ssu-nr99-classifier.qza

echo "Done"
