#!/bin/sh
#$ -cwd
#$ -N taxonomic_classification
#$ -l h_rt=12:00:00
#$ -l h_vmem=32G
#$ -pe sharedmem 4

echo "loading qiime2"

source /exports/applications/apps/RL9/anaconda/2024.02/etc/profile.d/conda.sh
conda activate rachis-qiime2-2026.4

echo "Differential abundance analysis with ANCOM-BC2"
qiime composition ancombc2 \
	--i-table table-filtered.qza \
	--m-metadata-file metadata.tsv \
	--p-fixed-effects-formula 'Malaria.status' \
	--o-ancombc2-output ancombc_malaria_status.qza

echo "barplot visualisation"
qiime composition ancombc2-visualizer \
	--i-data ancombc_malaria_status.qza \
	--i-taxonomy results_taxonomy.qza \
	--o-visualization ancombc-barplot.qzv

echo "Done! visualisation file; ancombc-barplot.qzv"

