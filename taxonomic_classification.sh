#!/bin/sh
#$ -cwd
#$ -N taxonomic_classification
#$ -l h_rt=12:00:00
#$ -l h_vmem=32G
#$ -pe sharedmem 4

echo "loading qiime2"

source /exports/applications/apps/RL9/anaconda/2024.02/etc/profile.d/conda.sh
conda activate rachis-qiime2-2026.4

echo "beginning the classification"
qiime feature-classifier classify-sklearn \
	--i-classifier silva-138.2-ssu-nr99-classifier.qza \
	--i-reads results_rep-seqs.qza \
	--o-classification results_taxonomy.qza

echo "converting to a visualisation file"
qiime metadata tabulate \
	--m-input-file results_taxonomy.qza \
	--o-visualization results_taxonomy.qzv

echo "creating a barplot for taxonomy and visualisations"
qiime taxa barplot \
	--i-table table-filtered.qza \
	--i-taxonomy results_taxonomy.qza \
	--metadata-file results_metadata.tsv \
	--o-visualization taxa-bar-plots.qzv

echo "creating a rooted tree for downstream diversity analysis"
qiime phylogeny align-to-tree-mafft-fasttree \
	--i-sequences results_rep-seqs.qza \
	--o-rooted-tree rooted-tree.qza \
	--o-alignment aligned-rep-seqs.qza \
	--o-masked-alignment masked-aligned-rep-seqs.qza \
	--o-tree unrooted-tree.qza

echo -e "Classification Done! Visualisation file: results_taxonomy.qzv, barplot: taxa-bar-plots.qzv, rooted tree: rooted-tree.qza."



