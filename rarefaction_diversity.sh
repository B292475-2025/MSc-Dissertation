#!/bin/sh
#$ -cwd
#$ -N rarefaction
#$ -l h_rt=12:00:00
#$ -l h_vmem=32G
#$ -pe sharedmem 4

echo "loading qiime2"

source /exports/applications/apps/RL9/anaconda/2024.02/etc/profile.d/conda.sh
conda activate rachis-qiime2-2026.4

echo -e "Creating a feature table for rarefaction"
qiime feature-table summarize \
        --i-table table-filtered.qza \
        --o-feature-frequencies feature-frequencies.qza \
        --o-sample-frequencies sample-frequencies.qza \
        --o-summary visual-summary.qzv
echo -e "Done!"

echo -e "starting alpha rarefaction"
qiime diversity alpha-rarefaction \
        --i-table table-filtered.qza \
        --m-metadata-file results_metadata.tsv \
        --o-visualization alpha-rarefaction-curves.qzv \
        --p-min-depth 8509 \
        --p-max-depth 24,564
echo -e "Done!"

echo -e "Diversity calculations"
qiime diversity core-metrics-phylogenetic \
        --i-table table-filtered.qza \
        --i-phylogeny rooted-tree.qza \
	--m-metadata-file results_metadata.tsv \
	--p-sampling-depth 8509 \
	--output-dir core-metrics-results
echo -e "Done!"

echo -e "visualising alpha diversity (shannon diveristy)"
qiime diversity alpha-group-significance \
	--i-alpha-diversity ./core-metrics-results/shannon_vector.qza \
	--m-metadata-file ./results_metadata.tsv \
	--o-visualization ./core-metrics-results/shannon_vector_statistics.qzv
echo -e "Done!"

echo -e "visualising beta diversity (bray curtis)"
qiime diversity beta-group-significance \
	--i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
	--m-metadata-file results_metadata.tsv \
	--m-metadata-column "Malaria status" \
	--o-visualization core-metrics-results/bray_curtis_distance_matrix.qzv
echo -e "Done!"

echo -e "Running the PERMANOVA"
qiime diversity adonis \
	--i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
	--m-metadata-file results_metadata.tsv \
	--o-visualization bray_curtis_adonis_age.qzv \
	--p-formula Malaria_status+Age

qiime diversity adonis \
	--i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
	--m-metadata-file results_metadata.tsv \
	--o-visualization bray_curtis_adonis_sex.qzv \
	--p-formula Malaria_status+Sex

qiime diversity adonis \
        --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
        --m-metadata-file results_metadata.tsv \
        --o-visualization bray_curtis_adonis_house.qzv \
        --p-formula Malaria_status+Household

qiime diversity adonis \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file results_metadata.tsv \
  --p-formula "Malaria_status + Age + Sex + Household" \
  --o-visualization bray_curtis_adonis_household.qzv

echo -e "Done!" 
