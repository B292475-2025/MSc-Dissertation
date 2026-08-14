#!/bin/sh
#$ -cwd
#$ -N denoising
#$ -l h_rt=12:00:00
#$ -l h_vmem=32G
#$ -pe sharedmem 4

echo "loading qiime2"

source /exports/applications/apps/RL9/anaconda/2024.02/etc/profile.d/conda.sh
conda activate rachis-qiime2-2026.4

 echo -e "running the denoising now."

qiime dada2 denoise-paired \
	--i-demultiplexed-seqs trimmed_results_demux.qza \
	--p-trunc-len-f 0 \
	--p-trunc-len-r 0 \
	--o-table results_table.qza \
	--o-representative-sequences results_rep-seqs.qza \
	--o-denoising-stats results_stats.qza \
	--o-base-transition-stats results_transitions.qza

echo -e "filtering out samples with only a few reads (this was actually done after looking at the visualisation file)"
qiime feature-table filter-samples \
	--i-table results_table.qza \
	--p-min-frequency 2000 \
	--o-filtered-table table-filtered.qza 

echo -e "converting to visualisation file"

qiime feature-table summarize \
	--i-table results_table.qza \
	--o-visualization results_table.qzv
	--m-sample-metadata-file results_metadata.tsv

echo -e "Done! Visualisation file is results_metadata.tsv "
