#!/bin/sh
#$ -cwd
#$ -N taxonomic_classification
#$ -l h_rt=12:00:00
#$ -l h_vmem=32G
#$ -pe sharedmem 4

echo "loading qiime2"

source /exports/applications/apps/RL9/anaconda/2024.02/etc/profile.d/conda.sh
conda activate rachis-qiime2-2026.4

qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path results_manifest.tsv \
	--output-path results_demux.qza \
	--input-format SingleEndFastqManifestPhred33V2

echo "creating visualisation file"
qiime demux summarize \
	--i-data results_demux.qza \
	--o-visualization results_demux.qzv

echo "Done"
