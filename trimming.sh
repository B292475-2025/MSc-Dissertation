#!/bin/bash

echo -e "This is the trimming script! The commands in this script will trim the adapter sequences from\n
 your fastq files using cutadapt! (currently trims polg and universal illumina adapter).\n 
This is step 2 after quality checking."

for p in *_R1_001.fastq.gz; do
	p2="${p/_R1_001.fastq.gz/_R2_001.fastq.gz}"
	cutadapt \
		-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
		-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
		-o trimmed_fastq/"$p" \
		-p trimmed_fastq/"$p2" \
		"$p" "$p2"
done

echo -e "Trimming finished!"
	

 
