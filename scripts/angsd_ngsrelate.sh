#!/bin/bash

prefix=$1

mkdir -p ngsrelate

# Calculate frequencies
zcat angsd/${prefix}_*.mafs.gz | cut -f5 | sed 1d > ngsrelate/${prefix}_freq.txt

# Make merged glf.gz file
zcat angsd/${prefix}_*.glf.gz | gzip > angsd/${prefix}_autosomes.glf.gz

tail -n+2 data/metadata_${prefix}.txt | cut -f 1 > ngsrelate/${prefix}_individuals.txt

ngsRelate -g angsd/${prefix}_autosomes.glf.gz \
	-p ${SLURM_CPUS_ON_NODE} -z ngsrelate/${prefix}_individuals.txt \
	-n $(wc -l ngsrelate/${prefix}_individuals.txt | tr -s ' ' | cut -d ' ' -f 2) \
	-f ngsrelate/${prefix}_freq.txt \
	-O ngsrelate/${prefix}_relatedness.txt

# Copy to data folder
cp ngsrelate/${prefix}_relatedness.txt data/${prefix}_ngsrelate_results.txt