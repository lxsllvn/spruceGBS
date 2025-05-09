#!/bin/bash

# Creates a new reference comprising only scaffolds with mapped reads;
# this step makes downstream analyses much smoother

# Set reference directory path;
# must contain Pabies1.0-genome.fa and scaffolds_with_coverage.txt
REF="${SPRUCE_PROJECT}/ref"

# Load required modules
ml GCC/13.2.0
ml SAMtools/1.19.2

echo "Extracting scaffolds at $(date)"

# Extract scaffolds from the full reference to obtain the new reference
xargs samtools faidx "${SPRUCE_PROJECT}/ref/Pabies1.0-genome.fa" \
	< "${SPRUCE_PROJECT}/ref/scaffolds_with_coverage.txt" \
	> "${SPRUCE_PROJECT}/ref/picea_newref.fa"

# ensure .fa is at least 1 s older than the .fai
sleep 1

# Index the new reference
samtools faidx "${SPRUCE_PROJECT}/ref/picea_newref.fa"

echo "Creating picea_newref.bed from FASTA index at $(date)"

# Create a BED file of all regions in the reduced reference
awk '{print $1 "\t0\t" $2}' "${REF}/picea_newref.fa.fai" \
	| sort -k1,1 -k2,2n > "${REF}/picea_newref.bed"