#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J remove_repeats
#SBATCH --output=remove_repeats.out
#SBATCH --error=remove_repeats.err
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00

set -euo pipefail
IFS=$'\n\t'

# Creates a BED file of regions that are:
# 1) On a scaffold with mapped reads
# 2) In or within ±500 bp of an annotated repeat
# 3) At least 1,000 bp long

# Set reference directory path
# Must contain picea_newref.fa.fai and spruce_repeats.bed
REF="${SPRUCE_PROJECT}/ref"

# Load required modules
ml GCC/12.3.0
ml BEDTools/2.31.0

echo "Preparing repeat BED at $(date)"

# Expand each region in spruce_repeats.bed by 500 bp on each side and sort
awk '{
    start = ($2 - 500 >= 0) ? $2 - 500 : 0
    end   = $3 + 500
    print $1 "\t" start "\t" end
}' \
    "${REF}/spruce_repeats.bed" \
	| sort -k1,1 -k2,2n \
	>  "${REF}/spruce_repeats_buffered.bed"

# Merge overlapping buffered regions
bedtools merge \
    -d 1000 \
    -i "${REF}/spruce_repeats_buffered.bed" \
    > "${REF}/spruce_repeats_buffered_merged.bed"

echo "Creating picea_newref.bed from FASTA index at $(date)"

# Create a BED file of all regions in the reduced reference
awk '{print $1 "\t0\t" $2}' "${REF}/picea_newref.fa.fai" \
	| sort -k1,1 -k2,2n > "${REF}/picea_newref.bed"

echo "Subtracting repeats and filtering regions at $(date)"

# Subtract buffered repeat regions and filter out regions < 1,000 bp
bedtools subtract \
	-a "${REF}/picea_newref.bed" \
	-b "${REF}/spruce_repeats_buffered_merged.bed" \
	| awk '{ if (($3 - $2) >= 1000) print }' \
	> "${REF}/picea_newref_target_regions.bed"

echo "Finished at $(date)"