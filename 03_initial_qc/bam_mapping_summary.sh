#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

# Takes a folder of BAM alignments and uses samtools idxstats to
# calculate the number of scaffolds with >= 1 mapped read and
# the total number of mapped reads per sample.

# Check if correct number of arguments is given
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <OUTPUT> <path/to/BAMS>"
    exit 1
fi

# Load required modules
ml GCC/13.2.0
ml SAMtools/1.19.2

# Get the output base name from the first command-line argument
OUTPUT="$1"

# Get the path to the BAM files from the second command-line argument
# Path should not have a trailing slash
BAMS="$2"

# Output header
echo -e "bam_code\tn_scaffolds_with_mapped_reads\ttotal_mapped_reads" > "${OUTPUT}.tsv"

# Loop through each BAM file
for BAM in "${BAMS}"/*.bam; do
    BAM_FILE=$(basename "$BAM")
    BAM_CODE="${BAM_FILE%%.*}"  # Extract text before first dot
    echo "Processing ${BAM_CODE} at $(date)"

    # Use samtools idxstats to count non-zero scaffolds and total mapped reads
    MAPPED_SCAFFOLDS=$(samtools idxstats "$BAM" | awk '$3 > 0' | wc -l)
    TOTAL_MAPPED=$(samtools idxstats "$BAM" | awk '{sum += $3} END {print sum}')

    echo -e "${BAM_CODE}\t${MAPPED_SCAFFOLDS}\t${TOTAL_MAPPED}" >> "${OUTPUT}.tsv"
done

echo "Finished at $(date)"

