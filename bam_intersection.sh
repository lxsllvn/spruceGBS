#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J intersect
#SBATCH -n 1
#SBATCH -c 6
#SBATCH -t 0-02:00:00

set -euo pipefail
IFS=$'\n\t'

# Intersects BAM files with scaffolds that had mapped reads
# Reduces the number of reference scaffolds and simplifies downstream analysis
# NB: Warnings are suppressed

# Get the sample name from the first command-line argument
SAMPLE="$1"
# Exit if no sample was provided
if [ -z "$SAMPLE" ]; then
  echo "Usage: $0 SAMPLE_NAME"
  exit 1
fi

# Set reference directory; must contain picea_newref.fa and picea_newref_target_regions.bed
REF="${SPRUCE_PROJECT}/ref"

# Set input BAM directory
INPUT="${SPRUCE_PROJECT}/bams/full_alignments"

# Set output directory
OUTDIR="${SPRUCE_PROJECT}/bams/intersected"
# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Load required modules
ml GCC/13.2.0
ml SAMtools/1.19.2

echo "Intersecting $SAMPLE at $(date)"

# Main processing pipeline (with suppressed warnings)
samtools view -h -F 4 -L "${REF}/picea_newref_target_regions.bed" \
	"${INPUT}/${SAMPLE}.bam" 2>/dev/null 
	| grep -v "^@SQ" \
	| samtools view -S -b -T "${REF}/picea_newref.fa" - 2>/dev/null \
	| samtools sort - \
	| samtools addreplacerg \
		-r "ID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
		-o "${OUTDIR}/${SAMPLE}.intersect.sorted.bam" -

# Index the BAM
samtools index "${OUTDIR}/${SAMPLE}.intersect.sorted.bam"

echo "Finished $SAMPLE at $(date)"
