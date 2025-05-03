#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J bwa
#SBATCH -n 1
#SBATCH -c 3
#SBATCH --constraint=skylake
#SBATCH -t 0-10:00:00

set -euo pipefail
IFS=$'\n\t'

# Aligns QC'd reads to the P. abies reference genome
# using bwa and calculates read depth for each 
# sample using samtools

# Get the sample name from the first command-line argument
SAMPLE="$1"
# Exit if no sample was provided
if [ -z "$SAMPLE" ]; then
  echo "Usage: $0 SAMPLE_NAME"
  exit 1
fi

# Set reference and read paths
REF="${SPRUCE_PROJECT}/ref/Pabies1.0-genome.fa"
READS="${SPRUCE_PROJECT}/reads/qcd"

# Set output directory
OUTDIR="${SPRUCE_PROJECT}/bams/full_alignments"
# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Load required modules
ml BWA/0.7.17
ml SAMtools/1.19.2

echo "Starting $SAMPLE at $(date)"

# Align reads and sort the output BAM
bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" "$REF" \
  "${READS}/fastp_${SAMPLE}_1.fq.gz" "${READS}/fastp_${SAMPLE}_2.fq.gz" | \
  samtools sort -o "${OUTDIR}/${SAMPLE}.bam"

# Index the BAM 
samtools index "${OUTDIR}/${SAMPLE}.bam"

# Calculate read depth
echo "Calculating depth for $SAMPLE at $(date)"
samtools depth -q 20 -Q 30 -J -o "${OUTDIR}/${SAMPLE}.depth" "${OUTDIR}/${SAMPLE}.bam"

echo "Finished $SAMPLE at $(date)"