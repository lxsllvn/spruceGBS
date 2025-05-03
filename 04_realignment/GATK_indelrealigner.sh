#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J indel_realignment
#SBATCH --constraint=skylake|zen3
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00

set -euo pipefail
IFS=$'\n\t'

# Runs IndelRealigner and indexes the realigned BAM

# Get the sample name from the first command-line argument
SAMPLE="$1"
# Exit if no sample was provided
if [ -z "$SAMPLE" ]; then
  echo "Usage: $0 SAMPLE_NAME"
  exit 1
fi

# Set paths
GATK="/pfs/proj/nobackup/fs/projnb10/wanglab-2020/Alex/alignments"
REF="${SPRUCE_PROJECT}/ref/picea_newref.fa"
INPUT="${SPRUCE_PROJECT}/bams/intersected"
TARGETS="${SPRUCE_PROJECT}/bams/realigned/realign_intervals"

# Set output directory
OUTDIR="${SPRUCE_PROJECT}/bams/realigned"

# Use local Java 8 install for compatibility with GATK 3.7
export JAVA_HOME="/pfs/proj/nobackup/fs/projnb10/wanglab-2020/Alex/java/jdk8u392-b08"
export PATH="$JAVA_HOME/bin:$PATH"

# Load required modules
ml GCC/13.2.0
ml SAMtools/1.19.2

echo "Starting IndelRealigner for $SAMPLE at $(date)"

java -jar "${GATK}/GenomeAnalysisTK.jar" \
  -T IndelRealigner \
  -R "$REF" \
  -entropy 0.05 \
  -LOD 3 \
  -I "${INPUT}/${SAMPLE}.intersect.sorted.bam" \
  -targetIntervals "${TARGETS}/${SAMPLE}.realign.intervals" \
  -o "${OUTDIR}/${SAMPLE}.realigned.bam"

echo "Indexing realigned BAM for $SAMPLE at $(date)"
samtools index "${OUTDIR}/${SAMPLE}.realigned.bam"

echo "Finished $SAMPLE at $(date)"
