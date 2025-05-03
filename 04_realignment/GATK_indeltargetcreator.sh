#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J indel_target
#SBATCH --constraint=skylake|zen3
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00

set -euo pipefail
IFS=$'\n\t'

# Runs RealignerTargetCreator

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

# Set output directory
OUTDIR="${SPRUCE_PROJECT}/bams/realigned/realign_intervals"
# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Use local Java 8 install for compatibility with GATK 3.7
export JAVA_HOME="/pfs/proj/nobackup/fs/projnb10/wanglab-2020/Alex/java/jdk8u392-b08"
export PATH="$JAVA_HOME/bin:$PATH"

# Load required modules
ml GCC/13.2.0
ml SAMtools/1.19.2

echo "Starting RealignerTargetCreator for $SAMPLE at $(date)"

java -jar "${GATK}/GenomeAnalysisTK.jar" \
  -T RealignerTargetCreator \
  -nt 1 \
  -R "$REF" \
  -L ${SPRUCE_PROJECT}/ref/picea_newref_target_regions.bed \
  -I "${INPUT}/${SAMPLE}.intersect.sorted.bam" \
  -o "${OUTDIR}/${SAMPLE}.realign.intervals" \
  --logging_level INFO

echo "Finished RealignerTargetCreator for $SAMPLE at $(date)"
