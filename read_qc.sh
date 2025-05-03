#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J bwa
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --constraint=skylake
#SBATCH -t 0-10:00:00

set -euo pipefail
IFS=$'\n\t'

# Trims low-quality bases with fastp

# Get the sample name from the first command-line argument
SAMPLE="$1"
# Exit if no sample was provided
if [ -z "$SAMPLE" ]; then
  echo "Usage: $0 SAMPLE_NAME"
  exit 1
fi

# Set input read path
INPUT="${SPRUCE_PROJECT}/reads"

# Set output directory
OUTDIR="${INPUT}/qcd"
# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Load required modules
ml GCC/13.2.0
ml fastp/0.23.4

echo "Running fastp for $SAMPLE at $(date)"

fastp \
  -i "${INPUT}/${SAMPLE}_1.fq.gz" \
  -I "${INPUT}/${SAMPLE}_2.fq.gz" \
  -o "${OUTDIR}/fastp_${SAMPLE}_1.fq.gz" \
  -O "${OUTDIR}/fastp_${SAMPLE}_2.fq.gz" \
  -h "${OUTDIR}/${SAMPLE}.html" \
  -j "${OUTDIR}/${SAMPLE}.json" \
  --cut_front --cut_right \
  --dont_eval_duplication \
  -p -y -x

echo "Finished fastp for $SAMPLE at $(date)"
