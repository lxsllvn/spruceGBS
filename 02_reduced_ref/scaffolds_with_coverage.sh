#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J scaff_cov
#SBATCH --output=scaff_cov.out
#SBATCH --error=scaff_cov.err
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 0-10:00:00

set -euo pipefail
IFS=$'\n\t'

# Returns a list of unique scaffolds with >= 5 mapped reads in any sample.

# Define paths
DEPTHS="${SPRUCE_PROJECT}/bams/full_alignments/read_depths"
SCRATCH="${SPRUCE_PROJECT}/scaff_cov_tmp_$$"

# Create scratch directory and register trap for cleanup
mkdir -p "$SCRATCH"
trap "rm -rf \"$SCRATCH\"" EXIT

echo "Starting scaffold coverage summary at $(date)"

# Extract scaffold names with at least 5 reads and ensure uniqueness
awk '$3 >= 5 {print $1}' "${DEPTHS}"/*.depth | \
  sort -u --temporary-directory="$SCRATCH" \
  > "${SPRUCE_PROJECT}/ref/scaffolds_with_coverage.txt"

echo "Finished at $(date)"