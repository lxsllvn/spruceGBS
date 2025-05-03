#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J scaff_cov
#SBATCH --output=scaff_cov.out
#SBATCH --error=scaff_cov.err
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 0-10:00:00

# Returns a list of scaffolds with >= 5 mapped reads in at least 100 samples.

# Define paths
DEPTHS="${SPRUCE_PROJECT}/bams/full_alignments/read_depths"
SCRATCH="${SPRUCE_PROJECT}/scaff_cov_tmp_$$"

# Create scratch directory and register trap for cleanup
mkdir -p "$SCRATCH"
trap "rm -rf \"$SCRATCH\"" EXIT

echo "Starting scaffold coverage summary at $(date)"

# For each depth file, extract scaffolds with depth >=5
for file in "${DEPTHS}"/*.depth; do
  awk '$3 >= 5 {print $1}' "$file"
done | \
  sort --temporary-directory="$SCRATCH" | \
  uniq -c | \
  awk '$1 >= 100 {print $2}' > scaffolds_with_coverage_min10samples.txt

echo "Finished at $(date)"
