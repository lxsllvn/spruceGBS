#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J depth_summary
#SBATCH --output=depth_summary.out
#SBATCH --error=depth_summary.err
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 0-10:00:00

set -euo pipefail
IFS=$'\n\t'

# Creates a sample-level summary of the output of samtools depth.
# Returns a tab-delimited file with the sample name, number of
# unique scaffolds with >= 5 reads, and the total number of mapped reads.

# Check if correct number of arguments is given
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <OUTPUT> <path/to/readdepths>"
    exit 1
fi

# Get the output base name from the first command-line argument
OUTPUT="$1"

# Get the path to the *.depth files from the second command-line argument
DEPTHS="$2"

# Print header
echo -e "bam_code\tn_scaff\tn_reads" > "${OUTPUT}.tsv"

# Process each sample
for FILE in "${DEPTHS}"/*.depth; do
  SAMPLE=$(basename "$FILE" .depth)

  awk -v sample="$SAMPLE" '
    $3 >= 5 {scaffolds[$1]++}
    {total += $3}
    END {
      printf "%s\t%d\t%d\n", sample, length(scaffolds), total
    }
  ' "$FILE"

done >> "${OUTPUT}.tsv"

