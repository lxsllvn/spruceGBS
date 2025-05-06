#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

# Returns a list of unique scaffolds with >= 5 mapped reads in any sample.

# Require path to read depth folder. Optionally take path for scratch directory and outfile name

# Require 1–3 arguments
if [ "$#" -lt 1 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 /path/to/read_depths [scratch_dir] [output_file]" >&2
    exit 1
fi

DEPTHS="$1"
SCRATCH="${2:-${SPRUCE_PROJECT}/scaff_cov_tmp_$$}"
OUTFILE="${3:-${SPRUCE_PROJECT}/ref/scaffolds_with_coverage.txt}"

echo "Using depth directory: $DEPTHS"
echo "Using scratch directory: $SCRATCH"
echo "Writing output to: $OUTFILE"

# Create scratch dir if needed and register trap for cleanup
mkdir -p "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT

echo "Starting scaffold coverage summary at $(date)"

# Fail if there are no .depth files
shopt -s nullglob
files=("$DEPTHS"/*.depth)
if [ "${#files[@]}" -eq 0 ]; then
    echo "No .depth files found in $DEPTHS" >&2
    exit 1
fi

# Extract scaffold names with at least 5 reads and ensure uniqueness
awk '$3 >= 5 {print $1}' "${DEPTHS}"/*.depth | \
  sort -u --temporary-directory="$SCRATCH" \
  > "${SPRUCE_PROJECT}/ref/scaffolds_with_coverage.txt"

echo "Finished at $(date)"