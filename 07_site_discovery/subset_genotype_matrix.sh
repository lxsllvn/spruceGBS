#!/usr/bin/env bash
set -euo pipefail

# Subset reheadered ANGSD counts (*.counts.gz) or BEAGLE GL (*.beagle.gz) files
# by SNPs and/or samples.
# - "Reheadered" means INPUT has informative snpcode/marker and sample names.
# - Subsetting requires providing --snps and/or --samples lists.
#
# Usage:
#   subset_by_list.sh INPUT OUTPUT [--snps SNPS.txt] [--samples SAMPLES.txt]
#
#   INPUT   : reheadered *.counts.gz or *.beagle.gz
#   OUTPUT  : output .gz file
#   --snps  : (optional) list of snpcodes/markers to keep (one per line)
#   --samples : (optional) list of sample names to keep (one per line)
#
# Example:
#   subset_by_list.sh \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered.beagle.gz \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/FennN.beagle.gz \
#     --samples $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/FennN_sample_list.txt \
#     --snps    $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered_snpcodes.txt

if (( $# < 2 )); then
  echo "Usage: $0 INPUT OUTPUT [--snps SNPS.txt] [--samples SAMPLES.txt]" >&2
  exit 1
fi

INPUT="$1"
OUTPUT="$2"
shift 2  # remaining args forwarded directly to geno-utils

# --- HPC modules / environment ---
module purge
module load GCC/13.3.0
module load Python/3.12.3
module load SciPy-bundle/2024.05

# --- Activate project venv ---
export REPO_DIR="${REPO_DIR:-$SCRIPTS/beagle-utils}"
if [[ ! -f "$REPO_DIR/pyproject.toml" ]]; then
  echo "[ERROR] REPO_DIR does not point to beagle-utils repo: $REPO_DIR" >&2
  exit 2
fi
# shellcheck source=/dev/null
source "$REPO_DIR/.venv/bin/activate"

# --- Input validation ---
if [[ ! -f "$INPUT" ]]; then
  echo "[ERROR] INPUT file not found: $INPUT" >&2
  exit 3
fi

# --- Run geno-utils subset ---
echo "[INFO] Subsetting $INPUT -> $OUTPUT ..."
geno-utils subset "$INPUT" "$OUTPUT" "$@"

echo "[OK] Wrote: $OUTPUT"
