#!/usr/bin/env bash
set -euo pipefail

# Reheader ANGSD counts (*.counts.gz) or BEAGLE GL (*.beagle.gz) using geno-utils.
#
# Usage:
#   reheader.sh INPUT OUTPUT SAMPLES [POS]
#
#   INPUT    : path to *.beagle.gz or *.counts.gz
#   OUTPUT   : output .gz path
#   SAMPLES  : text file with sample names (one per line)
#   POS      : (required for *.counts.gz) matching *.pos.gz to build snpcode
#
# Examples:
#   # BEAGLE
#   reheader.sh \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered.beagle.gz \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_reheadered.beagle.gz \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_sample_list.txt
#
#   # COUNTS (+ POS)
#   reheader.sh \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered.counts.gz \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_reheadered.counts.gz \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_sample_list.txt \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered.pos.gz

if (( $# < 3 )); then
  echo "Usage: $0 INPUT OUTPUT SAMPLES [POS]" >&2
  exit 1
fi

INPUT="$1"
OUTPUT="$2"
SAMPLES="$3"
POS="${4:-}"

# --- HPC modules / environment ---
module purge
module load GCC/13.3.0
module load Python/3.12.3
module load SciPy-bundle/2024.05

# --- Activate project venv ---
export REPO_DIR=$SCRIPTS/beagle-utils
source $REPO_DIR/.venv/bin/activate

if [[ ! -f "$REPO_DIR/pyproject.toml" ]]; then
  echo "[ERROR] REPO_DIR does not point to your beagle-utils repo: $REPO_DIR" >&2
  exit 2
fi

# --- Basic input validation ---
if [[ ! -f "$INPUT" ]]; then
  echo "[ERROR] INPUT not found: $INPUT" >&2
  exit 3
fi
if [[ ! -f "$SAMPLES" ]]; then
  echo "[ERROR] SAMPLES file not found: $SAMPLES" >&2
  exit 4
fi

is_counts=false
is_beagle=false
if [[ "$INPUT" == *.counts.gz ]]; then
  is_counts=true
elif [[ "$INPUT" == *.beagle.gz ]]; then
  is_beagle=true
else
  echo "[ERROR] INPUT must end with .counts.gz or .beagle.gz: $INPUT" >&2
  exit 5
fi

if $is_counts && [[ -z "$POS" ]]; then
  echo "[ERROR] POS (.pos.gz) is required when INPUT is a .counts.gz file." >&2
  exit 6
fi
if $is_counts && [[ ! -f "$POS" ]]; then
  echo "[ERROR] POS file not found: $POS" >&2
  exit 7
fi

# --- Run geno-utils reheader ---
if $is_counts; then
  echo "[INFO] Reheadering COUNTS with snpcode from POS... and column names from SAMPLES"
  geno-utils reheader "$INPUT" "$SAMPLES" "$OUTPUT" --pos "$POS"
else
  echo "[INFO] Reheadering BEAGLE with column names from SAMPLES..."
  geno-utils reheader "$INPUT" "$SAMPLES" "$OUTPUT"
fi

echo "[OK] Wrote: $OUTPUT"
