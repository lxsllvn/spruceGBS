#!/usr/bin/env bash
set -euo pipefail

# Split a *reheadered* ANGSD counts matrix (*.counts.gz) by genotype class (homRef/het/homAlt)
# using a matching *reheadered* BEAGLE GL file (*.beagle.gz).
# - Argmax over (AA,AB,BB) assigns the genotype per sample/site.
# - Exact ties => treated as uncallable (depth 0 in all three outputs).
# - Sample names and snpcodes/markers must match between BEAGLE and counts.
#
# Usage:
#   split_by_genotype.sh BEAGLE COUNTS OUT_PREFIX [SITE_STATS]
#
#   BEAGLE     : path to reheadered *.beagle.gz   (marker,allele1,allele2 + triplets per sample)
#   COUNTS     : path to reheadered *.counts.gz   (snpcode + one col per sample)
#   OUT_PREFIX : output prefix; writes:
#                  OUT_PREFIX.homRef.tsv.gz
#                  OUT_PREFIX.het.tsv.gz
#                  OUT_PREFIX.homAlt.tsv.gz
#                If SITE_STATS is provided, also writes *.maf05.tsv.gz versions.
#   SITE_STATS : (optional) *_site_summary.tsv.gz containing snpcode and MAF columns
#
# Example:
#   split_by_genotype.sh \
#     $SPRUCE_PROJECT/site_discovery/southern/southern.beagle.reheader.gz \
#     $SPRUCE_PROJECT/site_discovery/southern/southern.counts.reheader.gz \
#     southern \
#     $SPRUCE_PROJECT/site_discovery/southern/southern_site_summary.tsv.gz

if (( $# < 3 )); then
  echo "Usage: $0 BEAGLE COUNTS OUT_PREFIX [SITE_STATS]" >&2
  exit 1
fi

BEAGLE="$1"
COUNTS="$2"
OUT_PREFIX="$3"
SITE_STATS="${4:-}"

# --- HPC modules / environment ---
module purge
module load GCC/13.3.0
module load Python/3.12.3
module load SciPy-bundle/2024.05

# --- Activate project venv ---
# Set this once in your environment or edit the default below.
export REPO_DIR="${REPO_DIR:-$SCRIPTS/beagle-utils}"
if [[ ! -f "$REPO_DIR/pyproject.toml" ]]; then
  echo "[ERROR] REPO_DIR does not point to your beagle-utils repo: $REPO_DIR" >&2
  exit 2
fi
# shellcheck source=/dev/null
source "$REPO_DIR/.venv/bin/activate"

# --- Input validation ---
for f in "$BEAGLE" "$COUNTS"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] File not found: $f" >&2
    exit 3
  fi
done

if [[ "$BEAGLE" != *.beagle.gz ]]; then
  echo "[ERROR] BEAGLE must end with .beagle.gz: $BEAGLE" >&2
  exit 4
fi
if [[ "$COUNTS" != *.counts.gz ]]; then
  echo "[ERROR] COUNTS must end with .counts.gz: $COUNTS" >&2
  exit 5
fi

if [[ -n "$SITE_STATS" && ! -f "$SITE_STATS" ]]; then
  echo "[ERROR] SITE_STATS not found: $SITE_STATS" >&2
  exit 6
fi

# --- Run geno-utils split-genotype ---
echo "[INFO] Splitting counts by genotype using BEAGLE GLs…"
if [[ -n "$SITE_STATS" ]]; then
  # Default MAF threshold is 0.05 in your Python; override with --maf-threshold if needed.
  geno-utils split-genotype "$BEAGLE" "$COUNTS" "$OUT_PREFIX" \
    --site-stats "$SITE_STATS" \
    --maf-threshold 0.05
else
  geno-utils split-genotype "$BEAGLE" "$COUNTS" "$OUT_PREFIX"
fi

echo "[OK] Wrote:"
echo "  ${OUT_PREFIX}.homRef.tsv.gz"
echo "  ${OUT_PREFIX}.het.tsv.gz"
echo "  ${OUT_PREFIX}.homAlt.tsv.gz"
[[ -n "${SITE_STATS:-}" ]] && {
  echo "  ${OUT_PREFIX}.homRef.maf05.tsv.gz"
  echo "  ${OUT_PREFIX}.het.maf05.tsv.gz"
  echo "  ${OUT_PREFIX}.homAlt.maf05.tsv.gz"
}
