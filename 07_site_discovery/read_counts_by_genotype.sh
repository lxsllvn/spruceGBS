#!/bin/bash
set -euo pipefail

# Splits a reheadered ANGSD DP matrix (*.counts.gz) by genotype class (homRef, homAlt, and het) using 
# the beagle genotype likelihood file (*.beagle.gz). Counts are assigned to the genotype with the highest
# likelihood. If all three genotypes are equally likely for a sample, all three are assigned 0 depth.
# Note! "reheadered" means the input has informative column and row names. This is not the default
# output from ANGSD; run reheader_genotype_matrix.sh first to add SNP codes and sample names.
# Sample names and snpcodes in *.beagle.gz and *.counts.gz must match.
#
# Usage: $0 BEAGLE COUNTS OUT SITE_STATS
#   INPUT:       *.beagle.gz
#   COUNTS:      *.counts.gz
#   OUTPUT:      output file prefix, e.g. $OUPUT.homoRef.tsv.gz
#   SITE_STATS:  (optional) path to *_site_summary.tsv.gz; if provided, the script will write 
#                additional outputs for MAF > 0.05 sites
# Example usage
# $0 ./read_counts_by_genotype.sh \
#     $SPRUCE_PROJECT/site_discovery/southern/southern.beagle.reheader.gz \
#     $SPRUCE_PROJECT/site_discovery/southern/southern.counts.reheader.gz \
#     southern \
#     $SPRUCE_PROJECT/site_discovery/southern/southern_site_summary_maf05.tsv.gz

if [ $# -lt 3 ]; then
  echo "Usage: $0 BEAGLE COUNTS OUT [SITE_STATS]" >&2
  exit 1
fi

ml GCC/13.3.0
ml Python/3.12.3
ml SciPy-bundle/2024.05

BEAGLE="$1"
COUNTS="$2"
OUT="$3" 
SITE_STATS="${4:-}"

if [ -n "$SITE_STATS" ]; then
  python3 $SCRIPTS/07_site_discovery/read_counts_by_genotype.py "${BEAGLE}" "${COUNTS}" "${OUT}" "${SITE_STATS}"
else
  python3 $SCRIPTS/07_site_discovery/read_counts_by_genotype.py "${BEAGLE}" "${COUNTS}" "${OUT}" 
fi
