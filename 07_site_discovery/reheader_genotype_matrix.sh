#!/bin/bash
set -euo pipefail

# Adds SNP codes and sample names to ANGSD DP matrices (*.counts.gz) and beagle genotype 
# likelihood files (*.beagle.gz). By default, the DP matrix lacks a marker name column
# and meaningful sample names and the the beagle file has a marker name column but lacks 
# meaningful sample names.
#
# Usage: $0 INPUT SAMPLES ARGS
#   INPUT:     *.beagle.gz or *.counts.gz to subset 
#   OUTPUT:     output name
#   SAMPLES:   list of sample names to add
#   POS:       if the input is *.counts.gz, the corresponding *.pos.gz file must be provided
# Example usage for beagle input:
# $0 ./reheader_genotype_matrix.sh \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered.beagle.gz \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_reheadered.beagle.gz \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_sample_list.txt 
# Example usage for *.counts.gz input:
# $0 ./reheader_genotype_matrix.sh \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered.counts.gz \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_reheadered.counts.gz \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_sample_list.txt \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered.pos.gz
        
        
if [ $# -lt 3 ]; then
  echo "Usage: $0 INPUT OUTPUT SAMPLES [POS]" >&2
  exit 1
fi


ml GCC/13.3.0
ml Python/3.12.3
ml SciPy-bundle/2024.05

INPUT="$1"   
OUTPUT="$2" 
SAMPLES="$3"
POS="${4:-}"

if [ -n "$POS" ]; then
  python3 $SCRIPTS/07_site_discovery/reheader_genotype_matrix.py "${INPUT}" "${OUTPUT}" "${SAMPLES}" "${POS}"
else
  python3 $SCRIPTS/07_site_discovery/reheader_genotype_matrix.py "${INPUT}" "${OUTPUT}" "${SAMPLES}" 
fi
