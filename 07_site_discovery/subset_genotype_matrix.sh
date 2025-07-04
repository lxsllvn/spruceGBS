#!/bin/bash
set -euo pipefail

# Subsets reheadered ANGSD DP matrices (*.counts.gz) and beagle genotype likelihood files (*.beagle.gz)
# to a specified list of SNPs and/or samples. Note! "reheadered" means the input has informative
# column and row names. This is not the default output from ANGSD; run reheader_genotype_matrix.sh
# first to add SNP codes and sample names.
#
# Usage: $0 INPUT OUTPUT ARGS
#   INPUT:     *.beagle.gz or *.counts.gz to subset 
#   OUTPUT:    output name
#   ARGS:      --snps snps_to_keep.txt and/or --samples samples_to_keep.txt
# Example usage
# $0 ./subset_genotype_matrix.sh \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered.beagle.gz \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/FennN.beagle.gz \
#     --samples $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/FennN_sample_list.txt \ 
#     --snps $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered_snpcodes.txt

if [ $# -lt 3 ]; then
  echo "Usage: $0 INPUT OUTPUT ARGS" >&2
  exit 1
fi

ml GCC/13.3.0
ml Python/3.12.3
ml SciPy-bundle/2024.05

INPUT="$1"   
OUTPUT="$2" 
shift 2  # $@ all the ARGS

python3 $SCRIPTS/07_site_discovery/subset_genotype_matrix.py "${INPUT}" "${OUTPUT}" "$@"
