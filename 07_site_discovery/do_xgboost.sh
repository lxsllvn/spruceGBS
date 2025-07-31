#!/bin/bash
set -euo pipefail

# Calculate genotype-stratified read depth summaries (the input DP matrices must be created
# in advance with read_counts_by_genotype.py), prepares feature data for xgboost, 
# designates test and training datasets, conducts Bayesian hyperparameter optimization,
# trains the xgboost model, and finally makes predictions for the test and training data.

# Usage: $0 SNPSTATS_PATH ALT_PATH ref_path het_path out_path out_name
#  SNPSTATS_PATH:  path to *_site_summary_maf05.tsv.gz
#  ALT_PATH:  path to *.homAlt.tsv.gz; read depths for homozygous alternate genotypes
#  REF_PATH:  path to *.homRef.tsv.gz; read depths for homozygous reference genotypes
#  HET_PATH:  path *.het.tsv.gz; read depths for heterozygous genotypes
#  OUT_PATH:  path to output folder; will create if it doesn't exist
#  OUT_NAME:  basename for result files


if [ $# -lt 6 ]; then
  echo "Usage: $0 SNPSTATS_PATH ALT_PATH REF_PATH HET_PATH OUT_PATH OUT_NAME" >&2
  exit 1
fi

SNPSTATS_PATH="$1"
ALT_PATH="$2"
REF_PATH="$3"
HET_PATH="$4"
OUT_PATH="$5"
OUT_NAME="$6"

# Path to R script with the `depth_by_genotype` function
SOURCE_PATH="$SCRIPTS/07_site_discovery/discovery_and_filtering.R"

# Create folder for results if it doesn't exist
mkdir -p "${OUT_PATH}"

ml GCC/10.2.0  CUDA/11.1.1  OpenMPI/4.0.5 R/4.0.4

Rscript "$SCRIPTS/07_site_discovery/do_xgboost.R" \
"${SNPSTATS_PATH}" \
"${ALT_PATH}"  \
"${REF_PATH}" \
"${HET_PATH}" \
"${OUT_PATH}" \
"${OUT_NAME}" \
"${SOURCE_PATH}"
