#!/bin/bash
set -euo pipefail


# Identifies sites with >50/60/70/80% library call rates, creates *.beagle.gz subsets for each call rate, 
# and then runs PCAngsd on each *.beagle.gz subset
 
# Usage: $0 COUNTS BEAGLE META OUTPATH OUTNAME
#   COUNTS:    path to *.counts.gz; Note! Must have a snpcode column and sample names
#   BEAGLE:    path to *.beagle.gz; Note! Must have sample names
#   META:      path to metadata file; must have a bam_code column with same names and a library column with library membership 
#   OUTPATH:   path to save results
#   OUTNAME:   basename for result files; the output generated are:
#                ${OUTPATH}/${OUTNAME}_call_thresh0.${5..8};  snpcode lists with >50/60/70/80% call rate 
#                ${OUTPATH}/${OUTNAME}.beagle.ct${5..8}.gz;   beagle subsets for each call rate
#                ${OUTPATH}/${OUTNAME}.ct${5..8}.Pcangsd;     covariance matrices, selection statistics, etc from pcangsd
# Example usage
# $0 ./library_call_thresholds.sh \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered.counts.gz \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered.beagle.gz \
#     $SPRUCE_PROJECT/site_discovery/sequenced_samples_metadata.csv \
#     $SPRUCE_PROJECT/site_discovery/northern/domain_filtered/call_thresh
#     northern

if [ $# -lt 5 ]; then
  echo "Usage: $0 COUNTS BEAGLE META OUTPATH OUTNAME" >&2
  exit 1
fi

COUNTS="$1"
BEAGLE="$2"
META="$3"
OUTPATH="$4"
OUTNAME="$5"

# Step 1: identify sites with a > 50/60/70/80% call rate from $COUNTS
# and save snpcode lists to ${OUTPATH}/${OUTNAME}_call_thresh0.${5..8}
ml -r
ml GCC/10.2.0  CUDA/11.1.1  OpenMPI/4.0.5
ml R/4.0.4

# Create folder for results if it doesn't exist
mkdir -p "${OUTPATH}"

Rscript "$SCRIPTS/07_site_discovery/library_call_thresh.R" \
"${COUNTS}" \
"${META}" \
"${OUTPATH}/${OUTNAME}"


# Step 2: create subsets of $BEAGLE for sites with a >50/60/70/80% call rate
ml -r
ml GCC/13.3.0
ml Python/3.12.3
ml SciPy-bundle/2024.05

for i in 5 6 7 8
do
python3 $SCRIPTS/07_site_discovery/subset_genotype_matrix.py \
"${BEAGLE}" \
"${OUTPATH}/${OUTNAME}.beagle.ct${i}.gz" \
--snps "${OUTPATH}/${OUTNAME}_call_thresh0.${i}"
done

# Step 3: run Pcangsd on the beagle subsets
ml -r
ml GCC/13.2.0
ml SciPy-bundle/2023.11

for i in 5 6 7 8
do
pcangsd \
  -b "${OUTPATH}/${OUTNAME}.beagle.ct${i}.gz" \
  -o "${OUTPATH}/${OUTNAME}.ct${i}.Pcangsd" \
  --sites-save --snp-weights --pcadapt --selection --iter 500 --maf-iter 1000 -e 4
done
