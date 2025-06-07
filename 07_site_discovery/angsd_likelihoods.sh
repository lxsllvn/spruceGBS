#!/bin/bash
set -euo pipefail

# Load modules
ml GCC/10.2.0
ml angsd/0.935

# ANGSD script for genotype likelihoods and associated functions
# Usage: $0 REF REGION SITES BAMLIST OUTNAME [OUTDIR]
#   REF:      path to domain reference genome
#   REGION:   path to ANGSD region file
#   SITES:    path to ANGSD sites file
#   BAMLIST:  path to bamlist file
#   OUTNAME:  base name for output files
#   OUTDIR:   directory to save outputs (optional); current working directory if unspecified

if [ $# -lt 5 ]; then
  echo "Usage: $0 REF REGION SITES BAMLIST OUTNAME [OUTDIR]" >&2
  exit 1
fi

REF="$1"
REGION="$2"
SITES="$3"
BAM="$4"
OUTNAME="$5"
OUTDIR="${6:-.}"

mkdir -p "$OUTDIR"

OUTPATH="${OUTDIR%/}/$OUTNAME"

echo "Running ANGSD with:"
echo "REF=$REF"
echo "REGION=$REGION"
echo "SITES=$SITES"
echo "BAMLIST=$BAM"
echo "OUTNAME=$OUTNAME"
echo "OUTDIR=$OUTDIR"
echo "Output path prefix: $OUTPATH"

P="2" # number of threads

# Filters
MAP_FILTERS="-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -C 100"
QUAL_FILTERS="-minMapQ 50 -minQ 20 -baq 0"
DEPTH_FILTERS="-setMinDepthInd 3"
SITE_FILTERS="-skipTriallelic 1"

G="1"
TODO="-doCounts 1 -dumpCounts 2 -doMajorMinor 4 -doMaf 1 -doHWE 1 -dosnpstat 1 -doIBS 1 -doCov 1 -makeMatrix 1 -doGlf 2"

angsd -P $P -b "$BAM" -ref "$REF" -anc "$REF" -rf "$REGION" -sites "$SITES" -out "$OUTPATH" -GL $G $TODO $MAP_FILTERS $QUAL_FILTERS $DEPTH_FILTERS
