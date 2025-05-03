#!/bin/bash

# Load modules
ml GCC/10.2.0
ml angsd/0.935

# ------------------------
# SCRIPT ARGUMENTS
# ------------------------
if [ $# -lt 4 ]; then
  echo "Usage: $0 INDEX DOMAIN BAMLIST OUTDIR" >&2
  exit 1
fi

INDEX="$1" 
DOMAIN="$2"
BAMLIST="$3"
OUTDIR="$4"

mkdir -p "$OUTDIR"

# ------------------------
# SET PATHS
# ------------------------

REF="$SPRUCE_PROJECT/parameter_testing/exp_ref/experiment_ref_pt_a${INDEX}.fa"
REGION="$SPRUCE_PROJECT/parameter_testing/exp_ref/target_scaff_pt_a${INDEX}_regions"
SITES="$SPRUCE_PROJECT/parameter_testing/exp_ref/target_scaff_pt_a${INDEX}_sites"
OUTPUT="${OUTDIR}/${DOMAIN}_target_scaff_pt_a${INDEX}"

# ------------------------
# SET FILTERS
# ------------------------

n=$(wc -l < "$BAMLIST") #Count number of individuals in bam list
MININD=$(((${n} * 1/2))) #Minimum number of individuals; 1/2 for 50% etc.
DEFAULT_FILTERS="-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -setMinDepthInd 3 -minInd $MININD"

# Calculate domain-wide read count matrix 
angsd -bam "$BAMLIST" \
      -ref "$REF" \
      -anc "$REF" \
      -rf "$REGION" \
      -sites "$SITES" \
      -GL 1 \
      -doCounts 1 \
      -dumpCounts 2 \
      -baq 0 \
      -C 0 \
      -minQ 20 \
      -minMapQ 20 \
      $DEFAULT_FILTERS \
      -out "$OUTPUT"
