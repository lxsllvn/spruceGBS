#!/bin/bash
set -euo pipefail

# Load modules
ml GCC/10.2.0
ml angsd/0.935

# ------------------------
# Site Discovery per Domain with ANGSD
# ------------------------
# Usage: $0 INDEX DOMAIN BAMLIST REFDIR OUTDIR
#   INDEX:    subset index (e.g., 01)
#   DOMAIN:   sample domain (e.g., southern)
#   BAMLIST:  path to bamlist file
#   REFDIR:   directory where reference subsets are located
#   OUTDIR:   directory to save outputs


if [ $# -lt 5 ]; then
  echo "Usage: $0 INDEX DOMAIN BAMLIST REFDIR OUTDIR" >&2
  exit 1
fi

INDEX="$1"
DOMAIN="$2"
BAMLIST="$3"
REFDIR="$4"
OUTDIR="$5"

mkdir -p "$OUTDIR"

# ------------------------
# SET PATHS
# ------------------------
REF="$REFDIR/target_scaff_pt_${INDEX}.fa"
REGION="$REFDIR/target_scaff_pt_${INDEX}_regions"
SITES="$REFDIR/target_scaff_pt_${INDEX}_sites"
OUTNAME="${OUTDIR}/${DOMAIN}_pt_${INDEX}" #save directory/filename 

# ------------------------
# SET FILTERS
# ------------------------
n=$(wc -l < "$BAMLIST")
MININD=$(awk -v n="$n" 'BEGIN { printf("%d", n * 0.4 + 0.5) }') #Minimum number of individuals; n * 0.4 for 40% etc.
DEFAULT_FILTERS="-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -setMinDepthInd 3 -minInd $MININD"

echo "Running ANGSD for $DOMAIN, index $INDEX with minInd $MININD"
echo "REF=$REF"
echo "REGION=$REGION"
echo "SITES=$SITES"
echo "BAMLIST=$BAMLIST"
echo "OUTNAME=$OUTNAME"

angsd -bam "$BAMLIST" \
      -ref "$REF" \
      -anc "$REF" \
      -rf "$REGION" \
      -sites "$SITES" \
      -GL 1 \
      -doCounts 1 \
      -dumpCounts 2 \
      -baq 0 \
      -C 100 \
      -minQ 20 \
      -minMapQ 50 \
      $DEFAULT_FILTERS \
      -out "$OUTNAME"
