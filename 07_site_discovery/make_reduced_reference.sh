#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

# -----------------------------------------------------
# Reduced Reference Creation for ANGSD
# -----------------------------------------------------
# Usage:
#   ./make_reduced_reference.sh OUTDIR REF_BED REF_FASTA
# Arguments:
#   OUTDIR     - Output directory for all results
#   REF_BED    - BED file with target regions
#   REF_FASTA  - Reference genome FASTA file (indexed for samtools)

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 OUTDIR REF_BED REF_FASTA"
    exit 1
fi

OUTDIR="$1"
REF_BED="$2"
REF_FASTA="$3"

mkdir -p "$OUTDIR"

# Load required modules for samtools (edit/remove as needed for your system)
ml GCC/13.2.0
ml SAMtools/1.19.2

# 1) Get unique scaffold names from target region bed
cut -f1 "$REF_BED" | sort -u > "$OUTDIR/target_scaffs.txt"

# 2) Split into 20 subsets
split -n l/20 --numeric-suffixes=1 --suffix-length=2 --additional-suffix .txt \
    "$OUTDIR/target_scaffs.txt" \
    "$OUTDIR/target_scaff_pt_"

# 3) Extract FASTA for each subset and index
for idx in {01..20}; do
    list="$OUTDIR/target_scaff_pt_${idx}.txt"
    fa="$OUTDIR/target_scaff_pt_${idx}.fa"

    xargs samtools faidx "$REF_FASTA" < "$list" > "$fa"

    # ensure index files are at least 1 s older than sites
    sleep 1
    
    samtools faidx "$fa"
done

# 4) Load ANGSD modules
ml -r
ml GCC/10.2.0 angsd/0.935

# 5) Create region & site files for ANGSD, then index
for idx in {01..20}; do
    regions="$OUTDIR/target_scaff_pt_${idx}_regions"
    bedin="$REF_BED"
    bedout="$OUTDIR/target_scaff_pt_${idx}.bed"
    sites="$OUTDIR/target_scaff_pt_${idx}_sites"

    sort -k1,1 "$OUTDIR/target_scaff_pt_${idx}.txt" > "$regions"

    grep -Fwf "$regions" "$bedin" > "$bedout"

    # make 1-indexed sites for ANGSD
    awk 'BEGIN { OFS = "\t" } { print $1, $2+1, $3 }' "$bedout" > "$sites"

    # ensure index files are at least 1 s older than sites
    sleep 1

    angsd sites index "$sites"
done
