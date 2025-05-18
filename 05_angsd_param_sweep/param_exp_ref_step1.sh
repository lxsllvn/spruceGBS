#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

# -----------------------------------------------------
# Reduced Reference Creation for ANGSD Parameter Experiments
# -----------------------------------------------------
# Usage:
#   ./make_reduced_reference.sh OUTDIR REF_BED CONTAM_BED REF_FASTA
# Arguments:
#   OUTDIR     - Output directory for all results
#   REF_BED    - BED file with non-contaminant target regions
#   CONTAM_BED - BED file listing contaminant scaffolds
#   REF_FASTA  - Reference genome FASTA file (indexed for samtools)

# Check arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 OUTDIR REF_BED CONTAM_BED REF_FASTA"
    exit 1
fi

OUTDIR="$1"
REF_BED="$2"
CONTAM_BED="$3"
REF_FASTA="$4"

mkdir -p "$OUTDIR"

# Load required modules for samtools (edit/remove as needed for your system)
ml GCC/13.2.0
ml SAMtools/1.19.2

# 1) Select 100 Mbp from longest non-contaminant scaffolds
awk '{
    len = $3 - $2
    sum[$1] += len
}
END {
    for (scaff in sum) print scaff, sum[scaff]
}' "$REF_BED" \
    | sort -k2,2nr \
    | grep -Fvwf <(cut -f1 "$CONTAM_BED") \
    | grep -v 'chloroplast' \
    | awk 'BEGIN { OFS = "\t"; total = 0 }
        {
            total += $2
            print $1
            if (total > 100000000) exit
        }' \
    > "$OUTDIR/target_scaffs.txt"

# 2) Split into 10 subsets
split -n l/10 --additional-suffix .txt \
    "$OUTDIR/target_scaffs.txt" \
    "$OUTDIR/target_scaff_pt_"

# 3) Extract FASTA for each subset and index
for idx in {a..j}; do
    list="$OUTDIR/target_scaff_pt_a${idx}.txt"
    fa="$OUTDIR/experiment_ref_pt_a${idx}.fa"

    xargs samtools faidx "$REF_FASTA" < "$list" > "$fa"
    samtools faidx "$fa"
done

# 4) Load ANGSD modules
ml -r
ml GCC/10.2.0 angsd/0.935

# 5) Create region & site files for ANGSD, then index
for idx in {a..j}; do
    regions="$OUTDIR/target_scaff_pt_a${idx}_regions"
    bedin="$REF_BED"
    bedout="$OUTDIR/target_scaff_pt_a${idx}.bed"
    sites="$OUTDIR/target_scaff_pt_a${idx}_sites"

    sort -k1,1 "$OUTDIR/target_scaff_pt_a${idx}.txt" > "$regions"

    grep -Fwf "$regions" "$bedin" > "$bedout"

    # make 1-indexed sites for ANGSD
    awk 'BEGIN { OFS = "\t" } { print $1, $2+1, $3 }' "$bedout" > "$sites"

    # ensure index files are at least 1 s older than sites
    sleep 1

    angsd sites index "$sites"
done