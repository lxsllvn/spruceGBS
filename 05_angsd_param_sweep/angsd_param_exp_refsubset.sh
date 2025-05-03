#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

# Creates a reduced reference for the parameter experiments and divides it into 10 subsets for ANGSD

# ------------------------
# CONFIGURATION
# ------------------------
OUTDIR="${SPRUCE_PROJECT}/parameter_testing/exp_ref"
mkdir -p "$OUTDIR"

# Load required modules for samtools
ml GCC/13.2.0
ml SAMtools/1.19.2

# ------------------------
# 1) Select 100 Mbp from longest non-contaminant scaffolds
# ------------------------

awk '{
    len = $3 - $2
    sum[$1] += len
}
END {
    for (scaff in sum) print scaff, sum[scaff]
}' "${SPRUCE_PROJECT}/ref/picea_newref_target_regions.bed" \
    | sort -k2,2nr \
    | grep -Fvwf <(cut -f1 "${SPRUCE_PROJECT}/ref/ref_putative_contaminants.bed") \
    | grep -v 'chloroplast' \
    | awk 'BEGIN { OFS = "\t"; total = 0 }
        {
            total += $2
            print $1
            if (total > 100000000) exit
        }' \
    > "${OUTDIR}/target_scaffs.txt"

# ------------------------
# 2) Split into 10 subsets
# ------------------------
split -n l/10 --additional-suffix .txt \
    "${OUTDIR}/target_scaffs.txt" \
    "${OUTDIR}/target_scaff_pt_"
    
# ------------------------
# 3) Extract FASTA for each subset and index
# ------------------------
for idx in {a..j}; do
    list="${OUTDIR}/target_scaff_pt_a${idx}.txt"
    fa  ="${OUTDIR}/experiment_ref_pt_a${idx}.fa"

    xargs samtools faidx \
        "${SPRUCE_PROJECT}/ref/picea_newref.fa" \
        < "$list" \
        > "$fa"

    samtools faidx "$fa"
done

# ------------------------
# 4) Load ANGSD modules
# ------------------------
ml -r 
ml GCC/10.2.0 angsd/0.935

# ------------------------
# 5) Create region & site files for ANGSD, then index
# ------------------------
for idx in {a..j}; do
    regions="${OUTDIR}/target_scaff_pt_a${idx}_regions"
    bedin  ="${SPRUCE_PROJECT}/ref/picea_newref_target_regions.bed"
    bedout ="${OUTDIR}/target_scaff_pt_a${idx}.bed"
    sites  ="${OUTDIR}/target_scaff_pt_a${idx}_sites"

    sort -k1,1 "$OUTDIR/target_scaff_pt_a${idx}.txt" \
        > "$regions"

    grep -Fwf "$regions" \
        "$bedin" \
        > "$bedout"

    # make 1-indexed sites for ANGSD
    awk 'BEGIN { OFS = "\t" } { print $1, $2+1, $3 }' \
        "$bedout" \
        > "$sites"

    # ensure index files are at least 1 s older than sites
    sleep 1

    angsd sites index "$sites"
done
