#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
ANGSD reference preparation utility

Usage:
  $0 --splits DOMAIN INPUT REF OUTDIR
      DOMAIN:   sample domain (e.g., southern)
      INPUT:    directory containing *.pos.gz lists (e.g., southern_pt_01.pos.gz, etc)
      REF:      path to reference genome FASTA (indexed for samtools)
      OUTDIR:   directory to save merged .bed, .sites, .regions, .fa files

  $0 --sites SITES OUTDIR
      SITES:    1-indexed ANGSD sites file (tab-separated: scaffold, position, end)
      OUTDIR:   directory to save region and index files
EOF
  exit 1
}

if [ "$#" -lt 1 ]; then
  usage
fi

mode="$1"
shift

case "$mode" in
  --sites)
    if [ "$#" -ne 2 ]; then usage; fi
    SITES="$1"
    OUTDIR="$2"
    mkdir -p "$OUTDIR"

    # Create ANGSD region file from sites
    cut -f1 "$SITES" | sort -u > "${OUTDIR}/ref_regions"

    # ANGSD sites index
    ml -r
    ml GCC/10.2.0 angsd/0.935
    angsd sites index "$SITES"

    echo "Done (ANGSD sites/region files, --sites mode)"
    ;;

  --splits)
    if [ "$#" -ne 4 ]; then usage; fi
    DOMAIN="$1"
    INPUT="$2"
    REF="$3"
    OUTDIR="$4"
    mkdir -p "$OUTDIR"

    # Step 1: Convert each .pos.gz to BED
    for file in "${INPUT}/${DOMAIN}_pt_"*.pos.gz; do
      idx=$(basename "$file" | sed -E "s/^${DOMAIN}_pt_(.*)\.pos\.gz$/\1/")
      output="${INPUT}/${DOMAIN}_pt_${idx}.bed"
      zcat "$file" | tail -n +2 | awk '
        BEGIN { OFS = "\t" }
        NR == 1 { chr=$1; start=$2; prev=$2; next }
        {
          if ($1 == chr && $2 == prev + 1) { prev=$2 }
          else { print chr, start-1, prev; chr=$1; start=$2; prev=$2 }
      }
        END { print chr, start-1, prev }
        ' | sort -k1,1 -k2,2n > "$output"
    done

    # Step 2: Merge all BEDs into one sorted master BED
    cat "${INPUT}/${DOMAIN}"_pt_*.bed | sort -k1,1 -k2,2n > "${OUTDIR}/${DOMAIN}_ref.bed"

    # Step 3a: Create 1-indexed sites file
    awk 'BEGIN{OFS="\t"} {print $1, $2+1, $3}' "${OUTDIR}/${DOMAIN}_ref.bed" > "${OUTDIR}/${DOMAIN}_ref_sites"

    # Step 3b: Unique scaffold list for region
    cut -f1 "${OUTDIR}/${DOMAIN}_ref_sites" | sort -u > "${OUTDIR}/${DOMAIN}_ref_regions"

    # Step 4: Build & index reduced FASTA
    ml GCC/13.2.0 SAMtools/1.19.2
    xargs samtools faidx "$REF" < "${OUTDIR}/${DOMAIN}_ref_regions" > "${OUTDIR}/${DOMAIN}_ref.fa"
    sleep 1
    samtools faidx "${OUTDIR}/${DOMAIN}_ref.fa"

    # Step 5: ANGSD sites index
    ml -r
    ml GCC/10.2.0 angsd/0.935
    angsd sites index "${OUTDIR}/${DOMAIN}_ref_sites"

    echo "Done (full workflow, --splits mode)"
    ;;

  *)
    usage
    ;;
esac
