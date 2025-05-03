#!/usr/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -J refprep
#SBATCH --constraint=skylake
#SBATCH -t 0-06:00:00

set -euo pipefail
IFS=$'\n\t'

# ------------------------
# SCRIPT ARGUMENTS
# ------------------------
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 DOMAIN FOLDER OUTDIR" >&2
    exit 1
fi

DOMAIN="$1"   # e.g. “siberia”
FOLDER="$2"   # directory containing *.pos.gz lists
OUTDIR="$3"   # where to write *.bed, *.sites, *.regions, *.fa

mkdir -p "$OUTDIR"

# ------------------------
# 1) Convert each .pos.gz to a BED (group contiguous bases)
# ------------------------
for idx in {a..j}; do
    input="${FOLDER}/${DOMAIN}_target_scaff_pt_a${idx}.pos.gz"
    output="${FOLDER}/${DOMAIN}_pt_a${idx}.bed"

    awk '
        BEGIN { OFS = "\t" }
        NR == 1 {
            chr   = $1
            start = $2
            prev  = $2
            next
        }
        {
            if ($1 == chr && $2 == prev + 1) {
                prev = $2
            } else {
                print chr, start-1, prev
                chr   = $1
                start = $2
                prev  = $2
            }
        }
        END {
            print chr, start-1, prev
        }
    ' < <(zcat "$input") \
      | sort -k1,1 -k2,2n \
      > "$output"
done

# ------------------------
# 2) Merge all BEDs into one sorted experiment_ref.bed
# ------------------------
cat "${FOLDER}/${DOMAIN}"_pt_*.bed \
  | sort -k1,1 -k2,2n \
  > "${OUTDIR}/${DOMAIN}_experiment_ref.bed"

# ------------------------
# 3) Generate ANGSD sites & region files
# ------------------------
# 3a) 1-indexed sites
awk 'BEGIN{OFS="\t"} {print $1, $2+1, $3}' \
    "${OUTDIR}/${DOMAIN}_experiment_ref.bed" \
  > "${OUTDIR}/${DOMAIN}_experiment_ref_sites"

# 3b) Unique scaffold list for region
cut -f1 "${OUTDIR}/${DOMAIN}_experiment_ref_sites" \
  | sort -u \
  > "${OUTDIR}/${DOMAIN}_experiment_ref_regions"

# ------------------------
# 4) Build & index reduced FASTA
# ------------------------
ml GCC/13.2.0 SAMtools/1.19.2

xargs samtools faidx "$SPRUCE_PROJECT/ref/picea_newref.fa" \
    < "${OUTDIR}/${DOMAIN}_experiment_ref_regions" \
  > "${OUTDIR}/${DOMAIN}_experiment_ref.fa"

# ensure .fa is at least 1 s older than the .fai
sleep 1

samtools faidx "${OUTDIR}/${DOMAIN}_experiment_ref.fa"

# ------------------------
# 5) Index sites for ANGSD
# ------------------------
ml -r           # reset modules
ml GCC/10.2.0 angsd/0.935

angsd sites index "${OUTDIR}/${DOMAIN}_experiment_ref_sites"
