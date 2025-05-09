#!/bin/bash

# Exit on error, undefined var, or pipeline error
set -euo pipefail

# ------------------------
# CHECK REQUIRED ARGUMENTS
# ------------------------
if [ "$#" -ne 4 ]; then
  cat << EOF
Usage: $0 <REF> <DOMAIN> <POP_LIST> <RESULTS_DIR>
  REF         : Path to reference genome (FASTA)
  DOMAIN      : Domain name (e.g., southern, siberia)
  POP_LIST    : File listing populations (one per line)
  RESULTS_DIR : Top-level directory of parameter sweep results
EOF
  exit 1
fi

# ------------------------
# THIS SCRIPT...
# ------------------------
# 1) Runs ANGSD + realSFS to compute individual heterozygosity (*.sfs)
# 2) Collects all .sfs files into a summary TSV

# ------------------------
# USER-DEFINED PARAMETERS
# ------------------------
BAQ_LIST=(0 1 2)
C_LIST=(0 50)
MINQ_LIST=(20 30)
MINMAPQ_LIST=(20 30 40)
DEFAULT_FILTERS="-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -skipTriallelic 1 -setMinDepthInd 3"

# Load modules
ml GCC/10.2.0
ml angsd/0.935

# ------------------------
# SCRIPT ARGUMENTS
# ------------------------
REF="$1"
DOMAIN="$2"
POP_LIST="$3"
RESULTS_DIR="$(realpath "$4")"

# Directory containing per-population BAM lists:
POP_BAM_DIR="${DOMAIN}_populations"

# Validate POP_LIST and BAM list presence
while IFS= read -r POP; do
  bamlist="${POP_BAM_DIR}/${POP}.txt"
  if [[ ! -f "$bamlist" ]]; then
    echo "ERROR: BAM list for population '$POP' not found at $bamlist" >&2
    exit 1
  fi
done < "$POP_LIST"

# ------------------------
# MAIN PARAMETER SWEEP
# ------------------------
for baq in "${BAQ_LIST[@]}"; do
  for C in "${C_LIST[@]}"; do
    for minQ in "${MINQ_LIST[@]}"; do
      for minMapQ in "${MINMAPQ_LIST[@]}"; do

        PARAM_ID="baq${baq}_C${C}_q${minQ}_mq${minMapQ}"
        PREFIX_BASE="${RESULTS_DIR}/${DOMAIN}_${PARAM_ID}"
        mkdir -p "$PREFIX_BASE"
        echo "=== Processing parameter set: $PARAM_ID ==="

        # 1) Build region lists
        for ct in 4 5 6; do
          sites_file="${PREFIX_BASE}/domain_sfs_maf_sites_ct${ct}"
          regions_file="${PREFIX_BASE}/domain_sfs_maf_regions_ct${ct}"
          cut -f1 "$sites_file" | sort -u > "$regions_file"
        done

        # 2) Index site files
        angsd sites index "${PREFIX_BASE}/domain_sfs_maf_sites_ct4"
        angsd sites index "${PREFIX_BASE}/domain_sfs_maf_sites_ct5"
        angsd sites index "${PREFIX_BASE}/domain_sfs_maf_sites_ct6"

        # 3) Create per-pop directories
        while IFS= read -r POP; do
          for ct in 4 5 6; do
            mkdir -p "${PREFIX_BASE}/${POP}/${POP}_indvhet_ct${ct}"
          done
        done < "$POP_LIST"

        # 4) Run ANGSD + realSFS
        while IFS= read -r POP; do
          bamlist="${POP_BAM_DIR}/${POP}.txt"
          for ct in 4 5 6; do
            REGION="${PREFIX_BASE}/domain_sfs_maf_regions_ct${ct}"
            SITES="${PREFIX_BASE}/domain_sfs_maf_sites_ct${ct}"
            outdir="${PREFIX_BASE}/${POP}/${POP}_indvhet_ct${ct}"
            while IFS= read -r BAM_PATH; do
              BAM_FILE="$(basename "$BAM_PATH")"
              OUT="${BAM_FILE%%.*}"
              echo "--- $POP ct${ct} on $OUT ---"
              angsd \
                -i "$BAM_PATH" \
                -ref "$REF" \
                -anc "$REF" \
                -rf "$REGION" \
                -sites "$SITES" \
                -GL 1 \
                -doSaf 1 \
                -doCounts 1 \
                -baq $baq \
                -C $C \
                -minQ $minQ \
                -minMapQ $minMapQ \
                $DEFAULT_FILTERS \
                -out "${outdir}/${OUT}"
              realSFS "${outdir}/${OUT}.saf.idx" -tole 1e-7 -fold 1 > "${outdir}/${OUT}.sfs"
            done < "$bamlist"
          done
        done < "$POP_LIST"

      done
    done
  done
done

# ------------------------
# COLLECTION STEP
# ------------------------
echo "All parameter sweeps completed. Now collecting results…"
summary_file="${RESULTS_DIR}/indvhet_summary.tsv"
domain="$(basename "$RESULTS_DIR")"
domain="${domain%_results}"
{
  printf "domain\tbaq\tC\tq\tmq\tct\tpop_code\tbam_code\tref_sites\talt_sites\n"
  find "$RESULTS_DIR" -type f -name '*.sfs' | while IFS= read -r file; do
    rel="${file#$RESULTS_DIR/}"
    IFS='/' read -r param_id pop subdir bamfile <<< "$rel"
    pop_code="$pop"
    bam_code="${bamfile%.sfs}"
    ct="${subdir##*_ct}"
    params="${param_id#${domain}_}"
    IFS='_' read -r baq_val C_val q_val mq_val <<< "$params"
    baq="${baq_val#baq}"
    C="${C_val#C}"
    q="${q_val#q}"
    mq="${mq_val#mq}"
    read -r ref_sites alt_sites _ < "$file"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$domain" "$baq" "$C" "$q" "$mq" \
      "$ct" "$pop_code" "$bam_code" "$ref_sites" "$alt_sites"
  done
} > "$summary_file"
echo "Summary written to $summary_file"
