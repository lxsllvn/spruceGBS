#!/usr/bin/env bash
set -euo pipefail

# ------------------------
# Site Discovery per Domain with ANGSD
# ------------------------
# Usage: $0 DOMAIN BAMLIST REFDIR OUTDIR REF
#   DOMAIN:   sample domain name (e.g., "southern")
#   BAMLIST:  path to bamlist file (e.g., "southern_bam_list.txt"
#   REFDIR:   directory where reference subsets (from split_reference.sh) are located (e.g., "${SPRUCE_PROJECT}/ref/subsets")
#   OUTDIR:   directory name for outputs
#   REF:      path to complete reference FASTA (indexed for samtools)
#             e.g., "${SPRUCE_PROJECT}/ref/picea_newref.fa"

if [ "$#" -lt 5 ]; then
  echo "Usage: $0 DOMAIN BAMLIST REFDIR OUTDIR REF" >&2
  exit 1
fi

DOMAIN="$1"
BAMLIST="$2"
REFDIR="$3"
OUTDIR="$4"
REF="$5"

mkdir -p "$OUTDIR"

# Submit and return only the JobID (works with or without --parsable)
submit_and_get_id() {
  # shellcheck disable=SC2068
  "$SCRIPTS/submit.sh" $@ | awk '{print $NF}'
}

jobids=()

# 1) Submit 23 discovery jobs
for i in $(seq -w 01 23); do
  jid=$(submit_and_get_id \
    -J "site_discovery_${i}" \
    "$SCRIPTS/07_site_discovery/site_discovery.sh" \
      "${i}" \
      "${DOMAIN}" \
      "${BAMLIST}" \
      "${REFDIR}" \
      "${OUTDIR}")
  echo "Submitted discovery split ${i} as JobID ${jid}"
  jobids+=("$jid")
done

# Ensure we actually have 23 job IDs
if [[ "${#jobids[@]}" -ne 23 ]]; then
  echo "ERROR: expected 23 job IDs, got ${#jobids[@]}" >&2
  exit 1
fi

# 2) Build dependency string: afterok:job1:job2:...:job23
dep="afterok:$(IFS=:; echo "${jobids[*]}")"

# 3) Submit collator with dependency
# $OUTDIR from step one is the path to the discovery subsets
# we use the same $OUTDIR for this step
collate_jid=$(submit_and_get_id \
  -J "prepare_angsd_ref" \
  -d "$dep" \
  "$SCRIPTS/07_site_discovery/prepare_angsd_ref.sh" \
    --splits \
    "${DOMAIN}" \
    "${OUTDIR}" \
    "${REF}" \
    "${OUTDIR}")

echo "Submitted collator as JobID ${collate_jid} with dependency ${dep}"