#!/usr/bin/env bash
set -euo pipefail

ml GCC/12.3.0
ml HTSlib/1.18

if [ $# -lt 2 ]; then
  echo "Usage: $0 INPUT_GLOB_OR_FILE OUTNAME [OUTDIR]" >&2
  exit 1
fi

pattern="$1"                   # e.g. /path/.../mq_forensics/*.advanced_MQs.tsv
OUTNAME="$2"                   # e.g. southern_mqforensics_summary.tsv
OUTDIR="${3:-.}"

# Paths
mkdir -p "$OUTDIR"
OUTPATH="${OUTDIR%/}/$OUTNAME"
SORTPATH="${OUTDIR%/}/${OUTNAME%.tsv}_sorted.tsv"

# Temp dir for sort spill files — unique to this job
jobtag="${SLURM_JOB_ID:-$$}"
SCRATCH="${SPRUCE_PROJECT}/temp/sorttmp.${jobtag}"
mkdir -p "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT

# Expand the input pattern inside the job
shopt -s nullglob
# shellcheck disable=SC2206  # intentional glob expansion to array
files=( $pattern )
if ((${#files[@]} == 0)); then
  echo "No inputs matched: $pattern" >&2
  exit 3
fi

echo "Inputs (${#files[@]}):"
printf '  %s\n' "${files[@]}"
echo "Saving sorted input to: $SORTPATH"
echo "Saving summary stats to: $OUTPATH"

# -------------------------------
# Emit header once, unsorted,
# then sort only the data rows
# -------------------------------
{
  # 1) Header from the first file (this line contains hist_mq_b0, etc.)
  head -n1 "${files[0]}"

  # 2) All data rows (no headers) → sort
  {
    for f in "${files[@]}"; do
      tail -n +2 "$f"
    done
  } | LC_ALL=C sort --stable -t $'\t' -S 3G --parallel=1 \
                    --temporary-directory="$SCRATCH" \
                    -k1,1 -k2,2n
} \
| tee "$SORTPATH" \
| mq_forensics summarize -i - -o "$OUTPATH"

echo "Done."
