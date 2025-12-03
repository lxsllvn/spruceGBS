#!/bin/bash
set -euo pipefail

# Load modules
ml GCC/12.3.0
ml HTSlib/1.18

# Submits an mq_forensics analyze job
usage() {
  cat >&2 <<'USAGE'
Usage:
  do_mqforensics.sh [options]

Options:
  -B, --bam PATH         BAM alignment (required)
  -b, --bed PATH         BED intervals (required)
  -r, --ref PATH         Reference genome (required)
  -C, --cap INT          samtools-style -c to use (required)
  -d, --minDP INT        Minimum read depth (required) 
  -o, --outname STR      Output file prefix (required)
  -D, --outdir DIR       Output directory (default: current directory)
  -h, --help             Show this help
  
Example:
  do_mqforensics.sh --bam sample.bam --bed regions.bed --ref ref.fa --cap 100 --minDP 3 --outname sample --outdir results
USAGE
}

# Default: current working directory
OUTDIR="."

# Ensure GNU getopt is available
command -v getopt >/dev/null || { echo "[error] GNU getopt not found"; exit 1; }

# Parse options with GNU getopt
if ! opts=$(getopt -o B:b:r:C:d:o:D:h \
  --long bam:,bed:,ref:,cap:,minDP:,outname:,outdir:,help \
  -- "$@"
); then
  usage
  exit 1
fi

eval set -- "$opts"

while true; do
  case "$1" in
    -B|--bam)         BAM="$2"; shift 2 ;;
    -b|--bed)         BED="$2"; shift 2 ;;
    -r|--ref)         REF="$2"; shift 2 ;;
    -C|--cap)         C="$2"; shift 2 ;;
    -d|--minDP)       MIN_DP="$2"; shift 2 ;;
    -o|--outname)     OUTNAME="$2"; shift 2 ;;
    -D|--outdir)      OUTDIR="$2"; shift 2 ;;
    -h|--help)        usage; exit 0 ;;
    --)               shift; break ;;
    *)                echo "Internal parsing error: $1" >&2; exit 1 ;;
  esac
done

# ----- Required checks -----
missing=()
for v in BAM BED REF C MIN_DP OUTNAME; do
  [[ -n "${!v-}" ]] || missing+=("$v")
done
if (( ${#missing[@]} )); then
  echo "Missing required argument(s): ${missing[*]}" >&2
  usage
  exit 1
fi

# Ensure output dir exists
mkdir -p "$OUTDIR"

OUTDIR="${OUTDIR%/}"
OUTPATH="${OUTDIR}/${OUTNAME}"

echo "Starting mq_forensics analyze at $(date) with:"
echo "  BAM     = $BAM"
echo "  BED     = $BED"
echo "  REF     = $REF"
echo "  CAP (-C)= $C"
echo "  minDP   = $MIN_DP"
echo "  OUTNAME = $OUTNAME"
echo "  OUTDIR  = $OUTDIR"
echo "  OUTPATH = $OUTPATH"


# -----------------------
# mq_forensics configuration 
# -----------------------
mq_forensics analyze \
  -b "$BAM" \
  -r "$BED" \
  -C "$C" \
  -d "$MIN_DP" \
  --emit-suffstats \
  --emit-hist \
  -f "$REF" \
  --flank 10 \
  --ref-only \
  -o "${OUTPATH}_advanced_MQs.tsv" \
  -O "${OUTPATH}_intervals.tsv"

echo "mq_forensics analyze completed at $(date)"

