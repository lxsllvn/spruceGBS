#!/bin/bash
set -euo pipefail

# Load modules
ml GCC/10.2.0
ml angsd/0.935

# Submits an ANGSD genotype likelihood job and collects site-level SNP quality metrics 
usage() {
  cat >&2 <<'USAGE'
Usage:
  angsd_likelihoods.sh REF REGION SITES BAMLIST OUTNAME [OUTDIR] [-- ANGSD_EXTRAS...]

Options:
  -r, --ref PATH         Reference genome (required)
  -g, --region PATH      ANGSD region file (required)
  -s, --sites PATH       ANGSD sites file (required)
  -b, --bamlist PATH     Bamlist file (required)
  -o, --outname NAME     Base name for outputs (required)
  -d, --outdir DIR       Output directory (default: current directory)
  -x, --angsd-args STR   Optional ANGSD args as a single string (can repeat)
  -n, --no-summary       Skip building site summary tables
  -h, --help             Show this help
  
Examples:
  # With flags
  angsd_likelihoods.sh --ref ref.fa --region reg.txt --sites sites.txt \
            --bamlist bamlist.txt --outname run1 --outdir results \
            --angsd-args "-doIBS 1 -doCov 1 -makeMatrix 1 -doSaf 1"

  # Legacy positional for backward compatibility
  angsd_likelihoods.sh ref.fa reg.txt sites.txt bamlist.txt run1 results \
    -- -doGlf 2 -doSaf 1
USAGE
}

# -----------------------
# Parse input options
# -----------------------
OUTDIR="."        # Default: current working directory
EXTRA_OPTS=()     # Default: no extra ANGSD args
MAKE_SUMMARY=1    # Default: build site summary tables
MAF_CUTOFF="0.05" # Default: build MAF > "$MAF_CUTOFF" site summary table

# Ensure GNU getopt is available
command -v getopt >/dev/null || { echo "[error] GNU getopt not found"; exit 1; }

opts=$(getopt -o r:g:s:b:o:d:x:hn \
  --long ref:,region:,sites:,bamlist:,outname:,outdir:,angsd-args:,help,no-summary \
  -n "$0" -- "$@") || { usage; exit 1; }
eval set -- "$opts"

while true; do
  case "$1" in
    -r|--ref)         REF="$2"; shift 2 ;;
    -g|--region)      REGION="$2"; shift 2 ;;
    -s|--sites)       SITES="$2"; shift 2 ;;
    -b|--bamlist)     BAM="$2"; shift 2 ;;
    -o|--outname)     OUTNAME="$2"; shift 2 ;;
    -d|--outdir)      OUTDIR="$2"; shift 2 ;;
    -x|--angsd-args)  read -r -a _tmp <<< "$2"; EXTRA_OPTS+=("${_tmp[@]}"); shift 2 ;;
    -n|--no-summary)  MAKE_SUMMARY=0; shift ;;
    -h|--help)        usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Internal parsing error: $1" >&2; exit 1 ;;
  esac
done

# ----- Backward-compatible positionals -----
# Order: REF REGION SITES BAMLIST OUTNAME [OUTDIR]
if [[ $# -gt 0 && -z "${REF-}"     ]]; then REF="$1";     shift; fi
if [[ $# -gt 0 && -z "${REGION-}"  ]]; then REGION="$1";  shift; fi
if [[ $# -gt 0 && -z "${SITES-}"   ]]; then SITES="$1";   shift; fi
if [[ $# -gt 0 && -z "${BAM-}"     ]]; then BAM="$1";     shift; fi
if [[ $# -gt 0 && -z "${OUTNAME-}" ]]; then OUTNAME="$1"; shift; fi
if [[ $# -gt 0 && "${1-}" != "--"  ]]; then OUTDIR="${1:-$OUTDIR}"; shift || true; fi

# Drop a leading `--` (if user separated extras) then capture remaining as ANGSD extras
if [[ "${1-}" == "--" ]]; then shift; fi
if [[ $# -gt 0 ]]; then EXTRA_OPTS+=("$@"); fi

# ----- Required checks -----
missing=()
for v in REF REGION SITES BAM OUTNAME; do
  [[ -n "${!v-}" ]] || missing+=("$v")
done
if (( ${#missing[@]} )); then
  echo "Missing required argument(s): ${missing[*]}" >&2
  usage; exit 1
fi

# Ensure output dir exists
mkdir -p "$OUTDIR"

OUTPATH="${OUTDIR%/}/$OUTNAME"

echo "Starting ANGSD at $(date) with:"
echo "REF=$REF"
echo "REGION=$REGION"
echo "SITES=$SITES"
echo "BAMLIST=$BAM"
echo "OUTNAME=$OUTNAME"
echo "OUTDIR=$OUTDIR"
echo "Output path prefix: $OUTPATH"

# -----------------------
# ANGSD configuration 
# -----------------------
P="2" # number of threads

# Filters
MAP_FILTERS="-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -C 100"
QUAL_FILTERS="-minMapQ 50 -minQ 20 -baq 0"
DEPTH_FILTERS="-setMinDepthInd 3"
SITE_FILTERS="-skipTriallelic 1"

# Genotype likelihood model (1 = SAMtools)
G="1"

# Default ANGSD operations
TODO="-doCounts 1 -dumpCounts 2 -doMajorMinor 4 -doMaf 1 -doHWE 1 -dosnpstat 1 -doGlf 2"

# -----------------------
# Run ANGSD
# -----------------------
angsd -P $P -b "$BAM" -ref "$REF" -anc "$REF" -rf "$REGION" -sites "$SITES" -out "$OUTPATH" -GL "$G" $TODO $MAP_FILTERS $QUAL_FILTERS $DEPTH_FILTERS $SITE_FILTERS "${EXTRA_OPTS[@]}"

{ printf "ANGSD command used: angsd -P %s -b %q -ref %q -anc %q -rf %q -sites %q -out %q -GL %q %s %s %s %s" \
    "$P" "$BAM" "$REF" "$REF" "$REGION" "$SITES" "$OUTPATH" "$G" \
    "$TODO" "$MAP_FILTERS" "$QUAL_FILTERS" "$DEPTH_FILTERS";
  printf " %s" $SITE_FILTERS; printf " "; printf "%q " "${EXTRA_OPTS[@]}"; echo;
} >&2


echo "ANGSD completed at $(date)"

# -----------------------
# Site summary tables
# -----------------------
if [[ "$MAKE_SUMMARY" -eq 1 ]]; then
  echo "Building site summaries for ${OUTNAME} in ${OUTDIR} at $(date)"
  
# Set paths
POS_GZ="${OUTPATH}.pos.gz"
COUNTS_GZ="${OUTPATH}.counts.gz"
MAFS_GZ="${OUTPATH}.mafs.gz"
HWE_GZ="${OUTPATH}.hwe.gz"
SNPSTAT_GZ="${OUTPATH}.snpStat.gz"   
OUT_ALL_GZ="${OUTPATH}_site_summary.tsv.gz"
OUT_MAF05_GZ="${OUTPATH}_site_summary_maf05.tsv.gz"

# Create temp files
tpos=$(mktemp)
tcr=$(mktemp)
thwe=$(mktemp)
tsnp=$(mktemp)          
theader=$(mktemp)
trap 'rm -f "$tpos" "$tcr" "$thwe" "$tsnp" "$theader"' EXIT

# Extract snpcode from .pos.gz  (also captures chr and pos for safety)
#   Expect header with columns including 'chr' and 'pos'
#   Output columns: snpcode \t chr \t pos
gzip -cd "${POS_GZ}" \
| awk -F'\t' '
  BEGIN{OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){
      if($i=="chr" || $i=="Chromo") c=i
      if($i=="pos" || $i=="Position") p=i
    }
    if(!c||!p){print "[error] chr/pos header not found in .pos.gz" > "/dev/stderr"; exit 1}
    next
  }
  {
    printf "%s_%s\t%s\t%s\n", $c, $p, $c, $p
  }' > "${tpos}"

npos=$(wc -l < "${tpos}")
echo "[info] pos rows: ${npos}"

# call_rate from .mafs.gz (nInd / n_samples)
if [[ ! -s "${MAFS_GZ}" ]]; then
  echo "[error] ${MAFS_GZ} not found or empty." >&2; exit 1
fi
# Count non-empty lines in BAM as n_samples
n_samp=$(grep -cve '^[[:space:]]*$' "${BAM}")
if [[ "${n_samp}" -le 0 ]]; then
  echo "[error] BAM file appears empty (no non-empty lines)." >&2; exit 1
fi

gzip -cd "${MAFS_GZ}" \
| awk -v n="${n_samp}" -F'\t' 'BEGIN{OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){ if($i=="nInd") ncol=i }
    if(!ncol){print "[error] nInd column not found in .mafs.gz" > "/dev/stderr"; exit 1}
    next
  }
  {
    # Coerce to number; clamp to [0,n] to be safe
    val = $ncol + 0.0
    if(val < 0) val=0
    if(val > n) val=n
    printf "%.6f\n", (n>0 ? val/n : 0.0)
  }' > "${tcr}"
ncr=$(wc -l < "${tcr}"); echo "[info] call_rate rows: ${ncr}"

# MAF, F, Hexp, Hobs from .hwe.gz (row order assumed to match)
gzip -cd "${HWE_GZ}" \
| awk -F'\t' '
  BEGIN{OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){
      hn[$i]=i
    }
    req="hweFreq F"
    split(req,arr," ")
    for(k in arr){
      if(!(arr[k] in hn)){
        print "[error] Missing column in .hwe.gz: " arr[k] > "/dev/stderr"
        exit 1
      }
    }
    maf_i=hn["hweFreq"]; f_i=hn["F"]
    next
  }
  {
    maf = $maf_i+0.0
    F   = $f_i+0.0
    he  = 2.0*maf*(1.0-maf)
    ho  = he*(1.0 - F)
    printf "%.6f\t%.6f\t%.6f\t%.6f\n", maf, he, ho, F
  }' > "${thwe}"

nhwe=$(wc -l < "${thwe}")
echo "[info] hwe rows: ${nhwe}"

# snpStat columns
#  Note: the snpStat format uses tabs, spaces, and : as delimiters. The tab-delimited fields are:
#    1: Chromo
#    2: Position
#    3: "+Major +Minor -Major -Minor"  
#    4: "SB1:SB2:SB3"
#    5: "HWE_LRT:HWE_pval"
#    6: "baseQ_Z:baseQ_pval"
#    7: "mapQ_Z:mapQ_pval"
#    8: "edge_z:edge_pval"

gzip -cd "${SNPSTAT_GZ}" \
| awk -F'\t' '
  BEGIN{OFS="\t"}
  NR==1{next}  # skip header
  {
    if (NF < 8) {
      print "[error] Malformed line in snpStat file (expected ≥8 fields, got " NF "): " $0 > "/dev/stderr"
      exit 1
    }
    split($4, sb, ":")
    split($5, hw, ":")
    split($6, bq, ":")
    split($7, mq, ":")
    split($8, ed, ":")
    if (length(sb) < 3 || length(hw) < 2 || length(bq) < 2 || length(mq) < 2 || length(ed) < 2) {
      print "[error] Malformed tuple in snpStat file at line " NR ": " $0 > "/dev/stderr"
      exit 1
    }
    print sb[1], sb[2], sb[3], hw[1], hw[2], bq[1], bq[2], mq[1], mq[2], ed[1], ed[2]
  }' > "${tsnp}"

nsnp=$(wc -l < "${tsnp}")
echo "[info] snpStat rows: ${nsnp}"

# Sanity checks on row counts
if [[ "${npos}" -ne "${ncr}" || "${npos}" -ne "${nhwe}" || "${npos}" -ne "${nsnp}" ]]; then
  echo "[error] Row-count mismatch among inputs:
  pos=${npos}, call_rate=${ncr}, hwe=${nhwe}, snpstat=${nsnp}
  Aborting to avoid misalignment." >&2
  exit 1
fi

# Write headers
#   Final columns:
#     snpcode, call_rate, MAF, Hexp, Hobs, F, SB1, SB2, SB3, HWE_LRT, HWE_pval, baseQ_Z, baseQ_pval, mapQ_Z, mapQ_pval, edge_z, edge_pval
printf "snpcode\tcall_rate\tMAF\tHexp\tHobs\tF\tSB1\tSB2\tSB3\tHWE_LRT\tHWE_pval\tbaseQ_Z\tbaseQ_pval\tmapQ_Z\tmapQ_pval\tedge_z\tedge_pval\n" > "${theader}"

# Assemble and write outputs (all sites + MAF>0.05 & <0.95)
#   We only need snpcode from tpos (first field). Use paste to combine columns line-wise.
#   paste: tpos(1: snpcode) + tcr(call_rate) + thwe(MAF..F) + tsnp(11 cols)
{
  cat "${theader}"
  paste \
    <(cut -f1 "${tpos}") \
    "${tcr}" \
    "${thwe}" \
    "${tsnp}"
} | gzip -c > "${OUT_ALL_GZ}"

echo "Wrote ${OUT_ALL_GZ} at $(date)"

# Filter for MAF > 0.05 and MAF < 0.95 (MAF is column 3 in final table)
gzip -cd "${OUT_ALL_GZ}" \
| awk -F'\t' -v OFS="\t" -v m="${MAF_CUTOFF}" '
    NR==1{print; next}
    ($3+0.0) > m && ($3+0.0) < (1.0 - m)
  ' | gzip -c > "${OUT_MAF05_GZ}"

echo "Wrote ${OUT_MAF05_GZ} at $(date)"

# Quick diagnostics
n_all=$(gzip -cd "${OUT_ALL_GZ}" | wc -l | awk '{print $1-1}')
n_maf=$(gzip -cd "${OUT_MAF05_GZ}" | wc -l | awk '{print $1-1}')
ONE_MINUS=$(echo "1 - $MAF_CUTOFF" | bc -l)
echo "Final row counts: all=${n_all}, maf>${MAF_CUTOFF} & <${ONE_MINUS} = ${n_maf}"

# Optional: compare first 3 snpcodes across pos/hwe/snpstat alignment
echo "First 3 snpcodes (sanity):"
head -n 3 "${tpos}" | cut -f1

else
  echo "Skipping site summary tables (--no-summary enabled)"
fi
