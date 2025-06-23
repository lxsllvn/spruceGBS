#!/bin/bash
set -euo pipefail

# Creates site summary tables from *.pos.gz, *.counts.gz, *.hwe.gz, and *.snpStat.gz files

# Usage: $0 NAME INPATH [BASENAME]
#   INPATH:     path to folder containing *.pos.gz, *.counts.gz, *.hwe.gz, and *.snpStat.gz files
#   OUTNAME:    outname for site summary tables       
#   BASENAME:   the basename for the *.pos.gz, *.counts.gz etc. files; 
#               optional if INPATH only contains one set of input files. 

# Example usage
# Folder contains only one .pos.gz, one .counts.gz, etc:
# $0 ./summarize_site_stats.sh /path/to/folder output 

# Folder contains several (e.g., both siberia and southern files):
# $0 ./summarize_site_stats.sh /path/to/folder output siberia

if [ $# -lt 2 ]; then
  echo "Usage: $0 INPATH OUTNAME [BASENAME]" >&2
  exit 1
fi

INPATH="$1"
OUTNAME="$2"
BASENAME="${3:-}"

# Strip any trailing .tsv or .txt (case-insensitive) from OUTNAME
OUTNAME="${OUTNAME%.tsv}"
OUTNAME="${OUTNAME%.TSV}"
OUTNAME="${OUTNAME%.txt}"
OUTNAME="${OUTNAME%.TXT}"

ml GCC/13.3.0
ml Python/3.12.3
ml SciPy-bundle/2024.05

if [ -n "$BASENAME" ]; then
  python3 "$SCRIPTS/07_site_discovery/summarize_site_stats.py" "$INPATH" "${OUTNAME}.tsv" "$BASENAME"
else
  python3 "$SCRIPTS/07_site_discovery/summarize_site_stats.py" "$INPATH" "${OUTNAME}.tsv"
fi
