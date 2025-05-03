#!/usr/bin/bash
set -euo pipefail
IFS=$'\n\t'


# Script to compute per-sample call rates and mean depths from ANGSD counts.gz output
# Usage: summarize_counts.sh <input_counts_gz> <min_depth_threshold> <sample_names_file> <output_file>

# ------------------------ Arguments ------------------------
# Check if correct number of arguments is given
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_counts_gz_file> <minimumDp> <sample_names_file> <output_file>"
    exit 1
fi

input_file="$1"
threshold="$2"
column_names_file="$3"
output_file="$4"

# ------------------------ Read Sample Names ------------------------
# Read and clean column names from the file into an array
mapfile -t column_names < <(awk -F/ '{name=$NF; sub(/\..*/, "", name); print name}' "$column_names_file")

# ------------------------ Process Counts ------------------------
# Process the gzip file, removing the extra tab and skipping the first (header) row
zcat "$input_file" | awk -F'\t' -v threshold="$threshold" -v names_file="$column_names_file" -v outfile="$output_file" '
BEGIN {
    OFS = "\t";
    name_idx = 0;

    # Read cleaned column names into an array
    while ((getline name < names_file) > 0) {
        # Strip path and extension just in case (redundant if done already)
        gsub(/^.*\//, "", name);
        sub(/\..*/, "", name);
        column_names[++name_idx] = name;
    }
    close(names_file);
}
NR == 1 { next }  # Skip header row
{
    if ($NF == "") {
        NF--;
    }

    if (NR == 2) {
        num_cols = NF;
        print "DEBUG: Columns in counts.gz (after fixing extra tab):", num_cols > "/dev/stderr";
        print "DEBUG: Number of names in sample_names:", name_idx > "/dev/stderr";

        if (num_cols != name_idx) {
            print "Error: Number of column names does not match number of columns in data" > "/dev/stderr";
            exit 1;
        }
    }

    total_rows++;
    for (i = 1; i <= NF; i++) {
        if ($i > threshold) {
            count[i]++;
            sum[i] += $i;
        }
    }
}
END {
    print "bam_code", "call_rate", "mean_depth" > outfile;
    for (i = 1; i <= NF; i++) {
        call_rate = count[i] / total_rows;
        mean_above = (count[i] > 0) ? sum[i] / count[i] : "NA";
        print column_names[i], call_rate, mean_above > outfile;
    }
}'
