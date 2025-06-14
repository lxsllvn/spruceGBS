#!/bin/bash
set -euo pipefail

# Helper for error messages
die() { echo "$@" >&2; exit 1; }

# -------------- Parse arguments ----------------
counts_file=""
minDP=""
sample_list=""
output_file=""

# Parse long options
while [[ $# -gt 0 ]]; do
    key="$1"
    case "$key" in
        --counts)
            counts_file="$2"
            shift 2
            ;;
        --minDP)
            minDP="$2"
            shift 2
            ;;
        --sample_list)
            sample_list="$2"
            shift 2
            ;;
        --out)
            output_file="$2"
            shift 2
            ;;
        -*)
            die "Unknown option: $1"
            ;;
        *)
            die "Unexpected argument: $1"
            ;;
    esac
done

# ----------- Check required arguments ----------
[[ -z "$counts_file" ]] && die "Missing required argument: --counts counts.gz"
[[ -z "$minDP" ]] && die "Missing required argument: --minDP <minDP>"

[[ ! -f "$counts_file" ]] && die "Counts file '$counts_file' not found!"

# ------------- Process counts file ------------
tmpfile=$(mktemp)
zcat "$counts_file" > "$tmpfile"

header=$(head -n1 "$tmpfile")
if [[ "$header" =~ ^snpcode[[:space:]] ]]; then
    has_snpcode=1
else
    has_snpcode=0
fi

if [[ $has_snpcode -eq 1 ]]; then
    # snpcode column present; sample_list is optional
    :
else
    # snpcode column absent; sample_list required
    [[ -z "$sample_list" ]] && die "Missing required argument: --sample_list <samples.txt> (needed when no snpcode column in counts file!)"
    [[ ! -f "$sample_list" ]] && die "Sample list file '$sample_list' not found!"
fi

awk_script='
BEGIN {
    OFS="\t";
    skip_col = 0;
    num_samples = 0;
    if (has_snpcode == 1) {
        # Read header line for sample names
        getline header < tmpfile;
        split(header, h, "\t");
        skip_col = 1;
        for (i = 2; i <= length(h); i++) {
            column_names[i-1] = h[i];
        }
        num_samples = length(h) - 1;
    } else {
        # Read sample names from file
        name_idx = 0;
        while ((getline name < names_file) > 0) {
            column_names[++name_idx] = name;
        }
        num_samples = name_idx;
    }
}
NR == 1 { if (has_snpcode == 1) next }
{
    if ($NF == "") NF--;
    total_rows++;
    for (i = 1+skip_col; i <= NF; i++) {
        if ($i > threshold) {
            count[i-skip_col]++;
            sum_depth[i-skip_col] += $i;
        }
    }
}
END {
    print "sample", "call_rate", "mean_depth";
    for (i = 1; i <= num_samples; i++) {
        call_rate = (count[i] ? count[i] / total_rows : 0);
        mean_depth = (count[i] ? sum_depth[i] / count[i] : "NA");
        print column_names[i], call_rate, mean_depth;
    }
}'

# ------------- Run AWK --------------------------
if [[ $has_snpcode -eq 1 ]]; then
    if [[ -n "$output_file" ]]; then
        awk -F'\t' -v threshold="$minDP" -v has_snpcode=1 -v tmpfile="$tmpfile" "$awk_script" "$tmpfile" > "$output_file"
    else
        awk -F'\t' -v threshold="$minDP" -v has_snpcode=1 -v tmpfile="$tmpfile" "$awk_script" "$tmpfile"
    fi
else
    if [[ -n "$output_file" ]]; then
        awk -F'\t' -v threshold="$minDP" -v has_snpcode=0 -v names_file="$sample_list" "$awk_script" "$tmpfile" > "$output_file"
    else
        awk -F'\t' -v threshold="$minDP" -v has_snpcode=0 -v names_file="$sample_list" "$awk_script" "$tmpfile"
    fi
fi

rm "$tmpfile"
