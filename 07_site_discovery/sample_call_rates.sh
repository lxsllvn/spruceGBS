#!/bin/bash
set -euo pipefail

# -----------------------------------------------------
# Get per-sample call rates & mean depth from ANGSD *.count.gz files
# -----------------------------------------------------
# Usage:
#   ./sample_call_rates.sh <counts_gz_file(s)> <minimumDP> <sample_names_file> [output_file]
# Arguments:
#   <counts_gz_file(s)>  - e.g., *.counts.gz, mysamples.counts.gz, PATH/TO/*.counts.gz. 
#                          Don't quote the glob. 
#   <minimumDp>          - minimum DP for a call
#   [sample_names_file]  - list of sample codes, same order as ANGSD bamlist; optional if *.counts.gz has informative headers
#   [output_file]        - (optional) path to save output

#!/bin/bash
set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <counts_gz_file(s)> <minimumDp> [sample_names_file] [output_file]"
    exit 1
fi

# Get positional parameters from the end, but be flexible about presence/absence of sample_names_file
# We'll check for snpcode in the header to decide what is required
all_args=("$@")
num_args=$#
minDP="${all_args[$((num_args-2))]}"

# Combine all counts files into one temp file
tmpfile=$(mktemp)
first=1
for f in "${all_args[@]:0:$((num_args-2))}"; do
    if [ $first -eq 1 ]; then
        zcat "$f" > "$tmpfile"
        first=0
    else
        zcat "$f" | tail -n +2 >> "$tmpfile"
    fi
done

# Read header line and check for snpcode
header=$(head -n1 "$tmpfile")
if [[ "$header" =~ ^snpcode[[:space:]] ]]; then
    has_snpcode=1
else
    has_snpcode=0
fi

# Decide which args are sample_names_file and output_file
if [[ $has_snpcode -eq 1 ]]; then
    # snpcode present: sample_names_file is optional
    if [ "$num_args" -ge 3 ]; then
        possible_out="${all_args[$((num_args-1))]}"
        if [[ "$possible_out" =~ \.tsv$|\.txt$ ]]; then
            output_file="$possible_out"
        else
            output_file=""
        fi
    else
        output_file=""
    fi
    names_file=""
else
    # snpcode NOT present: sample_names_file is required
    if [ "$num_args" -lt 3 ]; then
        echo "Error: sample_names_file required if snpcode column is not present in input file(s)." >&2
        exit 1
    fi
    possible_out="${all_args[$((num_args-1))]}"
    possible_names="${all_args[$((num_args-2))]}"
    if [[ "$possible_out" =~ \.tsv$|\.txt$ ]]; then
        output_file="$possible_out"
        names_file="$possible_names"
    else
        output_file=""
        names_file="$possible_out"
    fi
    if [ ! -f "$names_file" ]; then
        echo "Error: sample_names_file ($names_file) does not exist or is not a file." >&2
        exit 1
    fi
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
# Export variables for awk
if [[ $has_snpcode -eq 1 ]]; then
    if [ -n "$output_file" ]; then
        awk -F'\t' -v threshold="$minDP" -v has_snpcode=1 -v tmpfile="$tmpfile" "$awk_script" "$tmpfile" > "$output_file"
    else
        awk -F'\t' -v threshold="$minDP" -v has_snpcode=1 -v tmpfile="$tmpfile" "$awk_script" "$tmpfile"
    fi
else
    if [ -n "$output_file" ]; then
        awk -F'\t' -v threshold="$minDP" -v has_snpcode=0 -v names_file="$names_file" "$awk_script" "$tmpfile" > "$output_file"
    else
        awk -F'\t' -v threshold="$minDP" -v has_snpcode=0 -v names_file="$names_file" "$awk_script" "$tmpfile"
    fi
fi

rm "$tmpfile"
