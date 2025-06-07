#!/usr/bin/awk -f
# Usage: ./filter_by_col.awk input.tsv SB1 '>0.2' call_rate '>0.72' ...

BEGIN {
    FS = OFS = "\t"
    # Parse command line arguments after the file name
    n_cond = (ARGC - 2) / 2
    for (i = 2; i < ARGC; i += 2) {
        cond_col[(i-2)/2] = ARGV[i]
        cond_expr[(i-2)/2] = ARGV[i+1]
        delete ARGV[i]
        delete ARGV[i+1]
    }
}

NR == 1 {
    # Build a mapping from column names to column numbers
    for (i = 1; i <= NF; i++) {
        colname[$i] = i
    }
    print
    next
}

{
    pass = 1
    for (j = 0; j < n_cond; j++) {
        cname = cond_col[j]
        expr  = cond_expr[j]
        cidx  = colname[cname]
        val   = $cidx

        # Parse the operator and threshold
        if (expr ~ /^ *>=/)      { op = ">="; thresh = substr(expr, 3) }
        else if (expr ~ /^ *<=/) { op = "<="; thresh = substr(expr, 3) }
        else if (expr ~ /^ *!=/) { op = "!="; thresh = substr(expr, 3) }
        else if (expr ~ /^ *==/) { op = "=="; thresh = substr(expr, 3) }
        else if (expr ~ /^ *</)  { op = "<";  thresh = substr(expr, 2) }
        else if (expr ~ /^ *>/)  { op = ">";  thresh = substr(expr, 2) }
        else { print "Error: unknown operator in condition: " expr > "/dev/stderr"; exit 1 }

        # Remove whitespace from threshold
        gsub(/^ +| +$/, "", thresh)

        # Numeric comparison
        if      (op == ">=") { pass = pass && (val+0 >= thresh+0) }
        else if (op == "<=") { pass = pass && (val+0 <= thresh+0) }
        else if (op == "!=") { pass = pass && (val+0 != thresh+0) }
        else if (op == "==") { pass = pass && (val+0 == thresh+0) }
        else if (op == "<")  { pass = pass && (val+0 <  thresh+0) }
        else if (op == ">")  { pass = pass && (val+0 >  thresh+0) }
    }
    if (pass) print
}
