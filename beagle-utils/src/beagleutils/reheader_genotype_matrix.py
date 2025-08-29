#!/usr/bin/env python3
import gzip
import sys

def get_sample_names(names_file):
    with open(names_file) as nf:
        sample_names = [line.strip() for line in nf if line.strip()]
    return sample_names

def get_snp_codes_from_pos(pos_gz):
    snpcodes = []
    with gzip.open(pos_gz, 'rt') as pf:
        pf.readline()  # skip header
        for line in pf:
            fields = line.rstrip('\n').split('\t')
            snpcodes.append(fields[0] + '_' + fields[1])
    return snpcodes

def reheader_counts(counts_gz, pos_gz, names_file, outfile):
    sample_names = get_sample_names(names_file)
    snpcodes = get_snp_codes_from_pos(pos_gz)
    n_samples = len(sample_names)

    with gzip.open(counts_gz, 'rt') as cf, gzip.open(outfile, 'wt') as out:
        cf.readline()  # skip original header
        out.write("snpcode\t" + "\t".join(sample_names) + "\n")
        for i, (snpcode, line) in enumerate(zip(snpcodes, cf)):
            vals = [x for x in line.rstrip('\n').split('\t') if x.strip() != ""]
            if len(vals) > n_samples:
                vals = vals[:n_samples]
            if len(vals) < n_samples:
                vals += [''] * (n_samples - len(vals))
            if i == 0:
                print(f"First data row: {len(vals)} columns; header: {n_samples} names", flush=True)
            out.write(f"{snpcode}\t" + "\t".join(vals) + "\n")

def reheader_beagle(beagle_gz, names_file, outfile):
    sample_names = get_sample_names(names_file)
    n_samples = len(sample_names)
    with gzip.open(beagle_gz, 'rt') as cf, gzip.open(outfile, 'wt') as out:
        header = cf.readline().rstrip('\n').split('\t')
        expected_cols = 3 + 3 * n_samples
        actual_cols = len(header)
        if actual_cols != expected_cols:
            raise ValueError(f"Header columns: {actual_cols}, expected: {expected_cols} (3 + 3 * n_samples). Check your sample_names.txt and beagle file format.")

        # Write new header: marker, allele1, allele2, then sample names x3
        new_header = header[:3] + [name for name in sample_names for _ in range(3)]
        out.write('\t'.join(new_header) + '\n')

        for i, line in enumerate(cf):
            vals = line.rstrip('\n').split('\t')
            # Pad/truncate to correct column count
            if len(vals) > expected_cols:
                vals = vals[:expected_cols]
            elif len(vals) < expected_cols:
                vals += [''] * (expected_cols - len(vals))
            if i == 0:
                print(f"First data row: {len(vals)} columns; header: {n_samples} names (expect {expected_cols} cols)", flush=True)
            out.write('\t'.join(vals) + '\n')

if __name__ == "__main__":
    if len(sys.argv) not in (4, 5):
        print(f"Usage: {sys.argv[0]} <data.gz> <sample_names.txt> <output.tsv.gz> [pos.gz]", file=sys.stderr)
        print("  If [pos.gz] is given, adds snpcode column from pos.gz. Otherwise, assumes first column is marker/BEAGLE.", file=sys.stderr)
        sys.exit(1)
    if len(sys.argv) == 5:
        # counts.gz mode (add snpcode)
        reheader_counts(sys.argv[1], sys.argv[4], sys.argv[2], sys.argv[3])
    else:
        # beagle.gz mode (header only, BEAGLE format with marker, allele1, allele2)
        reheader_beagle(sys.argv[1], sys.argv[2], sys.argv[3])
