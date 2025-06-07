#!/usr/bin/env python3
import sys
import gzip

def load_list(filename):
    if filename is None:
        return None
    with open(filename) as f:
        return set(line.strip() for line in f if line.strip())

def open_any(filename, mode="rt"):
    return gzip.open(filename, mode) if filename.endswith(".gz") else open(filename, mode)

def subset_matrix(input_file, output_file, snps_file=None, samples_file=None):
    snps_to_keep = load_list(snps_file)
    samples_to_keep = load_list(samples_file)

    with open_any(input_file, "rt") as infile, gzip.open(output_file, "wt") as out:
        header = infile.readline().rstrip('\n').split('\t')
        if header[0] not in ('snpcode', 'marker'):
            print(f"[ERROR] The input file does not have 'snpcode' or 'marker' as its first column!", file=sys.stderr)
            sys.exit(1)

        # Detect BEAGLE format (marker, allele1, allele2, then 3*sample columns)
        is_beagle = (len(header) > 3 and header[1] == 'allele1' and header[2] == 'allele2')

        if is_beagle:
            meta_cols = 3  # marker, allele1, allele2
            sample_names = header[meta_cols:]
            n_samples = len(sample_names) // 3
            true_sample_names = [sample_names[i*3] for i in range(n_samples)]

            # Decide which samples to keep
            if samples_to_keep is not None:
                wanted_sample_idx = [i for i, s in enumerate(true_sample_names) if s in samples_to_keep]
                keep_col_idx = list(range(meta_cols))  # always keep marker, allele1, allele2
                for idx in wanted_sample_idx:
                    keep_col_idx += [meta_cols + idx*3, meta_cols + idx*3 + 1, meta_cols + idx*3 + 2]
                # Build new header
                new_header = header[:meta_cols]
                for idx in wanted_sample_idx:
                    new_header += [sample_names[idx*3], sample_names[idx*3+1], sample_names[idx*3+2]]
            else:
                keep_col_idx = list(range(len(header)))
                new_header = header
        else:
            # counts or normal matrix: marker/snpcode + one col per sample
            meta_cols = 1
            sample_names = header[meta_cols:]
            if samples_to_keep is not None:
                wanted_sample_idx = [i for i, s in enumerate(sample_names) if s in samples_to_keep]
                keep_col_idx = [0] + [meta_cols + i for i in wanted_sample_idx]
                new_header = [header[0]] + [sample_names[i] for i in wanted_sample_idx]
            else:
                keep_col_idx = list(range(len(header)))
                new_header = header

        out.write('\t'.join(new_header) + '\n')

        for line in infile:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < len(header):
                # pad if needed (trailing missing columns)
                fields += [''] * (len(header) - len(fields))
            # Row filter: check SNPs/marker
            if snps_to_keep is not None and fields[0] not in snps_to_keep:
                continue
            # Column filter: subset samples (plus meta columns)
            row_out = [fields[i] for i in keep_col_idx]
            out.write('\t'.join(row_out) + '\n')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Subset reheadered counts/beagle file by snpcode/marker and/or sample names. Output is always gzipped."
    )
    parser.add_argument("input_file", help="Reheadered counts/beagle file (TSV, can be .gz)")
    parser.add_argument("output_file", help="Output file (will be gzipped)")
    parser.add_argument("--snps", help="List of snpcodes/markers to keep (one per line, optional)", default=None)
    parser.add_argument("--samples", help="List of sample names to keep (one per line, optional)", default=None)
    args = parser.parse_args()

    subset_matrix(
        input_file=args.input_file,
        output_file=args.output_file,
        snps_file=args.snps,
        samples_file=args.samples
    )
