#!/usr/bin/env python3
import pandas as pd
import numpy as np
import gzip
import sys
import os
import glob

def parse_counts_line(line):
    values = [float(x) for x in line.rstrip('\n').split('\t') if x.strip() != ""]
    return np.array(values)

def parse_snpstat_line(line):
    fields = line.rstrip('\n').split('\t')
    SB1, SB2, SB3 = fields[3].split(':')
    HWE_LRT, HWE_pval = fields[4].split(':')
    baseQ_Z, baseQ_pval = fields[5].split(':')
    mapQ_Z, mapQ_pval = fields[6].split(':')
    edge_z, edge_pval = fields[7].split(':')
    return [SB1, SB2, SB3, HWE_LRT, HWE_pval, baseQ_Z, baseQ_pval, mapQ_Z, mapQ_pval, edge_z, edge_pval]

def find_unique_file(folder, pattern, label, basename=None):
    if basename:
        files = glob.glob(os.path.join(folder, f"{basename}{pattern}"))
        if len(files) == 1:
            return files[0]
        elif len(files) == 0:
            sys.exit(f"[Error] No {label} file matching '{basename}{pattern}' in {folder}")
        else:
            sys.exit(f"[Error] Multiple {label} files matching '{basename}{pattern}' in {folder}: {files}")
    else:
        files = glob.glob(os.path.join(folder, f"*{pattern}"))
        if len(files) == 1:
            return files[0]
        elif len(files) == 0:
            sys.exit(f"[Error] No {label} file (*{pattern}) found in {folder}")
        else:
            sys.exit(f"[Error] More than one {label} file (*{pattern}) found in {folder}: {files}\nSpecify a basename to disambiguate.")

def main():
    if len(sys.argv) not in [3, 4]:
        print(f"Usage: {sys.argv[0]} FOLDER OUTFILE.tsv [BASENAME]", file=sys.stderr)
        sys.exit(1)
    folder = sys.argv[1]
    outfile = sys.argv[2]
    basename = sys.argv[3] if len(sys.argv) == 4 else None

    pos_path = find_unique_file(folder, '.pos.gz', "positions", basename)
    counts_path = find_unique_file(folder, '.counts.gz', "counts", basename)
    hwe_path = find_unique_file(folder, '.hwe.gz', "HWE", basename)
    snpstat_path = find_unique_file(folder, '.snpStat.gz', "snpStat", basename)
    use_snpstat = os.path.isfile(snpstat_path)

    if not use_snpstat:
        print(f"[Warning] {snpstat_path} not found, output columns will be 'NA'", file=sys.stderr)

    print("Reading all .pos.gz into memory...")
    pos_df = pd.read_csv(pos_path, sep='\t', compression='gzip', dtype={'chr': str, 'pos': int, 'totDepth': float})
    pos_df['snpcode'] = pos_df['chr'].astype(str) + '_' + pos_df['pos'].astype(str)
    poscodes = pos_df['snpcode'].tolist()
    total_depths = pos_df['totDepth'].tolist()
    del pos_df

    print("Reading all .hwe.gz into memory...")
    hwe_df = pd.read_csv(hwe_path, sep='\t', compression='gzip', dtype={'Chromo': str, 'Position': int})
    hwe_df['snpcode'] = hwe_df['Chromo'].astype(str) + '_' + hwe_df['Position'].astype(str)
    MAF = hwe_df['hweFreq'].tolist()
    F = hwe_df['F'].tolist()
    Hexp = [2 * maf * (1 - maf) for maf in MAF]
    Hobs = [he - f * he for he, f in zip(Hexp, F)]
    del hwe_df

    if use_snpstat:
        print("Preparing to stream .snpStat.gz...")
        snpstat_fh = gzip.open(snpstat_path, 'rt')
        snpstat_header = snpstat_fh.readline()  # skip
    else:
        snpstat_fh = None

    maf05_outfile = outfile.replace('.tsv', '_maf05.tsv')
    print("Streaming through .counts.gz row by row and writing output...")
    with gzip.open(counts_path, 'rt') as counts, \
         open(outfile, 'w') as out, \
         open(maf05_outfile, 'w') as out_maf05:

        summary_cols = ['snpcode','total_depth','mean_depth','median_depth','call_rate','cv_depth','rel_IQR_depth','MAF','Hexp','Hobs','F']
        snpstat_cols = ['SB1','SB2','SB3','HWE_LRT','HWE_pval','baseQ_Z','baseQ_pval','mapQ_Z','mapQ_pval','edge_z','edge_pval']
        header = counts.readline() # skip
        out.write('\t'.join(summary_cols + snpstat_cols) + '\n')
        out_maf05.write('\t'.join(summary_cols + snpstat_cols) + '\n')

        for i, line in enumerate(counts):
            depth = parse_counts_line(line)
            if len(depth) == 0:
                mean_depth = median_depth = call_rate = cv_depth = rel_IQR_depth = 0
            else:
                nz = depth[depth > 0]
                mean_depth = np.mean(nz) if nz.size > 0 else 0
                median_depth = np.median(nz) if nz.size > 0 else 0
                call_rate = np.sum(depth > 0) / len(depth)
                cv_depth = (np.std(nz) / mean_depth) if mean_depth > 0 and nz.size > 1 else 0
                q75, q25 = np.percentile(nz, [75, 25]) if nz.size > 0 else (0, 0)
                rel_IQR_depth = ((q75 - q25) / median_depth) if median_depth > 0 else 0

            row = [
                poscodes[i], total_depths[i], mean_depth, median_depth, call_rate,
                cv_depth, rel_IQR_depth, MAF[i], Hexp[i], Hobs[i], F[i]
            ]

            if use_snpstat:
                snpstat_row = parse_snpstat_line(snpstat_fh.readline())
            else:
                snpstat_row = ['NA'] * 11

            out.write('\t'.join([str(x) for x in row + snpstat_row]) + '\n')
            if MAF[i] > 0.05:
                out_maf05.write('\t'.join([str(x) for x in row + snpstat_row]) + '\n')

            if (i+1) % 50000 == 0:
                print(f"  Processed {i+1} sites...")

    print("Done! Output written to", outfile)
    print("Filtered MAF > 0.05 output written to", maf05_outfile)

if __name__ == "__main__":
    main()
