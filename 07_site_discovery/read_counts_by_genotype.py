import gzip
import sys

def read_counts_by_genotype(beagle_file, counts_file, out_prefix):
    """
    Splits an annotated read depth matrix (must have snpcode column added) by genotype class (homRef, homAlt, and het). 
    Sample names in .beagle.gz and .counts.gz must match.

    Parameters
    ----------
    beagle_file : str
        Path to gzipped BEAGLE file (.beagle.gz, with sample triplets per locus).
    counts_file : str
        Path to gzipped counts file (tab-separated, 1st col SNP, rest are samples).
    out_prefix : str
        Prefix for the output files (.homRef.tsv.gz, .het.tsv.gz, .homAlt.tsv.gz).

    Outputs
    -------
    Three gzipped tab-delimited files:
      - {out_prefix}.homRef.tsv.gz
      - {out_prefix}.het.tsv.gz
      - {out_prefix}.homAlt.tsv.gz
    Each contains: snpcode, genotype, sample1, sample2, ...
    """
    # Output files for each genotype
    outfiles = {
        'homRef': gzip.open(f'{out_prefix}.homRef.tsv.gz', 'wt'),
        'het': gzip.open(f'{out_prefix}.het.tsv.gz', 'wt'),
        'homAlt': gzip.open(f'{out_prefix}.homAlt.tsv.gz', 'wt')
    }

    with gzip.open(beagle_file, 'rt') as bfile, gzip.open(counts_file, 'rt') as cfile:
        # Read headers
        bheader = next(bfile).rstrip('\n').split('\t')
        cheader = next(cfile).rstrip('\n').split('\t')

        # Get sample list (should match between files)
        samples = cheader[1:]
        if len(samples) != (len(bheader) - 3) // 3:
            raise ValueError(f"Sample number mismatch between files: {len(samples)} vs {(len(bheader) - 3)//3}")

        # Write header to each output
        for fh in outfiles.values():
            fh.write('snpcode\tgenotype\t' + '\t'.join(samples) + '\n')

        # Process each line in parallel
        for bl, cl in zip(bfile, cfile):
            bparts = bl.rstrip('\n').split('\t')
            cparts = cl.rstrip('\n').split('\t')
            snpcode = bparts[0]
            if cparts[0] != snpcode:
                raise ValueError(f"SNP mismatch: {bparts[0]} vs {cparts[0]}")

            # Build genotype assignments: 0=homRef, 1=het, 2=homAlt, -1=NA
            gt_assign = []
            for i in range(len(samples)):
                start = 3 + i*3
                gts = [float(bparts[start]), float(bparts[start+1]), float(bparts[start+2])]
                maxval = max(gts)
                if gts.count(maxval) > 1:
                    gt_assign.append(-1)  # NA (tie)
                else:
                    gt_assign.append(gts.index(maxval))

            # Prepare output lines for each genotype
            for gt, label in zip([0, 1, 2], ['homRef', 'het', 'homAlt']):
                outrow = [snpcode, label]
                for i, assign in enumerate(gt_assign):
                    count = cparts[i+1] if assign == gt else '0'
                    outrow.append(count)
                outfiles[label].write('\t'.join(outrow) + '\n')

    # Close output files
    for fh in outfiles.values():
        fh.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: python3 {sys.argv[0]} beagle_file counts_file out_prefix", file=sys.stderr)
        sys.exit(1)
    beagle_file = sys.argv[1]
    counts_file = sys.argv[2]
    out_prefix = sys.argv[3]
    read_counts_by_genotype(beagle_file, counts_file, out_prefix)
