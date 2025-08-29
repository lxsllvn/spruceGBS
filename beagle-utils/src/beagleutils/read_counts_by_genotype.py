import gzip
import sys

def read_counts_by_genotype(beagle_file, counts_file, out_prefix, site_stats_file=None, maf_threshold=0.05):
    """
    Splits an annotated read depth matrix (must have snpcode column added) by genotype class (homRef, homAlt, and het).
    Optionally, outputs filtered files by MAF using site_stats.

    Parameters
    ----------
    beagle_file : str
        Path to gzipped BEAGLE file (.beagle.gz, with sample triplets per locus).
    counts_file : str
        Path to gzipped counts file (tab-separated, 1st col SNP, rest are samples).
    out_prefix : str
        Prefix for the output files (.homRef.tsv.gz, .het.tsv.gz, .homAlt.tsv.gz).
    site_stats_file : str or None
        Path to gzipped site stats file (must have snpcode and MAF columns).
    maf_threshold : float
        MAF threshold (default 0.05)

    Outputs
    -------
    Three gzipped tab-delimited files:
      - {out_prefix}.homRef.tsv.gz
      - {out_prefix}.het.tsv.gz
      - {out_prefix}.homAlt.tsv.gz
    Each contains: snpcode, genotype, sample1, sample2, ...

    If site_stats_file is provided, also outputs MAF-filtered files:
      - {out_prefix}.homRef.maf05.tsv.gz
      - {out_prefix}.het.maf05.tsv.gz
      - {out_prefix}.homAlt.maf05.tsv.gz
    """
    outfiles = {
        'homRef': gzip.open(f'{out_prefix}.homRef.tsv.gz', 'wt'),
        'het': gzip.open(f'{out_prefix}.het.tsv.gz', 'wt'),
        'homAlt': gzip.open(f'{out_prefix}.homAlt.tsv.gz', 'wt')
    }

    maf_snp_set = None
    maf_outfiles = None
    if site_stats_file:
        maf_snp_set = set()
        with gzip.open(site_stats_file, 'rt') as statfile:
            stat_header = next(statfile).rstrip('\n').split('\t')
            snpcode_idx = stat_header.index('snpcode')
            maf_idx = stat_header.index('MAF')
            for line in statfile:
                parts = line.rstrip('\n').split('\t')
                try:
                    maf = float(parts[maf_idx])
                except ValueError:
                    continue
                if maf > maf_threshold:
                    maf_snp_set.add(parts[snpcode_idx])
        maf_outfiles = {
            'homRef': gzip.open(f'{out_prefix}.homRef.maf05.tsv.gz', 'wt'),
            'het': gzip.open(f'{out_prefix}.het.maf05.tsv.gz', 'wt'),
            'homAlt': gzip.open(f'{out_prefix}.homAlt.maf05.tsv.gz', 'wt')
        }

    with gzip.open(beagle_file, 'rt') as bfile, gzip.open(counts_file, 'rt') as cfile:
        bheader = next(bfile).rstrip('\n').split('\t')
        cheader = next(cfile).rstrip('\n').split('\t')

        samples = cheader[1:]
        beagle_samples = [bheader[i] for i in range(3, len(bheader), 3)]
        if samples != beagle_samples:
            raise ValueError(
                f"Sample names/order mismatch between files!\n"
                f"Counts: {samples}\nBeagle: {beagle_samples}"
            )

        for fh in outfiles.values():
            fh.write('snpcode\tgenotype\t' + '\t'.join(samples) + '\n')
        if maf_outfiles:
            for fh in maf_outfiles.values():
                fh.write('snpcode\tgenotype\t' + '\t'.join(samples) + '\n')

        for bl, cl in zip(bfile, cfile):
            bparts = bl.rstrip('\n').split('\t')
            cparts = cl.rstrip('\n').split('\t')
            snpcode = bparts[0]
            if cparts[0] != snpcode:
                raise ValueError(f"SNP mismatch: {bparts[0]} vs {cparts[0]}")

            gt_assign = []
            for i in range(len(samples)):
                start = 3 + i*3
                gts = [float(bparts[start]), float(bparts[start+1]), float(bparts[start+2])]
                maxval = max(gts)
                if gts.count(maxval) > 1:
                    gt_assign.append(-1)
                else:
                    gt_assign.append(gts.index(maxval))

            for gt, label in zip([0, 1, 2], ['homRef', 'het', 'homAlt']):
                outrow = [snpcode, label]
                for i, assign in enumerate(gt_assign):
                    count = cparts[i+1] if assign == gt else '0'
                    outrow.append(count)
                outfiles[label].write('\t'.join(outrow) + '\n')
                if maf_outfiles and snpcode in maf_snp_set:
                    maf_outfiles[label].write('\t'.join(outrow) + '\n')

    for fh in outfiles.values():
        fh.close()
    if maf_outfiles:
        for fh in maf_outfiles.values():
            fh.close()

if __name__ == "__main__":
    if len(sys.argv) not in [4, 5]:
        print(f"Usage: python3 {sys.argv[0]} beagle_file counts_file out_prefix [site_stats.gz]", file=sys.stderr)
        sys.exit(1)
    beagle_file = sys.argv[1]
    counts_file = sys.argv[2]
    out_prefix = sys.argv[3]
    site_stats_file = sys.argv[4] if len(sys.argv) == 5 else None
    read_counts_by_genotype(beagle_file, counts_file, out_prefix, site_stats_file)
