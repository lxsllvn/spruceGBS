#!/usr/bin/env python3
import argparse, sys

# Import your existing modules (moved under src/beagleutils/)
from . import read_counts_by_genotype as rcg
from . import reheader_genotype_matrix as rh
from . import subset_genotype_matrix as sub

def build_parser():
    p = argparse.ArgumentParser(
        prog="geno-utils",
        description="Unified CLI for BEAGLE GL and counts matrix utilities."
    )
    subp = p.add_subparsers(dest="cmd", required=True)

    # reheader
    pr = subp.add_parser("reheader", help="Reheader BEAGLE or counts matrix")
    pr.add_argument("data_gz", help="Input (.beagle.gz or counts.tsv.gz)")
    pr.add_argument("names_file", help="sample_names.txt")
    pr.add_argument("outfile", help="Output .gz")
    pr.add_argument("--pos", help="sites.pos.gz (required for counts mode)", default=None)

    # subset
    ps = subp.add_parser("subset", help="Subset by SNPs and/or samples")
    ps.add_argument("input_file")
    ps.add_argument("output_file")
    ps.add_argument("--snps", default=None, help="list of snpcodes/markers")
    ps.add_argument("--samples", default=None, help="list of sample names")

    # split-genotype
    pg = subp.add_parser("split-genotype", help="Split counts by genotype using BEAGLE GLs")
    pg.add_argument("beagle_file")
    pg.add_argument("counts_file")
    pg.add_argument("out_prefix")
    pg.add_argument("--site-stats", dest="site_stats_file", default=None)
    pg.add_argument("--maf-threshold", type=float, default=0.05)

    return p

def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "reheader":
        if args.pos:
            rh.reheader_counts(args.data_gz, args.pos, args.names_file, args.outfile)
        else:
            rh.reheader_beagle(args.data_gz, args.names_file, args.outfile)

    elif args.cmd == "subset":
        sub.subset_matrix(args.input_file, args.output_file, args.snps, args.samples)

    elif args.cmd == "split-genotype":
        rcg.read_counts_by_genotype(
            beagle_file=args.beagle_file,
            counts_file=args.counts_file,
            out_prefix=args.out_prefix,
            site_stats_file=args.site_stats_file,
            maf_threshold=args.maf_threshold
        )

if __name__ == "__main__":
    sys.exit(main())
