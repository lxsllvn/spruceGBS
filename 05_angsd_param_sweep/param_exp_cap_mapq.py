#!/usr/bin/env python3
import argparse
import numpy as np
import pysam

def recalc_mq_for_read(read, fasta, C):
    """
    Given a pysam AlignedSegment and reference FASTA, compute:
      - AS1    : alignment score tag
      - SubQ   : sum of quality of high-Q mismatches (capped at Q33)
      - clipQ  : sum of soft-clip Q + 13*hard-clip
      - M      : count of non-N matched/mismatched, Q>=13 bases
      - T      : error penalty
      - NewMQ  : capped MQ per -C logic
      - RL     : read length
    Returns: (OrigMQ, NewMQ, AS1, SubQ, clipQ, T, RL)
    """
    # original mapping quality
    orig_mq = read.mapping_quality

    # alignment score tag
    try:
        AS1 = read.get_tag("AS")
    except KeyError:
        AS1 = None

    seq      = read.query_sequence
    quals    = read.query_qualities
    RL       = len(seq)
    ref_seq  = fasta.fetch(read.reference_name,
                            read.reference_start,
                            read.reference_end)

    mm = SubQ = clipQ = M = 0
    ref_idx = 0
    qpos    = 0

    for op, length in read.cigartuples:
        if op == 0:  # M/EQ/X
            for _ in range(length):
                q = quals[qpos]
                r = ref_seq[ref_idx]
                b = seq[qpos]
                if q >= 13 and r != 'N' and b != 'N':
                    M += 1
                    if b != r:
                        mm   += 1
                        SubQ += min(q, 33)
                qpos   += 1
                ref_idx += 1

        elif op == 1:  # insertion
            qpos += length

        elif op == 2:  # deletion
            ref_idx += length

        elif op == 4:  # soft-clip
            for _ in range(length):
                clipQ += quals[qpos]
                qpos   += 1

        elif op == 5:  # hard-clip
            clipQ += 13 * length

        else:
            # N, P, =, X etc: advance both
            ref_idx += length
            qpos    += length

    # Compute product = M^mm / mm!
    prod = 1.0
    for i in range(mm):
        prod *= M / (i + 1)

    # Compute T
    T = SubQ - 4.343 * np.log(prod) + clipQ / 5.0

    # Compute NewMQ
    if T > C:
        new_mq = 0
    else:
        if T < 0:
            T = 0
        new_mq = int(C * np.sqrt((C - T) / C) + 0.499)

    return orig_mq, new_mq, AS1, SubQ, clipQ, T, RL


def main():
    p = argparse.ArgumentParser(
        description="Recalculate MQ per -C caps from a BAM (with AS1, SubQ, clipQ, RL)"
    )
    p.add_argument("bam", help="Input BAM (must be indexed)")
    p.add_argument("fasta", help="Reference FASTA (must be indexed)")
    p.add_argument(
        "-C", "--caps", nargs="+", type=int,
        default=[0,50,60,75,100],
        help="List of cap values to test"
    )
    p.add_argument(
        "-m", "--max-reads", type=int, default=None,
        help="Max reads to process (for speed)"
    )
    p.add_argument(
        "-o", "--out", default="mq_recalc",
        help="CSV prefix (creates <out>_C<X>.csv files)"
    )
    args = p.parse_args()

    bam   = pysam.AlignmentFile(args.bam, "rb")
    fasta = pysam.FastaFile(args.fasta)

    # Open one CSV per C
    files = {}
    for C in args.caps:
        fn = f"{args.out}_C{C}.csv"
        f  = open(fn, "w")
        f.write("ReadName,OrigMQ,NewMQ,AS1,SubQ,ClipQ,RL,C,T\n")
        files[C] = f

    for i, read in enumerate(bam.fetch()):
        if args.max_reads and i >= args.max_reads:
            break
        if read.is_unmapped:
            continue

        for C, f in files.items():
            orig, new, AS1, subq, clipq, T, RL = recalc_mq_for_read(read, fasta, C)
            f.write(f"{read.query_name},{orig},{new},{AS1},{subq},{clipq},{RL},{C},{T:.2f}\n")

    bam.close()
    fasta.close()
    for f in files.values():
        f.close()
    print("Done. Outputs written with prefix:", args.out)

if __name__ == "__main__":
    main()
