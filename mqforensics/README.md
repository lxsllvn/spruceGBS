# mq_forensics

Fast C/htslib tool to extract per-site (within a BED) and per-interval mapping statistics from BAMs. Designed for GBS/amplicon panels (short loci), scalable to thousands of BAMs.

- Whole-read (WR) attribution for SubQ, ClipQ, capped MQ (`-C`) — each read contributes its single value to every site it overlaps.
- Local site features (mismatches, indel lengths at start position) derived from CIGAR + MD.
- Interval metric: aligned_bp_in_interval (area under coverage).
- Two analysis modes:
  - direct mode (per-site means/medians/SDs from raw read vectors)
  - sufficient statistics mode (`--emit-suffstats`) — compact, streaming-friendly per-site triplets (n,sum,sumsq) for each metric.
- Optional histograms (`--emit-hist`) allow divergences (KS, Wasserstein-1, Jensen–Shannon) and histogram-based medians to be computed later in summarize.
- Optional flanking-region context (`--flank N`, `--ref-only`) computes per-site windowed means of coverage and clip-fraction in the surrounding region, optionally restricted to ref-matching bases only.
- Summarizer reduces hundreds of per-sample TSVs (concatenated/sorted) into pooled per-site summary statistics across samples.
- Phase 2 adds per-strand base counts, base-composition entropy, GC fractions, and strand-bias metrics at every site.

# Build

Requires:
- GCC or Clang with C11 support
- htslib v1.18+ installed and discoverable via pkg-config

```bash
make
```

or manual:
```bash
gcc -O3 -std=c11 -o mq_forensics src/*.c \
  $(pkg-config --cflags --libs htslib) -lm
```

# Usage

## 1. Extract metrics from one BAM

```bash 
mq_forensics analyze -b <in.bam> -r <regions.bed> -C 50 -d 5 \
  -o per_site.tsv -O per_interval.tsv [options]
```

Required:

- `b` : input BAM (indexed .bai required)
- `r` : BED (0-based, half-open)
- `C` : cap threshold (samtools-style -C), e.g. 50
- `d` : minimum depth for site to emit numeric outputs

Optional:

`--emit-suffstats` : emit n,sum,sumsq triplets instead of direct means/medians (compact + required for summarize)
`--emit-hist`      : also emit histograms for MQ/effMQ/clipfrac
`-f/--fasta ref.fa`: FASTA reference (needed if using flank/ref-only)
`--flank N`        : include flanking context over ±N bp for per-site windows
`--ref-only`       : restrict flanking coverage/clip stats to reads matching the reference at that base

Example:

```bash
mq_forensics analyze -b sample.bam -r regions.bed -C 50 -d 5 \
  -o sample.site.tsv -O sample.iv.tsv --emit-suffstats --emit-hist \
  -f ref.fa --flank 5 --ref-only
```

Parallel many BAMs:

```bash
parallel -j 16 \
  'mq_forensics analyze -b {} -r regions.bed -C 50 -d 5 \
     -o {.}.site.tsv -O {.}.iv.tsv --emit-suffstats --emit-hist' \
  :::: bam.list
```

## 2. Summarize across samples
Concatenate per-site TSVs, sort, then summarize:

```bash
{ head -n1 *.site.tsv | head -n1
  tail -q -n +2 *.site.tsv \
  | LC_ALL=C sort -S 3G --parallel=1 -k1,1 -k2,2n
} \
| mq_forensics summarize -i - -o per_site_summary.tsv
```

- `-i` : sorted suffstats TSVs (chrom,pos key required)
- `-o` : summarized pooled TSV

# Outputs

**Per-site (direct mode)**

```bash
chrom  pos  depth  mismatch_bases  ins_len_sum  del_len_sum  clip_bases_sum
mq_mean  mq_median  mq_sd
capmq60_mean  capmq60_median  capmq60_sd
effmq_mean  effmq_median  effmq_sd
subQ_mean  subQ_median  subQ_sd
clipQ_mean  clipQ_median  clipQ_sd
[flank_cov_mean  flank_clipfrac_mean]   # if --flank enabled
nA_fwd  nC_fwd  nG_fwd  nT_fwd  nN_fwd  nA_rev  nC_rev  nG_rev  nT_rev  nN_rev  depth_fwd  depth_rev
entropy_pooled  alph_eff_pooled  entropy_fwd  alph_eff_fwd  entropy_rev  alph_eff_rev
gc_frac_pooled  gc_frac_fwd  gc_frac_rev
strand_bias_z
```

**Per-site (suffstats mode)**

```bash
chrom  pos  depth  mismatch_bases  ins_len_sum  del_len_sum  clip_bases_sum
n_mq  sum_mq  sumsq_mq
n_cap sum_cap sumsq_cap
n_eff sum_eff sumsq_eff
n_subQ sum_subQ sumsq_subQ
n_clipQ sum_clipQ sumsq_clipQ
n_clipfrac sum_clipfrac sumsq_clipfrac
n_capped sum_delta sumsq_delta
[n_flank_cov sum_flank_cov sumsq_flank_cov
 n_flank_cf  sum_flank_cf  sumsq_flank_cf]   # if --flank enabled
nA_fwd nC_fwd nG_fwd nT_fwd nN_fwd nA_rev nC_rev nG_rev nT_rev nN_rev
[hists …]   # if --emit-hist
```

**Per-interval**

```bash
chrom  start  end  aligned_bp_in_interval
```

**Summarized per-site (across samples)**

From suffstats input, `summarize` produces:
- depth-normalized global rates: mismatch/ins/del/clip
- pooled mean & SD for: MQ, capMQ, effMQ, SubQ, ClipQ, clipfrac
- capping load: frac_capped, delta mean/sd
- flanking context pooled stats: flank_cov_mean/sd, flank_cf_mean/sd (if present)
- per-strand base composition: entropy (bits and effective alphabet), GC fractions (pooled/fwd/rev), strand-bias z-score
- optional divergences & histogram-based medians if `--emit-hist`


# Requirements & assumptions

- BAM contains MD (and usually NM). If missing:

```bash
samtools calmd -bAr ref.fa in.bam > out.bam
```

- BAM must be indexed.
- BED is BED3+.



# Performance notes

- `analyze`: Single-pass per BAM; per-site accumulators reset each interval.
- `summarize`: Streaming reducer, memory ~ O(depth × #metrics). Sorting dominates runtime.
- Histograms make KS/W1/JS comparisons cheap in summarize.
- Flank stats add lightweight prefix-sum passes per interval.

