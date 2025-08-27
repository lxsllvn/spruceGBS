````markdown
# mq_forensics

Fast C/htslib tool to extract per-site (within a BED) and per-interval mapping statistics from BAMs.  Designed for GBS/amplicon panels (short loci), scalable to thousands of BAMs.

- Whole-read (WR) attribution for SubQ, ClipQ, capped MQ (`-C`) — each read contributes its single value to every site it overlaps.
- Local site features (mismatches, indel lengths at start position) derived from CIGAR + MD.
- Interval metric: `aligned_bp_in_interval` (area under coverage).

---

## Build

Requires:
- GCC or Clang with C11 support
- [htslib](https://github.com/samtools/htslib) installed and discoverable via `pkg-config`

Quick build:
```bash
make
```
Manual build (if you prefer no Makefile):

```bash
gcc -O3 -std=c11 mq_forensics.c -o mq_forensics \
  $(pkg-config --cflags --libs htslib) -lm
```

If `pkg-config` isn’t set up, specify include/lib paths:

```bash
gcc -O3 -std=c11 mq_forensics.c \
  -I/path/to/htslib/include -L/path/to/htslib/lib \
  -o mq_forensics -lhts -lpthread -lz -lbz2 -llzma -lm
```

## Usage

```
mq_forensics -b <in.bam> -r <regions.bed> -C <cap_threshold> -d <min_depth> \
  -o <per_site.tsv> -O <per_interval.tsv>
```

`-b` : input BAM (indexed `.bai` required)
- `-r` : BED (0-based, half-open)
- `-C` : cap threshold (samtools-style `-C`), e.g. `50`
- `-d` : **min depth**; if site depth `< d`, all numeric outputs for that site print as `NA` (default `0`)
- `-o` : per-site output TSV
- `-O` : per-interval output TSV

Example:

```bash
mq_forensics -b sample.bam -r regions.bed -C 50 -d 5 \
  -o sample.per_site.tsv \
  -O sample.per_interval.tsv
```
Many BAMs in parallel:

```bash
parallel -j 16 \
  'mq_forensics -b {} -r regions.bed -C 50 -d 5 -o {.}.site.tsv -O {.}.iv.tsv' \
  :::: bam.list
```

## Outputs

### Per-site TSV columns

```
chrom  pos  depth  mismatch_bases  ins_len_sum  del_len_sum  clip_bases_sum
mq_mean  mq_median  mq_sd
capmq60_mean  capmq60_median  capmq60_sd
effmq_mean  effmq_median  effmq_sd
subQ_mean  subQ_median  subQ_sd
clipQ_mean  clipQ_median  clipQ_sd
```

- `depth`: read count overlapping the site (reads with a deletion at the site are included — mpileup semantics).
- `mismatch_bases`: # of substitutions at the site (from MD).
- `ins_len_sum` / `del_len_sum`: summed lengths of insertions/deletions starting at this site.
- `clip_bases_sum`: sum of soft+hard clipped bases from all reads overlapping the site.
- `mq_*`: raw BAM-MEM MAPQ statistics per site (0–60, as reported by the aligner).
- `capmq60_*`: mapping quality caps computed by the -C formula (sam_cap_mapq), rescaled from [0..C] to the familiar [0..60] range. These are upper bounds; the cap can only lower mapping quality, never raise it.
- `effmq_*`: effective mapping qualities actually used downstream:
- `subQ_*`: sum of capped substitution base qualities per read (BQ≥13; each capped at 33), aggregated per site.
- `clipQ_*`: soft-clip quality sum + 13×hard-clip length per read, aggregated per site.

When `depth < min_depth` (or `depth == 0`), all numeric columns are `NA`.

### Per-interval TSV columns

```
chrom  start  end  aligned_bp_in_interval
```

- `aligned_bp_in_interval`: sum of M/= /X overlap with the interval (area under coverage).

## Requirements & assumptions

BAM contains `MD` (and usually `NM`) tags. If missing:

```bash
samtools calmd -bAr ref.fa in.bam > out.bam
```

- BAM must be indexed (`samtools index in.bam`).
- BED is BED3+ (extra columns ignored).


## Performance notes

- Single pass per BAM; per-site stats kept in memory only for current interval.
- Read vectors (for medians/SDs) scale with depth per site; tested with BAMs with average depths of ~30-80x.  For ultra-deep data, open an issue and we (realistically, ChatGPT) can switch medians to streaming/selection.


## Troubleshooting

- **Linker errors about `log10/sqrt`** → add `-lm` (Makefile already does).
- **`libhts.so` not found at runtime** → set `LD_LIBRARY_PATH` to your htslib lib dir or use a module.
- **`MD` missing** → run `samtools calmd` as above.
