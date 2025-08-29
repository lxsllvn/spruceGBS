# beagleutils

Utilities for BEAGLE genotype-likelihood files and ANGSD depth/counts matrices.
One installable CLI with subcommands:

- `reheader`: add sample names (and snpcode for counts) to raw outputs.
- `subset`: subset matrices by SNP list and/or sample list, works for counts or BEAGLE.
- `split-genotype`: split a counts matrix into homRef / het / homAlt using the BEAGLE GLs (argmax assignment), with optional MAF>threshold filtered outputs.

Designed for HPC pipelines: supports gzipped I/O, consistent TSV formats, and explicit erroring on format mismatches.


## Install

```bash
# from the repo root
pip install -e .
```

This provides the console entrypoint:

```
geno-utils --help
```

## Quick start

``` bash
# 1) Reheader (counts): add snpcode from *.pos.gz and sample names from a text file

geno-utils reheader domain_filtered.counts.gz sample_names.txt reheadered.counts.gz --pos domain_filtered.pos.gz

# 2) Reheader (beagle): rewrite header with sample names (3 GL columns per sample)
geno-utils reheader domain_filtered.beagle.gz sample_names.txt reheadered.beagle.gz

# 3) Subset to a region/sample set
geno-utils subset reheadered.beagle.gz region_subset.beagle.gz --snps keep_snpcodes.txt --samples keep_samples.txt

# 4) Split by genotype (writes three files: .homRef.tsv.gz, .het.tsv.gz, .homAlt.tsv.gz)
geno-utils split-genotype reheadered.beagle.gz reheadered.counts.gz OUTPREFIX \
  --site-stats site_summary.tsv.gz \
  --maf-threshold 0.05
```

## Commands

### `reheader` 
Adds informative headers so downstream tools don’t rely on position or anonymous sample columns.

- Counts mode (requires `--pos`):
  - Output columns: `snpcode  <sample_1> <sample_2>` ...
  - snpcode is created as `chrom_pos` from the provided `.pos.gz`.

- BEAGLE mode (no `--pos`):
  - Output columns: `marker allele1 allele2 <sample_1> <sample_1> <sample_1> <sample_2> <sample_2> <sample_2> ...`
  - Assumes triplets per sample (AA, AB, BB) and rewrites the header to repeat each sample name 3× in order. Errors if the column count doesn’t match `3 + 3 * n_samples`.


#### Usage

```
geno-utils reheader <data.gz> <sample_names.txt> <output.gz> [--pos sites.pos.gz]
```

#### Potential problems

- Counts input missing `--pos`: you’ll get the usage error; counts mode needs `.pos.gz` to build `snpcode`.
- Sample count vs columns (BEAGLE): CLI raises a clear error if `3 + 3*n_samples != header_cols`.


### `subset`

Subsets either a reheadered counts matrix or a reheadered BEAGLE file by SNPs and/or samples.

- Auto-detects BEAGLE layout from header (marker allele1 allele2 ...), otherwise treats as counts (snpcode  <samples...>).
- Keeps meta columns (marker/alleles for BEAGLE, snpcode for counts) and filters sample columns accordingly.

#### Usage

```
geno-utils subset <input_file> <output_file> [--snps snps.txt] [--samples samples.txt]
```

#### Potential problems

- Input not reheadered: tool expects first column to be `snpcode` (counts) or `marker` (BEAGLE). It will exit with an explicit error if missing.
- Sample names must match header exactly (case-sensitive).


### `split-genotype`

Splits a reheadered counts matrix into three files according to the argmax genotype from the matching reheadered BEAGLE file:
- For each site & sample, take the max of the three GLs (AA/AB/BB).
- If ties (e.g., flat likelihoods), that sample contributes 0 to all outputs for that site. This will mostly affect genotypes with no read data.

#### Outputs

Each file has: `snpcode genotype <sample_1> <sample_2>` ... with non-target genotypes zeroed per row.

- `OUTPREFIX.homRef.tsv.gz`: read counts for homozygous reference genotypes
- `OUTPREFIX.het.tsv.gz`: read counts for heterozygous genotypes
- `OUTPREFIX.homAlt.tsv.gz`: read counts for homozygous alternate genotypes

Optional MAF filter: if `--site-stats` is provided, also writes `*.maf05.tsv.gz` (or whatever threshold you set), keeping only rows where `MAF > threshold`. The site stats file must contain columns `snpcode` and `MAF`.

#### Usage

```
geno-utils split-genotype <beagle.gz> <counts.gz> <out_prefix> \
  [--site-stats site_summary.tsv.gz] \
  [--maf-threshold 0.05]
```

#### Potential problems

Sample order must match between BEAGLE and counts—a mismatch raises a loud error that prints the detected orders. Reheader first, then subset both files with the same lists if needed.
SNP codes must match site-by-site between both inputs; the script compares `snpcode` per row and errors on any mismatch.


## Troubleshooting

“Header columns … expected 3 + 3*n_samples” (BEAGLE):
 - Your BEAGLE file doesn’t align with the provided sample list; verify the sample count and that beagle indeed stores 3 GL columns per sample in the expected order. Re-generate the sample list or reheader.

“Sample names/order mismatch between files!” (split-genotype):
 - Reheader both files from the same sample_names.txt first. If subsetting, apply the same sample list to both files to preserve alignment.

“Input file does not have ‘snpcode’ or ‘marker’ as first column!” (subset):
 - You fed an unreheadered file; run reheader first.
