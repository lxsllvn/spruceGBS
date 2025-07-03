# Overview

Implements **Step 2: reduced reference preparation** of the spruceGBS pipeline. It uses read depths from the [full alignments](https://github.com/lxsllvn/spruceGBS/blob/main/01_read_alignment/README.md) and repeat annotations from the P.abies v. 1.0 [assembly](https://plantgenie.org/FTP) to:
1) identify scaffolds with mapped reads,
2) create a reduced reference comprising only those scaffolds,
3) make a buffered repeat mask and identify target regions for down-stream analysis, and
4) intersect the full alignments with the target regions.

---

# Contents

* **`scaffolds_with_coverage.sh`**: find scaffolds with ≥5 mapped reads in any sample
* **`reduced_reference_prep.sh`**: create new reference with only scaffolds with mapped reads
* **`picard_dictionary.sh`**: make sequence dictionary for GATK (used in [04_realignment](https://github.com/lxsllvn/spruceGBS/tree/main/04_realignment))
* **`find_targets.sh`**: find potential sites for downstream analyses; this removes regions near/within/book-ended by annotated repeats
* **`bam_intersection.sh`**: intersect [full alignments](https://github.com/lxsllvn/spruceGBS/tree/main/01_read_alignment) with the output of **`find_targets.sh`**

---

# Scaffold search

The *P. abies* reference genome is 12.4 Gb and has 10,253,694 scaffolds. This makes most downstream analyses extremely slow, have massive memory requirements, or just fail entirely (e.g. GATK, ANGSD). 

Fortunately, most scaffolds have zero mapped reads in any sample. **`scaffolds_with_coverage.sh`** takes the per-sample read depths produced in [01_read_alignment](https://github.com/lxsllvn/spruceGBS/tree/main/01_read_alignment) and identifies scaffolds with \>= 5 mapped reads (with Q >= 20 and MQ >= 30) in any sample. A single sample is permissive, but I tested requiring 100 samples and the resulting number of scaffolds were not really that different (218,545 vs. 162,766). 

## **`scaffolds_with_coverage.sh` usage**

```bash
#!/bin/bash
sbatch "$SCRIPTS/scaffolds_with_coverage.sh" <depth_dir> [scratch_dir] [output_file]
```

* `<depth_dir>`   (required): directory holding *.depth files
* `[scratch_dir]` (optional): temp dir for sorting (default: `${SPRUCE_PROJECT}/scaff_cov_tmp_$$`)
* `[output_file]` (optional): desired output name (default: `${SPRUCE_PROJECT}/ref/scaffolds_with_coverage.txt`)
  
---

# Reduced reference preparation

Next, we created a reduced reference comprising only the 218,545 scaffolds with mapped reads with **`reduced_reference_prep.sh`**. This script uses `samtools` to extract the scaffolds from the full reference and index the reduced reference, and creates a bed file from `picea_newref.fa.fai`. 

**`picard_dictionary.sh`** creates the sequence dictionary needed by GATK for [indel realignment](https://github.com/lxsllvn/spruceGBS/tree/main/04_realignment).

## **`reduced_reference_prep.sh` usage**

```bash
bash reduced_reference_prep.sh 
```

## **`picard_dictionary.sh` usage**

```bash
bash picard_dictionary.sh
```

---

# Identify target regions for analysis 

The *Picea* nuclear genome is ca. 70% transposable elements. Many of these are collapsed in the aging *P. abies* v. 1.0 reference assembly, but Illumina libraries are unlikely to offer much repeat resolution anyway. In a previous experiment, I found that sites in annotated repeats have a consistent deficit of heterozygotes, regardless of the stringency of the genotype calls/likelihoods. Because these sites are unlikely to survive to the final set of genotype calls/likelihoods, I decided to remove them at this stage to reduce the computational demand of the indel realignment. 

**`find_targets.sh`** takes a bed file of the annotated repeats, expands them by 500 bp on both sizes, subtracts these from `picea_newref.bed` and removes any short intervening regions (< 1000 bp) or scaffolds. This produces `picea_newref_target_regions.bed`, which are the initial set of sites for the subsequent steps. 

## **`find_targets.sh` usage**

```bash
bash find_targets.sh 
```

---

# BAM intersections

Finally, we intersected the [full alignments](https://github.com/lxsllvn/spruceGBS/blob/main/01_read_alignment/README.md) with `picea_newref_target_regions.bed` using the `-L` argument of `samtools view`. Reads mapping outside of these intervals are treated as unmapped and are not included in the output. These are the input for the next step, [indel realignment](https://github.com/lxsllvn/spruceGBS/tree/main/04_realignment).

## **`bam_intersection.sh` usage**

```bash
#!/bin/bash
while IFS= read -r sample; do
    sbatch "$SCRIPTS/bam_intersection.sh" "$sample"
done < sample.list
```

---
# Inputs & Outputs

**Inputs**:
* `${SPRUCE_PROJECT}/bams/full_alignments/`: alignments to the full P.abies reference, created [here](https://github.com/lxsllvn/spruceGBS/blob/main/01_read_alignment/)
* `${SPRUCE_PROJECT}/bams/full_alignments/read_depths`: read depths per sample, created during [alignment](https://github.com/lxsllvn/spruceGBS/blob/main/01_read_alignment/)
* `${SPRUCE_PROJECT}/ref/Pabies1.0-genome.fa`: the [reference genome](https://plantgenie.org/FTP)
* `${SPRUCE_PROJECT}/ref/spruce_repeats.bed`: annotated repeat bed file, converted from the gffs available [here](https://plantgenie.org/FTP)
    
**Outputs**:
* `${SPRUCE_PROJECT}/ref/picea_newref.fa`: reduced reference genome with fasta index (\*\.fai) and sequence dictionary (\*\.dict)
*  `${SPRUCE_PROJECT}/ref/picea_newref_target_regions.bed`: scaffolds with coverage, minus annotated repeats (+/- 500) bp and short (< 1,000 bp) intervening regions
* `${SPRUCE_PROJECT}/bams/intersected/`: alignments intersected by `picea_newref_target_regions.bed`
    
---

# Dependencies
* [samtools](https://www.htslib.org/) 1.19.2
* [bedtools](https://github.com/arq5x/bedtools2) 2.31.0
* [picard](https://github.com/broadinstitute/picard)  3.3.0
* Java 17.0.6

