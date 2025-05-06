## Overview

Briefly describe the purpose of this directory (one or two sentences).

---

## Contents

* **`scaffolds_with_coverage.sh`**: Short description of what this script does.
* **`remove_repeats.sh`**: Short description of what this analysis or step does.
* **`reduced_reference_prep.sh`**:
* **`picard_dictionary.sh`**:
* **`bam_intersection`**:
---

## Scaffold search
Why? Because the P. abies reference genome is 12.4 Gb and has 10,253,694 scaffolds.  Downstream analysis are either extremely slow or can't run at all (ANGSD, GATK).

Find scaffolds with with \>= 5 mapped reads in any sample **scaffolds_with_coverage.sh**. I tested requiring 100 samples, but the resulting number of scaffolds were not really that different (218,545 vs. 162,766). 

## **`scaffolds_with_coverage.sh` useage**

```bash
sbatch "$SCRIPTS/scaffolds_with_coverage.sh" <depth_dir> [scratch_dir] [output_file]
```
* `<depth_dir>` (required): directory holding your *.depth files
* `[scratch_dir]` (optional): temp dir for sorting (default: `${SPRUCE_PROJECT}/scaff_cov_tmp_$$`)
* `[output_file]` (optional): desired output name (default: `${SPRUCE_PROJECT}/ref/scaffolds_with_coverage.txt`)

---
## Reduce reference preperation
We created a reduced reference genome comprising only these scaffolds using **reduced_reference_prep.sh**.

use samtools faidx to extract these scaffolds
parse fasta index to bed file
make sequence dictionary for later 

## **`reduced_reference_prep.sh` useage**
```bash
bash <script1>.sh arg1 arg2
```
## **`picard_dictionary.sh` useage**
```bash
bash <script1>.sh arg1 arg2
```

---

## Identify target regions for analysis 
scaffolds in the reduced ref, minus annotated repeats (+/- 500 bp) and short intervening regions (< 1000 bp)


## **`remove_repeats.sh` useage**
```bash
bash <script1>.sh arg1 arg2
```

---
## BAM intersections
intersect bams with target regions

## **`bam_intersection.sh` useage**

```bash
#!/bin/bash
while IFS= read -r sample; do
    sbatch "$SCRIPTS/bam_intersection.sh" "$sample"
done < sample.list
```
---


## Inputs & Outputs

**Inputs**:
  * `${SPRUCE_PROJECT}/bams/full_alignments/`: alignments to the full P.abies reference, created [here](https://github.com/lxsllvn/spruceGBS/blob/main/01_read_alignment/).
  * `${SPRUCE_PROJECT}/bams/full_alignments/read_depths`: read depths per sample, created during [alignment](https://github.com/lxsllvn/spruceGBS/blob/main/01_read_alignment/).
  * `${SPRUCE_PROJECT}/ref/Pabies1.0-genome.fa`: the [reference genome](https://plantgenie.org/FTP)
  * `${SPRUCE_PROJECT}/ref/spruce_repeats.bed`: annotated repeat bed file, converted from the gffs available [here](https://plantgenie.org/FTP)
    
  **Outputs**:
  * `${SPRUCE_PROJECT}/ref/picea_newref.fa`: reduced reference genome with fasta index (\*\.fai) and sequence dictionary (\*\.dict).
  * `${SPRUCE_PROJECT}/ref/picea_newref_target_regions.bed`: scaffolds with coverage, minus annotated repeats (+/- 500) bp and short (< 1,000 bp) intervening regions
  * `${SPRUCE_PROJECT}/bams/intersected/`: alignments intersected by picea_newref_target_regions.bed

---

## Dependencies

List required modules, software, or packages:

* [samtools](https://www.htslib.org/) v. 1.19.2
* [bedtools](https://github.com/arq5x/bedtools2) v. 2.31.0
* [picard](https://github.com/broadinstitute/picard) v. 3.3.0
* Java v. 17.0.6


---

## Notes & Gotchas

* Any special instructions, known issues, or tips.
* For example: ensure you run this after the reduced reference is built.
* ...
