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
Why? Because the P. abies reference genome is 12.4 Gb and has 10,253,694 scaffolds. Downstream analysis are either extremely slow or can't run at all (ANGSD, GATK).

Find scaffolds with with \>= 5 mapped reads in any sample **scaffolds_with_coverage.sh**. I tested requiring 100 samples, but the resulting number of scaffolds were not really that different (218,545 vs. 162,766). 

## **`scaffolds_with_coverage.sh` useage**

---
## Reduce reference preperation


## **`reduced_reference_prep.sh` useage**
## **`picard_dictionary.sh` useage**


---

## Identify target regions for analysis 

## **`remove_repeats.sh` useage**

---
## BAM intersections

## **`bam_intersection.sh` useage**
---

## Usage

Provide example commands demonstrating a typical invocation:

```bash
bash <script1>.sh arg1 arg2
Rscript <script2>.R input_file output_file
```

---

## Inputs & Outputs

**Inputs**:
  *`${SPRUCE_PROJECT}/bams/full_alignments/`: alignments to the full P.abies reference, created [here](https://github.com/lxsllvn/spruceGBS/blob/main/01_read_alignment/).
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
