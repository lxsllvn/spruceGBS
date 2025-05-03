## Overview

---
## Contents

* **`read_qc.sh`**: Short description of what this script does.
* **`bwa.sh**: Short description of what this analysis or step does.
---

## Pre-processing and quality control

Reads were demultiplexed using the process_radtags module of Stacks v.2.0. Adapter sequences and low-quality bases were removed with fastp, using **read_qc.sh**.

Most libraries were sequenced on the HiSeq 2500, but two (*ca*. 600 samples) were sequenced on the Illumina NovaSeq. This is a bit unfortunate because the two platforms differ substantially. Unlike earlier platforms, Novoseq uses a two-channel system, and thus low-quality bases potentially result in spurious G calls, and quality scores are binned into Q10, Q20, Q30, and Q40. These differences can bias down-stream genetic analyses (e.g. [Lou and Therkildsen, 2021](https://doi.org/10.1111/1755-0998.13559)).

To help mitigate the potential effects of the differing platforms, I trimmed poly-X tails (`-x`), low entropy sequences (``-y``) and more aggressively filtered low-quality bases using the `--cut_front` and `--cut_right` options in fastp, which trim trailing and leading bases if the mean quality in a 4 bp window drops below 20. I analyzed overrepresented k-mers (`-p`) and visually assessed the reduction in poly-X, in addition to sanity-checking the filtered HiSeq *vs*. Novoseq reads.

## **read_qc.sh** usage

``` bash
#!/bin/bash
while IFS= read -r sample; do
    sbatch "$SCRIPTS/read_qc.sh" "$sample"
done < sample.list
```
---
## Alignment
We used mapped surviving paired-end reads to the *P. abies* reference genome with BWA-MEM v. 0.7.19, as implemented in **bwa.sh**. This script produces sorted, indexed BAMs. This script also calculates per-sample read depths, which are used to identify scaffolds to keep in the [reduced reference genome](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref)

## **bwa.sh** usage
```bash
#!/bin/bash
while IFS= read -r sample; do
    sbatch "$SCRIPTS/bwa.sh" "$sample"
done < sample.list
```
---
## Inputs & Outputs

* **Inputs**:
  * `\<${SPRUCE_PROJECT}/ref/Pabies1.0-genome.fa\>`: reference genome and BWA index files
  * `\<${SPRUCE_PROJECT}/reads>`: de-multiplexed reads
    
* **Outputs**:
  * `\<${SPRUCE_PROJECT}/reads/qc\>`: quality-controlled reads
  * `\<${SPRUCE_PROJECT}/bams/full_alignments\>`: sorted and indexed BAMs
  * `\<${SPRUCE_PROJECT}/bams/full_alignments/read_depths\>`: read depths per sample
---

## Dependencies
* [stacks](https://catchenlab.life.illinois.edu/stacks/) v. 2.0
* [fastp](https://github.com/OpenGene/fastp) v. 0.23.4
* [bwa](https://github.com/lh3/bwa) v. 0.7.18
* [samtools](https://www.htslib.org/) v. 1.19.2

---


