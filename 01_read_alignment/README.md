## Overview

## Contents

* **`<script1>.sh`**: Short description of what this script does.
* **`<script2>.R`**: Short description of what this analysis or step does.
* **`...`**

---

## Pre-processing and quality control

Reads were demultiplexed using the process_radtags module of Stacks v.2.0. Adapter sequences and low-quality bases were removed with fastp, using **read_qc.sh**.

Most libraries were sequenced on the HiSeq 2500, but two (*ca*. 600 samples) were sequenced on the Illumina NovaSeq. This is a bit unfortunate because the two platforms differ substantially. Unlike earlier platforms, Novoseq uses a two-channel system, and thus low-quality bases potentially result in spurious G calls, and quality scores are binned into Q10, Q20, Q30, and Q40. These differences can bias down-stream genetic analyses (e.g. [Lou and Therkildsen, 2021](https://doi.org/10.1111/1755-0998.13559)).

To help mitigate the potential effects of the differing platforms, I trimmed poly-X tails (`-x`), low entropy sequences (``-y``) and more aggressively filtered low-quality bases using the `--cut_front` and `--cut_right` options in fastp, which trim trailing and leading bases if the mean quality in a 4 bp window drops below 20. I analyzed overrepresented k-mers (`-p`) and visually assessed the reduction in poly-X, in addition to sanity-checking the filtered HiSeq *vs*. Novoseq reads.

---


## Usage

Provide example commands demonstrating a typical invocation:

```bash
bash <script1>.sh arg1 arg2
Rscript <script2>.R input_file output_file
```

---

## Inputs & Outputs

* **Inputs**:

  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...
* **Outputs**:

  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

---

## Dependencies

List required modules, software, or packages:

* Bash (with `set -euo pipefail` recommended)
* [samtools](https://www.htslib.org/) >= 1.9
* R (packages: tidyverse, data.table, ...)
* ...

---

## Notes & Gotchas



