## Overview

Briefly describe the purpose of this directory (one or two sentences).

---

## Contents

* **`<script1>.sh`**: Short description of what this script does.
* **`<script2>.R`**: Short description of what this analysis or step does.
* **`...`**

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

* Any special instructions, known issues, or tips.
* For example: ensure you run this after the reduced reference is built.
* ...
