# Overview

Implements Step 7: site discovery and filtering of the spruceGBS pipeline.

Helpful summary goes here. 

---

# Contents

* [Create ANGSD reference assembly](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* * [Split assembly](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
  * [Site discovery](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
  * [Create final ANGSD reference](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* [Calculate genotype likelihoods](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* [Site and sample quality filters](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* * [Site-level filters](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
  * [Sample filtering](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
  * [Apply filters](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* [Assess and remove batch effects](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* [Stratified MAF filters](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)

---

# Scripts

* **`make_reduced_reference.sh`**: does stuff
* **`domain_site_discovery.sh`**: does stuff
* **`prepare_angsd_ref.sh`**: does stuff
* **`angsd_likelihoods.sh`**: does stuff
* **`summarize_site_stats.py`**: does stuff
* **`summarize_site_stats.sh`**: does stuff
* **`auc_roc_filter.R`**: does stuff
* **`snp_stats_filter.awk`**: does stuff
* **`codeconvert`**: does stuff
* **`sample_call_rates.sh`**: does stuff
* **`reheader_genotype_matrix.py`**: does stuff
* **`subset_genotype_matrix.py`**: does stuff
* **`reheader_genotype_matrix.sh`**: does stuff
* **`subset_genotype_matrix.sh`**: does stuff
* **`library_call_thresh.R`**: does stuff
* **`library_call_thresholds.sh`**: does stuff
* **`pcangsd_batch_effects.R`**: does stuff



---

# Create ANGSD reference assembly

problem, solution, implementation. 

## Split assembly

### `make_reduced_reference.sh` usage

```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```

**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

## Site discovery

### `domain_site_discovery.sh` usage

```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

## Create final ANGSD reference

### `prepare_angsd_ref.sh` usage

```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

---

# Calculate genotype likelihoods

## `angsd_likelihoods.sh` usage

```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

---

# Site and sample quality filters

## Site-level filters
### `summarize_site_stats.sh` usage

```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...


### `auc_roc_filter.R` usage


```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

### `snp_stats_filter.awk` usage

```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

### `codeconvert` usage


```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

## Sample filtering

### `sample_call_rates.sh` usage


```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

## Apply filters 

### `reheader_genotype_matrix.sh` usage


```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

### `subset_genotype_matrix.sh` usage

```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

---

# Assess and remove batch effects

## `library_call_thresholds.sh` usage

```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...


## `pcangsd_batch_effects.R`

```bash
#!/bin/bash
sbatch solve_all_my_problems.sh

```
**Inputs**
  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...


# Dependencies

probably these 

* [samtools](https://www.htslib.org/) 1.19.2
* [angsd](https://github.com/ANGSD/angsd) 0.935
* [Pcangsd](https://github.com/Rosemeis/pcangsd) 1.36.4
* [Python 3.12.3](https://www.python.org/downloads/release/python-3123/) 
* [SciPy-bundle](https://docs.hpc2n.umu.se/software/libs/SciPy-bundle/) 2023.11 for pcangsd
* [R](https://www.r-project.org/) 4.2.1
  - dplyr 1.1.4
  - tidyr 1.3.0
  - data.table 1.14.8
  - tidyverse 2.0.0

Note! I am using an old version of R because I don't want to refactor my tidyverse-related code. Some functions may be deprecated in the current version.

---

# Notes & Gotchas

* Any special instructions, known issues, or tips.
