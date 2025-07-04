# Overview

Implements Step 7: site discovery and filtering of the spruceGBS pipeline.

Helpful summary goes here. 

---

# Contents

* [Create ANGSD reference assembly](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
  * [Split assembly](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
  * [Site discovery](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
  * [Create final ANGSD reference](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* [Calculate genotype likelihoods](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* [Site and sample quality filters](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
  * [Site-level filters](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
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

Preparing the reference is again a three-step procedure. 
First make the reduced reference subsets as previously, but this time with all target regions. Then, do site discovery to find sites that meet the filtering criteria in each domain. Next, collect all the sites over the subsets and create a single merged reference for each domain. 

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

With the ANGSD references created, we can now calculate genotype likelihoods for each domain, and output a bunch of useful things.

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

We summarize the site stats and use AUC-ROC to find filter parameters that can distinguish quite dubious from less dubious SNPs and their optimal thresholds. 

## Site-level filters

site-level summary data from angsd is written to four different, potentially very large files. python script collect the results, calculates a few depth-related statistics, and saves the results to two summary files (all sites/maf > 0.05 sites)

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

can any of the quality/depth metrics reliably distinguish egregiously bad from 'at least not awful' sites?

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

outputs a list of snpcodes passing the filters. not sure that a big awk script is really needed here, but this script accepts column names and correctly returns sites meeting all conditions and ignores -nan/nan/inf/-inf 

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

by default, the beagle/counts.gz files do not have useful column names and .counts.gz doesn't have marker names, either. reheader_genotype_matrix.py adds this data to the files, which then allows subset_genotype_matrix.py to select specified samples and/or sites. 

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

finds sites meeting minimum library call rate thresholds, creates beagle subsets, and then runs pcangsd

note! need pcangsd 1.34.6; the pcadapt-zscores in at least some previous versions are in the wrong order; i.e. column 1 doesn't have the z-scores for PC1. 

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
