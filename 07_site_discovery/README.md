# Overview

Implements Step 7: site discovery and filtering of the spruceGBS pipeline.

Helpful summary goes here. 

There is a lot of code, but most of the steps/code are just to get ANGSD to play nice with my dataset, get the results into more user-friendly formats, and filter the DP matrices/beagle likelihoods by sample and site. 

interesting analytical choices are really only made in seleciton of site-level filters/thresholds, the batch effect analysis, and the population-stratified MAF selection 

---

# Contents

* [Create ANGSD reference assemblies](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* [Calculate genotype likelihoods](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* [Site and sample quality filters](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
  * [Site-level filters](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
  * [Sample filtering](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
  * [Apply filters](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* [Assess and remove batch effects](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* [Stratified MAF filters](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)

---

# Scripts

* **`split_reference.sh`**: does stuff
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

# Create ANGSD reference assemblies

In [Step 2: Reduced reference genome preparation](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref), we identified regions in the *P. abies* assembly that are on scaffolds with mapped reads and outside of annotated repeats (+/- 500 bp), resulting in ~519 Mb across more than 100,000 scaffolds. ANGSD isn't designed to handle a dataset of this size, but we can further reduce the computational overhead by limiting analyses to sites passing some basic quality filters.

Creating the reduced references follows the same three-step process used in [the scaffold selection for the parameter sweep](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref):

1. `split_reference.sh`: Divides `picea_newref_target_regions.bed` into 23 subsets, extracts and indexes their FASTA records, and prepares their corresponding ANGSD site and region files. There’s nothing special about 23; based on trial runs before running out of memory, ANGSD could analyze ~5,000 scaffolds (+/- 15%) and 300 samples using ~36 Gb of memory, which is convenient for our cluster.

2. `domain_site_discovery.sh`: Runs ANGSD using the quality filters identified in the parameter sweep (`-minQ 20 -minMapQ 50 -C 100 -baq 0`) and finds sites with <60% missing data. The missing data cutoff is somewhat arbitrary; the goal is simply to make the reference small enough for analysis while leaving some leeway to optimize sample vs. site-level missing data.

3. `prepare_angsd_ref.sh`: Merges the passing sites and produces a single indexed FASTA and ANGSD site and region file.


## `split_reference.sh` usage

```bash
#!/bin/bash
$SCRIPTS/07_site_discovery/split_reference.sh \
 <outdir> \
 <reference.bed> \
 <reference.fa>
```
**Inputs**
  * `<outdir>`: Output directory for results
  * `<reference.bed>`: BED file with target regions
  * `<reference.fa>`: Reference genome FASTA file

**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

## `domain_site_discovery.sh` usage

```bash
#!/bin/bash
$SCRIPTS/07_site_discovery/domain_site_discovery.sh \
 <index> \
 <domain> \
 <bamlist> \
 <path/to/reference/subsets> \
 <outdir>
```

**Inputs**
  * `<index>`: numerical index denoting the 01...23 reference subset
  * `<domain>`: name of genetic domain, e.g. southern, northern, siberia. 
  * `<bamlist>`: file containing the paths to bam files to include in the analysis
  * `<path/to/reference/subsets>`:  path to \*.fa, \*\_region and \*\_site files for the reference subsets
  * `<outdir>`: desired output directory; will create if it doesn't exist
    
**Outputs**
  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...


## `prepare_angsd_ref.sh` usage

Has --splits/--sites modes. 

```bash
#!/bin/bash
$SCRIPTS/07_site_discovery/prepare_angsd_ref.sh \
 --splits \
 <domain>
 <input>
 <ref>
 <outdir>
```

**Inputs**
  * `<domain>`:   sample domain (e.g., southern)
  * `input>`:    directory containing *.pos.gz lists (e.g., southern_pt_01.pos.gz, etc)
  * `<ref>`:  path to reference genome FASTA (indexed for samtools
  * `<outdir>`:   directory to save merged .bed, .sites, .regions, .fa files

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
