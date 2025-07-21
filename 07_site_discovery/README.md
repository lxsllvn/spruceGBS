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

* **`split_reference.sh`**: divides reference genome into subsets, creating indexed FASTA files and corresponding ANGSD region/site files
* **`site_discovery.sh`**: runs ANGSD on each reference subset and finds sites passing minimal quality filters
* **`prepare_angsd_ref.sh`**: concatenates filtered position files and prepares final region, site, and indexed fasta files
* **`angsd_likelihoods.sh`**: does stuff
* **`summarize_site_stats.py`**: does stuff
* **`summarize_site_stats.sh`**: does stuff
* **`codeconvert`**: does stuff
* **`sample_call_rates.sh`**: does stuff
* **`reheader_genotype_matrix.py`**: does stuff
* **`subset_genotype_matrix.py`**: does stuff
* **`reheader_genotype_matrix.sh`**: does stuff
* **`subset_genotype_matrix.sh`**: does stuff
* **`library_call_thresh.R`**: does stuff
* **`batch_effects.sh`**: does stuff
* **`discovery_and_filtering.R`**: does stuff

---

# Create ANGSD reference assemblies

In [Step 2: Reduced reference genome preparation](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref), we identified regions in the *P. abies* assembly that are on scaffolds with mapped reads and outside of annotated repeats (+/- 500 bp), resulting in ~519 Mb across more than 100,000 scaffolds. ANGSD isn't designed to handle a dataset of this size, but we can further reduce the computational overhead by limiting analyses to sites passing some basic quality filters.

Creating the reduced references follows the same three-step process used in [the scaffold selection for the parameter sweep](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref):

1. `split_reference.sh`: Divides `picea_newref_target_regions.bed` into 23 subsets, extracts and indexes their FASTA records, and prepares their corresponding ANGSD site and region files. There’s nothing special about 23; based on trial runs before running out of memory, ANGSD could analyze ~5,000 scaffolds (+/- 15%) and 300 samples using ~36 Gb of memory, which is convenient for our cluster.

2. `site_discovery.sh`: Runs ANGSD using the quality filters identified in the parameter sweep (`-minQ 20 -minMapQ 50 -C 100 -baq 0`) and finds sites with <60% missing data. The missing data cutoff is somewhat arbitrary; the goal is simply to make the reference small enough for analysis while leaving some leeway to optimize sample vs. site-level missing data.

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

## `site_discovery.sh` usage

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

```R
source("discovery_and_filtering.R")

# Define group definitions using quoted logical expressions
group_defs <- list(
  posF = list(bad = "F > 0.5", good = "F < 0.5 & F > -0.5"),
  negF = list(bad = "F < -0.5", good = "F < 0.5 & F > -0.5"),
  Hobs = list(bad = "Hobs > 0.5", good = "Hobs < 0.5")
)

# Statistics to evaluate
stats <- c("total_depth", "mean_depth", "median_depth", "cv_depth", "rel_IQR_depth", "call_rate")

# Domain names to iterate over
domains <- c("northern_pt_aa", "northern_pt_ab", "northern_pt_ac", "southern", "siberia")

# List to collect summaries from all domains
all_summaries <- list()

for (dom in domains) {
  dat <- read.table(paste0(dom, "_site_summary_maf05.tsv"), header = TRUE)
  
  # Construct output filenames for each group_def
  output_file <- setNames(
    paste0(dom, "_", names(group_defs), "_site_filters_ROCs.png"),
    names(group_defs)
  )
  
  res <- do_auc(
    dat,
    group_defs, 
    stats,
    output_dir = "rocs",
    output_file = output_file
  )
  
  # Add domain info and store summary
  if (!is.null(res$summary) && nrow(res$summary) > 0) {
    res$summary$domain <- dom
    all_summaries[[dom]] <- res$summary
  }
}

# Combine all summaries into a single data.frame
summary_all_domains <- do.call(rbind, all_summaries)

# Write output to CSV (no row names, no quotes)
write.csv(
  summary_all_domains, 
  file = "site_filters_AUC_results.csv", 
  row.names = FALSE, 
  quote = FALSE
)
```

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

## `batch_effects.sh` usage

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


```R
source("discovery_and_filtering.R")
# PCA biplots
# Note: "sequenced_samples_metadata.csv" must be in the working directory

domains <- c("northern", "siberia", "southern")

# Generate PCA biplots for all domains and call thresholds
for (domain in domains) {
  base.names <- paste0(domain, ".ct", 5:8)
  for (base.name in base.names) {
    plot_pcas(
      cov_path = base.name,
      sample_list_path = paste0(domain, "_batch_samples.txt"),
      PCs = 1:10,
      fill_var = c("library", "region"),
      fill_pal = list("Set3", "Set1"),
      point_args = list(size = 1.5, alpha = 0.8),
      output_file = base.name
    )
  }
}

# Analyze pairwise PC dissimilarities vs. geographic distance and 
# binary library indicator ("same" or "different" library)
meta <- read.csv("sequenced_samples_metadata.csv")

all_domains_results <- list()

for (domain in domains) {
  base.names <- paste0(domain, ".ct", 5:8)
  res <- bind_rows(
    lapply(base.names, function(bn) {
      analyze_pca_distance_predictors(
        bn,
        paste0(domain, "_batch_samples.txt"),
        meta,
        n_pcs = 10
      )
    })
  ) %>% mutate(domain = domain)
  all_domains_results[[domain]] <- res
}

# Combine, sort, and format results
summary_table <- bind_rows(all_domains_results) %>%
  mutate(PC_num = as.integer(gsub("PC", "", PC))) %>%
  arrange(domain, PC_num, base.name, t_value) %>%
  select(-PC_num) %>%
  mutate(
    t_value = round(t_value, 2),
    p_value = formatC(p_value, format = "e", digits = 2)
  )

write.csv(
  summary_table,
  "PCA_batch_effects.csv",
  row.names = FALSE
)

# Identify SNPs contributing most to PCs dominated by batch effects
k_values <- c(siberia = 3, northern = 3, southern = 3)

for (domain in c("siberia", "southern", "northern")) {
  k <- k_values[[domain]]
  for (i in 5:8) {
    res <- pca_selection(
      basename = paste0(domain, ".ct", i),
      K = k,
      mode = "univariate",
      files = list(
        sites = paste0(domain, "_call_thresh0.", i),
        pca_sites = paste0(domain, ".ct", i, ".Pcangsd.sites"),
        pcadapt_zscores = paste0(domain, ".ct", i, ".Pcangsd.pcadapt.zscores"),
        selection = paste0(domain, ".ct", i, ".Pcangsd.selection")
      )
    )
    
    pc_name <- paste0("pcadapt", k, ".BH")
    sel_name <- paste0("selection", k, ".BH")
    
    write.unix(
      res$snpcode[res[[pc_name]] < 0.05 | res[[sel_name]] < 0.05],
      paste0(domain, "_ct", i, "k", k, "_Pcangsd_blacklist.txt")
    )
  }
}

```


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
