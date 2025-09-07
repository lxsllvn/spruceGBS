# Overview

Implements Step 7: site discovery and filtering of the spruceGBS pipeline. Helpful summary goes here. 

---

# Contents

* [Create ANGSD reference assemblies](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#create-angsd-reference-assemblies)
* [Calculate genotype likelihoods](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#create-angsd-reference-assemblies)
* [Site filtering](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#site-filtering)
  * [Interpretable machine learning](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#interpretable-machine-learning)
  * [Feature engineering](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#feature-engineering)
* [Assess and remove batch effects](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)
* [Stratified MAF filters](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section)

---

# Scripts

* **`split_reference.sh`**: divides reference genome into subsets, creating indexed FASTA files and corresponding ANGSD region/site files
* **`site_discovery.sh`**: runs ANGSD on each reference subset and finds sites passing minimal quality filters
* **`prepare_angsd_ref.sh`**: concatenates filtered position files and prepares final region, site, and indexed fasta files
* **`do_domain_discovery.sh`**: driver to launch `site_discovery.sh` and `prepare_angsd_ref.sh` jobs simultanously 
* **`angsd_likelihoods.sh`**: runs ANGSD on the prepared reference and outputs genotype likelihoods and site-level summary statistics
* **`summarize_site_stats.py`**: collects SNP diagnostic metrics from ANGSD outputs and site-level depth summaries
* **`summarize_site_stats.sh`**: helper for submitting `summarize_site_stats.py` jobs
* **`reheader_genotype_matrix.sh`**: helper for [`beagle-utils`](https://github.com/lxsllvn/spruceGBS/tree/main/beagle-utils) CLI to add informative headers to ANGSD `*.counts.gz` and `*.beagle.gz` files
* **`read_counts_by_genotype.sh`**: helper for [`beagle-utils`](https://github.com/lxsllvn/spruceGBS/tree/main/beagle-utils) CLI to split `*.counts.gz` into homRef/het/ homAlt using by most likely genotype from `*.beagle.gz`
* **`subset_genotype_matrix.sh`**:  helper for [`beagle-utils`](https://github.com/lxsllvn/spruceGBS/tree/main/beagle-utils) CLI to select samples and/or sites from `*.counts.gz` and `*.beagle.gz`
* **`do_xgboost.R`**: finalizes features, prepares test/train sets, runs Bayesian parameter optimization, and trains XGBoost model
* **`do_xgboost.sh`**: helper script to submit `do_xgboost.R` 
* **`GBDT_path_analysis.R`**: functions for gradient boosted decision tree (GBDT) ensemble path analyses
* **`discovery_and_filtering.R`**: does stuff
* **`codeconvert`**: does stuff
* **`sample_call_rates.sh`**: does stuff
* **`library_call_thresh.R`**: does stuff
* **`batch_effects.sh`**: does stuff

---

# Create ANGSD reference assemblies

In [Step 2: Reduced reference genome preparation](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref), we identified scaffolds with reads and, from these, designated target regions that fall outside of annotated repeats. These target regions comprise ~519 Mb across 100,000 scaffolds, a substantial reduction from the 12.4 Gb and 10,253,694 scaffolds of the original assembly, but still far too large for ANGSD to manage. 

Here, we repeat the methods used to create an ANGSD-friendly reference for the [the parameter sweep](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref), this time using all target regions and all samples that passed [initial quality control](https://github.com/lxsllvn/spruceGBS/tree/main/03_initial_qc). Our goal is simply to remove sites that fail some minimal call rate and quality filters, with some leeway for further parameter optimization. 

For the southern and Siberian domains, all samples (~300 each) were analyzed simultaneously to produce a single ANGSD reference per domain. To reduce the memory requirements for the northern domain, we divided the 775 samples in three groups and processed each separately. The resulting northern_pt_aa, northern_pt_ab, and northern_pt_ac references were analyzed seperately until after the site filtering stage. 

```bash
#!/bin/bash
shuf northern_bamlist.txt > tmp_shuf.txt
split -n l/3 --additional-suffix .txt tmp_shuf.txt northern_bamlist_pt_
rm tmp_shuf.txt
```

Creating ANGSD-ready references is a three-step process:

1. `split_reference.sh`: Divides `picea_newref_target_regions.bed` into 23 subsets, extracts and indexes their FASTA records, and prepares their corresponding ANGSD site and region files. This is is done once. 

There’s nothing special about 23; based on trial runs before running out of memory, ANGSD could analyze ~5,000 scaffolds (+/- 15%) and 300 samples using ~36 Gb of memory, which is convenient for our cluster. 

Then, site discovery on the reference subsets followed by the final reference preparation are done for each domain seperately: 

2. `site_discovery.sh`: Runs ANGSD using the quality filters identified in the parameter sweep (`-minQ 20 -minMapQ 50 -C 100 -baq 0`) and finds sites with <60% missing data. The missing data cutoff is somewhat arbitrary; the goal is simply to make the reference small enough for analysis while leaving some leeway to optimize sample vs. site-level missing data.

3. `prepare_angsd_ref.sh`: Merges the passing sites and produces a single indexed fasta and ANGSD site and region file.

You can submit both steps together using the driver `do_domain_discovery.sh`, which calls the SLURM helper [`submit.sh`](https://github.com/lxsllvn/spruceGBS/blob/main/submit.sh) to reduce job babysitting. The driver launches the 23 discovery jobs and then schedules the collation step with `--dependency=afterok`, so collation runs only if all splits succeed.
 
## `split_reference.sh` usage

```bash
#!/bin/bash
$SCRIPTS/07_site_discovery/split_reference.sh \
    "$SPRUCE_PROJECT/ref/subsets" \
    "$SPRUCE_PROJECT/ref/picea_newref_target_regions.bed" \
    "$SPRUCE_PROJECT/ref/picea_newref.fa"
```

**Inputs**
- `$1` - output directory for results
- `$2` - [target region BED file](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref)
- `$3` - [indexed reference fasta](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref)

**Outputs**

For each chunk:
- `target_scaffs_pt_{01..23}.bed` - BED coordinates
- `target_scaffs_pt_{01..23}.fa` - fasta records
- `target_scaffs_pt_{01..23}.fa.fai` - fasta index
- `target_scaffs_pt_{01..23}_regions` - ANGSD region file
- `target_scaffs_pt_{01..23}_sites` -  ANGSD sites file
- `target_scaffs_pt_{01..23}_sites.bin` - ANGSD sites index
- `target_scaffs_pt_{01..23}_sites.idx` - ANGSD sites index

## `do_domain_discovery.sh` usage

```bash
#!/bin/bash
$SCRIPTS/07_site_discovery/do_domain_discovery.sh \
"southern" \
"${SPRUCE_PROJECT}/site_discovery/southern_bamlist.txt" \
"${SPRUCE_PROJECT}/ref/subsets" \
"${SPRUCE_PROJECT}/ref/southern" \
"${SPRUCE_PROJECT}/ref/picea_newref.fa"
```

**Inputs**

- `$1` - domain or sample group (e.g. `southern`, `northern_pt_aa`)
- `$2` - path to the BAM list file (list of input BAMs)
- `$3` - path to directiory with the `*.fa`, `*.fai`, `*_region` and `*_site` files created by `split_reference.sh`
- `$4` - output directory; will create if it doesn't exist
- `$5` - [indexed reference fasta](https://github.com/lxsllvn/spruceGBS/tree/main/02_reduced_ref)
    
**Outputs**

An ANGSD-ready reference assembly:
- `${DOMAIN}_ref.bed`: merged and sorted BED 
- `${DOMAIN}_ref_sites`: ANGSD sites file
- `${DOMAIN}_ref_sites.bin`: ANGSD sites binary index 
- `${DOMAIN}_ref_sites.idx`: ANGSD sites binary index 
- `${DOMAIN}_ref_regions`: ANGSD region file
- `${DOMAIN}_ref.fa`: fasta records for unique scaffolds 
- `${DOMAIN}_ref.fa.fai`: fasta index file

Intermediate results: 
- `${DOMAIN}/${DOMAIN}_pt_{01..23}.counts.gz` - read count matrices for each chunk
- `${DOMAIN}/${DOMAIN}_pt_{01..23}.pos.gz` - scaffold and position names for each chunk

---

# Calculate genotype likelihoods

At this point, we have identified sites in each domain  that have > 3 uniquely-mapped, properly-paired reads with `-minQ 20 -C 100 -minMapQ 50` in ≥ 40% of the samples. Now, we can actually use `ANGSD` to produce genotype likelihoods in beagle format (`-doGlf 2`), along with allele frequencies (`-doMaf 1`), allele-based SNP quality metrics (`-dosnpstat 1`), and F<sub>is</sub> per site based on genotype likelihoods (`-doHWE 1`). As in the parameter sweep, we designated the major allele based on the reference state (`-doMaf 4`).

## `angsd_likelihoods.sh` usage

```bash
$SCRIPTS/07_site_discovery/angsd_likelihoods.sh \
 ${SPRUCE_PROJECT}/ref/siberia/siberia_ref.fa \
 ${SPRUCE_PROJECT}/ref/siberia/siberia_ref_regions \
 ${SPRUCE_PROJECT}/ref/siberia/siberia_ref_sites \
 siberia_bamlist.txt \
 siberia \
 $SPRUCE_PROJECT/site_discovery/siberia
```

**Inputs**
- `$1` - [domain reference genome](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#prepare_angsd_refsh-usage)
- `$2` - [domain ANGSD region file](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#prepare_angsd_refsh-usage)
- `$3` - [domain  ANGSD sites file](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#prepare_angsd_refsh-usage)
- `$4` - list of BAM files in the domain
- `$5` - base name for output files
- `$6` - save directory (optional); current working directory if unspecified

**Outputs**
- `*.beagle.gz` - genotype likelihoods, beagle format
- `*.counts.gz` - read count matrix
- `*.hwe.gz` - per-site inbreeding coefficients 
- `*.mafs.gz` - minor allele frequencies 
- `*.pos.gz` - scaffold names and positions
- `*.snpStat.gz` - base quality bias, mapping quality bias, strand bias, and edge bias between reference and alternate allele

---

# Site filtering

The filters evaluated during our parameter sweep primarily target base and mapping qualities of individual reads. These filters prevent the lowest-quality data from informing genotype likelihoods: 

```math
\ell(g_L \mid d_L) \;\propto\; P(d_L \mid g_L, L),
```

where $d_L$ is the observed sequencing data and $g_l$ is a candidate genotype at locus $L$.  While read-level filters should improve genotype likelihood estimates, they cannot guarantee that $L$ is a biologically real locus, or that all genotypes $g_l$ produce equally observable sequencing data.

Real datasets violate these assumptions. In large and repetitive genomes, reads mapped to the same coordinates may originate from multiple genomic regions due to collapsed paralogs or repeats longer than the read length. Conversely, mutations in restriction enzyme recognition sites can reduce or eliminate the detection of certain alleles, resulting in undetected homozygotes and undercalled heterozygotes.

To mitigate such artefacts, most variant calling pipelines apply additional site-level filters. These typically involve simple, univariate thresholds on summary statistics—such as mean mapping quality, depth, strand bias, or the results of rank-sum tests comparing the distribution of these metrics between reference and alternate alleles. While convenient, these filters are often borrowed from other studies, with little validation in the specific context of the dataset at hand. The filtering metrics reported by tools like `samtools` or `GATK` were primarily designed for human datasets and may not capture the error profiles associated with different genome architectures, references, or sequencing protocols.

Even in their original context, these filters often show low discriminatory power when used in isolation. Their continued use reflects practical constraints rather than demonstrated efficacy. In particular, designing new filters tailored to a specific dataset is non-trivial. Without access to gold-standard variant calls, new heuristics inherit many of the same limitations as existing ones—along with additional uncertainty. Filters that appear theoretically sound may behave unpredictably in practice. For example, a threshold on call rate could unintentionally enrich for collapsed paralogs if coverage artefacts co-vary with local mappability. Interactions among filters may further complicate interpretation, and ad hoc thresholds risk introducing systematic, yet unrecognized, biases.

## Interpretable machine learning 

Gradient boosting decision trees (GBDTs) provide a flexible, data-driven framework for modeling site-level quality in sequencing data. These models are particularly well-suited to this context: they capture nonlinear relationships, tolerate noisy or correlated predictors, and learn higher-order interactions among features. By iteratively fitting shallow decision trees, GBDTs can transform a collection of weakly informative metrics into a stronger composite model of site behavior.

Unlike manual filtering, GBDTs automatically optimize the contribution of each feature, accounting not only for its individual predictive power but also for how it interacts with other variables. They are robust to multicollinearity, naturally down-weighting redundant or uninformative inputs without requiring explicit feature selection. Critically, GBDTs remain interpretable. Tools such as feature importance scores and SHAP (SHapley Additive exPlanations) values allow us to quantify which features most strongly influence model outputs and how they interact—making the logic of the model transparent rather than opaque. This combination of flexibility and interpretability makes GBDTs well-suited for diagnosing how site-level properties affect downstream genetic inference, even in the absence of ground truth labels.

Machine learning applications are commonly prediction-oriented and rely on labeled training data, for example, a list sites known to be either true variants or artefacts. In our case, no such gold-standard labels exist, and we cannot directly estimate the probability that a given site reflects a real biological polymorphism. Instead, we use GBDTs as exploratory diagnostic tools, modeling how site-level quality metrics relate to downstream population genetic statistics such as observed heterozygosity (H<sub>o</sub>), expected heterozygosity (H<sub>e</sub>​), and inbreeding coefficients (F<sub>is</sub>​). This reframes the task from prediction to inference: rather than classifying sites, the goal is to understand how variation in sequencing quality influences the estimates that underpin biological interpretation.

This approach captures nonlinear, multivariate relationships between quality control features and genetic summary statistics. It can reveal interactions among features that shape diversity estimates, and help identify vulnerable regions of parameter space where inference may be distorted by systematic measurement error. In particular, it highlights cases where genetic statistics are strongly influenced by technical factors, providing a means to diagnose potential artefacts in downstream analyses.

While this method does not produce definitive classifications of site quality, it enables the identification of conditions under which diversity estimates are especially sensitive to the measurement process. When combined with domain knowledge, these insights can inform the design of more effective, interpretable, and biologically grounded filtering strategies.

## Feature engineering 

Feature engineering refers to the process of creating and transforming variables from raw data, analogous to selecting predictor variables for a statistical model.  While some ML jargon is just re-packaged statistics, the expectations from a well-designed set of 'features' are quite different from a predictor matrix, both philosophically and practically.

In regression, predictors are typically chosen for their _a priori_ mechanistic relevance and transformed to satisfy model assumptions about error structure, linearity, or collinearity. By contrast, GBDTs can flexibly approximate nonlinearities and interactions, and issues such as correlated variables, weak predictors, or differences in scale do not undermine the validity of the model itself. The focus shifts instead to extracting as much complementary information as possible about the processes that could generate the observed data.

This shift introduces new risks: features can capture patterns that are highly predictive but irrelevant to the underlying biology and leak information about the target variable. For instance, mapping-quality distributions differ systematically between reference and alternate alleles because mismatches against the reference are penalized during alignment. Features based on MQs must therefore be designed carefully to avoid dependence on true genetic diversity. Indirect leakage can arise when sample-size–sensitive statistics such as skewness or kurtosis are calculated by genotype class, allowing the model to “predict” population-level parameters (e.g. F<sub>is</sub>​) by inferring genotype frequencies, a discovery that is technically correct but wholly unhelpful.


To get started, `ANGSD` provides several commonly used diagnostic metrics based on allele states, including:
- **`baseQ_Z`**: Wilcoxon Mann–Whitney test comparing base qualities of major vs. minor alleles.
- **`mapQ_Z`**: analogous statistic for mapping qualities.
- **`edge_Z`**: Wilcoxon test comparing distances from read edges between alleles.
- **`SB1`, `SB2`, `SB3`**: three measures of strand bias.

At present, these are collected by `summarize_site_stats.py`, which also calculates some site-level depth summaries from the `*.counts.gz` matrix.  This setup is awkward because we calculate more (and partially overlapping) statistics from `*.counts.gz` later on. I hope to streamline this later on. 

### `summarize_site_stats.sh` usage

```bash
#!/bin/bash
$SCRIPTS/07_site_discovery/summarize_site_stats.sh \
$SPRUCE_PROJECT/site_discovery/siberia \
siberia 
```

**Inputs**
* `$1` - path to folder containing the `*.pos.gz`, `*.counts.gz`, `*.hwe.gz`, and `*.snpStat.gz` files
* `$2` - output basename for site summary tables
* `$3` - basename of the `*.pos.gz`, `*.counts.gz`, etc; optional if folder only contains one set of input files. 

**Outputs**
  * `*_site_summary.tsv` - summary table with columns `snpcode`, `total_depth`, `mean_depth`, `median_depth`, `call_rate`, `cv_depth`, `rel_IQR_depth`, `MAF`, `Hexp`, `Hobs`, `F` , `SB1` ,`SB2`, `SB3`, `HWE_LRT`, `HWE_pval`, `baseQ_Z`, `baseQ_pval`, `mapQ_Z`, `mapQ_pval`, `edge_z`, `edge_pval`
  * `*_site_summary_maf05.tsv` - as above, but only for MAF > 0.05 sites

### Depth-based features

Beyond these built-in metrics, we also examined read depth patterns. Unusually high, low, or variable coverage can signal problematic sites, especially when depth differs systematically by genotype. For example, null alleles can result in reduced read depth for apparent homozygotes compared to heterozygotes at the same locus ([McCouch et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4734769/)).


#### Prerequisites

Before we can define more sophisticated depth features, we need add informative column names to `*.counts.gz` and `*.beagle.gz`. For `*.counts.gz`, we also need a column with a `chrom_pos` code that matches the `marker` column  the corresponding `*.beagle.gz`. 

This can be done with the `reheader` command from [`beagle-utils`](https://github.com/lxsllvn/spruceGBS/tree/main/beagle-utils), a small Python tool for manipulating the `*.counts.gz`/`*.beagle.gz` files from ANGSD,  using the `reheader_genotype_matrix.sh` helper script. 

##### `reheader_genotype_matrix.sh` usage

Both `*.counts.gz` and `*.beagle.gz` require a text file of sample names (one per line). The number of lines must equal the number of columns in `*.counts.gz`. For a `*.beagle.gz` file, the keeps the first three columns untouched (`marker` `allele1` `allele2`) and then repeats each sample name in triplicate. 

Replacing the `*.beagle.gz` header only requires `sample_list.txt`:

```bash
$SCRIPTS/07_site_discovery/reheader_genotype_matrix.sh \
$SPRUCE_PROJECT/site_discovery/siberia/siberia.beagle.gz \
$SPRUCE_PROJECT/site_discovery/siberia/siberia.reheader.beagle.gz \
$SPRUCE_PROJECT/site_discovery/siberia_sample_list.txt
```

**Inputs**
- `$1` - path to `*.beagle.gz`
- `$2` - output path
- `$3` - text file with sample names (one per line) 

**Outputs**
- a `*.beagle.gz` file with the same number of columns as the original, with informative sample names.

Re-headering `*.counts.gz` additionally requires the path to the corresponding `*.pos.gz`file generated by `ANGSD`:

```bash
$SCRIPTS/07_site_discovery/reheader_genotype_matrix.sh \
$SPRUCE_PROJECT/site_discovery/siberia/siberia.beagle.gz \
$SPRUCE_PROJECT/site_discovery/siberia/siberia.reheader.beagle.gz \
$SPRUCE_PROJECT/site_discovery/siberia_sample_list.txt \
$SPRUCE_PROJECT/site_discovery/northern/domain_filtered/northern_domain_filtered.pos.gz

```
**Inputs**
- `$1` - path to `*.counts.gz`
- `$2` - output path
- `$3` - text file with sample names (one per line) 
- `$4` - path to corresponding `*.pos.gz`

**Outputs**
- a `*.counts.gz` file with an additional `snpcode` column (`chrom_pos` from `*.poz.gz`) and with informative sample names.

##### `read_counts_by_genotype.sh`

ANGSD produces genotype likelihoods, not genotype calls. To summarize read depths by genotype, we assigned the reads indexed by $counts[i,j]$ to the genotype with the highest likelihood for sample $j$. If all three genotypes had the same likelihood (as with missing data), all three were assigned zero depth.

This can be done with the `split-genotype` command in [`beagle-utils`](https://github.com/lxsllvn/spruceGBS/tree/main/beagle-utils). This requires `*.beagle.gz` and `*.counts.gz` to have matching sample names and snpcodes/markers columns (i.e., the output of `reheader`)

```bash 
$SCRIPTS/07_site_discovery/split_by_genotype.sh \
$SPRUCE_PROJECT/site_discovery/southern/southern.beagle.reheader.gz \
$SPRUCE_PROJECT/site_discovery/southern/southern.counts.reheader.gz \
southern \
$SPRUCE_PROJECT/site_discovery/southern/southern_site_summary.tsv.gz
```

**Inputs**
- `$1` - path to `*.reheader.beagle.gz`
- `$2` - path to `*.reheader.counts.gz`
- `$3` - output file prefix
- `$4` - (optional) path to `*_site_summary.tsv.gz` (from [`summarize_site_stats.sh`](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#summarize_site_statssh-usage)) with `snpcode` and `MAF` columns; if provided, creates MAF > 0.05 subsets of the count matrices. 

**Outputs**
- `OUT_PREFIX.homRef.tsv.gz`, `OUT_PREFIX.het.tsv.gz` and `OUT_PREFIX.homAlt.tsv.gz` - read counts assigned to genotypes
- `OUT_PREFIX.homRef.maf05.tsv.gz`, `OUT_PREFIX.het.maf05.tsv.gz` and `OUT_PREFIX.homAlt.maf05.tsv.gz` if `*_site_summary.tsv.gz` was provided. 

#### Transformations 
Because such patterns are often noisy and confounded by differences in sequencing effort among individuals and libraries, we adopted a feature engineering approach. Specifically, we derived a set of non-collinear variables that capture multiple aspects of read depth variation, using several normalization schemes to quantify the “unusualness” of a site’s depth distribution:

- **No normalization**
- **Locus-standardized** (depth relative to other loci within an individual)
- **Sample-standardized** (depth relative to other individuals at the same locus)
- **Double standardizations** (locus–sample and sample–locus scaling).

For each normalization, we applied the summaries both to total read counts per site and to genotype-specific read counts per site. We then calculated a comprehensive panel of distributional statistics, including: mean, median, standard deviation, median absolute deviation (MAD), coefficient of variation (CV), interquartile range (IQR), kurtosis, skewness, Shannon entropy, proportion of extreme values, and dip test for unimodality.


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
