# Overview

Implements **Step 6: Analyze parameter sweep of the spruceGBS pipeline**. Selected results are also shown and discussed here.

---

# Contents

* **`analysis_functions.R`**: Parse, summarize, model, and plot results of the [individual heterozygosity](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#individual-heterozygosity), [PCA/RDA](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#rda-on-principal-components), and [site and population-level statistics](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#site-and-population-level-statistics) experiments.
* **`example_usage.R`**: Workflow used to create the plots summarized here. 

---

# `analysis_functions.R` Index

## I. Data Import and Utility

- **`import_summary(filepath)`**  
    Import parameter summary table, parse columns, and annotate with metadata.
    
- **`import_indvhet(filepath)`**  
    Import and summarize individual heterozygosity by parameter combination and population.
    
- **`read_emms(path)`**  
    Read EMMs tables from CSV, automatically detecting numeric vs character columns.

- **`write.unix(x, out.name, ...)`**  
    Write a data frame/tibble as a UNIX-style TSV (tab-separated, no quotes).
    
- **`get.info(sample.list, metadata)`**  
    Load a sample list and join with sample metadata.
    
- **`filter_by_maf(df, maf_col, maf_threshold)`**  
    Replace F/absF with NA for rows below a MAF threshold.
    
- **`%||% (a, b)`**  
    Null coalescing operator: returns `a` if not null, else `b`.

## II. Statistical Summaries

- **`summarize_stats(df, stat_cols, param_cols, mode)`**  
    Summarize statistics (mean, median, SD, IQR, etc.) by parameter sets for loci or populations.
    
## III. Model Fitting and Selection

- **`get_factor_vars(model)`**  
    Extract factor variables from a mixed-effects model object.
    
- **`model_selection(data, response, fixed_formula, random_formula, ...)`**  
    Fit maximal mixed models, perform model selection (MuMIn), and extract marginal means.
    
- **`write_models(result_list, base_name, outdir)`**  
    Write model results (formulas, variance components, marginal means) to CSV files.

- **`pcangsd_rda(domain, cov_dir_suffix , n_pc, output_dir)`**  
    Partition variance in principal components explained by region or library for all parameter combos.
  
### IV. Plotting Functions

- **`plot_stat_distributions(data, stat_name, value_col, ..., ncol)`**  
    Ridgeline density plots of statistics by parameter levels, multi-panel by parameter.
    
- **`plot_indvhet(df, fill_var, plot_type, output_dir, output_file)`**  
    Heatmaps of individual heterozygosity, faceted by parameter and population/grouping.
    
- **`plot_manova(df, fill_type, plot_type, output_dir, output_file)`**  
    Heatmaps of MANOVA variance explained by library/region/difference.
    
- **`plot_pcangsd_sweep(domains, var_list, group_var, ...)`**  
    Batch PCA biplots for each domain/parameter combination.
    
- **`plot_popvar(df, statistic, fill_var, output_dir, output_file, ...)`**  
    Heatmaps of within-region, among-population variation in pogen statistics.
    
- **`plot_pop_heatmaps(pop_summary, statistic, fill_var, output_file, output_dir, ...)`**  
    Heatmaps of statistic means per population, faceted by filter parameters.
    
- **`plot_RMSE(pop_summary, output_dir, output_file, ...)`**  
    RMSE heatmaps for F, by filter parameters.
    
- **`plot_marginal_means(df, target_stat, x_var, ..., output_file, output_dir)`**  
    Scatterplots of estimated marginal means for a statistic, with error bars and auto-extraction of grouping variables.

- **`plot_rda_heatmap(df, fill_type, plot_type, output_dir, output_file)`**
    Heatmaps of variance components from RDA by parameter setting
---

# Usage

See `example_usage.R` for the specific workflow used to create the plots summarized here, among many others.

---

# Inputs & Outputs

* **Inputs**:

  * `\<path/to/input1\>`: Description of the expected input file or directory.
  * `\<path/to/input2\>`: ...
* **Outputs**:

  * `\<path/to/output1\>`: Description of the generated output.
  * `\<path/to/output2\>`: ...

---

# Dependencies

List required modules, software, or packages:

* Bash (with `set -euo pipefail` recommended)
* [samtools](https://www.htslib.org/) >= 1.9
* R (packages: tidyverse, data.table, ...)
* ...

---

# Notes & Gotchas

* I am using R v. 4.2.1 and some tidyverse-related functions may be depreceated in their current releases.
