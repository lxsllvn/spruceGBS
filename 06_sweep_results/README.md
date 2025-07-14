# Overview

Implements **Step 6: Analyze parameter sweep** of the spruceGBS pipeline. Selected results are also shown and discussed here.

---

# Contents
* Summary
* Core parameter combination results
  * Individual heterozygosity
    * Mixed-effect models
    * Within-population variation
  * RDA and PCA
* Extended parameter results
  * RDA and PCA
  * Individual heterozygosity
    * Mixed-effect models
    * Within-population variation
  * Genetic diversity results
* Scripts 
     
---

---

# Summary
We first examined a coarse grid of parameter settings to identify the most promising set for more detailed investigation. The 'core' set of combinations comprised:
- **Base Quality (`-minQ`)**: 20, 30
- **Mapping Quality (`-minMapQ`)**: 20, 30, 40
- **BAQ Model (`-baq`)**: 0 (off), 1 (simple), 2 (extended)
- **Mapping Quality Capping Coefficient (`-C`)**: 0, 50
- **Per-Library Call Rate Filters**: 40%, 50%, and 60%

In summary, we found that combinations with `-baq 1` were highly sensitive to minimum base quality and showed inconsistent performance across tests and domains. Our dataset includes reads with both binned (NovaSeq) and continuous (HiSeq) quality scores, so shifting alignment uncertainty to this metric is risky. No recalibration (`-baq 0`) and the extended model (`-baq 2`) produced very similar or indistinguishable results. Uncapped (`-C 0`) and recalibrated (`-C 50`) mapping qualities also produced inconsistent results across tests and domains, with `-C 50` reducing library effects in individual heterozygosity for the Siberian but not the northern or southern populations. Capped qualities did not reduce the library effects inferred from RDA, but biplots revealed that individual variance within a NovaSeq library dominated the first two principal components for the northern domain. Minimum base qualities (`-minQ`) and mapping qualities (`-minMapQ`) did not strongly influence the results unless combined with `-baq 1` and `-C 50`, respectively.

Next, we focused on a wider range of `-C` and `-minMapQ` values and omitted BAQ recalibration (`-baq 0`). This 'extended' set of parameter combinations included:
- **Base Quality (`-minQ`)**: 20, 30
- **Mapping Quality (`-minMapQ`)**: 20, 30, 40, 50
- **BAQ Model (`-baq`)**: 0 (off)
- **Mapping Quality Capping Coefficient (`-C`)**: 0, 60, 75, 100
- **Per-Library Call Rate Filters**: 40%, 50%, and 60%

While no parameter combination emerged as a clear winner, `-C 0` and `60` produced inconsistent results depending on the domain and test. In particular, if one combination with `-C 60` led to the *weakest* library effect for a given domain and test, a different `-C 60` often produced the *strongest* effect. This pattern is attributable to the increased importance of the minimum mapping qualities setting, with the difference between `-minMapQ 40` and `50` already enough to lead to contrasting results. While `-C 60` seems to downgrade mapping qualities too harshly, the two more relaxed settings (`-C 75` and `100`) generally performed better than uncapped mapping qualities (`-C 0`). This is most striking in the northern domain PCA biplots, where `-C 0 -minMapQ 50` did not ameliorate the signal from the NovaSeq library, but even `-C 100 -minMapQ 20` reduced it substantially.

The differences between `-C 75` and `100` were small, and the choice between them ultimately comes down to personal preference. For this dataset, I chose `-C 100` and `-minMapQ 50` because this combination minimized spurious structure in the PCAs most effectively, and it does not appear too stringent for non-reference populations based on the marginal means of diversity statistics from the mixed effect models. However, the small quantitative and qualitative differences between `-C 75` and `100` over a range of minimum mapping and base qualities give me confidence that this region of parameter space is reasonable. I also selected the more permissive minimum base quality (`-minQ 20`), given the lack of difference with `-minQ 30`.

More detailed results are summarized below.

# Core parameter combinations

## Individual heterozygosity

### Mixed-effect models 

### Within-population variation

## RDA and PCA 

# Extended parameter combinations

## RDA and PCA

## Individual heterozygosity

### Mixed-effect models 

### Within-population variation

### Genetic diversity 

# Scripts

* **`analysis_functions.R`**: Parse, summarize, model, and plot results of the [individual heterozygosity](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#individual-heterozygosity), [PCA/RDA](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#rda-on-principal-components), and [site and population-level statistics](https://github.com/lxsllvn/spruceGBS/tree/main/05_angsd_param_sweep#site-and-population-level-statistics) experiments.
* **`example_usage.R`**: Workflow used to create the plots summarized here. 

## `analysis_functions.R` Index

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
