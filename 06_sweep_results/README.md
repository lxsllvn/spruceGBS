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
- **Base Alignment Quality (BAQ) Model (`-baq`)**: 0 (off), 1 (simple), 2 (extended)
- **Mapping Quality Capping Coefficient (`-C`)**: 0, 50
- **Per-Library Call Rate Filters**: 40%, 50%, and 60%

We found that combinations with `-baq 1` were highly sensitive to minimum base quality and showed inconsistent performance across tests and domains. Our dataset includes reads with both binned (NovaSeq) and continuous (HiSeq) quality scores, so shifting alignment uncertainty to this metric is risky. No recalibration (`-baq 0`) and the extended model (`-baq 2`) produced very similar or indistinguishable results. Uncapped (`-C 0`) and recalibrated (`-C 50`) mapping qualities also produced inconsistent results across tests and domains, with `-C 50` reducing library effects in individual heterozygosity for the Siberian but not the northern or southern populations. Capped qualities did not reduce the library effects inferred from RDA, but biplots revealed that individual variance within a NovaSeq library dominated the first two principal components for the northern domain. Minimum base qualities (`-minQ`) and mapping qualities (`-minMapQ`) did not strongly influence the results unless combined with `-baq 1` and `-C 50`, respectively.

Next, we focused on a wider range of `-C` and `-minMapQ` values and omitted BAQ recalibration (`-baq 0`). This 'extended' set of parameter combinations included:
- **Base Quality (`-minQ`)**: 20, 30
- **Mapping Quality (`-minMapQ`)**: 20, 30, 40, 50
- **BAQ Model (`-baq`)**: 0 (off)
- **Mapping Quality Capping Coefficient (`-C`)**: 0, 60, 75, 100
- **Per-Library Call Rate Filters**: 40%, 50%, and 60%

While no parameter combination emerged as a clear winner, `-C 0` and `60` produced inconsistent results depending on the domain and test. In particular, if one combination with `-C 60` led to the *weakest* library effect for a given domain and test, a different `-C 60` often produced the *strongest* effect. This pattern is attributable to the increased importance of the minimum mapping qualities setting, with the difference between `-minMapQ 40` and `50` already enough to lead to contrasting results. While `-C 60` seems to downgrade mapping qualities too harshly, the two more relaxed settings (`-C 75` and `100`) generally performed better than uncapped mapping qualities (`-C 0`). This is most striking in the northern domain PCA biplots, where `-C 0 -minMapQ 50` did not ameliorate the signal from the NovaSeq library, but even `-C 100 -minMapQ 20` reduced it substantially.

The differences between `-C 75` and `100` were small, and the choice between them ultimately comes down to personal preference. For this dataset, I chose `-C 100` and `-minMapQ 50` because this combination minimized spurious structure in the PCAs most effectively, and it does not appear too stringent for non-reference populations based on the marginal means of diversity statistics from the mixed effect models. However, the small quantitative and qualitative differences between `-C 75` and `100` over a range of minimum mapping and base qualities give me confidence that this region of parameter space is reasonable. I also selected the more permissive minimum base quality (`-minQ 20`), given the lack of difference with `-minQ 30`.

More detailed results are summarized below.

Note to self: should summarize the analyses applied to the indvidual heterozygosity, PCA, and per-locus diversity metrics. Random effects models (indv_het ~ 1 + (1|library) and exploratory figures for individual heterozygosity, RDA and qualitative intepretation of PCA bipolots, model selection on maximal mixed effect models (y ~ ... + (1 | snpcode) + (1 | pop_code)) for diversity metrics plus an unholy number of exploratory plots. 

# Core parameter combinations

Differences among parameter combinations were dominated by the choice of BAQ model. Across all three domains,  `-baq 0` or `2` resulted in very similar or indistinguishable results, but `-baq 1` tended to produce the best and worst results for a given domain and evaluation test. Over all, the effects of `-baq 1` were overall highly idiosyncratic and varied by test, domain, and minimum base quality (`-minQ 20` and `30`).

Uncapped (`-C 0`) minimum mapping qualities had no clear effect on the evaluation results. When capped mapping qualities (`-C 50`) were invoked, minimum mapping qualities had a more noticeable influence.  However, whether `-C 50` reduced or exacerbated batch effects, and whether it should be paired with higher or lower MQ cutoffs, varied by domain and evaluation test. 

While support for `-C 50` was weak, the `-C 0` PCA biplots showed high variation in northern domain samples from one NovaSeq library. Variation among these samples dominated the second principal component and was not ameliorated by higher minimum MQ cutoffs. Unfortunately, this pattern was not detected by the dd-RDA. Subjectively, I still found the dd-RDA summaries useful, but this approach is not a substitute for visual examination of the PCA biplots, at least not in the current implementation.

## Individual heterozygosity

### Mixed-effect models 

In the northern (Fig. 1) domain, parameter combinations invoking `-baq 1` resulted in both the highest and lowest library effect on individual heterozygosity, depending on the minimum base quality setting.  Out of all combinations for this domain, the strongest library effects were found with the `-baq 1 -minQ 20` combination, whereas the weakest were found with `-baq 1 -minQ 30`. Other parameters had no clear effect: results were similar over the `-C`,  `-minMapQ` and `-ct` settings (Figs. 1). While the `-minQ` setting was very influential in combination with `-baq 1`, both values produced similar effects with `-baq 0` and `-baq 2`.

<img width="3000" height="1800" alt="northern_core_varcomp" src="https://github.com/user-attachments/assets/37c8caf7-a14b-4d07-b5ae-63e7169f7938" />

Figure 1.  Northern domain: proportion of variation in individual heterozygosity explained by library membership. Parameter combinations are faceted by base alignment quality model (`-baq`) in rows and the combined setting for the mapping quality capping coefficient (`-C`) and minimum per-library call threshold (`ct`) in columns. Within each facet, the minimum base quality (`-minQ`) is shown on the x-axis and minimum mapping quality (`-minMapQ`) on the y-axis. For example, the top left box shows results over varying `-minQ` and `-minMapQ` settings, with the remaining parameters fixed (`-baq 0 -C 0 -ct 0.40`). 

For Siberia, some `-baq 1` combinations led to the strongest library effects, as observed in the northern domain, but none were among the weakest (Fig. 2). In contrast to the northern domain (Fig. 1), `-baq 1 -minQ 30` resulted in stronger library effects than `-baq 1 -minQ 20`. Also unlike the northern domain, the choice of `-C` had a relatively large effect on the results (Fig. 2). Invoking `-C 50` resulted in less variation attributable to library membership and also made `-minMapQ` more impactful, with higher values resulting in a smaller library effects (Fig. 2).

<img width="3000" height="1800" alt="siberia_core_varcomp" src="https://github.com/user-attachments/assets/05a509d7-5bb9-4878-90a7-24922bfe18e9" />

Figure 2. Siberia domain: proportion of variation in individual heterozygosity explained by library membership. 

For this particular test, the overall response of the southern domain was very similar to the north so I will not give them more attention here. 

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
