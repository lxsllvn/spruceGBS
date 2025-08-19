# Overview

Implements **Step 6: Analyze parameter sweep** of the spruceGBS pipeline. Selected results are also shown and discussed here.

---

# Contents
* [Summary](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#summary)
* [Core parameter combination results](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#core-parameter-combinations)
  * [Individual heterozygosity](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#individual-heterozygosity)
    * [Mixed-effect models](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#mixed-effect-models)
    * [Within-population variation](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#within-population-variation)
  * [RDA and PCA](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#rda-and-pca)
* [Extended parameter results](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#extended-parameter-combinations)
  * [RDA and PCA](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#rda-and-pca-1)
  * [Individual heterozygosity](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#individual-heterozygosity-1)
    * [Mixed-effect models](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#mixed-effect-models-1)
    * [Within-population variation](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#within-population-variation-1)
  * [Genetic diversity results](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#genetic-diversity)
* [Scripts](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#scripts)
  * [Function index](https://github.com/lxsllvn/spruceGBS/tree/main/06_sweep_results#analysis_functionsr-index)
     
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

Differences among parameter combinations were dominated by the choice of BAQ model. Across all three domains, `-baq 0` or `2` resulted in very similar or indistinguishable results, but `-baq 1` tended to produce the best and worst results for a given domain and evaluation test. Over all, the effects of `-baq 1` were overall highly idiosyncratic and varied by test, domain, and minimum base quality (`-minQ 20` and `30`).

Uncapped (`-C 0`) minimum mapping qualities had no clear effect on the evaluation results. When capped mapping qualities (`-C 50`) were invoked, minimum mapping qualities had a more noticeable influence. However, whether `-C 50` reduced or exacerbated batch effects, and whether it should be paired with higher or lower MQ cutoffs, varied by domain and evaluation test. 

While support for `-C 50` was weak, the `-C 0` PCA biplots showed high variation in northern domain samples from one NovaSeq library. Variation among these samples dominated the second principal component and was not ameliorated by higher minimum MQ cutoffs. Unfortunately, this pattern was not detected by the dd-RDA. Subjectively, I still found the dd-RDA summaries useful, but this approach is not a substitute for visual examination of the PCA biplots, at least not in the current implementation.

## Individual heterozygosity

### Mixed-effect models 

In the northern (Fig. 1) domain, parameter combinations invoking `-baq 1` resulted in both the highest and lowest library effect on individual heterozygosity, depending on the minimum base quality setting. Out of all combinations for this domain, the strongest library effects were found with the `-baq 1 -minQ 20` combination, whereas the weakest were found with `-baq 1 -minQ 30`. Other parameters had no clear effect: results were similar over the `-C`, `-minMapQ` and `-ct` settings (Fig. 1). While the `-minQ` setting was very influential in combination with `-baq 1`, both values produced similar effects with `-baq 0` and `-baq 2`.

<img width="3000" height="1800" alt="northern_core_varcomp" src="https://github.com/user-attachments/assets/37c8caf7-a14b-4d07-b5ae-63e7169f7938" />

Figure 1. Northern domain: proportion of variation in individual heterozygosity explained by library membership. Parameter combinations are faceted by base alignment quality model (`-baq`) in rows and the combined setting for the mapping quality capping coefficient (`-C`) and minimum per-library call threshold (`ct`) in columns. Within each facet, the minimum base quality (`-minQ`) is shown on the x-axis and minimum mapping quality (`-minMapQ`) on the y-axis. For example, the top left box shows results over varying `-minQ` and `-minMapQ` settings, with the remaining parameters fixed (`-baq 0 -C 0 -ct 0.40`). 

For Siberia, some `-baq 1` combinations led to the strongest library effects, as observed in the northern domain, but none were among the weakest (Fig. 2). In contrast to the northern domain (Fig. 1), `-baq 1 -minQ 30` resulted in stronger library effects than `-baq 1 -minQ 20`. Also unlike the northern domain, the choice of `-C` had a relatively large effect on the results (Fig. 2). Invoking `-C 50` resulted in less variation attributable to library membership and also made `-minMapQ` more impactful, with higher values resulting in a smaller library effects (Fig. 2).

<img width="3000" height="1800" alt="siberia_core_varcomp" src="https://github.com/user-attachments/assets/05a509d7-5bb9-4878-90a7-24922bfe18e9" />

Figure 2. Siberia domain: proportion of variation in individual heterozygosity explained by library membership. See the Fig. 1 caption for more details.

For this particular test, the overall response of the southern domain was very similar to the north so I will not give them more attention here. 

### Within-population variation

Here, I focus on the summary statistics by parameter combination, which show the mean or median of the population coefficients of variation (CV) and relative interquartile ranges (IQR). In this case (*cf*. dd-RDA), I found that these summary plots supported the same conclusions as the more complex population plots (using`plot_pop_heatmaps()` in `analysis_functions.R`), so I omit these for simplicity. 

In the northern domain, the `-C` parameter had a strong influence on within-population variation in individual heterozygosity (Fig. 3). Whether quantified as the mean (Fig. 3) or median (not shown) CV over the five populations, `-C 50` resulted more heterogeneity among individuals, particularly paired with `-baq 0` and `-baq 2`. Combinations with `-baq 1 -C 0` resulted in higher mean CVs than either `-baq 0` or `2`, although this increase was apparently not attributable to library membership (*cf*. Figs. 1 & 3). With `-C 50`, differences among BAQ models were less pronounced, but `-minMapQ` had a larger effect in combination with `-baq 0` or `2`. 

<img width="3000" height="1800" alt="northern_core_mean_cv_het" src="https://github.com/user-attachments/assets/9e21f6ca-b7cc-4c24-9d30-5ee6750a6632" />

Figure 3. Northern domain: mean coefficient of variation (CV) of within-population individual heterozygosity. Parameter combinations are faceted by base alignment quality model (`-baq`) in rows and the combined setting for the mapping quality capping coefficient (`-C`) and minimum per-library call threshold (`ct`) in columns. Within each facet, the minimum base quality (`-minQ`) is shown on the x-axis and minimum mapping quality (`-minMapQ`) on the y-axis. For example, the top left box shows results over varying `-minQ` and `-minMapQ` settings, with the remaining parameters fixed (`-baq 0 -C 0 -ct 0.40`). 

Median CVs showed similar patterns, but with an even more difference between results obtained with `-baq 1` versus either `-baq 0` or `2`. The higher heterogeneity generally found with `-C 50`, however, seems largely due to a few individuals. Mean within-population relative IQRs (IQR/median) showed no clear difference between `-C 0` and `-C 50`, except when combined with `-minMapQ 40` (Fig. 4). Because IQRs are less sensitive to within-population extremes, this suggests that heterozygosity estimates for some individuals are highly sensitive to `-C`, while most are not.

More anecdotally, I identified the individuals responsible for this this pattern but could not find any unifying characteristics. The majority were from populations with relatively large sample sizes (ca. 20), which were randomized between two HiSeq libraries. These sensitive individuals came from both libraries and did not differ in sequencing depth or call rate from the other individuals. In addition, the starting DNA for these individuals (and the rest of their populations) was isolated from cotyledons and should be the highest quality out of the entire dataset. 

<img width="3000" height="1800" alt="northern_core_mean_rel_iqr_het" src="https://github.com/user-attachments/assets/1bb6059b-5d5e-4824-a5c4-bcd4aa6a4d50" />

Figure 4. Northern domain: mean relative IQR of within-population individual heterozygosity. See the Fig. 3 caption for more details.

In the Siberian domain, the BAQ model had the strongest influence on mean within-population CVs (Fig. 5). Combinations with `-baq 1` yielded in the highest and lowest mean CV, depending on the minimum base quality threshold. Median CVs, along with mean and median standardized IQRs, responded similarly, but with `-C 50` more clearly resulting in less heterogeneity than `-C 0`. 

<img width="3000" height="1800" alt="siberia_core_mean_cv_het" src="https://github.com/user-attachments/assets/3d042abf-2c7b-40ca-ac54-d3cee4b2f56a" />

Figure 5. Siberian domain: mean coefficient of variation (CV) of within-population individual heterozygosity. See the Fig. 3 caption for more details.

Minimum mapping qualities had a large effect on mean CVs in the southern domain (Fig. 6), but their effect was modulated by the `-C` setting. Mean CVs were generally lower with `-C 0` but were the lowest with `-C 50` and `-minMapQ 20`. However, median CVs tended to be lower with `-C 50`, with the lowest values again with `-minMapQ 20` or `30`. 

<img width="3000" height="1800" alt="southern_core_mean_cv_het" src="https://github.com/user-attachments/assets/c00dddd8-e67b-439c-9d47-b6f4ac65ef12" />

Figure 6. Southern domain: mean coefficient of variation (CV) of within-population individual heterozygosity. See the Fig. 3 caption for more details.

## RDA and PCA 

In all three domains, parameter combinations with `-baq 1 -minQ 20` resulted in the smallest proportion of variation attributable to geography. (Figs. 7-9). Some RDAs on combinations with `-baq 1` had relatively large proportion of constrained variance explained by geography, particularly with higher call rate thresholds for the Siberian and northern domains (Figs. 7-8).

<img width="3000" height="1800" alt="northern_core_region_prop" src="https://github.com/user-attachments/assets/963e16c5-42df-4b51-a344-3ef3992bf5e5" />

Figure 7. Northern domain: proportion of constrained variation uniquely explained by region in dd-RDA. Parameter combinations are faceted by base alignment quality model (`-baq`) in rows and the combined setting for the mapping quality capping coefficient (`-C`) and minimum per-library call threshold (`ct`) in columns. Within each facet, the minimum base quality (`-minQ`) is shown on the x-axis and minimum mapping quality (`-minMapQ`) on the y-axis. For example, the top left box shows results over varying `-minQ` and `-minMapQ` settings, with the remaining parameters fixed (`-baq 0 -C 0 -ct 0.40`). 

<img width="3000" height="1800" alt="siberia_core_region_prop" src="https://github.com/user-attachments/assets/b1bb78a4-9458-493d-aaae-c85b41de1269" />

Figure 8. Siberia domain: proportion of constrained variation uniquely explained by region in dd-RDA. See the Fig. 7 caption for more details.

As in the individual heterozygosity analysis, the effect of `-C` differed among domains. While results for either setting produced similar results in the Siberian and southern domains (Figs. 8-9), combinations with `-C 0` led to higher proportion of geographic variation in the northern domain (Fig. 7), particularly combined with `-baq 0`. This particular example is also the only case where the choice between `-baq 0` and `2` seems impactful. 

<img width="3000" height="1800" alt="southern_core_region_prop" src="https://github.com/user-attachments/assets/cffe6102-8366-4c2d-93f4-ae9827e723d2" />

Figure 9. Southern domain: proportion of constrained variation uniquely explained by region in dd-RDA. See the Fig. 7 caption for more details.

For the southern and Siberian domains, the dd-RDA summaries were consistent with my visual interpretation of their PCA score plots. However, for the northern domain, the PCA plots show an important pattern that the dd-RDA summaries were not sensitive enough to detect. While Fig. 7 indicates that parameter combinations with `-C 50` recovered less geographic structure, the score plots show that variation within a NovaSeq library (NO19) dominates the second principal component and that this pattern was not affected by the minimum mapping quality (Fig. 10). Note that these plots also show results with `-minMapQ 50`, which was included in the second parameter sweep, described below. A similar pattern was not evident in another NovaSeq library (ALEX1; Fig. 10), but this may due to the relatively few samples from this library that were included in the parameter sweep analyses. 

<img width="2700" height="1500" alt="northern_baq0_C0_pcangsd_library ct5" src="https://github.com/user-attachments/assets/0db79ee3-df78-4479-b447-3d4b080ab2ae" />

Figure 10. Northern domain: PCA score plots for `-baq 0 -C 0 -ct 0.5` parameter combinations. Plots differ in their `-minQ` and `-minMapQ` values. Samples are colored by library membership. 

The NO19 samples originated from two of the four categorical geographic regions, and both show the same spread along PC2 (Fig. 11). This is probably why my dd-RDA summaries failed to indicate an issue with the `-C 0` combinations. [In subsequent filtering steps](https://github.com/lxsllvn/spruceGBS/tree/main/07_site_discovery#section), I compared pairwise geographic differences to a binary same/different library membership matrix instead. That approach may have worked better in the parameter sweep as well, but I have not tested it yet. 

<img width="2700" height="1500" alt="northern_baq0_C0_pcangsd_region ct5" src="https://github.com/user-attachments/assets/14dda2a6-ce08-49d6-887f-d68bb62b99ab" />

Figure 11. Northern domain: PCA score plots for `-baq 0 -C 0 -ct 0.5` parameter combinations. Plots differ in their `-minQ` and `-minMapQ` values. Samples are colored by geographic region. 

For comparison, the library effects indicated by dd-RDA for the `-baq 1` combinations are consistent with visual inspection of the score plots (Figs. 12-13), which are colored by library and geographic region, respectively. This also illustrates why some quantitative summary of the score plots is desirable: the entire parameter sweep generated an 868 plots that required visual inspection and interpretation. Fortunately, the PCA results were not sensitive to every parameter, which simplified the process somewhat. In the future, I will try to find a more reliable, scalable summary approach that can be confidently applied to other datasets. 

<img width="2700" height="1500" alt="northern_baq1_C0_pcangsd_library ct5" src="https://github.com/user-attachments/assets/e58f2832-dd58-448a-9ce5-b64fab21c4fd" />

Figure 12. Northern domain: PCA score plots for `-baq 1 -C 0 -ct 0.5` parameter combinations. Plots differ in their `-minQ` and `-minMapQ` values. Samples are colored by library membership. 

<img width="2700" height="1500" alt="northern_baq1_C0_pcangsd_region ct5" src="https://github.com/user-attachments/assets/59a74374-f756-4ccd-a217-5dd7447535cc" />

Figure 13. Northern domain: PCA score plots for `-baq 1 -C 0 -ct 0.5` parameter combinations. Plots differ in their `-minQ` and `-minMapQ` values. Samples are colored by geographic region. 

# Extended parameter combinations

Based on the results of the initial parameter sweep, we excluded `-baq 1` and `-baq 2` from further consideration. Both capped mapping quality coefficients produced concerning results in some scenarios, so we focused on exploring a wide range of `-C` values in the second sweep: 

- **Base Quality (`-minQ`)**: 20, 30
- **Mapping Quality (`-minMapQ`)**: 20, 30, 40, 50
- **BAQ Model (`-baq`)**: 0 (off)
- **Mapping Quality Capping Coefficient (`-C`)**: 0, 60, 75, 100
- **Per-Library Call Rate Filters**: 40%, 50%, and 60%
 
Minimum base quality and per-library call rate cutoffs performed similarly in the initial sweep, so we neither expanded our search over their space nor eliminated any value from further consideration. The initial parameter sweep did not include `-minMapQ 50`, so we included it here to evaluate the `-C 0 -minMapQ 50` combination.

When `-C 60` was invoked, results were highly sensitive to the minimum mapping quality. Overall, `-C 60 -minMapQ 50` and `-C 60 -minMapQ40` tended to produce contrasting results: if `-minMapQ 50` resulted in strong batch effects, then `-minMapQ 40` showed weak evidence of batch effects, and vice-versa. Similarly to the pattern observed with `-baq 1` in the initial sweep, the performance of `-C 60` was idiosyncratic among evaluation tests and domains.

In the southern and Siberian domains, similar variation in individual heterozygosity was observed with `-C 75` and `-C 100`. Higher `-minMapQ` tended to produce less variation attributable to library membership, but the level of within-population variation overall was less sensitive. In the northern domain, however, higher `-minMapQ` apparently exacerbated library effects, and within-population variation was lowest with `-C 0`. 

In all three domains, `-C 100` resulted in the greatest proportion of variance attributable to geography in the dd-RDA analysis, with minimum mapping qualities ≥ 30 giving the best results. In the southern and Siberian domains, combinations with `-C 75` and `-C 0` were comparable. In the northern domain, increasing `-minMapQ` to `50` did not resolve the high individual variation observed from the PCA score plots of the northern domain (Figs. 8-9), but even the relatively lenient combination of `-C 100 -minMapQ 20` removed it entirely. 

We quantified the effect of `-C`, `-minQ`, and `-minMapQ` on the marginal means of H<sub>e</sub>, H<sub>o</sub>, F, π, Θ<sub>w</sub>, Tajima's D, and MAF using mixed effect models. Consistent with the results from the dd-RDA and individual heterozygosity analyses, we found that `-C` and `-minMapQ` were the most influential parameters overall, with their interaction term included in 13 of the 21 best-performing models (three domains X 7 statistics). Minimum base quality (`-minQ`) was not included in any model.

The largest differences among parameter combinations were observed in the northern domain, although this may be a result of the larger sample size. Overall, `-C 60 -minMapQ 50` resulted in the highest Fs and lowest H<sub>e</sub>, H<sub>o</sub>, π, and Tajima's D in all models. Lower minimum mapping qualities paired with `-C 60` also tended to result in distinctly different 95% CIs. As expected, `-C 0` and `-C 100` produced the most similar estimates and results from `-C 75` tended to be more similar to these two more lenient settings than `-C 60`.

## RDA and PCA

In the extended parameter sweep, the strongest geographic signal in the dd-RDA was consistently recovered using `-C 100`, particularly in combination with higher minimum mapping qualities (`-minMapQ 40` to `50`) and more stringent call thresholds (60%). This trend was clearest in the northern domain (Fig. 15, where increasing both `-minMapQ` and `-ct` improved the proportion of constrained variance attributable to geography. 

In the northern domain, increasing `-minMapQ` to `50` did not resolve the high individual variation observed from the PCA plots of the northern domain (Figs. 8-9), but even the relatively lenient combination of `-C 100 -minMapQ 20` removed it entirely. 

<img width="3000" height="1800" alt="northern_extended_region_prop" src="https://github.com/user-attachments/assets/b64b7855-8ef9-4b5a-a8f2-01af7128c2c9" />

Figure 14. Northern domain: proportion of constrained variation uniquely explained by region in dd-RDA. Parameter combinations are faceted the mapping quality capping coefficient (`-C`) in rows and minimum per-library call threshold (`ct`) in columns. Within each facet, the minimum base quality (`-minQ`) is shown on the x-axis and minimum mapping quality (`-minMapQ`) on the y-axis. For example, the top left box shows results over varying `-minQ` and `-minMapQ` settings with `-C 0` and `-ct 0.40`. 

In the southern domain (Fig. 15), geographic signal was again maximized with `-C 100`, followed by `-C 75`, and then `-C 0`. Minimum mapping quality around 40 appeared optimal. Minimum base quality only influenced a few specific combinations without a consistent pattern across tests or domains.

<img width="3000" height="1800" alt="southern_extended_region_prop" src="https://github.com/user-attachments/assets/f222b88f-e6e7-4daf-b50b-2f1a01caa1c7" />

Figure 15. Southern domain: proportion of constrained variation uniquely explained by region in dd-RDA. See the Fig. 14 caption for more details.

In the Siberian domain (Fig. 16), the highest geographic signal was also observed with `-C 100` and `-minMapQ` ≥ 40. However, a paradoxical pattern emerged with `-C 0`, where increasing call thresholds sometimes reduced geographic signal, particularly at the highest mapping quality settings.

<img width="3000" height="1800" alt="siberia_extended_region_prop" src="https://github.com/user-attachments/assets/8caac9e5-6571-48f5-9a07-f4b85cc2275b" />

Figure 16. Siberian domain: proportion of constrained variation uniquely explained by region in dd-RDA. See the Fig. 14 caption for more details.

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
