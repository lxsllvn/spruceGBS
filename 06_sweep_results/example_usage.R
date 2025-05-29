# -------------------------------------------------------
# Summary heatmaps of within-population variation in
# individual heterozygosity
# -------------------------------------------------------

domains   <- c("northern", "southern", "siberia")
plot_types <- c("core", "extended", "pop")
fill_vars  <- c("rel_iqr_het", "iqr_het", "cv_het")

for (dom in domains) {
  dat <- import_indvhet(paste0(dom, "_indvhet_summary.tsv"))
  for (type in plot_types) {
    for (fill in fill_vars) {
      plot_indvhet(
        dat,
        fill_var    = fill,
        plot_type   = type,
        output_dir  = "indvhet",
        output_file = paste0(dom, "_", type, "_", fill, ".png")
      )
    }
  }
}

# -------------------------------------------------------
# Heatmaps of MANOVA variance explained by library, region, and their difference
# -------------------------------------------------------

domains    <- c("northern", "southern", "siberia")
plot_types <- c("core", "extended")
fill_types <- c("library", "region", "difference")

for (dom in domains) {
  dat <- read.csv(paste0(dom, "_pca_manova_summary.csv"))
  for (type in plot_types) {
    for (fill in fill_types) {
      plot_manova(
        dat,
        fill_type   = fill,
        plot_type   = type,
        output_dir  = "manova",
        output_file = paste0(dom, "_", type, "_", fill, ".png")
      )
    }
  }
}

# -------------------------------------------------------
# PCA biplots of genetic structure, colored by library and region
# -------------------------------------------------------

plot_pcangsd_sweep(
  domains        = c("siberia", "southern", "northern"),
  cov_dir_suffix = "collected_covs",
  ct_values      = c(4, 5, 6),
  var_list       = c("library", "region"),
  group_var      = c("baq", "C"),
  output_dir     = "sweep_biplots"
)

# -------------------------------------------------------
# Mixed effect models: effect of filter parameters on population genetic statistics
# -------------------------------------------------------

domains <- c("northern", "southern", "siberia")

for (dom in domains) {
  # Load and annotate data; filter F/absF by MAF (replaces with NA if MAF < maf_threshold)
  # Models are fit only to parameter combinations where:
  #   baq = 0, clipC in c(0, 60, 75, 100), call_thresh = 0.60
  df <- import_summary(paste0(dom, "_maf05_angsd_param_summaries.csv")) %>%
    filter_by_maf(maf_col = "MAF", maf_threshold = 0.05) %>%
    filter(baq == "0") %>%
    filter(clipC != "50") %>%
    filter(call_thresh == "0.6") %>%
    filter(F != "NA")

  # Fit maximal mixed effect models, perform model selection (MuMIn), and calculate marginal means
  result <- setNames(
    lapply(VALUE_STATS, function(stat) {
      model_selection(
        data           = df,
        response       = stat,
        fixed_formula  = "clipC * minMapQ * minQ",
        random_formula = "(1 | snpcode) + (1 | pop_code)",
        delta_thresh   = 2,
        rank           = "AIC"
      )$models
    }),
    VALUE_STATS
  )

  # Save results: best fitting models to three tables per domain
  write_models(
    result_list = result,
    base_name   = dom,
    outdir      = "output_tables"
  )
}

# -------------------------------------------------------
# Plot marginal means (for each domain/statistic)
# -------------------------------------------------------

for (dom in domains) {
  marg <- read_emms(paste0("output_tables/", dom, "_marginal_means.csv"))
  for (stat in VALUE_STATS) {
    if (!any(marg$statistic == stat & !is.na(marg$emmean))) {
      message(sprintf("No data for %s in domain %s; skipping.", stat, dom))
      next
    }
    plot_marginal_means(
      marg,
      target_stat  = stat,
      output_dir   = "emmeans",
      output_file  = paste0(dom, "_", stat, "_marginals.png")
    )
  }
}

# -------------------------------------------------------
# Exploratory figures by locus and population
# -------------------------------------------------------

# Load metadata for experiment populations
pop_status <- read.csv("angsd_parameter_exp_pops.csv")

for (dom in domains) {
  df <- import_summary(paste0(dom, "_maf05_angsd_param_summaries.csv")) %>%
    filter_by_maf(maf_col = "MAF", maf_threshold = 0.05) %>%
    filter(baq == "0") %>%
    filter(clipC != "50") %>%
    filter(call_thresh == "0.6") %>%
    filter(F != "NA")

  # Summarize by population
  pop_summary <- summarize_stats(
    df,
    stat_cols  = VALUE_STATS,
    param_cols = FILTER_PARAMS,
    mode       = "pop"
  ) %>%
    left_join(
      pop_status[, c("pop_code", "majority_library", "total_samples", "region", "is_mixed", "is_paired")],
      relationship = "many-to-many"
    )

  # Estimate within-region, among-population variation
  within_region_var <- pop_summary %>%
    group_by(across(all_of(FILTER_PARAMS)), statistic, region) %>%
    summarize(
      among_pop_cv       = sd(mean_value, na.rm = TRUE) / abs(mean(mean_value, na.rm = TRUE)),
      among_pop_iqr      = IQR(median_value, na.rm = TRUE),
      among_pop_rel_iqr  = IQR(median_value, na.rm = TRUE) / abs(median(median_value, na.rm = TRUE)),
      mean_within_iqr    = mean(iqr_value, na.rm = TRUE),
      mean_within_rel_iqr= mean(rel_iqr_value, na.rm = TRUE),
      mean_within_cv     = mean(cv_value, na.rm = TRUE),
      n_pops             = n_distinct(pop_code),
      .groups = "drop"
    ) %>%
    filter(n_pops >= 3)

  # Summarize by locus
  loci_summary <- summarize_stats(
    df,
    stat_cols  = VALUE_STATS,
    param_cols = FILTER_PARAMS,
    mode       = "loci"
  )

  # RMSE of F by filter parameters (heatmap)
  plot_RMSE(
    pop_summary = pop_summary,
    output_dir  = "popgen_plots",
    output_file = paste0(dom, "_F_RMSE_heatmap.png")
  )

  # Plots of population genetic statistics
  for (stat in VALUE_STATS) {
    # Heatmap of population means (faceted by clipC, minQ, minMapQ)
    plot_pop_heatmaps(
      pop_summary = pop_summary,
      statistic   = stat,
      fill_var    = "mean_value",
      output_dir  = "popgen_plots",
      output_file = paste0(dom, "_", stat, "_pop_heatmap.png")
    )

    # Density plots of population means by filter parameter level
    plot_stat_distributions(
      pop_summary,
      stat_name   = stat,
      value_col   = "mean_value",
      param_list  = c("clipC", "minMapQ", "minQ"),
      scale       = 5,
      alpha       = 0.4,
      quantile_lines = TRUE,
      ncol        = 2,
      output_dir  = "popgen_plots",
      output_file = paste0(dom, "_", stat, "_pop_density.png")
    )

    # Density plots of locus means by filter parameter level
    plot_stat_distributions(
      loci_summary,
      stat_name   = stat,
      value_col   = "mean_value",
      param_list  = c("clipC", "minMapQ", "minQ"),
      scale       = 5,
      alpha       = 0.4,
      quantile_lines = TRUE,
      ncol        = 2,
      output_dir  = "popgen_plots",
      output_file = paste0(dom, "_", stat, "_loci_density.png")
    )

    # Plots of within-region variation among populations
    for (var in c("among_pop_cv", "among_pop_rel_iqr")) {
      plot_popvar(
        within_region_var,
        statistic   = stat,
        fill_var    = var,
        output_dir  = "popgen_plots",
        output_file = paste0(dom, "_", stat, "_region_var.png")
      )
    }
  }
}
