#' ---
#' Title: ANGSD parameter sweep analysis functions
#' Description: Import, summarize, model, and plot
#' results of the individual heterozygosity, 
#' PCA/MANOVA, RDA on PCs, and site and 
#' population-level statistics experiments.
#' ---

# ============================
#      Package Loading
# ============================
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(lme4)
library(MuMIn)
library(emmeans)
library(purrr)
library(ggridges)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(vegan)

# ============================
#      Global Constants
# ============================
FILTER_PARAMS <- c("baq", "clipC", "minQ", "minMapQ", "call_thresh")
VALUE_STATS   <- c("pi", "theta_W", "MAF", "Hexp", "Hobs", "F", "absF", "tajima_d")

# ============================
# 1. Import/Utility Functions
# ============================
#'
#' Null coalescing operator
#'
#' Returns the left-hand side if it is not \code{NULL}, otherwise returns the right-hand side.
#'
#' @param a First value to test for \code{NULL}.
#' @param b Value to return if \code{a} is \code{NULL}.
#' @return \code{a} if not \code{NULL}, otherwise \code{b}.
#' @examples
#' 1 %||% 2      # returns 1
#' NULL %||% 5   # returns 5
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Import and parse summary CSV file.
#'
#' @param filepath Path to CSV file.
#' @return Tibble with parsed parameter columns.
import_summary <- function(filepath) {
  df <- read_csv(filepath, show_col_types = FALSE)
  df <- df %>%
    mutate(
      pi       = exp(pi),
      theta_W  = exp(theta_W),
      absF     = abs(F),
      tajima_d = pi - theta_W,
      domain   = str_extract(param_id, "^[^_]+"),
      baq      = as.factor(str_extract(param_id, "(?<=_baq)\\d+")),
      clipC    = as.factor(str_extract(param_id, "(?<=_C)\\d+")),
      minQ     = as.factor(str_extract(param_id, "(?<=_q)\\d+")),
      minMapQ  = as.factor(str_extract(param_id, "(?<=_mq)\\d+")),
      call_thresh = as.factor(call_thresh)
    ) %>%
    select(snpcode, all_of(VALUE_STATS), pop_code, param_id, call_thresh, domain, all_of(FILTER_PARAMS))
  return(df)
}

#' Import individual heterozygosity summary table.
#'
#' @param filepath Path to table.
#' @return Tibble of summarized heterozygosity stats.
import_indvhet <- function(filepath) {
  df <- read.table(filepath, header = TRUE)
  df <- df %>%
    group_by(baq, C, q, mq, ct, pop_code) %>%
    summarize(
      mean_het    = mean(alt_sites/(alt_sites+ref_sites), na.rm = TRUE),
      sd_het      = sd(alt_sites/(alt_sites+ref_sites),   na.rm = TRUE),
      cv_het      = ifelse(mean_het == 0, NA, sd_het / abs(mean_het)),
      iqr_het     = IQR(alt_sites/(alt_sites+ref_sites), na.rm=TRUE),
      rel_iqr_het = IQR(alt_sites/(alt_sites+ref_sites), na.rm=TRUE) / abs(median(alt_sites/(alt_sites+ref_sites), na.rm=TRUE)),
      n_samples   = n(),
      .groups = "drop"
    ) %>%
    mutate(baq    = factor(baq, levels = 0:2, labels = paste0("baq: ", 0:2)),
           minQ   = factor(q),
           minMapQ= factor(mq),
           call_thresh = factor(ct),
           clipC  = factor(C))
  return(df)
}

#' Read marginal means table from CSV and flexibly type columns
#'
#' Reads a CSV of estimated marginal means (EMMs), coercing all columns to numeric
#' except for \code{statistic} and \code{formula}, which are always character.
#' Suppresses parsing warnings (e.g., from blank/incomplete lines) and filters out rows
#' with NA in \code{emmean}.
#'
#' @param path Path to the CSV file.
#'
#' @return A tibble where all columns are numeric except \code{statistic} and \code{formula}, which are character.
#'
#' @details
#' All columns except \code{statistic} and \code{formula} are coerced to numeric. 
#' If you have additional non-numeric columns, they will be read as numeric (becoming NA if not coercible).
#' Suppresses parsing warnings by default; remove \code{suppressWarnings()} if you want to see them.
#'
#' @examples
#' df <- read_emms("my_emms_table.csv")
#'
#' @export
read_emms <- function(path) {
  # Get column names from the first row
  col_names <- names(readr::read_csv(path, n_max = 0, show_col_types = FALSE))
  # Prepare col_types: everything is col_double() except for statistic and formula
  col_types <- rep(list(readr::col_double()), length(col_names))
  names(col_types) <- col_names
  col_types[["statistic"]] <- readr::col_character()
  col_types[["formula"]]   <- readr::col_character()
  
  df <- suppressWarnings(readr::read_csv(
    path,
    col_types = do.call(readr::cols, col_types),
    show_col_types = FALSE
  ))
  if ("emmean" %in% colnames(df)) {
    df <- dplyr::filter(df, !is.na(.data$emmean))
  }
  return(df)
}

#' Write a UNIX-style TSV file.
#'
#' @param x Data.frame or tibble to write.
#' @param out.name Output filename.
#' @param ... Passed to write.table().
#' @return Invisibly NULL.
write.unix <- function(x, out.name, row.names = FALSE, col.names = FALSE,
                       quote = FALSE, sep = "\t", ...) {
  con <- file(out.name, "wb")
  write.table(x, file = con, row.names = row.names, col.names = col.names,
              quote = quote, sep = sep, ...)
  close(con)
  invisible(NULL)
}

#' Load sample list and add metadata.
#'
#' @param sample.list Path to sample list.
#' @param metadata Metadata CSV file.
#' @return Data.frame with merged metadata.
get.info <- function(sample.list, metadata) {
  indv.info <- read.table(paste0(sample.list))
  colnames(indv.info)[1] <- "bam_code"
  meta <- read.csv(paste0(metadata))
  indv.info <- left_join(indv.info, 
                         meta[,c("bam_code", "pop_code", "library", "domain", "region", 
                                 "latitude", "longitude")])
  return(indv.info)
}

#' Filter F and absF estimates by MAF threshold.
#'
#' @param df Data.frame with F and MAF columns.
#' @param maf_col Name of MAF column.
#' @param maf_threshold MAF threshold (default 0.05).
#' @return Data.frame with F and absF set to NA if MAF < threshold.
filter_by_maf <- function(df, maf_col = "MAF", maf_threshold = 0.05) {
  out <- df
  if ("F" %in% names(out) && maf_col %in% names(out)) {
    out$F[out[[maf_col]] < maf_threshold] <- NA
  }
  if ("absF" %in% names(out) && maf_col %in% names(out)) {
    out$absF[out[[maf_col]] < maf_threshold] <- NA
  }
  return(out)
}

# ==============================
# 2. Data Summarization Functions
# ==============================

#' Summarize statistics by parameter set.
#'
#' @param df Data.frame to summarize.
#' @param stat_cols Columns to summarize (default: VALUE_STATS).
#' @param param_cols Parameter columns to group by (default: FILTER_PARAMS).
#' @param mode "loci" or "pop".
#' @return Tibble with summary statistics.
summarize_stats <- function(
    df,
    stat_cols   = VALUE_STATS,
    param_cols  = FILTER_PARAMS,
    mode        = c("loci", "pop")
) {
  mode <- match.arg(mode)
  group_cols <- if (mode == "pop") c("pop_code", param_cols) else param_cols
  bind_rows(lapply(stat_cols, function(stat) {
    if (!(stat %in% names(df))) return(NULL)
    df %>%
      group_by(across(all_of(group_cols))) %>%
      summarize(
        statistic     = stat,
        mean_value    = mean(.data[[stat]], na.rm=TRUE),
        median_value  = median(.data[[stat]], na.rm=TRUE),
        sd_value      = sd(.data[[stat]], na.rm=TRUE),
        iqr_value     = IQR(.data[[stat]], na.rm=TRUE),
        var_value     = var(.data[[stat]], na.rm=TRUE),
        cv_value      = ifelse(mean_value == 0, NA, sd_value / abs(mean_value)),
        q25           = quantile(.data[[stat]], 0.25, na.rm=TRUE),
        q75           = quantile(.data[[stat]], 0.75, na.rm=TRUE),
        rel_iqr_value = IQR(.data[[stat]], na.rm=TRUE) / abs(median_value),
        n_snps        = n_distinct(snpcode[!is.na(.data[[stat]])]), 
        .groups = "drop"
      )
  }))
}

# ====================
# 3. Modeling Functions
# ====================

#' Extract factor variables from a mixed effect model.
#'
#' @param model lmer model object.
#' @return Character vector of factor variable names.
get_factor_vars <- function(model) {
  mf <- tryCatch(model@frame, error = function(e) NULL)
  if (is.null(mf)) return(character(0))
  terms_obj <- terms(model)
  terms <- attr(terms_obj, "term.labels")
  if (is.null(terms) || length(terms) == 0) return(character(0))
  terms <- as.character(terms)
  terms[sapply(terms, function(x) {
    if (length(x) == 0 || !nzchar(x)) return(FALSE)
    all(all.vars(as.formula(paste("~", x))) %in% names(mf)) && is.factor(mf[[x]])
  })]
}

#' Mixed model selection with MuMIn dredge.
#'
#' @param data Data frame.
#' @param response Response variable name.
#' @param fixed_formula Fixed effect formula string.
#' @param random_formula Random effect formula string.
#' @param delta_thresh Model selection threshold (default 2).
#' @param rank Model selection rank metric (default "AICc").
#' @return List of selected models with formula, variance components, marginal means, and delta.
model_selection <- function(
    data,
    response,
    fixed_formula,
    random_formula,
    delta_thresh = 2,
    rank = "AICc"
) {
  full_formula <- as.formula(paste(response, "~", fixed_formula, "+", random_formula))
  global_model <- lmer(full_formula, data = data)
  options(na.action = "na.fail")
  model_set <- MuMIn::dredge(global_model, rank = rank, trace = FALSE)
  best_models <- MuMIn::get.models(model_set, subset = delta < delta_thresh)
  get_factor_vars <- function(model) {
    mf <- tryCatch(model@frame, error = function(e) NULL)
    if (is.null(mf)) return(character(0))
    terms_obj <- terms(model)
    terms <- attr(terms_obj, "term.labels")
    if (is.null(terms) || length(terms) == 0) return(character(0))
    terms <- as.character(terms)
    terms[sapply(terms, function(x) {
      if (length(x) == 0 || !nzchar(x)) return(FALSE)
      all(all.vars(as.formula(paste("~", x))) %in% names(mf)) && is.factor(mf[[x]])
    })]
  }
  extract_model_info <- function(mod, idx) {
    formula_chr <- format(formula(mod))
    model_delta <- model_set$delta[as.numeric(names(best_models))[idx]]
    varcomp <- tryCatch({
      as.data.frame(VarCorr(mod))[, c("grp", "var1", "var2", "vcov", "sdcor")]
    }, error = function(e) NULL)
    factor_vars <- get_factor_vars(mod)
    marginal_means <- NULL
    if (length(factor_vars) > 0) {
      em_formula <- as.formula(paste("~", paste(factor_vars, collapse = "*")))
      marginal_means <- tryCatch({
        as.data.frame(emmeans(mod, em_formula))
      }, error = function(e) NULL)
    }
    list(
      formula         = formula_chr,
      var_components  = varcomp,
      marginal_means  = marginal_means,
      delta           = model_delta
    )
  }
  results <- purrr::imap(best_models, extract_model_info)
  names(results) <- vapply(results, function(x) x$formula, character(1L))
  list(models = results)
}

#' Write model selection results to CSV.
#'
#' @param result_list Output from model_selection.
#' @param base_name Prefix for output files.
#' @param outdir Output directory.
#' @return List of tibbles invisibly.
write_models <- function(
    result_list,
    base_name = "results",
    outdir = "model_results"
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  # Formulas
  all_formulas <- purrr::map2_dfr(
    names(result_list), result_list,
    function(stat, stat_list) {
      purrr::map_dfr(names(stat_list), function(model_name) {
        tibble::tibble(
          statistic = stat,
          formula   = stat_list[[model_name]]$formula
        )
      })
    }
  )
  write.csv(
    all_formulas,
    file.path(outdir, paste0(base_name, "_all_formulas.csv")),
    row.names = FALSE
  )
  # Variance components
  all_varcomps <- purrr::map2_dfr(
    names(result_list), result_list,
    function(stat, stat_list) {
      purrr::map_dfr(names(stat_list), function(model_name) {
        vc <- stat_list[[model_name]]$var_components
        if (!is.null(vc)) dplyr::mutate(vc, statistic = stat, formula = stat_list[[model_name]]$formula)
      })
    }
  )
  write.csv(
    all_varcomps,
    file.path(outdir, paste0(base_name, "_var_components.csv")),
    row.names = FALSE
  )
  # Marginal means
  all_margmeans <- purrr::map2_dfr(
    names(result_list), result_list,
    function(stat, stat_list) {
      purrr::map_dfr(names(stat_list), function(model_name) {
        mm <- stat_list[[model_name]]$marginal_means
        if (!is.null(mm)) dplyr::mutate(mm, statistic = stat, formula = stat_list[[model_name]]$formula)
      })
    }
  )
  write.csv(
    all_margmeans,
    file.path(outdir, paste0(base_name, "_marginal_means.csv")),
    row.names = FALSE
  )
  message("All summary tables written to: ", normalizePath(outdir))
  invisible(list(
    all_formulas   = all_formulas,
    var_components = all_varcomps,
    marginal_means = all_margmeans
  ))
}

#' Run variance partitioning (varpart) or RDA for all parameter combos in a PCangsd covariance folder.
#'
#' @param domain         Character. E.g. "northern"
#' @param cov_dir_suffix Suffix for covariance directory (default: "collected_covs")
#' @param n_pc           Number of PCs to use (default: 2)
#' @param output_dir     Optional directory to save output CSV
#' @return               Data frame of results
#' @export
pcangsd_rda <- function(domain,
                        cov_dir_suffix = "collected_covs",
                        n_pc           = 2,
                        output_dir     = NULL) {
  
  library(dplyr); library(stringr); library(tibble); library(vegan)
  
  # Get sample metadata
  indv_info <- get.info(paste0("parameter_exp_", domain, "_bamlist.txt"),
                        "sequenced_samples_metadata.csv")
  # List all covariance files
  cov_folder <- file.path(paste0(domain, "_", cov_dir_suffix))
  files <- list.files(path = cov_folder, pattern = "\\.Pcangsd\\.cov$", full.names = TRUE)
  if (length(files) == 0) stop("No covariance files found.")
  
  results <- list()
  for (f in files) {
    # Parse parameters from filename
    param_id <- sub("\\.Pcangsd\\.cov$", "", basename(f))
    params <- tibble(
      param_id = param_id,
      baq      = str_extract(param_id, "(?<=_baq)\\d+"),
      C        = str_extract(param_id, "(?<=_C)\\d+"),
      minQ     = str_extract(param_id, "(?<=_q)\\d+"),
      minMapQ  = str_extract(param_id, "(?<=_mq)\\d+"),
      ct       = str_extract(param_id, "(?<=_ct)\\d+")
    )
    
    # Read covariance, get eigenvectors (PCs)
    mat <- as.matrix(read.table(f, header = FALSE))
    eig <- eigen(mat)
    pc_names <- paste0("PC", seq_len(n_pc))
    df_ev <- indv_info %>%
      mutate(!!!setNames(as.data.frame(eig$vectors[, seq_len(n_pc)]), pc_names)) %>%
      filter(complete.cases(.))
    Y <- as.matrix(df_ev[, pc_names])
    
    # Run variance partitioning (varpart)
    vp <- tryCatch({
      vegan::varpart(Y, ~region, ~library, data = df_ev)
    }, error = function(e) NULL)
    
    # Extract results
    if (!is.null(vp)) {
      adj <- vp$part$indfract$Adj.R.squared
      names(adj) <- c("region_unique", "library_unique", "shared", "residual")
      res <- bind_cols(params, as_tibble_row(adj))
      results[[length(results)+1]] <- res
    }
  }
  
  results_df <- bind_rows(results)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    out_file <- file.path(output_dir, paste0(domain, "_rda_varpart_summary.csv"))
    readr::write_csv(results_df, out_file)
    message("Results written to: ", out_file)
  }
  return(results_df)
}


# ===========================
# 4. Plotting Functions
# ===========================

# --- Ridgeline Distribution Plot ---

#' Plot distribution of a statistic by parameter levels.
#'
#' @param data Dataframe to plot.
#' @param stat_name Statistic name (e.g. "F").
#' @param value_col Column with values (default "value").
#' @param param_list List of parameters to facet by.
#' @param output_dir Directory to save plot.
#' @param output_file Filename for output.
#' @param ncol Number of columns for facets.
#' @param ... Additional arguments to geom_density_ridges.
#' @return Plot or saves image.
plot_stat_distributions <- function(
    data,
    stat_name,
    value_col     = "value",
    output_dir    = NULL,
    output_file   = NULL,
    param_list    = FILTER_PARAMS,
    filter_by_quantiles = FALSE,
    quantile_range      = c(0.025, 0.975),
    scale         = 0.9,
    alpha         = 0.7,
    maf_threshold = 0.05,
    plot_title    = NULL,
    ncol          = NULL,    
    ...
) {
    df <- data %>% filter(statistic == stat_name)

  # Compute x-limits if requested
  xlims <- NULL
  if (filter_by_quantiles) {
    xlims <- quantile(df[[value_col]], probs = quantile_range, na.rm = TRUE)
  }
  # Pivot parameters into long for ridge plotting
  long_data <- df %>%
    tidyr::pivot_longer(cols = all_of(param_list),
                        names_to  = "parameter",
                        values_to = "level")
  # Build title
  title <- if (is.null(plot_title)) {
    paste(stat_name, "distribution by parameter")
  } else plot_title
  # Facet setup
  facet_args <- list(~parameter, scales = "free_y")
  if (!is.null(ncol)) facet_args$ncol <- ncol
  # Create plot
  p <- ggplot(long_data, aes(x = .data[[value_col]],
                             y = factor(level),
                             fill = factor(level),
                             group = level)) +
    ggridges::geom_density_ridges(scale = scale, alpha = alpha,
                                  rel_min_height = 0.01,
                                  color = "gray40", ...) +
    do.call(facet_wrap, facet_args) +
    theme_minimal() +
    labs(title = title, x = "", y = "parameter level") +
    theme(legend.position = "none",
          strip.text = element_text(size = 8))
  # Apply x-limits
  if (!is.null(xlims)) {
    p <- p + coord_cartesian(xlim = xlims)
  }
  # Output
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (is.null(output_file)) output_file <- paste0(stat_name, "_all_params.png")
    save_path <- file.path(output_dir, output_file)
    ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
    message("Saved plot to ", save_path)
  } else {
    print(p)
  }
}

# --- Heatmaps of Individual Heterozygosity ---

#' Plot heatmap of individual heterozygosity.
#'
#' @param df Data.frame from import_indvhet.
#' @param fill_var Column to use as fill.
#' @param plot_type Type ("core", "extended", "pop").
#' @param output_dir Output directory.
#' @param output_file Output file name.
plot_indvhet <- function(
    df,
    fill_var    = "iqr_het",
    plot_type   = c("core", "extended", "pop"),
    output_dir  = NULL,
    output_file = NULL
) {
  plot_type <- match.arg(plot_type)
  if (plot_type == "core") {
    df_plot <- df %>% 
      filter(minMapQ != 50) %>%
      filter(clipC %in% c("0", "50")) %>%
      mutate(C_ct = paste0("c: ", clipC, " / ct: ", call_thresh))
    p <- ggplot(df_plot, aes(x = minQ, y = minMapQ, fill = .data[[fill_var]])) +
      geom_tile(linetype = "blank") +
      facet_grid(baq ~ C_ct, labeller = label_value)
  }
  if (plot_type == "extended") {
    df_plot <- df %>% filter(baq == "baq: 0")
    p <- ggplot(df_plot, aes(x = minQ, y = minMapQ, fill = .data[[fill_var]])) +
      geom_tile(linetype = "blank") +
      facet_grid(clipC ~ call_thresh, labeller = label_value)
  }
  if (plot_type == "pop") {
    df_plot <- df %>% filter(baq == "baq: 0", call_thresh == "6", clipC %in% c("0", "60", "75", "100"))
    p <- ggplot(df_plot, aes(x = minQ, y = pop_code, fill = .data[[fill_var]])) +
      geom_tile(linetype = "blank") +
      facet_grid(clipC ~ minMapQ, labeller = label_value)
  }
  p <- p +
    scale_fill_viridis_c() +
    labs(x = "minQ", y = ifelse(plot_type == "pop", "pop_code", "minMapQ")) +
    theme_minimal() +
    theme(
      panel.spacing = unit(0, "pt"),
      strip.placement = "outside",
      strip.text.x.top = element_text(angle = 45),
      strip.text.y.right = element_text(angle = 0)
    )
  # Output
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (is.null(output_file)) output_file <- paste0("indv_het_", fill_var, "_", plot_type, ".png")
    save_path <- file.path(output_dir, output_file)
    ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
    message("Saved plot to ", save_path)
  } else {
    print(p)
  }
  invisible(p)
}

# --- Heatmaps for MANOVA Effects ---

#' Plot MANOVA heatmap by parameter settings.
#'
#' @param df Data frame with MANOVA summary.
#' @param fill_type Fill variable: "library", "region", or "difference".
#' @param plot_type "core" or "extended".
#' @param output_dir Output directory.
#' @param output_file Output file name.
plot_manova <- function(
    df,
    fill_type   = c("library", "region", "difference"),
    plot_type   = c("core", "extended"),
    output_dir  = NULL,
    output_file = NULL
) {
  fill_type <- match.arg(fill_type)
  plot_type <- match.arg(plot_type)
  df <- df %>%
    dplyr::mutate(
      baq        = factor(baq, levels = 0:2, labels = paste0("baq: ", 0:2)),
      minQ       = factor(q),
      minMapQ    = factor(mq),
      call_thresh= factor(ct),
      clipC      = factor(C),
      C_ct       = paste0("c: ", clipC, " / ct: ", call_thresh)
    )
  fill_var <- switch(
    fill_type,
    "library"   = "R2_library",
    "region"    = "R2_region",
    "difference"= NULL
  )
  fill_label <- switch(
    fill_type,
    "library"   = expression(R^2[library]),
    "region"    = expression(R^2[region]),
    "difference"= expression(R^2[region] - R^2[library])
  )
  if (fill_type == "difference") {
    df <- df %>%
      dplyr::mutate(R2_diff = R2_region - R2_library)
    fill_var <- "R2_diff"
  }
  if (plot_type == "core") {
    df_plot <- df %>%
      dplyr::filter(minMapQ != 50, clipC %in% c("0", "50"))
    p <- ggplot(df_plot, aes(x = minQ, y = minMapQ, fill = .data[[fill_var]])) +
      geom_tile(linetype = "blank") +
      facet_grid(baq ~ C_ct, labeller = label_value)
  }
  if (plot_type == "extended") {
    df_plot <- df %>% dplyr::filter(baq == "baq: 0")
    p <- ggplot(df_plot, aes(x = minQ, y = minMapQ, fill = .data[[fill_var]])) +
      geom_tile(linetype = "blank") +
      facet_grid(clipC ~ call_thresh, labeller = label_value)
  }
  p <- p +
    scale_fill_viridis_c(name = fill_label) +
    labs(x = "minQ", y = "minMapQ") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      panel.spacing = unit(0, "pt"),
      strip.placement = "outside",
      strip.text.x.top = element_text(angle = 45),
      strip.text.y.right = element_text(angle = 0)
    )
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (is.null(output_file)) {
      output_file <- paste0("manova_heatmap_", fill_type, "_", plot_type, ".png")
    }
    save_path <- file.path(output_dir, output_file)
    ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
    message("Saved plot to ", save_path)
  } else {
    print(p)
  }
  invisible(p)
}

# --- PCA Plots for PCANGSD Covariances ---

#' Plot PCANGSD PCA plots for all domains/parameter sets.
#'
#' @param domains Character vector of domain names.
#' @param var_list Metadata columns to color by.
#' @param group_var Covariate names for PCA grouping.
#' @param cov_dir_suffix Covariance folder suffix.
#' @param ct_values Vector of call_thresh values.
#' @param output_dir Output directory.
plot_pcangsd_sweep <- function(domains,
                         var_list,
                         group_var,
                         cov_dir_suffix = "collected_covs",
                         ct_values = c(4,5,6),
                         output_dir = getwd()) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  message("Working directory: ", normalizePath(getwd()))
  message("Output directory: ", normalizePath(output_dir))
  pal_Set2 <- brewer.pal(8, "Set2")
  pal_Set3 <- brewer.pal(9, "Set1") %>% setdiff("#FFFF33")
  for (domain in domains) {
    indv_info <- get.info(paste0("parameter_exp_", domain, "_bamlist.txt"),
                          "sequenced_samples_metadata.csv")
    cov_folder <- file.path(paste0(domain, "_", cov_dir_suffix))
    for (ct in ct_values) {
      pattern <- paste0("_ct", ct, "\\.Pcangsd\\.cov$")
      files <- list.files(path = cov_folder,
                          pattern = pattern,
                          full.names = TRUE)
      if (length(files) == 0) {
        warning("No covariance files found in ", cov_folder, " matching pattern ", pattern)
        next
      }
      cov_list <- lapply(files, read.table, header = FALSE)
      names(cov_list) <- sub("\\.Pcangsd\\.cov$", "", basename(files))
      lookup <- tibble(param_id = names(cov_list)) %>%
        mutate(
          domain  = str_extract(param_id, "^[^_]+"),
          baq     = str_extract(param_id, "(?<=_baq)\\d+"),
          C       = str_extract(param_id, "(?<=_C)\\d+"),
          minQ    = str_extract(param_id, "(?<=_q)\\d+"),
          minMapQ = str_extract(param_id, "(?<=_mq)\\d+"),
          ct      = str_extract(param_id, "(?<=_ct)\\d+")
        ) %>%
        mutate(across(all_of(group_var), as.factor),
               row_id = row_number())
      group_rows <- lookup %>%
        group_by(across(all_of(group_var))) %>%
        summarise(rows = list(row_id), .groups = "drop")
      for (color_var in var_list) {
        for (i in seq_len(nrow(group_rows))) {
          idxs <- group_rows$rows[[i]]
          plots <- lapply(idxs, function(j) {
            mat <- as.matrix(cov_list[[j]])
            eig <- eigen(mat)
            df_ev <- indv_info %>%
              mutate(
                PC1 = eig$vectors[,1],
                PC2 = eig$vectors[,2]
              )
            ggplot(df_ev, aes_string(x = "PC1", y = "PC2", color = color_var)) +
              geom_point(size = 2, alpha = 0.5, shape = 16) +
              labs(
                title = names(cov_list)[j],
                color = color_var
              ) +
              {if(color_var=="region") {
                scale_color_brewer(palette = "Set2")
              } else {
                scale_color_manual(values = pal_Set3)
              }} +
              theme_minimal() +
              theme(plot.title = element_text(size = 8),
                    axis.title = element_text(size = 8))
          })
          grid <- wrap_plots(plotlist = plots,
                             ncol     = min(3, length(plots)),
                             nrow     = ceiling(length(plots)/3),
                             guides   = "collect") &
            theme(legend.position = "right",
                  legend.spacing = unit(1, "mm"),
                  legend.key.spacing.x =  unit(1, "mm"),
                  legend.key.spacing.y =  unit(1, "mm"))
          first_id <- names(cov_list)[idxs][1]
          parts <- strsplit(first_id, "_")[[1]][1:3]
          out_prefix <- paste(parts, collapse = "_")
          out_file <- file.path(output_dir,
                                paste0(out_prefix, "_pcangsd_", color_var, ".ct", ct, ".png"))
          message("Saving plot to: ", out_file)
          ggsave(filename = out_file,
                 plot     = grid,
                 width    = 9,
                 height   = 5,
                 dpi      = 300)
        }
      }
    }
  }
}

# --- Plot of marginal means from mixed-effect models---

#' Plot marginal means with automated y-axis from formula
#'
#' Reads y-axis variables from the formula column, groups and deduplicates appropriately,
#' and plots the marginal means for a given statistic.
#' Optionally saves the plot to a file.
#'
#' @param df Data frame with marginal means (must include columns: emmean, formula, statistic)
#' @param target_stat Statistic to plot (e.g. "pi")
#' @param x_var X-axis variable (default "emmean")
#' @param lower Name for lower CI column (optional, for error bars)
#' @param upper Name for upper CI column (optional, for error bars)
#' @param output_file Output file name (optional, e.g. "pi_marginal_means.png")
#' @param output_dir Directory to save plot (optional, e.g. "results/plots")
#' @param width Plot width for saving (default 8)
#' @param height Plot height for saving (default 6)
#' @param dpi Plot resolution for saving (default 300)
#' @param ... Other args to ggplot2::geom_point()
#'
#' @export
plot_marginal_means <- function(
    df,
    target_stat,
    x_var = "emmean",
    lower = "asymp.LCL",
    upper = "asymp.UCL",
    output_file = NULL,
    output_dir = NULL,
    width = 8,
    height = 6,
    dpi = 300,
    ...
) {
  # Helper to extract grouping vars from formula column
  extract_yvars_from_formula <- function(formula_str) {
    rhs <- sub(".*~", "", formula_str)
    rhs <- gsub("\\s+", "", rhs)
    rhs <- gsub("\\([^)]*\\)", "", rhs)  # remove (1|something)
    vars <- unlist(strsplit(rhs, "[+:]"))
    vars <- vars[vars != ""]
    unique(vars)
  }
  
  # 1. Filter to the right statistic
  df_stat <- df %>%
    dplyr::filter(statistic == target_stat, !is.na(.data[[x_var]]))
  
  # 2. Get y_vars from formula
  formula_str <- unique(df_stat$formula)[1]
  y_vars <- extract_yvars_from_formula(formula_str)
  
  # 3. Deduplicate if necessary
  df_stat <- df_stat %>%
    dplyr::distinct(dplyr::across(all_of(y_vars)), .keep_all = TRUE)
  
  # 4. Create y column (for interaction if multiple)
  y_is_interaction <- length(y_vars) > 1
  df_stat <- df_stat %>%
    dplyr::mutate(
      y = if (y_is_interaction) interaction(!!!rlang::syms(y_vars), sep = ":")
      else .data[[y_vars]]
    )
  
  # 5. Ensure numeric x and error bars (if needed)
  df_stat[[x_var]] <- as.numeric(df_stat[[x_var]])
  if (lower %in% names(df_stat) && upper %in% names(df_stat)) {
    df_stat[[lower]] <- as.numeric(df_stat[[lower]])
    df_stat[[upper]] <- as.numeric(df_stat[[upper]])
  }
  
  # 6. Plot!
  p <- ggplot2::ggplot(df_stat, ggplot2::aes(x = .data[[x_var]], y = y)) +
    ggplot2::geom_point(...) +
    ggplot2::labs(
      x = x_var,
      y = if (y_is_interaction) paste(y_vars, collapse = ":") else y_vars
    )
  
  # Add error bars if present
  if (lower %in% names(df_stat) && upper %in% names(df_stat)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data[[lower]], xmax = .data[[upper]]), width = 0.15)
  }
  
  # Clean up x axis a bit
  p <- p +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    ggplot2::theme_minimal()
  
  # 7. Save or print
  if (!is.null(output_file)) {
    save_path <- output_file
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      save_path <- file.path(output_dir, output_file)
    }
    ggplot2::ggsave(
      filename = save_path, plot = p, width = width, height = height, dpi = dpi, bg = "white"
    )
    message("Plot saved to: ", save_path)
    invisible(save_path)
  } else {
    print(p)
    invisible(p)
  }
}

# --- Heatmaps of Population Variation ---

#' Plot heatmap of within-region, among-population variation.
#'
#' @param df Data frame with pop summary.
#' @param statistic Statistic (e.g. "F").
#' @param fill_var Column for fill.
#' @param title Plot title.
#' @param output_dir Output directory.
#' @param output_file Output filename.
plot_popvar <- function(
    df,
    statistic   = "F",
    fill_var    = "among_pop_cv",
    title       = NULL,
    output_dir  = NULL,
    output_file = NULL,
    width       = 7,
    height      = 5,
    dpi         = 300
) {
  dat <- df %>% dplyr::filter(statistic == !!statistic)
  p <- ggplot(dat, aes(x = minQ, y = minMapQ, fill = .data[[fill_var]])) +
    geom_tile(linetype = "blank") +
    facet_grid(region~clipC, labeller = label_value) +
    scale_fill_viridis_c(name = NULL) +
    labs(
      title = title %||% paste0(fill_var, " for ", statistic),
      x = "minQ", y = "minMapQ"
    ) +
    theme_minimal()
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (is.null(output_file)) {
      output_file <- paste0("cv_heatmap_", statistic, "_", fill_var, ".png")
    }
    save_path <- file.path(output_dir, output_file)
    ggsave(save_path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
    message("Saved plot to ", save_path)
    invisible(save_path)
  } else {
    print(p)
    invisible(p)
  }
}

# --- Heatmaps for Population-Specific Plots ---

#' Plot per-population heatmaps.
#'
#' @param pop_summary Pop summary dataframe.
#' @param statistic Statistic (e.g. "F").
#' @param fill_var Fill variable.
#' @param output_file Output filename.
#' @param output_dir Output directory.
plot_pop_heatmaps <- function(
    pop_summary,
    statistic   = "F",
    fill_var    = "mean_value",
    output_file = NULL,
    output_dir  = NULL,
    width       = 8,
    height      = 6,
    dpi         = 300
) {
  dat <- pop_summary %>% dplyr::filter(statistic == !!statistic)
  if (nrow(dat) == 0) stop("No data for the specified statistic.")
  fill_range <- range(dat[[fill_var]], na.rm = TRUE)
  pop_codes <- unique(dat$pop_code)
  n_pops <- length(pop_codes)
  ncol  <- ceiling(sqrt(n_pops))
  nrow  <- ceiling(n_pops / ncol)
  fill_scale <- if (statistic == "F") {
    scale_fill_distiller(
      type = "div",
      palette = "RdBu",
      limits = fill_range, 
      name = fill_var,
      oob = scales::squish
    )
  } else {
    scale_fill_viridis_c(
      limits = fill_range, 
      name = fill_var,
      oob = scales::squish
    )
  }
  plots <- lapply(seq_along(pop_codes), function(i) {
    pop <- pop_codes[i]
    row_idx <- ceiling(i / ncol)
    col_idx <- ifelse(i %% ncol == 0, ncol, i %% ncol)
    dat_pop <- dat[dat$pop_code == pop,]
    ggplot(dat_pop, aes(y = minMapQ, x = minQ, fill = .data[[fill_var]])) +
      geom_tile(linetype = "blank") +
      facet_grid(~clipC, labeller = label_value) +
      fill_scale +
      labs(
        title = pop,
        x = if (row_idx == nrow) "min Q" else NULL,
        y = if (col_idx == 1) "min MQ" else NULL
      ) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.x = if (row_idx == nrow) element_text() else element_blank(),
        axis.text.y = if (col_idx == 1) element_text() else element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 8, hjust = 0.5),
        strip.text.y.right = element_text(angle = 0)
      )
  })
  legend_plot <- ggplot(dat, aes(y = minMapQ, x = minQ, fill = .data[[fill_var]])) +
    geom_tile(linetype = "blank") +
    fill_scale +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.key.height = unit(0.15, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7)
    )
  legend <- ggpubr::get_legend(legend_plot)
  final_plot <- (wrap_plots(plots, ncol = ncol) / patchwork::wrap_elements(full = legend)) +
    plot_layout(heights = c(10, 0.75)) &
    theme(plot.margin = margin(0,2,2,0, unit="mm"),
          panel.spacing.x = unit(0, "mm"),
          strip.text.x = element_text(margin = margin(b = 0)),
          strip.text.y = element_text(margin = margin(r = 0)))
  if (!is.null(output_file)) {
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      save_path <- file.path(output_dir, output_file)
    } else {
      save_path <- output_file
    }
    ggsave(filename = save_path, plot = final_plot, width = width, height = height, dpi = dpi, bg = "white")
    message("Plot saved to: ", save_path)
    invisible(save_path)
  } else {
    print(final_plot)
    invisible(final_plot)
  }
}

# --- RMSE Plot ---

#' Plot RMSE of mean_value by parameter settings.
#'
#' @param pop_summary Pop summary dataframe.
#' @param output_dir Output directory.
#' @param output_file Output filename.
#' @param width, height, dpi Figure parameters.
plot_RMSE <- function(pop_summary,
                      output_dir  = NULL,
                      output_file = NULL,
                      width       = 7,
                      height      = 5,
                      dpi         = 300) {
  p <- pop_summary %>%
    filter(statistic == "F") %>%
    group_by(across(all_of(FILTER_PARAMS))) %>%
    summarize(RMSE = sqrt(mean(mean_value^2, na.rm = TRUE)),
              .groups = "drop")  %>%
    ggplot(aes(y = minMapQ, x = minQ, fill = RMSE)) +
    geom_tile(linetype = "blank") +
    facet_grid(~clipC+baq, labeller = label_value) +
    scale_fill_viridis_c() +
    theme_minimal()
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (is.null(output_file)) {
      output_file <- "F_RMSE.png"
    }
    save_path <- file.path(output_dir, output_file)
    ggsave(save_path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
    message("Saved plot to ", save_path)
    invisible(save_path)
  } else {
    print(p)
    invisible(p)
  }
}

# # --- Raincloud plot of two-parameter interactions ---

#' Create a half-eye + boxplot "raincloud" for a given statistic split by two filter parameters.
#'
#' @param data Long-format tibble or summarized tibble with 'statistic', 'value', and parameter columns.
#' @param stat_name Character; name of the statistic to filter by.
#' @param param1 Character; first parameter column for the half-eye.
#' @param param2 Character; second parameter column for facetting.
#' @param maf_threshold Numeric; MAF threshold (only applied if raw 'value').
#' @param output_dir Optional path to save the plot; if NULL, prints to device.
#' @param output_file Optional filename for saving; defaults to '<stat_name>_<param1>_<param2>_raincloud.png'.
#' @return Invisibly returns the ggplot object after printing or saving.
#' @examples
#' raincloud_two_param(df_long, "F", "baq", "minQ")
#' raincloud_two_param(pop_summary, "absF", "baq", "minQ", maf_threshold = NA)
raincloud_two_param <- function(
  data,
  stat_name,
  param1,
  param2,
  maf_threshold = 0.05,
  output_dir    = NULL,
  output_file   = NULL
) {
  if (!is.na(maf_threshold)) {
    df <- shape_and_filter(data, maf_threshold) %>% 
    filter(statistic == stat_name)
  } else {
    df <- data %>% filter(statistic == stat_name)
  }
  # Build plot
  p <- ggplot(df, aes(
      x     = factor(.data[[param1]]),
      y     = value,
      fill  = factor(.data[[param1]]),
      group = .data[[param1]]
    )) +
    ggdist::stat_halfeye(
      adjust        = 0.5,
      width         = 0.6,
      .width        = 0,
      justification = -0.3,
      point_colour  = NA
    ) +
    geom_violin(
      width = 0.6,
      alpha = 0.5,
      trim  = FALSE
    ) +
    facet_wrap(as.formula(paste("~", param2)), scales = "free_y") +
    theme_minimal() +
    labs(
      title = paste("Effect of", param1, "and", param2, "on", stat_name),
      x     = param1,
      y     = stat_name
    ) +
    theme(
      legend.position   = "none",
      strip.background  = element_blank(),
      strip.text        = element_text(size = 12)
    )

  # Save or print
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (is.null(output_file)) {
      output_file <- paste0(stat_name, "_", param1, "_", param2, "_raincloud.png")
    }
    save_path <- file.path(output_dir, output_file)
    ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
    message("Saved plot to ", save_path)
  } else {
    print(p)
  }
}

#' Plot heatmaps of RDA/varpart results by filter parameters
#'
#' @param df         Data frame (must have columns: baq, C, minQ, minMapQ, ct, etc.)
#' @param fill_var   Column name (string) to use for fill
#' @param plot_type  "core" or "extended"
#' @param legend_title Optional string for the legend title (defaults to fill_var)
#' @param output_dir Optional directory to save plot
#' @param output_file Optional output filename
#' @param ...        Additional arguments passed to ggplot2::labs() or theme()
#' @return           Invisibly, the ggplot object
#' @export
plot_rda_heatmap <- function(
    df,
    fill_var,
    plot_type    = c("core", "extended"),
    legend_title = NULL,
    output_dir   = NULL,
    output_file  = NULL,
    ...
) {
  plot_type <- match.arg(plot_type)
  # Defensive: make sure relevant columns exist and are factors for plotting
  df <- df %>%
    dplyr::mutate(
      baq        = factor(baq, levels = 0:2, labels = paste0("baq: ", 0:2)),
      minQ       = factor(minQ),
      minMapQ    = factor(minMapQ),
      call_thresh= factor(ct),
      clipC      = factor(C),
      C_ct       = paste0("c: ", clipC, " / ct: ", call_thresh)
    )
  
  if (!fill_var %in% names(df)) stop(sprintf("Column '%s' not found in dataframe.", fill_var))
  
  if (plot_type == "core") {
    df_plot <- df %>% dplyr::filter(minMapQ != "50", clipC %in% c("0", "50"))
    p <- ggplot(df_plot, aes(x = minQ, y = minMapQ, fill = .data[[fill_var]])) +
      geom_tile(linetype = "blank") +
      facet_grid(baq ~ C_ct, labeller = label_value)
  } else if (plot_type == "extended") {
    df_plot <- df %>% dplyr::filter(baq == "baq: 0" & clipC %in% c(0, 100, 60, 75))
    p <- ggplot(df_plot, aes(x = minQ, y = minMapQ, fill = .data[[fill_var]])) +
      geom_tile(linetype = "blank") +
      facet_grid(clipC ~ call_thresh, labeller = label_value)
  }
  
  if (is.null(legend_title)) legend_title <- fill_var
  
  # Allow user to override or add labels/theme options via ...
  gg_extra <- list(...)
  labs_args <- gg_extra[names(gg_extra) %in% names(formals(ggplot2::labs))]
  theme_args <- gg_extra[names(gg_extra) %in% names(formals(ggplot2::theme))]
  
  p <- p +
    scale_fill_viridis_c(name = legend_title) +
    labs(x = "minQ", y = "minMapQ", !!!labs_args) +
    theme_minimal() +
    do.call(theme, c(
      list(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing = unit(0, "pt"),
        strip.placement = "outside",
        strip.text.x.top = element_text(angle = 45),
        strip.text.y.right = element_text(angle = 0)
      ),
      theme_args
    ))
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (is.null(output_file)) {
      output_file <- paste0("rda_heatmap_", fill_var, "_", plot_type, ".png")
    }
    save_path <- file.path(output_dir, output_file)
    ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
    message("Saved plot to ", save_path)
  } else {
    print(p)
  }
  invisible(p)
}
# --- End of Script ---
