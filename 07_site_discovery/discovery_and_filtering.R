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

#' Calculate and plot AUC-ROCs for multiple group/statistic sets
#'
#' This function calculates ROC curves and AUC values for multiple sets of group definitions and statistics.
#' For each group definition, it constructs positive ("good") and negative ("bad") subsets of the data using filter expressions,
#' computes ROC and AUC for each statistic, summarizes thresholds and performance, and optionally saves ROC plots to disk.
#'
#' @param dat Data frame. Input data containing all relevant variables for filtering and statistics.
#' @param group_defs Named list. Each element defines a group set (e.g., "set1", "set2"), with elements \code{bad} and \code{good}
#'   as filter expressions (strings) for the negative and positive classes, respectively.
#'   Example: \code{list(myset = list(bad = "X < 1", good = "X >= 1"))}
#' @param stats Character vector. Names of the columns in \code{dat} to use as test statistics for ROC/AUC.
#' @param output_dir Character (optional). If supplied, ROC plots are saved to this directory (created if needed). If NULL, plots are printed to screen.
#' @param output_file Character vector or list (optional). Named vector/list of output file names for plots, corresponding to names in \code{group_defs}.
#'   If not specified, defaults to \code{"ROC_<setname>.png"}.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{summary}{A data frame summarizing, for each group and statistic: AUC, optimal threshold, sensitivity, specificity, and group medians.}
#'     \item{plots}{A named list of \code{ggplot} ROC curve objects, one per group set.}
#'   }
#' @details
#' Requires the \pkg{pROC}, \pkg{dplyr}, \pkg{purrr}, and \pkg{ggplot2} packages. The function expects
#' that \code{dat} contains all columns used in \code{stats} and in the filter expressions.
#' "Bad" and "good" group filter expressions are evaluated with \code{rlang::parse_expr()}.
#'
#' If \code{output_dir} is supplied, plots are saved and not printed.
#' If \code{output_file} is provided, it should be a named vector/list mapping each set name to a file name.
#' 
#' @export
do_auc <- function(
    dat, 
    group_defs,   # named list: each element is list(bad=filter_expr, good=filter_expr)
    stats, 
    output_dir = NULL,       # if non-NULL, saves plots to this dir and doesn't print them
    output_file = NULL       # named character vector/list of output file names
) {
  library(pROC)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  
  all_results <- list()
  all_plots <- list()
  
  for (setname in names(group_defs)) {
    defs <- group_defs[[setname]]
    df_bad  <- dat %>% filter(!!rlang::parse_expr(defs$bad))
    df_good <- dat %>% filter(!!rlang::parse_expr(defs$good))
    
    df_bad$group <- 0 # bad (negative class)
    df_good$group <- 1 # good (positive class)
    snps <- bind_rows(df_bad, df_good)
    
    # Skip if not enough data
    if (nrow(snps) == 0 || any(sapply(stats, function(s) all(is.na(snps[[s]]))))) next
    
    # ROC objects, tidy ROC curves, and summary table
    roc_list <- lapply(stats, function(stat) roc(snps$group, snps[[stat]], quiet=TRUE))
    names(roc_list) <- stats
    
    roc_df <- imap_dfr(roc_list, function(rocobj, stat) {
      data.frame(
        fpr = 1 - rocobj$specificities,
        tpr = rocobj$sensitivities,
        stat = stat,
        auc = as.numeric(auc(rocobj)),
        set = setname
      )
    })
    
    summary_df <- map_dfr(stats, function(stat) {
      rocobj <- roc_list[[stat]]
      best <- coords(rocobj, "best", ret = c("threshold", "sensitivity", "specificity"))
      # Calculate medians for each group
      med_good <- median(snps[[stat]][snps$group == 1], na.rm = TRUE)
      med_bad  <- median(snps[[stat]][snps$group == 0], na.rm = TRUE)
      tibble(
        set = setname,
        stat = stat,
        auc = as.numeric(auc(rocobj)),
        threshold = as.numeric(best["threshold"]),
        sensitivity = as.numeric(best["sensitivity"]),
        specificity = as.numeric(best["specificity"]),
        median_good = med_good,
        median_bad  = med_bad
      )
    })
    
    # ROC plot (multi-panel)
    panel_titles <- summary_df %>% 
      mutate(stat_auc = paste0(stat, " (AUC = ", round(auc, 3), ")"))
    names(panel_titles$stat_auc) <- panel_titles$stat
    p <- ggplot(roc_df, aes(x = fpr, y = tpr, color = stat)) +
      geom_line(size = 1) +
      geom_abline(linetype = "dashed", color = "grey") +
      facet_wrap(~ stat, ncol = 2, labeller = as_labeller(panel_titles$stat_auc)) +
      labs(
        x = "False Positive Rate (1 - Specificity)",
        y = "True Positive Rate (Sensitivity)",
        title = paste0("ROC Curves for Each Statistic: ", setname)
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      # Choose the filename
      out_file <- if (!is.null(output_file) && !is.null(output_file[[setname]])) {
        file.path(output_dir, output_file[[setname]])
      } else {
        file.path(output_dir, paste0("ROC_", setname, ".png"))
      }
      ggsave(out_file, plot = p, width = 9, height = 6, dpi = 300)
      message("Saved: ", out_file)
      # Don't print!
    } else {
      print(p)
    }
    
    all_results[[setname]] <- summary_df
    all_plots[[setname]] <- p
  }
  
  # Combine all summary tables
  summary_all <- bind_rows(all_results)
  
  return(list(summary = summary_all, plots = all_plots))
}


#' Analyze pairwise PC dissimilarities vs. geography/library using one covariance file.
#'
#' @param base.name Character: Basename for analysis (e.g. "northern_domain_filtered.ct5")
#' @param sample_list Character: Path to sample metadata file to use for all analyses
#' @param meta Dataframe: Must include columns bam_code, pop_code, region, library, latitude, longitude
#' @param n_pcs Integer: Number of PCs to analyze (default 2)
#' @return Dataframe with base.name, PC, predictor, t_value, p_value
analyze_pca_distance_predictors <- function(
    base.name,
    sample_list,
    meta,
    n_pcs = 2
) {
  require(dplyr)
  require(geosphere)
  
  # Build covariance file path
  cov_file <- paste0(base.name, ".Pcangsd.cov")
  
  # If file doesn't exist, skip and return NULL
  if (!file.exists(cov_file)) return(NULL)
  if (!file.exists(sample_list)) stop("Sample list file not found: ", sample_list)
  
  samples <- read.table(sample_list)
  colnames(samples)[1] <- "bam_code"
  samples <- left_join(samples, meta[,c("bam_code", "pop_code", "region", "library", "longitude", "latitude")])
  
  cov <- read.table(cov_file)
  colnames(cov) <- samples$bam_code
  rownames(cov) <- samples$bam_code
  eig <- eigen(cov)
  pc_mat <- eig$vectors[, 1:n_pcs, drop = FALSE]
  pc_names <- paste0("PC", 1:n_pcs)
  pca_df <- cbind(
    samples[, c("bam_code", "pop_code", "region", "library", "latitude", "longitude")],
    as.data.frame(pc_mat)
  )
  colnames(pca_df)[(ncol(pca_df)-n_pcs+1):ncol(pca_df)] <- pc_names
  n <- nrow(pca_df)
  pairs <- expand.grid(i = 1:n, j = 1:n) %>%
    filter(i < j) %>%
    mutate(
      bam_code_i = pca_df$bam_code[i],
      bam_code_j = pca_df$bam_code[j],
      lib_i = pca_df$library[i],
      lib_j = pca_df$library[j],
      lon_i = pca_df$longitude[i],
      lat_i = pca_df$latitude[i],
      lon_j = pca_df$longitude[j],
      lat_j = pca_df$latitude[j],
      geo_dist_km = distGeo(matrix(c(lon_i, lat_i), ncol = 2),
                            matrix(c(lon_j, lat_j), ncol = 2)) / 1000,
      diff_library = as.integer(lib_i != lib_j)
    )
  # Add pairwise distances for each PC
  for (k in 1:n_pcs) {
    pairs[[paste0("PC", k, "_i")]] <- pca_df[[paste0("PC", k)]][pairs$i]
    pairs[[paste0("PC", k, "_j")]] <- pca_df[[paste0("PC", k)]][pairs$j]
    pairs[[paste0("PC", k, "_dist")]] <- abs(pairs[[paste0("PC", k, "_i")]] - pairs[[paste0("PC", k, "_j")]])
  }
  get_t_p <- function(model, predictor) {
    summ <- summary(model)$coefficients
    if (predictor %in% rownames(summ)) {
      tval <- as.numeric(summ[predictor, "t value"])
      pval <- as.numeric(summ[predictor, "Pr(>|t|)"])
    } else {
      tval <- NA; pval <- NA
    }
    c(t_value = tval, p_value = pval)
  }
  # Gather results for all PCs/predictors
  out <- data.frame()
  for (k in 1:n_pcs) {
    pc_dist_col <- paste0("PC", k, "_dist")
    model <- lm(pairs[[pc_dist_col]] ~ geo_dist_km + diff_library, data = pairs)
    for (pred in c("geo_dist_km", "diff_library")) {
      res <- get_t_p(model, pred)
      out <- rbind(out, data.frame(
        base.name = base.name,
        PC = paste0("PC", k),
        predictor = pred,
        t_value = as.numeric(res["t_value"]),
        p_value = as.numeric(res["p_value"])
      ))
    }
  }
  rownames(out) <- NULL
  out
}

#' Identify and summarize selection candidates from extended fastPCA and pcadapt models (PCAngsd)
#'
#' Loads SNP site and selection tables, calculates p-values, and returns a summary table
#' with p-values for each site (plus BH-adjusted p-values). The pcadapt z-scores can be summarized per-PC ("univariate" mode)
#' or as a multivariate outlier statistic ("multivariate" mode). The fastPCA selection statistics
#' are always processed per component.
#'
#' @param basename Character. Basename for file naming, used if \code{files} is not specified.
#' @param K Integer. Number of principal components (PCs) to analyze.
#' @param mode Character. Either \code{"univariate"} (default; analyzes each PC separately)
#'   or \code{"multivariate"} (analyzes all PCs jointly; requires \code{K} >= 2).
#' @param files List (optional). A named list allowing override of file paths. Entries may include:
#'   \describe{
#'     \item{sites}{Path to snpcodes file (default: \code{paste0(basename, "_snpcodes.txt")})}
#'     \item{pca_sites}{Path to PCAngsd .sites file (default: \code{paste0(basename, ".Pcangsd.sites")})}
#'     \item{pcadapt_zscores}{Path to pcadapt z-scores file (default: \code{paste0(basename, ".Pcangsd.pcadapt.zscores")})}
#'     \item{selection}{Path to PCAngsd selection file (default: \code{paste0(basename, ".Pcangsd.selection")})}
#'   }
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{snpcode}{SNP identifier (scaffold_position).}
#'     \item{pcadapt1, ..., pcadaptK}{(Univariate mode) P-value for each PC from the extended pcadapt test, computed from squared z-score deviation from the median for each component.}
#'     \item{pcadapt1}{(Multivariate mode) Single multivariate p-value per SNP from the robust Mahalanobis distance across PCs (requires \code{K} >= 2).}
#'     \item{selection1, ..., selectionK}{P-values from the extended fastPCA selection statistic, for each PC (always univariate).}
#'   }
#' @details
#' All input tables must contain the same SNPs in the same order.  
#' In \code{mode = "univariate"}, each principal component is processed separately;  
#' in \code{mode = "multivariate"}, a single robust Mahalanobis distance is computed across all \code{K} components.
#'
#' @importFrom bigutilsr dist_ogk
#' @export
pca_selection <- function(
    basename, K,
    mode = "univariate",
    files = list(
      sites = NULL,
      pca_sites = NULL,
      pcadapt_zscores = NULL,
      selection = NULL
    )
) {
  require(bigutilsr)

  # Set file paths, allowing override from files argument
  f_sites <- if (!is.null(files$sites)) files$sites else paste0(basename, "_snpcodes.txt")
  f_pca_sites <- if (!is.null(files$pca_sites)) files$pca_sites else paste0(basename, ".Pcangsd.sites")
  f_pcadapt_z <- if (!is.null(files$pcadapt_zscores)) files$pcadapt_zscores else paste0(basename, ".Pcangsd.pcadapt.zscores")
  f_selection <- if (!is.null(files$selection)) files$selection else paste0(basename, ".Pcangsd.selection")

  # Read input files
  sites <- read.table(f_sites, header = FALSE, stringsAsFactors = FALSE)
  colnames(sites) <- "snpcode"
  pca_sites <- read.table(f_pca_sites)
  colnames(pca_sites) <- "index"
  rownames(pca_sites) <- sites$snpcode
  site.codes <- rownames(subset(pca_sites, index == 1))

  pcadapt_zscores <- read.table(f_pcadapt_z, header = FALSE)
  rownames(pcadapt_zscores) <- site.codes

  # --- Process pcadapt z-scores according to mode ---
  if (mode == "univariate") {
    # For each PC, convert z to squared deviation from median, then to p-value
    d2 <- as.data.frame(apply(pcadapt_zscores, 2, function(z) (z - median(z))^2))
    pvals <- as.data.frame(sapply(seq_len(K), function(k) pchisq(d2[, k], 1, lower.tail = FALSE)))
    colnames(pvals) <- paste0("pcadapt", seq_len(K))
  } else if (mode == "multivariate") {
    if (K < 2) {
      stop("At least two principal components (K >= 2) are required for multivariate analysis.")
    }
    d2 <- as.vector(dist_ogk(as.matrix(pcadapt_zscores[, 1:K])))
    pvals <- data.frame(pcadapt1 = pchisq(d2, K, lower.tail = FALSE))
  } else {
    stop("mode must be 'univariate' or 'multivariate'")
  }

  # --- Handle fastPCA-style selection statistics (always univariate) ---
  selection <- read.table(f_selection, header = FALSE)
  rownames(selection) <- site.codes
  sel_pvals <- as.data.frame(sapply(seq_len(K), function(k) pchisq(selection[, k], 1, lower.tail = FALSE)))
  colnames(sel_pvals) <- paste0("selection", seq_len(K))

  # --- Combine results ---
  results <- data.frame(
    snpcode = site.codes,
    pvals,
    sel_pvals,
    stringsAsFactors = FALSE
  )

  # --- Add BH-adjusted p-values ---
  for (col in grep("^pcadapt|^selection", colnames(results), value = TRUE)) {
    results[[paste0(col, ".BH")]] <- p.adjust(results[[col]], method = "BH")
  }

  return(results)
}

#' PCA biplot with flexible palettes and multi-panel support
#'
#' Produces biplots from a covariance matrix and sample info, colored by sample metadata.
#' Handles flexible file naming, discrete or continuous fill variables, and multiple panels.
#' Accepts color palettes as character vectors (e.g., from \code{paletteer_d} or \code{paletteer_c}).
#' Auto-wraps palette vectors as lists if needed for single-panel plotting.
#'
#' @param cov_path Character. Path (or prefix) to covariance file (e.g. \code{"mysamples.Pcangsd.cov"}). Tries \code{.Pcangsd.cov} and \code{.cov} suffixes if needed.
#' @param sample_list_path Character. Path to sample list file (passed to \code{get.info}).
#' @param PCs Integer vector. Principal components to plot (default \code{c(1,2)}).
#' @param fill_var Character vector. Metadata columns in sample info to color points by (e.g. \code{"library"}, \code{"region"}, \code{"latitude"}).
#' @param fill_pal List or vector of palettes for coloring points (optional). Each entry can be a paletteer vector or custom color vector, e.g. \code{paletteer_d("ggsci::default_jama", 3)}. Auto-wrapped in a list if needed. Also accepts "Set1", "Set2", and "Set3" as presents; these are the Rcolorbrewer palettes of the same name but with yellow removed.
#' @param point_args List. Named arguments to \code{geom_point} (e.g. \code{list(size=2, alpha=0.8)}).
#' @param plot_title Character. Optional title for the plot.
#' @param output_file Character. Optional filename to save PNG (if provided).
#' @param output_dir Character. Output directory for PNGs (default current directory).
#'
#' @return A \code{ggplot} object, or a named list of ggplot objects if multiple PC pairs/panels.
#'
#' @details
#' - Color palette length must match the number of unique levels in the fill variable if discrete.
#' - Pass continuous palettes using \code{paletteer_c(..., n)}.
#' - For multi-panel plots, supply \code{fill_pal} as a list of palettes matching \code{fill_var}.
#' - If \code{fill_pal} is a vector, it is auto-wrapped in a list for you.
#'
#' @examples
#' # Discrete fill (library)
#' plot_pcas(
#'   "mysamples.Pcangsd.cov",
#'   "mysamples_sample_list.txt",
#'   PCs = c(1,2),
#'   fill_var = "library",
#'   fill_pal = paletteer_d("ggsci::default_jama", 3),
#'   point_args = list(size = 2, alpha = 0.85)
#' )
#'
#' # Continuous fill (latitude)
#' plot_pcas(
#'   "mysamples.Pcangsd.cov",
#'   "mysamples_sample_list.txt",
#'   PCs = c(1,2),
#'   fill_var = "latitude",
#'   fill_pal = paletteer_c("viridis::plasma", 30),
#'   point_args = list(size = 1.5, alpha = 0.8)
#' )
#'
#' # Multi-panel: color by library and region with custom palettes
#' plot_pcas(
#'   "mysamples.Pcangsd.cov",
#'   "mysamples_sample_list.txt",
#'   PCs = c(1,2),
#'   fill_var = c("library", "region"),
#'   fill_pal = list(
#'     paletteer_d("ggsci::default_jama", 3),
#'     paletteer_d("nord::aurora", 3)
#'   )
#' )
#'
#' @export
#' PCA biplot with flexible palettes and multi-panel support
#'
#' Produces biplots from a covariance matrix and sample info, colored by sample metadata.
#' Handles flexible file naming, discrete or continuous fill variables, and multiple panels.
#' Accepts color palettes as character vectors (e.g., from \code{paletteer_d} or \code{paletteer_c}).
#' Auto-wraps palette vectors as lists if needed for single-panel plotting.
#'
#' @param cov_path Character. Path (or prefix) to covariance file (e.g. \code{"mysamples.Pcangsd.cov"}). Tries \code{.Pcangsd.cov} and \code{.cov} suffixes if needed.
#' @param sample_list_path Character. Path to sample list file (passed to \code{get.info}).
#' @param PCs Integer vector. Principal components to plot (default \code{c(1,2)}).
#' @param fill_var Character vector. Metadata columns in sample info to color points by (e.g. \code{"library"}, \code{"region"}, \code{"latitude"}).
#' @param fill_pal List or vector of palettes for coloring points (optional). Each entry can be a paletteer vector or custom color vector, e.g. \code{paletteer_d("ggsci::default_jama", 3)}. Auto-wrapped in a list if needed. Also accepts "Set1", "Set2", and "Set3" as presents; these are the Rcolorbrewer palettes of the same name but with yellow removed.
#' @param point_args List. Named arguments to \code{geom_point} (e.g. \code{list(size=2, alpha=0.8)}).
#' @param plot_title Character. Optional title for the plot.
#' @param output_file Character. Optional filename to save PNG (if provided).
#' @param output_dir Character. Output directory for PNGs (default current directory).
#'
#' @return A \code{ggplot} object, or a named list of ggplot objects if multiple PC pairs/panels.
#'
#' @details
#' - Color palette length must match the number of unique levels in the fill variable if discrete.
#' - Pass continuous palettes using \code{paletteer_c(..., n)}.
#' - For multi-panel plots, supply \code{fill_pal} as a list of palettes matching \code{fill_var}.
#' - If \code{fill_pal} is a vector, it is auto-wrapped in a list for you.
#'
#' @examples
#' # Discrete fill (library)
#' plot_pcas(
#'   "mysamples.Pcangsd.cov",
#'   "mysamples_sample_list.txt",
#'   PCs = c(1,2),
#'   fill_var = "library",
#'   fill_pal = paletteer_d("ggsci::default_jama", 3),
#'   point_args = list(size = 2, alpha = 0.85)
#' )
#'
#' # Continuous fill (latitude)
#' plot_pcas(
#'   "mysamples.Pcangsd.cov",
#'   "mysamples_sample_list.txt",
#'   PCs = c(1,2),
#'   fill_var = "latitude",
#'   fill_pal = paletteer_c("viridis::plasma", 30),
#'   point_args = list(size = 1.5, alpha = 0.8)
#' )
#'
#' # Multi-panel: color by library and region with custom palettes
#' plot_pcas(
#'   "mysamples.Pcangsd.cov",
#'   "mysamples_sample_list.txt",
#'   PCs = c(1,2),
#'   fill_var = c("library", "region"),
#'   fill_pal = list(
#'     paletteer_d("ggsci::default_jama", 3),
#'     paletteer_d("nord::aurora", 3)
#'   )
#' )
#'
#' @export
plot_pcas <- function(
    cov_path,
    sample_list_path,
    PCs = c(1,2),
    fill_var = "library",
    fill_pal = NULL,
    point_args = list(size = 2, alpha = 0.9),
    plot_title = NULL,
    output_file = NULL,
    output_dir = "."
) {
  require(ggplot2)
  require(dplyr)
  require(RColorBrewer)
  require(paletteer)

  # Handle covariance file naming
  test_paths <- c(
    cov_path,
    paste0(cov_path, ".Pcangsd.cov"),
    paste0(cov_path, ".cov")
  )
  cov_file <- NULL
  for (p in test_paths) {
    if (file.exists(p)) { cov_file <- p; break }
  }
  if (is.null(cov_file)) stop("Covariance file not found with any naming convention: ", paste(test_paths, collapse=", "))

  # Read sample and covariance
  samples <- get.info(sample_list_path, "sequenced_samples_metadata.csv")
  cov <- read.table(cov_file)
  colnames(cov) <- samples$bam_code
  rownames(cov) <- samples$bam_code

  # PCA
  eig <- eigen(cov)
  n_samples <- nrow(cov)
  max_pc <- min(n_samples, ncol(cov))
  n_pc_to_plot <- max(PCs)
  if (n_pc_to_plot > max_pc) stop("Requested PCs exceed the number of possible PCs.")

  # Make PC matrix
  pca_df <- cbind(
    samples,
    as.data.frame(eig$vectors[, seq_len(n_pc_to_plot), drop=FALSE])
  )
  colnames(pca_df)[(ncol(pca_df)-n_pc_to_plot+1):ncol(pca_df)] <- paste0("PC", 1:n_pc_to_plot)

  # Chunk PCs into pairs
  pc_pairs <- split(PCs, ceiling(seq_along(PCs)/2))
  if (length(pc_pairs[[length(pc_pairs)]]) == 1) {
    warning("Last PC chunk only has 1 PC. Skipping it.")
    pc_pairs <- pc_pairs[sapply(pc_pairs, length) == 2]
  }
  if (length(pc_pairs) == 0) stop("No valid PC pairs to plot.")

  # Palette presets for when fill_pal is NULL only
  presets <- list(
    Set1 = setdiff(brewer.pal(9, "Set1"), "#FFFF33"),
    Set2 = brewer.pal(8, "Set2"),
    Set3 = setdiff(brewer.pal(12, "Set3"), c("#FFFFB3","#FFED6F"))
  )

  # --- Palette logic ---
  # Make fill_pal a list with one palette per fill_var if it isn't already
  if (!is.null(fill_pal) && !is.list(fill_pal)) {
    fill_pal <- rep(list(fill_pal), length(fill_var))
  }
  # For NULL, fill with NULLs
  if (is.null(fill_pal)) {
    fill_pal <- rep(list(NULL), length(fill_var))
  }

  # Helper for palette selection
  get_palette <- function(var, pal) {
    is_discrete <- !is.numeric(pca_df[[var]])
    # If pal is a preset name, use that
    if (!is.null(pal) && is.character(pal) && length(pal) == 1 && pal %in% names(presets)) {
      return(presets[[pal]])
    }
    # If pal is a character vector (multiple colors), use as is
    if (!is.null(pal) && (is.character(pal) || is.numeric(pal))) {
      return(as.character(pal))
    }
    # Defaults
    if (is_discrete) {
      return(presets$Set3)
    } else {
      return(as.character(paletteer_c("grDevices::Viridis", 30)))
    }
  }

  # Plotting
  all_plots <- list()
  for (pc_idx in seq_along(pc_pairs)) {
    pcs <- pc_pairs[[pc_idx]]
    x_pc <- paste0("PC", pcs[1])
    y_pc <- paste0("PC", pcs[2])
    plot_list <- list()
    for (j in seq_along(fill_var)) {
      var <- fill_var[j]
      pal_vec <- get_palette(var, fill_pal[[j]])
      is_discrete <- !is.numeric(pca_df[[var]])
      mapping <- aes_string(x = x_pc, y = y_pc, col = var)

      # Theme tweaks: only set legend key sizes for discrete
      legend_theme <- if (is_discrete) {
        theme(
          legend.key.height = unit(6, "pt"),
          legend.key.width = unit(6, "pt")
        )
      } else {
        theme()
      }

      g <- ggplot(pca_df, mapping) +
        do.call(geom_point, c(list(), point_args)) +
        theme_minimal(base_size = 8) +
        theme(
          axis.title = element_text(size=8),
          legend.title = element_text(size=8),
          axis.text = element_text(size=6),
          legend.position = "right",
          legend.justification = c("left", "center"),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.spacing.y = unit(0, "pt"),
          panel.spacing = margin(0, 3, 0, 3, unit = "mm"),
          legend.text = element_text(size = 6, margin = margin(l = 2, r = 2)),
          plot.title = element_text(size=10),
          plot.margin = margin(0, 1, 0, 1, unit = "mm")
        ) +
        legend_theme  # Add discrete-specific theme tweaks

      # Add discrete guide if needed
      if (is_discrete) {
        g <- g + guides(colour = guide_legend(override.aes = list(size=1.5)))
      }

      if (!is.null(plot_title) && j == 1) {
        g <- g + ggtitle(plot_title)
      }
      if (!is_discrete) {
        g <- g + scale_color_gradientn(colors = as.character(pal_vec))
      } else {
        # Always recalculate the palette to the number of unique groups
        group_levels <- sort(unique(pca_df[[var]]))
        n_groups <- length(group_levels)
        pal_final <- as.character(pal_vec)
        if (length(pal_final) < n_groups) {
          stop(sprintf("Palette has %d colors but %d are needed for '%s'", length(pal_final), n_groups, var))
        }
        # Ensure palette matches group levels order
        names(pal_final) <- group_levels
        g <- g + scale_color_manual(values = pal_final)
      }
      plot_list[[j]] <- g
    }
    # Multipanel arrangement if multiple fill_vars
    if (length(plot_list) == 1) {
      final_plot <- plot_list[[1]]
    } else {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        library(patchwork)
        final_plot <- wrap_plots(plot_list, ncol = min(3, length(plot_list)))
      } else if (requireNamespace("cowplot", quietly = TRUE)) {
        library(cowplot)
        final_plot <- plot_grid(plotlist = plot_list, ncol = min(3, length(plot_list)), labels = NULL)
      } else {
        warning("Install 'patchwork' or 'cowplot' for multi-panel support. Returning first plot only.")
        final_plot <- plot_list[[1]]
      }
    }
    all_plots[[paste0("PCs", pcs[1], "and", pcs[2])]] <- final_plot

    # Save output if requested
    if (!is.null(output_file)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      save_file <- if (length(pc_pairs) == 1) {
        file.path(output_dir, output_file)
      } else {
        file.path(output_dir, paste0(sub("\\.png$","",output_file), "_PCs", pcs[1], "and", pcs[2], ".png"))
      }
      ggsave(save_file, plot = final_plot, width = 5, height = 4, dpi = 300, bg = "white")
    }
  }

  if (!is.null(output_file)) {
    return(invisible(NULL))
  }
  if (length(all_plots) == 1) return(all_plots[[1]])
  return(all_plots)
}


#' Create and (optionally) save an UpSet plot
#'
#' Generates an UpSet plot (visualizing set intersections) from a named list of character vectors (each representing a group).
#' If \code{output_file} is supplied, the plot is saved as a PNG in \code{output_dir} and nothing is rendered to the screen.
#' If \code{output_file} is \code{NULL}, the ggplot object is returned for further display or modification.
#'
#' @param upset.list A named list of character vectors, each vector giving the members of a set/group.
#' @param output_dir Directory to save the plot if \code{output_file} is provided (default: current directory).
#' @param output_file Optional file name for PNG output (do not include directory; see \code{output_dir}).
#'
#' @return If \code{output_file} is \code{NULL}, returns a ggplot object representing the UpSet plot. Otherwise, saves the plot to file and returns \code{NULL} (invisibly).
#'
#' @import UpSetR
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggsave
#' @examples
#' \dontrun{
#' # Example list of sets
#' sets <- list(
#'   groupA = c("gene1", "gene2", "gene3"),
#'   groupB = c("gene2", "gene3", "gene4"),
#'   groupC = c("gene3", "gene4", "gene5")
#' )
#' # Display plot in RStudio
#' plot.Upset(sets)
#' # Save plot to file
#' plot.Upset(sets, output_dir = "plots", output_file = "example_upset")
#' }
plot.Upset <- function(upset.list,
                       output_dir = ".",
                       output_file = NULL) {
  require(UpSetR)
  f <- fromList(upset.list)
  u <- upset(f,
             sets = colnames(f),
             order.by = "freq",
             nintersects = 15,
             text.scale = 1,
             point.size = 2.0,
             keep.order = TRUE,
             decreasing = TRUE)
  
  uplot <- plot_grid(u$Main_bar,
                     NULL,
                     u$Matrix,
                     nrow = 3,
                     align = 'v',
                     rel_heights = c(3, 0.25, 2),
                     rel_widths = c(4, 0, 4))
  
  if (!is.null(output_file)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    save_file <- file.path(output_dir, paste0(sub("\\.png$","",output_file), ".png"))
    ggsave(filename = save_file, plot = uplot, width = 4, height = 5, dpi = 300, bg = "white")
    return(invisible(NULL))
  } else {
    return(uplot)
  }
}

#' Summarize Depth by Genotype Class with Optional Standardization
#'
#' Reads in per-genotype depth tables and site summary statistics, applies several standardization procedures,
#' and summarizes per-SNP, per-genotype statistics (mean, median, IQR, etc.), with optional differences
#' between genotype classes.
#'
#' @param snpstats_path Character. Path to SNP statistics table (must have "snpcode" and "MAF" columns).
#' @param alt_path      Character. Path to "homoAlt" genotype matrix (samples in columns, SNPs in rows).
#' @param ref_path      Character. Path to "homoRef" genotype matrix (same structure as above).
#' @param het_path      Character. Path to "hetero" genotype matrix (same structure as above).
#' @param standards     Character vector. Which standardizations to perform. 
#'                      Options: "none", "sample_z", "loci_z", "double_loci_sample", "double_sample_loci",
#'
#' @return
#' Tibble: one row per SNP. Columns include per-standardization summaries and difference columns
#'         for each genotype class and statistic.
#'
#' @examples
#' \dontrun{
#' summary_df <- depth_by_genotype(
#'   snpstats_path = "southern_site_summary_maf05.tsv.gz",
#'   alt_path = "batch_effect/southern.homoAlt.maf05.tsv",
#'   ref_path = "batch_effect/southern.homoRef.maf05.tsv",
#'   het_path = "batch_effect/southern.hetero.maf05.tsv"
#' )
#' }
#' @import dplyr data.table tidyr matrixStats
#' @export
depth_by_genotype <- function(
    snpstats_path,
    alt_path,
    ref_path,
    het_path,
    standards = c(
      "none", "sample_z", "loci_z", "double_loci_sample", "double_sample_loci", "aggregate"
    )
) {
  cat("=== depth_by_genotype() DEBUG RUN ===\n")

  # --- Load packages ---
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(matrixStats)
  library(diptest)

  # --- Standardization map ---
  std_map <- list(
    none               = "raw",
    sample_z           = "sz",
    loci_z             = "lz",
    double_loci_sample = "dlsz",
    double_sample_loci = "dslz",
    aggregate          = "aggregate"
  )
  allowed_standards <- names(std_map)
  if (!all(standards %in% allowed_standards)) {
    stop("Allowed standardizations: ", paste(allowed_standards, collapse=", "))
  }
  std_names <- unname(unlist(std_map[standards]))

  # --- 1. Read & join data ---
  cat("Reading input tables...\n")
  snpstats <- fread(snpstats_path, data.table = FALSE)
  alt <- fread(alt_path, data.table = FALSE)
  ref <- fread(ref_path, data.table = FALSE)
  het <- fread(het_path, data.table = FALSE)
  cat(sprintf("Read alt: %d ref: %d het: %d\n", nrow(alt), nrow(ref), nrow(het)))

  sample_cols <- setdiff(names(alt), c("genotype", "snpcode"))
  cat("Number of sample columns:", length(sample_cols), "\n")

  # --- 2. Nonzero filtering and intersection ---
  alt$nonzero <- rowSums(as.matrix(alt[, sample_cols]) != 0 & !is.na(alt[, sample_cols]))
  ref$nonzero <- rowSums(as.matrix(ref[, sample_cols]) != 0 & !is.na(ref[, sample_cols]))
  het$nonzero <- rowSums(as.matrix(het[, sample_cols]) != 0 & !is.na(het[, sample_cols]))
  cat(sprintf("Alt with >=5 nonzero: %d\n", sum(alt$nonzero >= 5)))
  cat(sprintf("Ref with >=5 nonzero: %d\n", sum(ref$nonzero >= 5)))
  cat(sprintf("Het with >=5 nonzero: %d\n", sum(het$nonzero >= 5)))

  alt <- alt[alt$nonzero >= 5, , drop = FALSE]
  ref <- ref[ref$nonzero >= 5, , drop = FALSE]
  het <- het[het$nonzero >= 5, , drop = FALSE]

  snps_intersect <- Reduce(intersect, list(alt$snpcode, ref$snpcode, het$snpcode))
  cat(sprintf("Snps in all three after filtering: %d\n", length(snps_intersect)))
  alt <- alt[alt$snpcode %in% snps_intersect, , drop = FALSE]
  ref <- ref[ref$snpcode %in% snps_intersect, , drop = FALSE]
  het <- het[het$snpcode %in% snps_intersect, , drop = FALSE]
  cat(sprintf("Filtered alt: %d ref: %d het: %d\n", nrow(alt), nrow(ref), nrow(het)))
  cat("Unique snpcodes after filtering: alt:", length(unique(alt$snpcode)),
      "ref:", length(unique(ref$snpcode)), "het:", length(unique(het$snpcode)), "\n")

  alt$genotype <- "homoAlt"
  ref$genotype <- "homoRef"
  het$genotype <- "hetero"

  # --- 3. Combine and annotate ---
  all <- bind_rows(alt, ref, het) %>%
    left_join(snpstats[, c("snpcode", "MAF", "call_rate")], by = "snpcode") %>%
    mutate(genotype_majmin = case_when(
      genotype == "homoRef" & MAF <= 0.5 ~ "homoMaj",
      genotype == "homoRef" & MAF > 0.5  ~ "homoMin",
      genotype == "homoAlt" & MAF <= 0.5 ~ "homoMin",
      genotype == "homoAlt" & MAF > 0.5  ~ "homoMaj",
      TRUE ~ genotype
    )) %>%
    subset(call_rate > 0.4 & MAF < 0.95)

  cat(sprintf("Rows after join and annotation: %d\n", nrow(all)))
  cat("Table of genotypes:\n")
  print(table(all$genotype))
  cat("Table of genotype_majmin:\n")
  print(table(all$genotype_majmin))
  cat("Snpcodes present after join:", length(unique(all$snpcode)), "\n")
  cat("Counts per snpcode (should be 3 for all):\n")
  print(table(table(all$snpcode)))

  all <- all[, !names(all) %in% c("call_rate", "nonzero")]
  sample_cols <- setdiff(names(all), c("snpcode", "genotype", "genotype_majmin", "MAF"))
  meta_cols <- c("snpcode", "genotype_majmin")

  # --- Z-score helpers ---
  zscore_col <- function(mat) {
    col_mean <- colMeans(mat, na.rm = TRUE)
    col_sd <- matrixStats::colSds(mat, na.rm = TRUE)
    z <- sweep(sweep(mat, 2, col_mean, "-"), 2, col_sd, "/")
    z[is.nan(z)] <- NA
    z
  }
  zscore_row <- function(mat) {
    row_mean <- rowMeans(mat, na.rm = TRUE)
    row_sd <- matrixStats::rowSds(mat, na.rm = TRUE)
    z <- sweep(sweep(mat, 1, row_mean, "-"), 1, row_sd, "/")
    z[is.nan(z)] <- NA
    z
  }

  # --- All dataframes for standardizations ---
  mat <- as.matrix(all[, sample_cols])
  mat[mat == 0] <- NA
  cat("mat dims for raw:", dim(mat), "\n")
  dfs <- list()
  dfs$raw <- cbind(mat, all[, meta_cols])

  if ("sz" %in% std_names) {
    mat_sample_z <- zscore_col(mat)
    cat("mat_sample_z dims:", dim(mat_sample_z), "\n")
    colnames(mat_sample_z) <- paste0(sample_cols, "_sample_z")
    dfs$sz <- cbind(mat_sample_z, all[, meta_cols])
    rm(mat_sample_z); gc()
  }
  if ("lz" %in% std_names) {
    mat_loci_z <- zscore_row(mat)
    cat("mat_loci_z dims:", dim(mat_loci_z), "\n")
    colnames(mat_loci_z) <- paste0(sample_cols, "_loci_z")
    dfs$lz <- cbind(mat_loci_z, all[, meta_cols])
    rm(mat_loci_z); gc()
  }
  if ("dlsz" %in% std_names) {
    mat_loci_z <- zscore_row(mat)
    mat_double_loci_sample_z <- zscore_col(mat_loci_z)
    cat("mat_double_loci_sample_z dims:", dim(mat_double_loci_sample_z), "\n")
    colnames(mat_double_loci_sample_z) <- paste0(sample_cols, "_double_loci_sample")
    dfs$dlsz <- cbind(mat_double_loci_sample_z, all[, meta_cols])
    rm(mat_loci_z, mat_double_loci_sample_z); gc()
  }
  if ("dslz" %in% std_names) {
    mat_sample_z <- zscore_col(mat)
    mat_double_sample_loci_z <- zscore_row(mat_sample_z)
    cat("mat_double_sample_loci_z dims:", dim(mat_double_sample_loci_z), "\n")
    colnames(mat_double_sample_loci_z) <- paste0(sample_cols, "_double_sample_loci")
    dfs$dslz <- cbind(mat_double_sample_loci_z, all[, meta_cols])
    rm(mat_sample_z, mat_double_sample_loci_z); gc()
  }
  if ("aggregate" %in% std_names) {
    mat_aggregate <- rowsum(as.matrix(all[, sample_cols, drop=FALSE]), group = all$snpcode, na.rm=TRUE)
    cat("mat_aggregate dims:", dim(mat_aggregate), "\n")
    mat_aggregate[mat_aggregate == 0] <- NA
    colnames(mat_aggregate) <- paste0(sample_cols, "_aggregate")
    mat_aggregate <- data.frame(
      snpcode = rownames(mat_aggregate),
      genotype_majmin = NA_character_,
      mat_aggregate,
      row.names = NULL,
      check.names = FALSE
    )
    dfs$aggregate <- mat_aggregate
    rm(mat_aggregate); gc()
  }

  # --- Skewness and kurtosis helpers ---
  rowSkewness_fast <- function(mat, na.rm = TRUE) { 
    m <- rowMeans(mat, na.rm = na.rm)
    s <- matrixStats::rowSds(mat, na.rm = na.rm)
    n <- rowSums(!is.na(mat))
    m_mat <- matrix(m, nrow = length(m), ncol = ncol(mat))
    s_mat <- matrix(s, nrow = length(s), ncol = ncol(mat))
    m3 <- rowMeans((mat - m_mat)^3, na.rm = na.rm)
    skew <- m3 / (s^3)
    result <- skew * sqrt(n * (n - 1)) / (n - 2)
    result[n < 5 | !is.finite(result)] <- NA_real_
    result
  }
  rowKurtosis_fast <- function(mat, na.rm = TRUE) {
    m <- rowMeans(mat, na.rm = na.rm)
    s <- matrixStats::rowSds(mat, na.rm = na.rm)
    n <- rowSums(!is.na(mat))
    m_mat <- matrix(m, nrow = length(m), ncol = ncol(mat))
    m4 <- rowMeans((mat - m_mat)^4, na.rm = na.rm)
    kurt <- m4 / (s^4)
    result <- kurt - 3
    result[n < 5 | !is.finite(result)] <- NA_real_
    result
  }
  rowEntropy <- function(mat) {
    rowsums <- rowSums(mat, na.rm=TRUE)
    p <- mat / rowsums
    p[is.na(p)] <- 0
    logp <- log(p)
    logp[p == 0] <- 0
    ent <- -rowSums(p * logp)
    N <- rowSums(!is.na(mat))
    ent_max <- log(N) 
    out <- ent / ent_max # standardized by maximum entropy
    out[rowsums == 0 | N <= 1] <- NA_real_
    out
  }
  row_prop_extreme <- function(mat) {
    q1 <- rowQuantiles(mat, probs = 0.25, na.rm = TRUE)
    q3 <- rowQuantiles(mat, probs = 0.75, na.rm = TRUE)
    iqr_vals <- q3 - q1
    n <- rowSums(!is.na(mat))
    out <- sapply(1:nrow(mat), function(i) {
      x <- mat[i, ]; ql <- q1[i]; qh <- q3[i]; iqr <- iqr_vals[i]
      if (is.na(ql) || is.na(qh) || is.na(iqr)) return(NA)
      extreme <- (x < (ql - 1.5*iqr)) | (x > (qh + 1.5*iqr))
      mean(extreme, na.rm = TRUE)
    })
    out[n < 5 | !is.finite(out)] <- NA_real_
    out
  }
  # --- Row summary function (by locus) ---
  row_summaries_fast <- function(mat, 
                                 prefix, 
                                 which_stats = c("mean", "median", "sd", "mad", "cv", "iqr", "kurtosis", "skewness", "entropy", "extreme", "dip")) {
    out <- list()
    for (stat in which_stats) {
      if (stat == "mean")     out$mean <- rowMeans(mat, na.rm = TRUE)
      if (stat == "median")   out$median <- rowMedians(mat, na.rm = TRUE)
      if (stat == "sd")       out$sd <- rowSds(mat, na.rm = TRUE)
      if (stat == "mad")      out$mad <- rowMads(mat, na.rm = TRUE)
      if (stat == "cv")       out$cv <- rowSds(mat, na.rm = TRUE)/rowMeans(mat, na.rm = TRUE)
      if (stat == "iqr")      out$iqr <- rowIQRs(mat, na.rm = TRUE)
      if (stat == "kurtosis") out$kurtosis <- rowKurtosis_fast(mat, na.rm = TRUE)
      if (stat == "skewness") out$skewness <- rowSkewness_fast(mat, na.rm = TRUE)
      if (stat == "entropy")  out$entropy <- rowEntropy(mat)
      if (stat == "extreme")  out$extreme <- row_prop_extreme(mat)
      if (stat == "dip") {
        out$dip <- apply(mat, 1, function(x) {
          x <- x[is.finite(x) & !is.na(x)]
          if (length(x) < 5) return(NA_real_)
          suppressWarnings(diptest::dip(x))
        })
      }
    }
    # Mask all non-finite values (Inf, -Inf, NaN) as NA in all columns
    out <- lapply(out, function(v) { v[!is.finite(v)] <- NA_real_; v })
    df <- as.data.frame(out)
    colnames(df) <- paste0(prefix, which_stats)
    df
  }

  # --- Summaries ---
  snpcode_vec <- all$snpcode
  geno_vec    <- all$genotype_majmin
  cat("Genotype_majmin counts:\n")
  print(table(geno_vec))

  summaries <- list()
  for (std in std_names) {
    cat("Summarizing", std, "...\n")
    if (std == "raw") {
      mat <- as.matrix(dfs$raw[, sample_cols])
      stats_to_calc <- c("mean", "median", "sd", "mad", "cv", "iqr", "extreme")
      summaries[[std]] <- row_summaries_fast(mat, "raw_", which_stats=stats_to_calc)
      cat("raw summary dims:", dim(summaries[[std]]), "\n")
    } else if (std == "sz") {
      mat <- as.matrix(dfs$sz[, paste0(sample_cols, "_sample_z")])
      stats_to_calc <- c("mean", "median", "sd", "mad", "cv",  "iqr")
      summaries[[std]] <- row_summaries_fast(mat, "sz_", which_stats=stats_to_calc)
      cat("sz summary dims:", dim(summaries[[std]]), "\n")
    } else if (std == "lz") {
      mat <- as.matrix(dfs$lz[, paste0(sample_cols, "_loci_z")])
      stats_to_calc <- c("median", "mad")
      summaries[[std]] <- row_summaries_fast(mat, "lz_", which_stats=stats_to_calc)
      cat("lz summary dims:", dim(summaries[[std]]), "\n")
    } else if (std == "dlsz") {
      mat <- as.matrix(dfs$dlsz[, paste0(sample_cols, "_double_loci_sample")])
      stats_to_calc <- c("mean", "median", "sd", "mad", "iqr")
      summaries[[std]] <- row_summaries_fast(mat, "dlsz_", which_stats=stats_to_calc)
      cat("dlsz summary dims:", dim(summaries[[std]]), "\n")
    } else if (std == "dslz") {
      mat <- as.matrix(dfs$dslz[, paste0(sample_cols, "_double_sample_loci")])
      stats_to_calc <- c("median", "mad")
      summaries[[std]] <- row_summaries_fast(mat, "dslz_", which_stats=stats_to_calc)
      cat("dslz summary dims:", dim(summaries[[std]]), "\n")
    } else if (std == "aggregate") {
      mat <- as.matrix(dfs$aggregate[, paste0(sample_cols, "_aggregate")])
      stats_to_calc <- c("mean", "median", "sd", "mad", "cv", "iqr", "kurtosis", "skewness", "entropy", "extreme", "dip")
      summaries[[std]] <- row_summaries_fast(mat, "aggregate_", which_stats=stats_to_calc)
      cat("aggregate summary dims:", dim(summaries[[std]]), "\n")
    }
  }

  # --- Output Construction and Pivot ---
  if (identical(std_names, "aggregate")) {
    # Only aggregate: return aggregate stats alone (one row per snpcode)
    agg_stats <- summaries[["aggregate"]]
    agg_out <- tibble(snpcode = dfs$aggregate$snpcode)
    agg_out <- bind_cols(agg_out, agg_stats)
    cat("Returning only aggregate summary stats; dims:", dim(agg_out), "\n")
    return(agg_out)
  }

  # If aggregate is present with other standards, join its summary stats later
  summary_df <- tibble(
    snpcode = snpcode_vec,
    genotype_majmin = geno_vec
  )
  for (std in std_names) {
    if (std != "aggregate") {
      summary_df <- bind_cols(summary_df, summaries[[std]])
    }
  }
  cat("summary_df dims:", dim(summary_df), "\n")
  print(head(summary_df))

  all_summary_long <- summary_df

  # 2. Pivot longer/wider to produce summary stats by genotype and stat
  all_summary <- all_summary_long %>%
    pivot_longer(
      cols = matches("_(mean|median|iqr|cv|sd|mad|kurtosis|skewness|entropy|extreme|dip)$"),
      names_to = "stat",
      values_to = "value"
    ) %>%
    mutate(
      stat_geno = if_else(
        is.na(genotype_majmin),
        stat,
        paste0(stat, "_", genotype_majmin)
      )
    ) %>%
    pivot_wider(
      id_cols = c(snpcode),    
      names_from = stat_geno,
      values_from = value
    )

  # 3. Attach aggregate summary stats, if present
  if ("aggregate" %in% std_names) {
    agg_stats <- summaries[["aggregate"]]
    agg_df <- tibble(snpcode = dfs$aggregate$snpcode)
    agg_df <- bind_cols(agg_df, agg_stats)
    all_summary <- left_join(all_summary, agg_df, by = "snpcode")
  }

  cat("all_summary dims after pivot_wider:", dim(all_summary), "\n")

  # --- Add differences and log ratios between genotype classes ---
  add_diff_cols <- function(df, stat_prefix) {
    if (stat_prefix %in% c("lz_", "dslz_")) {
      suffixes <- "median"
    } else if (stat_prefix %in% "raw_") {
      suffixes <- c("mean", "median", "iqr")
    } else {
      suffixes <- c("mean", "median", "iqr")
    }
    for (suf in suffixes) {
      maj <- paste0(stat_prefix, suf, "_homoMaj")
      min <- paste0(stat_prefix, suf, "_homoMin")
      het <- paste0(stat_prefix, suf, "_hetero")
      diff1 <- paste0(stat_prefix, suf, "_homoDiff")
      diff2 <- paste0(stat_prefix, suf, "_MajHetDiff")
      diff3 <- paste0(stat_prefix, suf, "_MinHetDiff")

      # Only calculate diffs if columns exist
      if (all(c(maj, min) %in% names(df))) df[[diff1]] <- df[[maj]] - df[[min]]
      if (all(c(maj, het) %in% names(df))) df[[diff2]] <- df[[maj]] - df[[het]]
      if (all(c(min, het) %in% names(df))) df[[diff3]] <- df[[min]] - df[[het]]

      # Log-ratios only for raw_
      if (stat_prefix == "raw_") {
        # Log(Maj/Min)
        logratio1 <- paste0(stat_prefix, suf, "_MajMinLogRatio")
        if (all(c(maj, min) %in% names(df))) {
          m <- df[[maj]]
          n <- df[[min]]
          out <- rep(NA_real_, length(m))
          idx <- is.finite(m) & is.finite(n) & m > 0 & n > 0
          out[idx] <- log(m[idx] / n[idx])
          df[[logratio1]] <- out
        }
        # Log(Maj/Het)
        logratio2 <- paste0(stat_prefix, suf, "_MajHetLogRatio")
        if (all(c(maj, het) %in% names(df))) {
          m <- df[[maj]]
          h <- df[[het]]
          out <- rep(NA_real_, length(m))
          idx <- is.finite(m) & is.finite(h) & m > 0 & h > 0
          out[idx] <- log(m[idx] / h[idx])
          df[[logratio2]] <- out
        }
        # Log(Min/Het)
        logratio3 <- paste0(stat_prefix, suf, "_MinHetLogRatio")
        if (all(c(min, het) %in% names(df))) {
          n <- df[[min]]
          h <- df[[het]]
          out <- rep(NA_real_, length(n))
          idx <- is.finite(n) & is.finite(h) & n > 0 & h > 0
          out[idx] <- log(n[idx] / h[idx])
          df[[logratio3]] <- out
        }
      }
    }
    df
  }

  for (std in std_names) {
    prefix <- switch(std,
                     raw   = "raw_",
                     sz    = "sz_",
                     lz    = "lz_",
                     dlsz  = "dlsz_",
                     dslz  = "dslz_",
                     NULL)
    if (!is.null(prefix)) {
      all_summary <- add_diff_cols(all_summary, prefix)
    }
  }

  cat("Final summary NA counts (post-diff):\n")
  print(colSums(is.na(all_summary)))
  cat("Rows in final output:", nrow(all_summary), "\n")
  cat("=== END OF depth_by_genotype() RUN ===\n")
  return(all_summary)
}



