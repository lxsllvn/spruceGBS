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
  
  # 1. Handle covariance file naming
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
  
  # 2. Read sample and covariance
  samples <- get.info(sample_list_path, "sequenced_samples_metadata.csv")
  cov <- read.table(cov_file)
  colnames(cov) <- samples$bam_code
  rownames(cov) <- samples$bam_code
  
  # 3. PCA
  eig <- eigen(cov)
  n_samples <- nrow(cov)
  max_pc <- min(n_samples, ncol(cov))
  n_pc_to_plot <- max(PCs)
  if (n_pc_to_plot > max_pc) stop("Requested PCs exceed the number of possible PCs.")
  
  # 4. Make PC matrix
  pca_df <- cbind(
    samples,
    as.data.frame(eig$vectors[, seq_len(n_pc_to_plot), drop=FALSE])
  )
  colnames(pca_df)[(ncol(pca_df)-n_pc_to_plot+1):ncol(pca_df)] <- paste0("PC", 1:n_pc_to_plot)
  
  # 5. Chunk PCs into pairs
  pc_pairs <- split(PCs, ceiling(seq_along(PCs)/2))
  if (length(pc_pairs[[length(pc_pairs)]]) == 1) {
    warning("Last PC chunk only has 1 PC. Skipping it.")
    pc_pairs <- pc_pairs[sapply(pc_pairs, length) == 2]
  }
  if (length(pc_pairs) == 0) stop("No valid PC pairs to plot.")
  
  # 6. Palette presets for when fill_pal is NULL only
  presets <- list(
    Set1 = setdiff(brewer.pal(9, "Set1"), "#FFFF33"),
    Set2 = brewer.pal(8, "Set2"),
    Set3 = setdiff(brewer.pal(12, "Set3"), c("#FFFFB3","#FFED6F"))
  )
  
  # Helper for palette selection
  get_palette <- function(var, pal) {
    is_discrete <- !is.numeric(pca_df[[var]])
    # If pal is one of the preset names, return preset colors
    if (is.character(pal) && pal %in% names(presets)) {
      return(presets[[pal]])
    }
    if (!is.null(pal)) {
      return(as.character(pal))
    }
    if (is_discrete) {
      return(presets$Set3)
    } else {
      return(as.character(paletteer_c("grDevices::Viridis", 30)))
    }
  }
  
  # 7. Plotting
  all_plots <- list()
  for (pc_idx in seq_along(pc_pairs)) {
    pcs <- pc_pairs[[pc_idx]]
    x_pc <- paste0("PC", pcs[1])
    y_pc <- paste0("PC", pcs[2])
    plot_list <- list()
    for (j in seq_along(fill_var)) {
      var <- fill_var[j]
      pal <- if (!is.null(fill_pal) && length(fill_pal) >= j) fill_pal[[j]] else NULL
      pal_vec <- get_palette(var, pal)
      is_discrete <- !is.numeric(pca_df[[var]])
      mapping <- aes_string(x = x_pc, y = y_pc, col = var)
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
          legend.key.height = unit(6, "pt"),
          legend.key.width = unit(6, "pt"),
          panel.spacing = margin(0, 3, 0, 3, unit = "mm"),
          legend.text = element_text(size = 6, margin = margin(l = 2, r = 2)),
          plot.title = element_text(size=10),
          plot.margin = margin(0, 1, 0, 1, unit = "mm")
        ) +
        guides(colour = guide_legend(override.aes = list(size=1.5)))
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
        # If length is 1 but n_groups > 1, assume user forgot to set n
        if (length(pal_final) == 1 && n_groups > 1) {
          # Try to "re-call" paletteer with n_groups if possible (parse string)
          try_palette <- try(eval(parse(text = pal_final)), silent = TRUE)
          if (inherits(try_palette, "colors")) {
            pal_final <- as.character(try_palette(n_groups))
          } else {
            # Otherwise, repeat color to needed length (fallback)
            pal_final <- rep(pal_final, n_groups)
          }
        }
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
    
    # 8. Save output if requested
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

