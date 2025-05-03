# angsd_param_summary.R
# Filters the results from ANGSD_param_sweep.sh by library call rate thresholds
# aggregates per-locus, population-level statistics
# and subsets domain_sfs.beagle.gz by callable sites per threshold

library(data.table)
library(tidyverse)

# ---------- INPUT ARGUMENTS ----------
args <- commandArgs(trailingOnly = TRUE)
project_dir  <- args[1]
meta_file    <- args[2]
bamlist_file <- args[3]
out_name      <- args[4]

# -------- Helpers --------

normalize_columns <- function(df) {
  df %>%
    rename_with(~ str_replace_all(., regex(".*chr.*", ignore_case = TRUE), "scaffold")) %>%
    rename_with(~ str_replace_all(., regex(".*pos.*", ignore_case = TRUE), "pos"))
}

get_callable_sites <- function(counts, meta, threshold = 0.5, min_lib = 0.1) {
  meta <- meta %>% filter(bam_code %in% colnames(counts))

  libs <- meta %>%
    count(library) %>%
    filter(n / sum(n) >= min_lib) %>%
    pull(library)

  valid <- vector("list", length(libs))
  names(valid) <- libs
  for (lib in libs) {
    samples <- meta %>% filter(library == lib) %>% pull(bam_code)
    present_counts <- counts[, samples, drop = FALSE] > 0
    keep <- rowSums(present_counts) >= threshold * length(samples)
    valid[[lib]] <- rownames(counts)[keep]
  }

  Reduce(intersect, valid)
}

get_all_callable_sites <- function(count_path, meta_path, bamlist_path, thresholds = c(0.4, 0.5, 0.6)) {
  counts_df <- read.table(gzfile(count_path), header = TRUE, check.names = FALSE)
  pos       <- read.table(gsub("counts.gz", "pos.gz", count_path), header = TRUE)
  meta      <- read.csv(meta_path)

  samples <- readLines(bamlist_path) %>% basename() %>% str_remove("\\.realigned\\.bam$")
  if (length(samples) != ncol(counts_df)) stop("Mismatch between number of samples in BAM list and columns in counts matrix")

  colnames(counts_df) <- samples
  snpcode             <- paste0(pos$chr, "_", pos$pos)
  rownames(counts_df) <- snpcode

  counts <- as.matrix(counts_df)
  rm(counts_df); gc()
  counts <- counts > 0

  maf_file <- gsub("counts.gz", "mafs.gz", count_path)
  mafs     <- read.table(gzfile(maf_file), header = TRUE, check.names = FALSE) %>% normalize_columns() %>% mutate(snpcode = paste0(scaffold, "_", pos))
  maf_pass <- mafs %>% filter(knownEM >= 0.05) %>% pull(snpcode)
  rm(mafs); gc()

  callable_sites        <- list()
  callable_maf_filtered <- list()
  for (thresh in thresholds) {
    message("Processing threshold ", thresh)
    snps <- get_callable_sites(counts, meta, threshold = thresh)
    callable_sites[[as.character(thresh)]]        <- snps
    callable_maf_filtered[[as.character(thresh)]] <- intersect(snps, maf_pass)
  }

  list(all = callable_sites, maf_filtered = callable_maf_filtered)
}

parse_population_stats <- function(pop_dir) {
  theta_file <- file.path(pop_dir, "thetas.persite.txt")
  hwe_file   <- file.path(pop_dir, "gl.hwe.gz")

  thetas <- fread(theta_file, header = TRUE) %>% normalize_columns() %>% transmute(
    snpcode  = paste0(scaffold, "_", pos),
    pi        = Pairwise,
    theta_W   = Watterson,
    theta_F   = thetaSingleton,
    theta_H   = thetaH,
    theta_L   = thetaL
  )

  hwe <- fread(hwe_file, header = TRUE) %>% normalize_columns() %>% transmute(
    snpcode = paste0(scaffold, "_", pos),
    MAF     = hweFreq,
    F       = F,
    Hexp    = 2 * hweFreq * (1 - hweFreq),
    Hobs    = Hexp - F * Hexp
  )

  left_join(thetas, hwe, by = "snpcode")
}

parse_and_filter_population <- function(pop_dir, pop_code, param_id, callable_sites, which = "all") {
  stats    <- parse_population_stats(pop_dir)
  site_set <- callable_sites[[which]]
  dfs <- lapply(names(site_set), function(thresh) {
    stats %>%
      filter(snpcode %in% site_set[[thresh]]) %>%
      mutate(pop_code = pop_code, param_id = param_id, call_thresh = as.numeric(thresh))
  })
  bind_rows(dfs)
}

# -------- the pipeline --------
main <- function(project_dir = ".", meta_file = "metadata.csv", bamlist_file = NULL, out_name = NULL) {
  param_folders <- list.dirs(project_dir, recursive = FALSE)
  out_all <- file.path(project_dir, paste0(out_name, "_angsd_param_summaries.csv"))
  out_maf <- file.path(project_dir, paste0(out_name, "_maf05_angsd_param_summaries.csv"))
  first_write <- TRUE

  for (param_path in param_folders) {
    param_id <- basename(param_path)
    cat("Processing param combo:", param_id, "\n")

    count_file <- file.path(param_path, "domain_sfs.counts.gz")
    callable   <- get_all_callable_sites(count_file, meta_file, bamlist_file)

    # subset beagle by maf‐filtered snpcodes per threshold
    beagle_in <- file.path(param_path, "domain_sfs.beagle.gz")
    beagle_dt <- fread(beagle_in, check.names = FALSE)
    for (thresh in names(callable$maf_filtered)) {
      snps  <- callable$maf_filtered[[thresh]]
      suffix <- sub("^0\\.(.*)$", "\\1", thresh)
      out_beagle <- file.path(param_path, paste0(param_id, "_ct", suffix, ".beagle.gz"))
      subset_dt <- beagle_dt[marker %in% snps]
      fwrite(subset_dt, out_beagle, sep = " ", quote = FALSE, compress = "gzip")
    }
    rm(beagle_dt); gc()

    pop_dirs <- list.dirs(param_path, recursive = FALSE)
    pop_dirs <- pop_dirs[basename(pop_dirs) != param_id]

    for (pop_path in pop_dirs) {
      pop_code <- basename(pop_path)
      cat("  -> Population:", pop_code, "\n")

      dat     <- parse_and_filter_population(pop_path, pop_code, param_id, callable, which = "all")
      maf_dat <- parse_and_filter_population(pop_path, pop_code, param_id, callable, which = "maf_filtered")

      if (first_write) {
        write_csv(dat,  out_all, append = FALSE)
        write_csv(maf_dat, out_maf, append = FALSE)
        first_write <- FALSE
      } else {
        write_csv(dat,  out_all, append = TRUE)
        write_csv(maf_dat, out_maf, append = TRUE)
      }
      rm(dat, maf_dat); gc()
    }
    rm(callable); gc()
  }
  cat("All done. Summaries written to:\n  ", out_all, "\n  ", out_maf, "\n")
}

# -------- run the pipeline --------
results <- main(project_dir, meta_file, bamlist_file, out_name)
