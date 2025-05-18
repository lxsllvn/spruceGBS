# ---------------------------------------------------------
# Sample Filtering & Population Pairing Pipeline
#
# This script filters sequencing samples, identifies mixed-library populations,
# and finds candidate population pairs for parameter experiments.
#
# Dependencies: dplyr, tidyr
#
# Usage: See example code at the end of this script.
# ---------------------------------------------------------

library(dplyr)
library(tidyr)

#' Load and filter sample depth summary
#'
#' @param path Path to sample depth summary table.
#' @param meta_path Path to metadata CSV file.
#' @param meta_cols Metadata columns to join.
#' @param read_quantile Lower quantile for n_reads filter.
#' @param scaff_quantile Lower quantile for n_scaffolds filter.
#' @return Filtered and annotated sample tibble.
load_and_filter_samples <- function(path, meta_path, meta_cols = c("bam_code", "pop_code", "domain", "region", "library", "latitude", "longitude"),
                                   read_quantile = 0.25, scaff_quantile = 0.25) {
  if (!file.exists(path)) stop("Sample depth file not found: ", path)
  if (!file.exists(meta_path)) stop("Metadata file not found: ", meta_path)
  df <- read.table(path, header = TRUE)
  meta <- read.csv(meta_path, header = TRUE)
  depths <- df %>%
    left_join(meta %>% select(all_of(meta_cols)), by = "bam_code") %>%
    mutate(
      log_reads = log(n_reads),
      log_scaffs = log(n_scaffolds)
    )
  filtered <- depths %>%
    filter(
      log_reads > quantile(log_reads, read_quantile),
      log_scaffs > quantile(log_scaffs, scaff_quantile)
    ) %>%
    select(-log_reads, -log_scaffs)
  return(filtered)
}

#' Summarize population/library structure
#'
#' @param samples Data frame with filtered/annotated samples.
#' @param min_samples Minimum number of samples per population.
#' @param max_prop Max proportion allowed from a single library (for mixed_only).
#' @param mixed_only Logical: only return mixed-library populations?
#' @param keep_cols Metadata columns to keep (e.g., "domain", "region").
#' @return Data frame summarizing populations and library composition.
summarize_populations <- function(samples, 
                                  min_samples = 5, 
                                  max_prop = 0.8, 
                                  mixed_only = TRUE,
                                  keep_cols = c("domain")) {
  out <- samples %>%
    count(across(all_of(c(keep_cols, "pop_code", "library"))), name = "n_samples") %>%
    group_by(across(all_of(c(keep_cols, "pop_code")))) %>%
    mutate(
      total_samples = sum(n_samples),
      prop = n_samples / total_samples
    ) %>%
    slice_max(n_samples, with_ties = FALSE) %>%
    ungroup() %>%
    filter(total_samples >= min_samples)
  if (mixed_only) {
    out <- out %>% filter(prop < max_prop)
  }
  summary <- out %>%
    transmute(
      pop_code,
      across(all_of(keep_cols)),
      majority_library = library,
      max_library_prop = prop,
      total_samples
    )
  return(summary)
}

#' Find population pairs by geodesic distance and cross-library status
#'
#' @param geo_matrix_path Path to distance matrix file (tab-delimited).
#' @param lib_lookup Output of summarize_populations(), must include 'domain'.
#' @param siberia_thresh Max km for Siberian pairs.
#' @param other_thresh Max km for other domains.
#' @return Data frame of cross-library population pairs within distance thresholds.
find_population_pairs <- function(geo_matrix_path, lib_lookup, 
                                  siberia_thresh = 250, other_thresh = 100) {
  if (!file.exists(geo_matrix_path)) stop("Geodesic distance file not found: ", geo_matrix_path)
  geo_matrix <- read.table(geo_matrix_path, header = TRUE, check.names = FALSE)
  pop_codes <- lib_lookup$pop_code
  pop_info <- lib_lookup %>% select(pop_code, domain)
  pairs <- geo_matrix %>%
    filter(pop_code %in% pop_codes) %>%
    pivot_longer(-pop_code, names_to = "pop2", values_to = "distance_km") %>%
    rename(pop1 = pop_code) %>%
    filter(pop1 != pop2, !is.na(distance_km)) %>%
    left_join(pop_info, by = c("pop1" = "pop_code")) %>% rename(domain1 = domain) %>%
    left_join(pop_info, by = c("pop2" = "pop_code")) %>% rename(domain2 = domain) %>%
    mutate(
      is_siberia = domain1 == "siberia" | domain2 == "siberia",
      keep = if_else(is_siberia, distance_km <= siberia_thresh, distance_km <= other_thresh),
      pair_id = paste0(pmin(pop1, pop2), "_", pmax(pop1, pop2))
    ) %>%
    filter(keep) %>%
    distinct(pair_id, .keep_all = TRUE)
  pairs <- pairs %>%
    left_join(lib_lookup %>% select(pop_code, majority_library), by = c("pop1" = "pop_code")) %>%
    rename(lib1 = majority_library) %>%
    left_join(lib_lookup %>% select(pop_code, majority_library), by = c("pop2" = "pop_code")) %>%
    rename(lib2 = majority_library) %>%
    filter(lib1 != lib2) %>%
    select(pop1, pop2, distance_km)
  return(pairs)
}

#' Build final population status table
#'
#' @param mixed_tbl Output of summarize_populations() with mixed_only=TRUE.
#' @param paired_tbl Output of find_population_pairs().
#' @param lib_lookup Output of summarize_populations() with mixed_only=FALSE.
#' @return Status table summarizing populations, library, and pairing info.
build_population_status <- function(mixed_tbl, paired_tbl, lib_lookup) {
  extra_cols <- setdiff(names(mixed_tbl), c("pop_code", "majority_library", "max_library_prop", "total_samples"))
  all_pops <- union(mixed_tbl$pop_code, unique(c(paired_tbl$pop1, paired_tbl$pop2)))
  status <- tibble(pop_code = all_pops)
  status <- status %>%
    left_join(
      mixed_tbl %>% select(pop_code, all_of(extra_cols), max_library_prop, total_samples),
      by = "pop_code"
    ) %>%
    mutate(is_mixed = !is.na(max_library_prop)) %>%
    left_join(lib_lookup %>% select(pop_code, majority_library, total_samples), by = "pop_code", suffix = c("", ".lookup")) %>%
    mutate(total_samples = coalesce(total_samples, total_samples.lookup)) %>%
    select(-total_samples.lookup)
  paired_long <- bind_rows(
    paired_tbl %>% transmute(pop_code = pop1, paired_pop = pop2, distance_km),
    paired_tbl %>% transmute(pop_code = pop2, paired_pop = pop1, distance_km)
  )
  status <- status %>%
    left_join(paired_long, by = "pop_code") %>%
    mutate(is_paired = !is.na(paired_pop))
  final <- status %>%
    select(pop_code, all_of(extra_cols), majority_library, max_library_prop, total_samples, is_mixed, is_paired, paired_pop, distance_km) %>%
    arrange(across(all_of(extra_cols)), pop_code)
  return(final)
}

# ==========================
# Example usage (uncomment to run)
# ==========================

# Step 1: Filter samples
# samples <- load_and_filter_samples("depth_summary.tsv", "metadata.csv")

# Step 2: Identify mixed-library populations and build lookup
# mixed_pops <- summarize_populations(samples, min_samples=5, max_prop=0.8, mixed_only=TRUE, keep_cols="domain")
# lib_lookup <- summarize_populations(samples, min_samples=5, mixed_only=FALSE, keep_cols="domain")

# Step 3: Find population pairs using distance matrix path
# pairs <- find_population_pairs("geo_matrix.tsv", lib_lookup)

# Step 4: Build final population status table
# status_table <- build_population_status(mixed_pops, pairs, lib_lookup)

# View output
# print(status_table)
