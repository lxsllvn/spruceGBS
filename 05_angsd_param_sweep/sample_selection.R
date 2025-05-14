# =========================================
# Sample Filtering for Parameter Experiment
# =========================================

library(dplyr)
library(tidyr)

#--------------------------------------
# 1. Load and threshold depth summary
#--------------------------------------
load_depth_summary <- function(path,
                               read_quantile = 0.25,
                               scaff_quantile = 0.25) {
  df <- read.table(path, header = TRUE)

  df %>%
    mutate(
      log_reads  = log(total_mapped_reads),
      log_scaffs = log(n_scaffolds_with_mapped_reads)
    ) -> depths

  cut_reads  <- quantile(depths$log_reads,  probs = read_quantile)
  cut_scaffs <- quantile(depths$log_scaffs, probs = scaff_quantile)

  depths %>%
    filter(
      log_reads  > cut_reads,
      log_scaffs > cut_scaffs
    ) %>%
    select(-log_reads, -log_scaffs)
}

#--------------------------------------
# 2. Annotate with metadata
#--------------------------------------
annotate_metadata <- function(filtered_depths,
                              metadata,
                              meta_cols = c(
                                "bam_code", "pop_code", "domain", "region",
                                "library", "latitude", "longitude"
                              )) {
  filtered_depths %>%
    left_join(metadata %>% select(all_of(meta_cols)), by = "bam_code")
}

#--------------------------------------
# 3. Identify candidate populations
#   (>= min_samples & <= max_prop from a single library)
#--------------------------------------
find_candidates <- function(samples,
                            min_samples = 5,
                            max_prop     = 0.8) {
  samples %>%
    count(domain, pop_code, library, name = "n_samples") %>%
    group_by(domain, pop_code) %>%
    mutate(
      total_samples = sum(n_samples),
      prop          = n_samples / total_samples
    ) %>%
    slice_max(n_samples, with_ties = FALSE) %>%
    ungroup() %>%
    filter(
      total_samples >= min_samples,
      prop < max_prop
    ) %>%
    transmute(
      pop_code,
      majority_library = library,
      max_library_prop = prop,
      total_samples
    )
}

#--------------------------------------
# 4. Build full-population library lookup
#   (>= min_samples, regardless of prop)
#--------------------------------------
build_library_lookup <- function(pass, min_samples = 5) {
  pass %>%
    count(pop_code, library, name = "n_samples") %>%
    group_by(pop_code) %>%
    mutate(total_samples = sum(n_samples)) %>%
    filter(total_samples >= min_samples) %>%
    slice_max(n_samples, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      pop_code,
      majority_library = library,
      total_samples
    )
}

#--------------------------------------
# 5. Geodesic filtering of population pairs
#--------------------------------------
filter_pairs_by_distance <- function(geo_matrix,
                                     pops,
                                     domain_map,
                                     siberia_thresh = 250,
                                     other_thresh   = 100) {
  geo_matrix %>%
    filter(pop_code %in% pops) %>%
    pivot_longer(-pop_code, names_to = "pop2", values_to = "distance_km") %>%
    rename(pop1 = pop_code) %>%
    filter(pop1 != pop2, !is.na(distance_km)) %>%
    left_join(domain_map, by = c("pop1" = "pop_code")) %>%
    rename(domain1 = domain) %>%
    left_join(domain_map, by = c("pop2" = "pop_code")) %>%
    rename(domain2 = domain) %>%
    mutate(
      is_siberia = domain1 == "siberia" | domain2 == "siberia",
      keep       = if_else(is_siberia, distance_km <= siberia_thresh, distance_km <= other_thresh),
      pair_id    = paste0(pmin(pop1, pop2), "_", pmax(pop1, pop2))
    ) %>%
    filter(keep) %>%
    distinct(pair_id, .keep_all = TRUE) %>%
    select(pop1, pop2, distance_km)
}

#--------------------------------------
# 6. Exclude same-library pairs
#--------------------------------------
filter_cross_library <- function(pairs, lib_lookup) {
  pairs %>%
    left_join(lib_lookup, by = c("pop1" = "pop_code")) %>%
    rename(lib1 = majority_library) %>%
    left_join(lib_lookup, by = c("pop2" = "pop_code")) %>%
    rename(lib2 = majority_library) %>%
    filter(lib1 != lib2) %>%
    select(pop1, pop2, distance_km)
}

#--------------------------------------
# 7. Build final population status table
#--------------------------------------
#' Build a summary of population-level experiment design
#'
#' @param candidates_tbl Tibble of candidate populations
#' @param paired_tbl Tibble of paired population relations
#' @param domain_map Tibble mapping pop_code -> domain
#' @param lib_lookup Tibble mapping every pop_code -> majority_library + total_samples
#' @return Tibble with pop_code, domain, is_mixed, is_paired,
#'   majority_library, max_library_prop, total_samples,
#'   paired_pop, distance_km
build_population_status_table <- function(candidates_tbl,
                                          paired_tbl,
                                          domain_map,
                                          lib_lookup) {
  # Mixed populations: candidates
  mixed_tbl <- candidates_tbl %>%
    mutate(is_mixed = TRUE, is_paired = FALSE)

  # Paired population endpoints
  paired_info <- bind_rows(
    paired_tbl %>% transmute(pop_code = pop1, paired_pop = pop2, distance_km),
    paired_tbl %>% transmute(pop_code = pop2, paired_pop = pop1, distance_km)
  ) %>%
    left_join(candidates_tbl, by = "pop_code") %>%
    transmute(
      pop_code,
      is_mixed    = FALSE,
      is_paired   = TRUE,
      max_library_prop,
      total_samples = NA_integer_,
      paired_pop,
      distance_km
    )

  # Combine mixed + paired
  combined <- bind_rows(
    mixed_tbl %>% select(pop_code, is_mixed, is_paired, max_library_prop, total_samples),
    paired_info
  )

    # Join lib_lookup to get majority_library and true total_samples
  combined <- combined %>%
    left_join(lib_lookup, by = "pop_code", suffix = c("", ".lookup")) %>%
    mutate(
      total_samples = coalesce(total_samples, total_samples.lookup)
    ) %>%
    select(
      pop_code,
      majority_library,
      max_library_prop,
      total_samples,
      is_mixed,
      is_paired,
      paired_pop,
      distance_km
    )

  # Attach domain and finalize
  combined %>%
    left_join(domain_map, by = "pop_code") %>%
    arrange(domain, pop_code)
  combined %>%
    left_join(domain_map, by = "pop_code") %>%
    arrange(domain, pop_code)
}
