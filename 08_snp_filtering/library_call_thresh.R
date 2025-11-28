# library_call_tresh.R
# Identifies and creates snpcode lists for sites with a 
# minimum 50/60/70/80% per-library call rate

library(tidyverse)

# ---------- INPUT ARGUMENTS ----------
args <- commandArgs(trailingOnly = TRUE)
count_path  <- args[1]
meta_path   <- args[2]
out_base    <- args[3]

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
    valid[[lib]] <- counts[keep, "snpcode"]
  }
  
  Reduce(intersect, valid)
}

meta      <- read.csv(meta_path)
counts_df <- read.table(gzfile(count_path), header = TRUE, check.names = FALSE)

for (i in c(0.5, 0.6, 0.7, 0.8)) {
  
sites <- get_callable_sites(counts_df, meta, i)
write.table(x = sites,
            file = paste0(out_base, "_call_thresh", i),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
}
