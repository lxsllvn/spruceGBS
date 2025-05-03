# run_pca_manova.R

# ---------- INPUT ARGUMENTS ----------
args        <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript run_pca_manova.R <domain> <results_dir> <bamlist_path> <meta_path>")
}
domain       <- args[1]
results_dir  <- args[2]
bamlist_path <- args[3]
meta_path    <- args[4]

# ---------- LIBRARIES ----------
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
})

# ---------- FUNCTION DEFINITION ----------
summarize_domain_covs <- function(domain,
                                  results_dir,
                                  bamlist_path,
                                  meta_path) {
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  bam_codes <- readLines(bamlist_path) %>%
    basename() %>%
    str_remove("\\.realigned\\.bam$")
  indv_info <- tibble(bam_code = bam_codes) %>%
    left_join(meta, by = "bam_code")
  combo_dirs <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)
  
  map_dfr(combo_dirs, function(combo_dir) {
    fname  <- basename(combo_dir)
    params <- str_match(fname,
      sprintf("^%s_baq(\\d+)_C(\\d+)_q(\\d+)_mq(\\d+)$", domain))
    baq    <- as.integer(params[2])
    C      <- as.integer(params[3])
    q      <- as.integer(params[4])
    mq     <- as.integer(params[5])
    
    cov_files <- list.files(combo_dir,
                            pattern = paste0("^", fname, "_ct\\d+\\.Pcangsd\\.cov$"),
                            full.names = TRUE)
    
    map_dfr(cov_files, function(cov_path) {
      ct <- as.numeric(str_extract(basename(cov_path), "(?<=_ct)\\d+")) / 10
      mat <- as.matrix(read.table(cov_path, header = FALSE))
      colnames(mat) <- indv_info$bam_code
      rownames(mat) <- indv_info$bam_code
      
      eig   <- eigen(mat)
      PC1   <- eig$vectors[,1]
      PC2   <- eig$vectors[,2]
      pca_df <- indv_info %>% mutate(PC1 = PC1, PC2 = PC2)
      
      m   <- manova(cbind(PC1, PC2) ~ library + region, data = pca_df)
      sm  <- summary(m, test = "Wilks")
      SSl <- sm$SS$library
      SSr <- sm$SS$region
      SSx <- sm$SS$Residuals
      SSt <- SSl + SSr + SSx
      trace <- function(x) sum(diag(x))
      
      tibble(
        domain      = domain,
        baq         = baq,
        C           = C,
        q           = q,
        mq          = mq,
        ct          = ct,
        R2_library  = trace(SSl) / trace(SSt),
        R2_region   = trace(SSr) / trace(SSt)
      )
    })
  })
}

# ---------- RUN & WRITE OUT ----------
res <- summarize_domain_covs(domain, results_dir, bamlist_path, meta_path)



# save to a CSV for downstream steps:
write.csv(res,
          file = paste0(domain, "_pca_manova_summary.csv"),
          row.names = FALSE)
