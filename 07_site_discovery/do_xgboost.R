library(ggplot2)
library(xgboost)
library(caret)
library(ParBayesianOptimization)
library(doParallel)

# -------------------------------------------------------
# Input arguments
# -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
snpstats_path <- args[1]
alt_path <- args[2]
ref_path <- args[3]
het_path <- args[4]
out_path <- args[5]
out_name <- args[6]
source_path <- args[7]

source(source_path) # discovery_and_filtering.R for depth_by_genotype()
out <- paste0(out_path, "/", out_name)

# -------------------------------------------------------
# Genotype-level site summaries
# -------------------------------------------------------

# Calculate genotype-level read depth summaries, e.g.,
# the mean depth of homoMaj, homoMin, and heterozygote 
# genotypes at a given locus.
# The coding for the minor allele is determined by the 
# reference state in the read count matrices (to ensure 
# consistent notation across domains), but we recode them 
# here based on their frequency.

all_summary <- depth_by_genotype(
  snpstats_path = snpstats_path,
  alt_path = alt_path,
  ref_path = ref_path,
  het_path = het_path,
  standards = c("none", 
                "sample_z", 
                "loci_z", 
                "double_loci_sample",
                "double_genotype_sample"))

# Join the genotype-level statistics with the 
# site-level statistics previously calculated with
# summarize_site_stats.py.

snpstats <- fread(snpstats_path,
                  data.table = FALSE)
snpstats <- left_join(all_summary, snpstats)

# An earlier version of summarize_site_stats.py did not
# remove sites where the frequency of the reference allele
# was < 5%. Double-check that those are removed here.
# Also remove sites without at least 2 observations of the 
# homoRef, homoAlt, and heterozygous genotypes, or with a
# < 40% call rate.

snpstats <- subset(snpstats, MAF < 0.95 & 
                   !is.na(raw_mean_hetero) &
                   !is.na(raw_mean_homoMin) & 
                   !is.na(raw_mean_homoMaj) & 
                   raw_iqr_homoMin > 0 & 
                   raw_iqr_homoMaj > 0 & 
                   raw_iqr_hetero > 0 & 
                   call_rate > 0.40)

rm(all_summary); gc()

write.csv(snpstats,
          paste0(out, "_all_site_summary_maf05.csv"), 
          row.names = FALSE, quote = FALSE)

# -------------------------------------------------------
# GBM preparation
# -------------------------------------------------------

# Check for multicollinearity. This is not a problem for 
# the GBM, but strongly correlated variables make the 
# interpretation more difficult.

drop_cols <- c("snpcode", "MAF", "Hexp", 
               "Hobs", "F", "HWE_LRT", 
               "HWE_pval", "baseQ_pval", 
               "mapQ_pval", "edge_pval")

# 'predictors' is a character vector of columns to check.
predictors <- colnames(snpstats)[!colnames(snpstats) %in% drop_cols]
dat <- snpstats[, predictors, drop = FALSE]
to_drop <- caret::findCorrelation(cor(dat, use = "pairwise.complete.obs"), cutoff = 0.7)
filtered_feats <- dat[, -to_drop, drop = FALSE]
# Add binary indicator for SB2 == NA or SB2 != NA.
filtered_feats$SBbin <- ifelse(is.na(filtered_feats[,"SB2"]), yes = 1, no = 0)

# Combine yvar columns and filtered features into gbmMat.
yvars <- c("snpcode", "MAF", "Hexp", "Hobs", "F")
gbmMat <- cbind(snpstats[,yvars], filtered_feats) 
rm(snpstats, filtered_feats); gc()

# Split into training and test sets, and add a set indicator to gbmMat.
set.seed(1232434234)
train.index <- sample(seq_len(nrow(gbmMat)), size = round(nrow(gbmMat) * 0.7))
gbmMat$set <- "test"
gbmMat$set[train.index] <- "train"

# Prepare matrices for xgboost.
# Remove yvar columns and set indicator for training and test data.
feature_cols <- setdiff(colnames(gbmMat), c(yvars, "set"))
train.dat <- as.matrix(gbmMat[gbmMat$set == "train", feature_cols, drop=FALSE])
test.dat  <- as.matrix(gbmMat[gbmMat$set == "test", feature_cols, drop=FALSE])

# Prepare outcome variables.
f.train <- gbmMat$F[gbmMat$set == "train"]
f.test  <- gbmMat$F[gbmMat$set == "test"]

# Save dataset.
write.csv(gbmMat, 
          paste0(out, "_GBM_alldat.csv"), 
          row.names = FALSE, quote = FALSE)

# -------------------------------------------------------
# Bayesian hyperparameter optimization
# -------------------------------------------------------

# Set up cross-validation folds.
folds <- list(
  fold1 = as.integer(seq(1,nrow(train.dat),by = 5)),
  fold2 = as.integer(seq(2,nrow(train.dat),by = 5)),
  fold3 = as.integer(seq(3,nrow(train.dat),by = 5)),
  fold4 = as.integer(seq(4,nrow(train.dat),by = 5)),
  fold5 = as.integer(seq(5,nrow(train.dat),by = 5))
)

# Define function for optimization. Must take the 
# hyperparameters as arguments.
obj_func <- function(max_depth, 
                     min_child_weight, 
                     subsample, 
                     colsample_bytree,
                     colsample_bynode) {
  # Hyperparameters. 
  param <- list(
    eta = 0.1, # Step size shrinkage used in updates to prevent overfitting.
    max_depth = max_depth, # Maximum depth of a tree.
    min_child_weight = min_child_weight, # Minimum sum of instance weights needed in a child.
    subsample = subsample, # Subsample ratio of the training data (rows).
    colsample_bytree =  colsample_bytree,  # Subsample ratio of columns when creating a new tree.
    colsample_bynode =  colsample_bynode # Subsample ratio of columns used in each split.
  )
  xgbcv <-  xgb.cv(
    params = param,
    data = train.dat,
    label = f.train,
    nrounds = 5000, # Largely depends on eta; lower eta needs more rounds. 5000 should be more than enough for eta = 0.1.
    folds = folds,
    objective = "reg:squarederror",  
    eval_metric = "rmse",
    verbose = 0,             
    early_stopping_rounds = 20, # Stop if no improvement.
    maximize = FALSE
  )
  lst <- list(
    # First argument must be named "Score".
    # Function finds maxima so inverting the output.
    Score = -1*min(xgbcv$evaluation_log$test_rmse_mean),
    # Get number of trees for the best performing model.
    nrounds = which.min(xgbcv$evaluation_log$test_rmse_mean)
  )
  return(lst)
}

# Define the search space boundaries.
bounds <- list(max_depth = c(1L, 7L),
               min_child_weight = c(1, 50),
               subsample = c(0.5, 1),
               colsample_bytree = c(0.5, 1),
               colsample_bynode = c(0.5, 1))

# Set up to run in parallel.
cl <- makeCluster(3)
registerDoParallel(cl)
clusterExport(cl, c("folds", "train.dat", "f.train"))
clusterEvalQ(cl, expr= {
  library(xgboost)
})

# Run the Bayesian optimization
# If running interactively, start with a small
# number of iterations (e.g., iters.n = initPoints + 10), inspect output, and then add more
# as needed with addIterations(). In my testing, gpuUtility converged
# around iters.n = 20, so 40 is to be extra-sure I won't need to 
# rerun the optimization.
# initPoints is the initial number of evenly spaced points in 
# the parameter space to sample. Set it too low and optimization
# will take forever to converge; too high is a waste of computational
# resources. scikit-learn's optimizer uses 10 as a default; ParBayesianOptimization
# uses n.parameters + 2 in examples, and mlrMBO uses n.parameters*4 by default.

bayes_out <- bayesOpt(FUN = obj_func,
                      bounds = bounds, 
                      initPoints = length(bounds) + 6, 
                      iters.n = 40,
                      iters.k = 3, 
                      parallel = TRUE)
stopCluster(cl)

# Temporary: save session in case more iterations are needed
save.image(paste0(out, "_ParOpt.RData")) 

opt_params <- append(list(eta = 0.1), 
                     getBestPars(bayes_out))
n.rounds <- bayes_out$scoreSummary$nrounds[which(bayes_out$scoreSummary$Score == max(bayes_out$scoreSummary$Score))]

# Save the parameter values.
write.csv(unlist(append(list(eta = 0.1, 
                             nrounds = n.rounds), 
                        getBestPars(bayes_out))),
          paste0(out, "_GBM_hyperparameters.csv"),
          row.names = FALSE, quote = FALSE)

# Save the score summary.
write.csv(bayes_out$scoreSummary,
          paste0(out, "_BayesOpt_score_summary.csv"),
          row.names = FALSE, quote = FALSE)

# -------------------------------------------------------
# Fit GBM with xgboost
# -------------------------------------------------------

model <- xgboost(
  data = train.dat,
  label = f.train,
  params = opt_params,
  nrounds = n.rounds,
  objective = "reg:squarederror",
  eval_metric = "rmse",
  verbose = 0 
)

# -------------------------------------------------------
# GBM predictions
# -------------------------------------------------------

train.pred <- predict(model, train.dat)
test.pred <- predict(model, test.dat)

results <- rbind(
  data.frame(set = "test",  actual = f.test,  pred = test.pred),
  data.frame(set = "train", actual = f.train, pred = train.pred)
)

# Add squared errors.
results$sqerr <- (results$pred - results$actual)^2

# Save predicted and actual values
write.csv(results,
          paste0(out, "_GBM_predvactual.csv"),
          row.names = FALSE, quote = FALSE)

# Plot predicted vs. actual F.
p <- ggplot() +
  geom_hex(data = subset(results, set == "test"), aes(x = actual, y = pred), bins = 30) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + 
  geom_hline(yintercept = mean(subset(results, set == "test")$pred),  linetype = "dashed") +
  geom_vline(xintercept = mean(subset(results, set == "test")$pred), linetype = "dashed") +
  scale_fill_viridis_c() +
  xlab("observed F") +
  ylab("predicted F") +
  ggtitle(paste0("RMSE = ", round(sqrt(mean(subset(results, set == "test")$sqerr)),2))) + 
  theme_minimal()

ggsave(paste0(out, "_test_predvactual.png"), 
       plot = p, width = 6, 
       height = 6, dpi = 300, 
       bg = "white")

p <- ggplot() +
  geom_hex(data = subset(results, set == "train"), aes(x = actual, y = pred), bins = 30) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + 
  geom_hline(yintercept = mean(subset(results, set == "train")$pred),  linetype = "dashed") +
  geom_vline(xintercept = mean(subset(results, set == "train")$pred), linetype = "dashed") +
  scale_fill_viridis_c() +
  xlab("observed F") +
  ylab("predicted F") +
  ggtitle(paste0("RMSE = ", round(sqrt(mean(subset(results, set == "train")$sqerr)),2))) + 
  theme_minimal()

ggsave(paste0(out, "_train_predvactual.png"), 
       plot = p, width = 6, 
       height = 6, dpi = 300, 
       bg = "white")
rm(p); gc()
