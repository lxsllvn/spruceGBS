suppressPackageStartupMessages({
  library(xgboost)
  library(data.table)
})

#============================#
# Small utilities
#============================#
timestamped_file <- function(prefix, ext) {
  sprintf("%s_%s%s", prefix, format(Sys.time(), "%Y%m%d_%H%M%S"), ext)
}
.info     <- function(...) cat(sprintf("[INFO] %s\n", sprintf(...)))
.warn_bad <- function(...) warning(sprintf("[WARN] %s", sprintf(...)), call. = FALSE)
.stop_bad <- function(...) stop(sprintf("[ERROR] %s", sprintf(...)))

#============================#
# Parse model → compact tree table (model → tdt)
#============================#
# Purpose: read xgboost booster into a tidy, consistent data.table
parse_xgb_tree <- function(model) {
  stopifnot(inherits(model, "xgb.Booster"))
  dt <- xgboost::xgb.model.dt.tree(model = model)
  data.table::setDT(dt)
  
  to_int_child <- function(x) {
    # children look like "0-123" → take the number after the hyphen
    y <- suppressWarnings(as.integer(sub(".*-", "", x)))
    y[is.na(x)] <- NA_integer_
    y
  }
  
  if (!"Node" %in% names(dt)) stop("[parse_xgb_tree] column 'Node' not found")
  
  # Basic normalization
  dt[, `:=`(
    Tree = as.integer(Tree),
    ID   = as.integer(Node)
  )]
  
  for (col in c("Yes", "No", "Missing")) {
    if (col %in% names(dt)) dt[, (col) := to_int_child(get(col))] else dt[, (col) := NA_integer_]
  }
  
  if (!"Split" %in% names(dt)) dt[, Split := NA_real_]
  if (!"Cover" %in% names(dt)) dt[, Cover := NA_real_]
  dt[, `:=`(Split = suppressWarnings(as.numeric(Split)),
            Cover = suppressWarnings(as.numeric(Cover)))]
  
  # Leaf flag
  dt[, Leaf := (Feature == "Leaf")]
  
  has_quality <- "Quality" %in% names(dt)
  has_gaincol <- "Gain" %in% names(dt)
  
  # Leaf values and gains
  if (has_quality) {
    q <- suppressWarnings(as.numeric(dt$Quality))
    # Quality = leaf value on leaves; gain on splits
    dt[, LeafVal := ifelse(Leaf, q, NA_real_)]
    if (has_gaincol) {
      g <- suppressWarnings(as.numeric(dt$Gain))
      dt[, Gain := ifelse(!Leaf, g, NA_real_)]
    } else {
      dt[, Gain := ifelse(!Leaf, q, NA_real_)]
    }
  } else {
    # no Quality; try Gain for splits, and Split for leaves as a last resort
    if (has_gaincol) {
      g <- suppressWarnings(as.numeric(dt$Gain))
      dt[, Gain := ifelse(!Leaf, g, NA_real_)]
    } else {
      dt[, Gain := NA_real_]
    }
    s <- suppressWarnings(as.numeric(dt$Split))
    dt[, LeafVal := ifelse(Leaf, s, NA_real_)]
  }
  
  # Clean up
  drop_cols <- intersect(c("Node", "Quality"), names(dt))
  if (length(drop_cols)) dt[, (drop_cols) := NULL]
  
  data.table::setkey(dt, Tree, ID)
  data.table::setorder(dt, Tree, ID)
  
  dt[, .(Tree, ID, Feature, Split, Yes, No, Missing, Leaf, LeafVal, Gain, Cover)]
}

#============================#
# Structure helpers (tdt → helpers)
#============================#
# Builds cached child→parent lookups, per-tree depth maps, and leaf paths.
build_structure_helpers <- function(tdt) {
  stopifnot(is.data.frame(tdt),
            all(c("Tree","ID","Yes","No","Missing","Leaf","Gain","Cover") %in% names(tdt)))
  tdt <- as.data.table(tdt)
  trees <- sort(unique(tdt$Tree))
  
  # child -> parent maps, per tree
  parent_env <- new.env(parent = emptyenv())
  for (tr in trees) {
    tt <- tdt[Tree == tr, .(ID, Yes, No, Missing)]
    max_id <- max(tt$ID, 0L)
    pmap <- rep(NA_integer_, max_id + 1L)
    for (col in c("Yes","No","Missing")) {
      ch <- tt[[col]]; ok <- !is.na(ch)
      if (any(ok)) pmap[ch[ok] + 1L] <- tt$ID[ok]
    }
    parent_env[[as.character(tr)]] <- pmap
  }
  
  split_gc <- tdt[Leaf == FALSE, .(Tree, ID, Gain, SplitCover = Cover)]
  setkey(split_gc, Tree, ID)
  
  depth_env  <- new.env(parent = emptyenv())
  maxd_env   <- new.env(parent = emptyenv())
  leaf_cache <- new.env(parent = emptyenv())
  
  get_depth_map <- function(tr) {
    key <- as.character(tr)
    dep <- depth_env[[key]]
    if (!is.null(dep)) return(dep)
    pmap <- parent_env[[key]]; if (is.null(pmap)) return(NULL)
    max_id <- length(pmap) - 1L
    dep <- rep(NA_integer_, max_id + 1L)
    dep[1] <- 0L
    
    # BFS queue (preallocated)
    q <- integer(max_id + 1L); head <- 1L; tail <- 1L; q[1] <- 0L
    while (head <= tail) {
      id <- q[head]; head <- head + 1L
      ch_idx <- which(pmap == id)
      if (length(ch_idx)) {
        new <- ch_idx - 1L
        unseen <- is.na(dep[new + 1L])
        if (any(unseen)) {
          dep[new[unseen] + 1L] <- dep[id + 1L] + 1L
          n_new <- sum(unseen)
          q[tail + seq_len(n_new)] <- new[unseen]
          tail <- tail + n_new
        }
      }
    }
    depth_env[[key]] <- dep
    maxd_env[[key]]  <- max(dep, na.rm = TRUE)
    dep
  }
  get_max_depth <- function(tr) {
    key <- as.character(tr)
    md <- maxd_env[[key]]; if (!is.null(md)) return(md)
    dep <- get_depth_map(tr); if (is.null(dep)) return(0L)
    md <- max(dep, na.rm = TRUE); maxd_env[[key]] <- md; md
  }
  # root→leaf node IDs + gain/cover along the path
  get_leaf_path <- function(tr, leaf_id) {
    tkey <- as.character(tr)
    env  <- leaf_cache[[tkey]]
    if (is.null(env)) { env <- new.env(parent = emptyenv()); leaf_cache[[tkey]] <- env }
    hit <- env[[as.character(leaf_id)]]
    if (!is.null(hit)) return(hit)
    
    pmap <- parent_env[[tkey]]
    if (is.null(pmap)) {
      out <- list(nodes=integer(0), gain=numeric(0), gcover=numeric(0))
      leaf_cache[[tkey]][[as.character(leaf_id)]] <- out; return(out)
    }
    # path length k
    cur <- leaf_id; k <- 0L
    while (!is.na(cur) && cur != 0L) { par <- pmap[cur + 1L]; if (is.na(par)) break; k <- k + 1L; cur <- par }
    if (k == 0L) {
      out <- list(nodes=integer(0), gain=numeric(0), gcover=numeric(0))
      leaf_cache[[tkey]][[as.character(leaf_id)]] <- out; return(out)
    }
    nodes <- integer(k); cur <- leaf_id
    for (pos in k:1) {
      par <- pmap[cur + 1L]; nodes[pos] <- par; cur <- par
      if (is.na(cur) || cur == 0L) break
    }
    sgc <- split_gc[.(tr, nodes), .(Gain, SplitCover)]
    out <- list(nodes = nodes,
                gain  = as.numeric(sgc$Gain),
                gcover= as.numeric(sgc$SplitCover))
    leaf_cache[[tkey]][[as.character(leaf_id)]] <- out
    out
  }
  
  list(get_depth_map = get_depth_map,
       get_max_depth = get_max_depth,
       get_leaf_path = get_leaf_path)
}

#============================#
# Modal extreme path per tree
#============================#
modal_extreme_paths <- function(leaves_mat, extreme_idx, helpers) {
  Tm <- ncol(leaves_mat)
  out <- vector("list", Tm)
  for (tt in seq_len(Tm)) {
    tr <- tt - 1L
    lt <- leaves_mat[, tt]
    tab <- sort(table(lt[extreme_idx]), decreasing = TRUE)
    if (!length(tab)) { out[[tt]] <- list(nodes=integer(0), gain=numeric(0), gcover=numeric(0)); next }
    modal_leaf <- as.integer(names(tab)[1])
    out[[tt]] <- helpers$get_leaf_path(tr, modal_leaf)
  }
  out
}

#============================#
# LCP helper
#============================#
.lcp_len <- function(a, b) {
  if (!length(a) || !length(b)) return(0L)
  k <- min(length(a), length(b)); i <- 0L
  while (i < k && a[i + 1L] == b[i + 1L]) i <- i + 1L
  i
}

#============================#
# Path similarity (per-SNP)
#============================#
# Returns per-SNP similarity + LCP metrics, plus reusable leaves/tdt/helpers.
compute_snp_similarity <- function(model, X, f_vec,
                                   lower_tail     = FALSE,
                                   extreme_k      = 1,
                                   progress_every = 100,
                                   show_progress  = TRUE) {
  .info("[compute_snp_similarity] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  stopifnot(inherits(model, "xgb.Booster"))
  X <- if (inherits(X, "dgCMatrix") || is.matrix(X)) X else as.matrix(X)
  stopifnot(is.numeric(f_vec), length(f_vec) == nrow(X))
  
  ## Extreme set --------------------------------------------------------------
  mu  <- mean(f_vec, na.rm = TRUE)
  sdv <- sd(f_vec,  na.rm = TRUE)
  if (!lower_tail) {
    thr <- mu + extreme_k * sdv; extreme_idx <- which(f_vec >= thr); tail_dir <- "high"
  } else {
    thr <- mu - extreme_k * sdv; extreme_idx <- which(f_vec <= thr); tail_dir <- "low"
  }
  if (!length(extreme_idx)) .stop_bad("No rows met the extreme-F criterion; adjust extreme_k/tail.")
  .info("Extreme-F (%s tail) threshold=%.6f; n_extreme=%d", tail_dir, thr, length(extreme_idx))
  
  ## Leaves (n x Tm) ----------------------------------------------------------
  leaves <- predict(model, X, predleaf = TRUE)
  if (is.null(dim(leaves))) leaves <- matrix(leaves, ncol = 1L)
  storage.mode(leaves) <- "integer"
  n  <- nrow(leaves); Tm <- ncol(leaves)
  .info("Leaves dim: %d x %d (rows x trees)", n, Tm)
  
  ## Trees + helpers ----------------------------------------------------------
  tdt <- fast_xgb_tree_dt(model)  
  if ((max(tdt$Tree) + 1L) != Tm) {
    .warn_bad("Tree count mismatch: dump=%d vs leaves=%d.", max(tdt$Tree)+1L, Tm)
  }
  .info("Building structure helpers...")
  helpers <- build_structure_helpers(tdt)
  
  .info("Computing modal extreme paths...")
  modal <- modal_extreme_paths(leaves, extreme_idx, helpers)
  
  ## Accumulators -------------------------------------------------------------
  score_raw    <- numeric(n)
  score_depthW <- numeric(n)
  lcp_sum      <- numeric(n)
  lcp_g_sum    <- numeric(n)
  lcp_gc_sum   <- numeric(n)
  lcp_n        <- integer(n)
  
  # Distribution-aware similarity & context
  ratio_sum      <- numeric(n)  # avg(prop_ext / p_leaf)
  loglr_sum      <- numeric(n)  # avg(log(prop_ext / p_leaf))
  zmean_sum      <- numeric(n)  # avg standardized residual
  leafdepth_sum  <- numeric(n)  # avg raw leaf depth (depth of leaf node)
  leafcover_sum  <- numeric(n)  # avg leaf cover (training count)
  leafval_sum    <- numeric(n)  # avg leaf value (xgb leaf prediction)
  
  eps <- 1e-12
  m   <- length(extreme_idx)
  
  ## Per-tree loop ---------------------------------------
  for (tt in seq_len(Tm)) {
    tr <- tt - 1L
    lt <- leaves[, tt]
    
    # Count extreme rows per leaf in this tree
    ft_tab <- table(lt[extreme_idx])
    if (!length(ft_tab)) {
      if (show_progress && (tt %% progress_every == 0L || tt == Tm))
        .info("processed %d/%d trees (%.1f%%)", tt, Tm, 100*tt/Tm)
      next
    }
    
    leaf_ids_ext <- as.integer(names(ft_tab))
    idx_ext <- match(lt, leaf_ids_ext); hit_ext <- !is.na(idx_ext)
    cnt_ext <- integer(n); if (any(hit_ext)) cnt_ext[hit_ext] <- as.integer(ft_tab[idx_ext[hit_ext]])
    prop_ext <- cnt_ext / m
    
    # similarity (raw)
    score_raw <- score_raw + prop_ext
    
    # depth weights
    dep_vals <- helpers$get_depth_map(tr)
    dep_max  <- helpers$get_max_depth(tr)
    depthW   <- if (!is.null(dep_vals) && dep_max > 0L) (dep_vals[lt + 1L] / dep_max) else rep(0, n)
    depthW[!hit_ext] <- 0
    score_depthW <- score_depthW + (depthW * prop_ext)
    
    # baseline global leaf frequency in this tree (for ratio/logLR/z)
    tab_all <- table(lt)
    leaf_ids_all <- as.integer(names(tab_all))
    idx_all <- match(lt, leaf_ids_all); hit_all <- !is.na(idx_all)
    cnt_all <- integer(n); if (any(hit_all)) cnt_all[hit_all] <- as.integer(tab_all[idx_all[hit_all]])
    p_leaf  <- cnt_all / n
    
    # standardized residual (binomial null), ratio, and log-likelihood
    var0 <- p_leaf * (1 - p_leaf)
    zval <- (cnt_ext - m * p_leaf) / sqrt(pmax(eps, m * var0))
    ratio <- prop_ext / pmax(eps, p_leaf)
    llr   <- ifelse(prop_ext > 0 & p_leaf > 0, log(prop_ext / p_leaf), 0)
    
    ratio_sum  <- ratio_sum  + ratio
    loglr_sum  <- loglr_sum  + llr
    zmean_sum  <- zmean_sum  + zval
    
    # modal path pieces 
    modal_nodes <- modal[[tt]]$nodes
    modal_g     <- modal[[tt]]$gain
    modal_gc    <- modal[[tt]]$gcover
    sum_modal_g  <- if (length(modal_g))  sum(modal_g,  na.rm = TRUE) else 0
    sum_modal_gc <- if (length(modal_gc)) sum(modal_gc, na.rm = TRUE) else 0
    
    if (any(hit_ext)) {
      uniq_nids <- unique(lt[hit_ext])
      # precompute LCP per unique leaf used in this tree
      lcp_cache <- vector("list", length(uniq_nids)); names(lcp_cache) <- as.character(uniq_nids)
      for (i in seq_along(uniq_nids)) {
        nid <- uniq_nids[i]
        lp <- helpers$get_leaf_path(tr, nid)
        pn <- lp$nodes
        L  <- if (length(pn) && length(modal_nodes)) .lcp_len(pn, modal_nodes) else 0L
        lcp_norm <- if (dep_max > 0L) (L / dep_max) else 0
        if (L > 0L) {
          sum_g_common  <- sum(lp$gain[seq_len(L)],   na.rm = TRUE)
          sum_gc_common <- sum(lp$gcover[seq_len(L)], na.rm = TRUE)
        } else {
          sum_g_common <- 0; sum_gc_common <- 0
        }
        lcp_gain_norm   <- if (sum_modal_g  > 0) sum_g_common  / sum_modal_g  else 0
        lcp_gcover_norm <- if (sum_modal_gc > 0) sum_gc_common / sum_modal_gc else 0
        lcp_cache[[i]] <- c(lcp_norm, lcp_gain_norm, lcp_gcover_norm)
      }
      # apply to all rows by their leaf
      for (i in seq_along(uniq_nids)) {
        nid <- uniq_nids[i]
        rows_i <- which(lt == nid)
        vals <- lcp_cache[[i]]
        lcp_sum[rows_i]    <- lcp_sum[rows_i]    + vals[1]
        lcp_g_sum[rows_i]  <- lcp_g_sum[rows_i]  + vals[2]
        lcp_gc_sum[rows_i] <- lcp_gc_sum[rows_i] + vals[3]
        lcp_n[rows_i]      <- lcp_n[rows_i] + 1L
      }
    }
    
    # leaf depth/cover/value summaries (cheap)
    if (!is.null(dep_vals)) leafdepth_sum <- leafdepth_sum + dep_vals[lt + 1L] else leafdepth_sum <- leafdepth_sum + 0
    leaf_rows <- tdt[Tree == tr & Leaf == TRUE, .(ID, LeafVal, LeafCover = Cover)]
    if (nrow(leaf_rows)) {
      data.table::setkey(leaf_rows, ID)
      LR <- leaf_rows[.(as.integer(lt))]
      leafcover_sum <- leafcover_sum + as.numeric(LR$LeafCover)
      leafval_sum   <- leafval_sum   + as.numeric(LR$LeafVal)
    }
    
    if (show_progress && (tt %% progress_every == 0L || tt == Tm)) {
      .info("processed %d/%d trees (%.1f%%)", tt, Tm, 100*tt/Tm)
      if ((tt %% 50L) == 0L) gc(FALSE)
    }
  }
  
  ## Final per-SNP metrics ----------------------------------------------------
  sim_raw        <- score_raw / Tm
  sim_depthW     <- score_depthW / Tm
  path_lcp_mean  <- ifelse(lcp_n > 0L, lcp_sum    / lcp_n, NA_real_)
  path_lcpG_mean <- ifelse(lcp_n > 0L, lcp_g_sum  / lcp_n, NA_real_)
  path_lcpGC_mean<- ifelse(lcp_n > 0L, lcp_gc_sum / lcp_n, NA_real_)
  
  sim_ratio      <- ratio_sum     / Tm
  sim_logLR      <- loglr_sum     / Tm
  sim_zmean      <- zmean_sum     / Tm
  mean_leaf_depth<- leafdepth_sum / Tm
  mean_leaf_cover<- leafcover_sum / Tm
  mean_leaf_value<- leafval_sum   / Tm
  
  summaries <- data.table::data.table(
    row                         = seq_len(n),
    F                           = as.numeric(f_vec),
    similarity                  = sim_raw,
    similarity_depthW           = sim_depthW,
    path_overlap_lcp_mean       = path_lcp_mean,
    lcp_gain_weighted_mean      = path_lcpG_mean,
    lcp_gainCover_weighted_mean = path_lcpGC_mean,
    sim_ratio                   = sim_ratio,
    sim_logLR                   = sim_logLR,
    sim_zmean                   = sim_zmean,
    mean_leaf_depth             = mean_leaf_depth,
    mean_leaf_cover             = mean_leaf_cover,
    mean_leaf_value             = mean_leaf_value
  )
  
  list(
    summaries   = summaries,
    leaves      = leaves,   # n x Tm integer
    tdt         = tdt,
    helpers     = helpers,
    extreme_idx = extreme_idx,
    tail_dir    = if (lower_tail) "low" else "high"
  )
}

#============================#
# Build per-leaf steps 
#============================#
# Yields long table: (Tree, leaf_id, depth, feature, direction, split, thresh_bin)
build_leaf_steps <- function(leaves,
                             tdt,
                             helpers,
                             # Binning control
                             binning         = c("digits", "adaptive"),
                             bin_digits      = 3L,      # used if binning == "digits"
                             # Adaptive-binning knobs (used if binning == "adaptive")
                             target_bins     = 8L,      # aim ~this many per-feature bins
                             min_per_bin     = 50L,     # avoid tiny bins across all trees/uses
                             winsor_prob     = 0.01,    # trim tails for stability (two-sided)
                             method          = c("fd", "quantile"),
                             # Back-compat shim (ignored unless provided)
                             min_bin         = NULL,    # deprecated; mapped to min_per_bin if given
                             # Performance
                             trees_per_batch = 200L,
                             progress_every  = 500L) {
  
  .info("[build_leaf_steps] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  stopifnot(is.matrix(leaves) || inherits(leaves, "dgCMatrix"))
  Tm <- ncol(leaves)
  
  # ---- normalize args / back-compat ----
  binning <- match.arg(binning)
  method  <- match.arg(method)
  if (!is.null(min_bin)) {
    .warn_bad("Argument 'min_bin' is deprecated; using it as min_per_bin = %s", as.integer(min_bin))
    min_per_bin <- as.integer(min_bin)
  }
  
  # ---- index tdt for fast joins ----
  tdt <- data.table::as.data.table(tdt)
  data.table::setkey(tdt, Tree, ID)
  
  # ---- helpers: build raw steps (no binning yet) ----
  build_raw_steps_batch <- function(tt) {
    tr <- tt - 1L
    lt <- leaves[, tt]
    uleaf <- sort(unique(as.integer(lt)))
    if (!length(uleaf)) return(data.table::data.table())
    
    tt_dt <- tdt[.(tr)]
    step_rows <- vector("list", length(uleaf)); si <- 1L
    
    for (leaf_id in uleaf) {
      lp <- helpers$get_leaf_path(tr, leaf_id)
      pn <- lp$nodes
      if (!length(pn)) { step_rows[[si]] <- NULL; si <- si + 1L; next }
      
      parents <- tt_dt[.(tr, pn)]
      yes  <- parents$Yes
      no   <- parents$No
      mis  <- parents$Missing
      feat <- parents$Feature
      splt <- parents$Split
      
      k <- length(pn)
      child_along <- integer(k)
      for (i in seq_len(k)) {
        next_child <- if (i < k) pn[i + 1L] else leaf_id
        child_along[i] <- next_child
      }
      direction <- ifelse(yes == child_along, "<",
                          ifelse(no == child_along, ">=", "missing"))
      depth <- seq.int(0L, k - 1L)
      
      step_rows[[si]] <- data.table::data.table(
        Tree       = tr,
        leaf_id    = as.integer(leaf_id),
        depth      = as.integer(depth),
        feature    = as.character(feat),
        direction  = as.character(direction),
        split      = suppressWarnings(as.numeric(splt))
      )
      si <- si + 1L
    }
    
    data.table::rbindlist(step_rows, use.names = TRUE, fill = TRUE)
  }
  
  # ---- stream trees in batches to control RAM ----
  out_list <- vector("list", ceiling(Tm / trees_per_batch)); oi <- 1L
  for (batch_start in seq(1L, Tm, by = trees_per_batch)) {
    batch_end <- min(batch_start + trees_per_batch - 1L, Tm)
    rows_batch <- vector("list", length = batch_end - batch_start + 1L); bi <- 1L
    
    for (tt in batch_start:batch_end) {
      rows_batch[[bi]] <- build_raw_steps_batch(tt)
      bi <- bi + 1L
      if ((tt %% progress_every) == 0L || tt == Tm) {
        message(sprintf("[build_leaf_steps] built up to tree %d / %d (%.1f%%)", tt, Tm, 100*tt/Tm))
      }
    }
    
    out_list[[oi]] <- data.table::rbindlist(rows_batch, use.names = TRUE, fill = TRUE)
    oi <- oi + 1L
    gc(FALSE)
  }
  
  LS <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (!nrow(LS)) {
    data.table::setkey(LS, Tree, leaf_id)
    return(LS[])
  }
  
  # ---- normalize columns ----
  LS <- LS[!is.na(feature) & nzchar(feature)]
  LS[, `:=`(
    Tree      = as.integer(Tree),
    leaf_id   = as.integer(leaf_id),
    depth     = as.integer(depth),
    feature   = as.character(feature),
    direction = as.character(direction),
    split     = as.numeric(split)
  )]
  
  # ---- binning: digits vs adaptive ----
  if (binning == "digits") {
    # classic behavior
    LS[, thresh_bin := base::signif(split, as.integer(bin_digits))]
    data.table::setkey(LS, Tree, leaf_id)
    return(LS[])
  }
  
  # ---- adaptive binning per feature ----
  # Winsorize tails, derive breaks via FD (or quantile), enforce min_per_bin by merging adjacent sparse bins.
  .info("[build_leaf_steps] computing adaptive bins per feature (target_bins=%d, min_per_bin=%d, method=%s, winsor=%.3g)",
        as.integer(target_bins), as.integer(min_per_bin), method, as.numeric(winsor_prob))
  
  LS_split <- LS[is.finite(split)]
  feats <- sort(unique(LS_split$feature))
  
  # envs to store per-feature breaks and midpoints
  breaks_env <- new.env(parent = emptyenv())
  mids_env   <- new.env(parent = emptyenv())
  
  propose_breaks <- function(v, max_bins, method, winsor) {
    v <- sort(v[is.finite(v)])
    if (length(v) <= 1L) return(range(v, na.rm = TRUE))
    lo <- suppressWarnings(stats::quantile(v, probs = winsor, names = FALSE))
    hi <- suppressWarnings(stats::quantile(v, probs = 1 - winsor, names = FALSE))
    v  <- v[v >= lo & v <= hi]
    if (length(v) <= 1L) return(range(v, na.rm = TRUE))
    
    if (method == "fd") {
      h <- 2 * stats::IQR(v) / (length(v)^(1/3))
      if (!is.finite(h) || h <= 0) h <- (max(v) - min(v)) / max(2L, max_bins)
      nb <- ceiling((max(v) - min(v)) / max(h, .Machine$double.eps))
      nb <- min(max(nb, 2L), max_bins)
      br <- pretty(v, n = nb)
    } else {
      nb <- max(2L, as.integer(max_bins))
      probs <- seq(0, 1, length.out = nb + 1L)
      br <- unique(stats::quantile(v, probs = probs, names = FALSE, type = 7))
      if (length(br) < 3L) br <- pretty(v, n = min(max_bins, 3L))
    }
    br <- sort(unique(as.numeric(br)))
    if (length(br) < 3L) br <- unique(c(min(v), median(v), max(v)))
    br
  }
  
  for (f in feats) {
    v <- LS_split[feature == f, split]
    if (!length(v)) next
    br <- propose_breaks(v, max_bins = as.integer(target_bins), method = method, winsor = as.numeric(winsor_prob))
    
    # enforce min_per_bin by merging sparse neighbors
    bin_idx <- findInterval(v, br, all.inside = TRUE)
    tab <- tabulate(bin_idx, nbins = length(br) - 1L)
    
    while ((any(tab < min_per_bin)) && (length(tab) > 2L)) {
      k <- which.min(tab)
      if (k == 1L) {
        rm_pos <- 2L
      } else if (k == length(tab)) {
        rm_pos <- length(tab)    # drop last internal break
      } else {
        rm_pos <- if (tab[k - 1L] <= tab[k + 1L]) k + 1L else k + 1L  # boundary on right of the smaller neighbor
      }
      br <- br[-rm_pos]
      bin_idx <- findInterval(v, br, all.inside = TRUE)
      tab <- tabulate(bin_idx, nbins = length(br) - 1L)
      # hard cap: keep bins <= target_bins
      if ((length(tab) > target_bins)) {
        k2 <- which.min(tab)
        if (k2 == 1L) br <- br[-2L] else br <- br[-(k2 + 1L)]
        bin_idx <- findInterval(v, br, all.inside = TRUE)
        tab <- tabulate(bin_idx, nbins = length(br) - 1L)
      }
    }
    
    mids <- (br[-1L] + br[-length(br)]) / 2
    breaks_env[[f]] <- br
    mids_env[[f]]   <- mids
  }
  
  # stamp thresh_bin by feature
  LS[, thresh_bin := {
    br <- breaks_env[[feature[1L]]]
    md <- mids_env[[feature[1L]]]
    if (is.null(br) || is.null(md) || !is.finite(split[1L])) {
      NA_real_
    } else {
      idx <- findInterval(split, br, all.inside = TRUE)
      as.numeric(md[idx])
    }
  }, by = feature]
  
  data.table::setkey(LS, Tree, leaf_id)
  LS[]
}

#============================#
# Keys (feature, direction, binned threshold)
#============================#
# Assign stable integer key_id to each unique split condition (feature, direction, thresh_bin),
# and map each (Tree, leaf_id, depth, ...) step to that key_id.
build_keys_from_leaf_steps <- function(leaf_steps) {
  stopifnot(data.table::is.data.table(leaf_steps) || is.data.frame(leaf_steps))
  LS <- data.table::as.data.table(leaf_steps)
  
  # Must already have thresh_bin from build_leaf_steps(..., binning="digits"/"adaptive")
  need <- c("Tree","leaf_id","feature","direction","thresh_bin")
  miss <- setdiff(need, names(LS))
  if (length(miss)) {
    stop(sprintf("[keys] leaf_steps missing: %s. Make sure build_leaf_steps() produced 'thresh_bin'.",
                 paste(miss, collapse=", ")))
  }
  
  LS <- LS[!is.na(feature) & nzchar(feature)]
  LS[, `:=`(
    Tree       = as.integer(Tree),
    leaf_id    = as.integer(leaf_id),
    depth      = as.integer(if ("depth" %in% names(LS)) depth else NA_integer_),
    feature    = as.character(feature),
    direction  = as.character(direction),
    thresh_bin = as.numeric(thresh_bin)
  )]
  
  # Unique conditions → key_id
  keys_tbl <- unique(LS[, .(feature, direction, thresh_bin)])
  data.table::setorder(keys_tbl, feature, direction, thresh_bin)
  keys_tbl[, key_id := .I]
  data.table::setkey(keys_tbl, feature, direction, thresh_bin)
  
  # Map steps to key_id (keep only what drivers need)
  LS_keys <- keys_tbl[
    LS, on = .(feature, direction, thresh_bin), nomatch = 0L
  ][, .(Tree, leaf_id, depth, key_id)]
  
  data.table::setkey(LS_keys, Tree, leaf_id)
  list(keys_tbl = keys_tbl[], LS_keys = LS_keys[])
}

#============================#
# Feature concentration (per SNP)
#============================#
# Input can be long (feature/depth) or wide (path_features list)
per_leaf_topk_share <- function(leaf_steps,
                                top_k = 2L,
                                progress_every = 5000L) {
  stopifnot(data.table::is.data.table(leaf_steps) || is.data.frame(leaf_steps))
  LS <- data.table::as.data.table(leaf_steps)
  
  has_long  <- "feature" %in% names(LS)
  has_list  <- "path_features" %in% names(LS)
  if (!has_long && !has_list) stop("[per_leaf_topk_share] Need 'feature' (long) or 'path_features' (list).")
  
  # Count occurrences per (Tree, leaf_id, feature)
  if (has_long) {
    counts <- LS[!is.na(feature) & nzchar(feature),
                 .(count = .N),
                 by = .(Tree = as.integer(Tree),
                        leaf_id = as.integer(leaf_id),
                        feature = as.character(feature))]
  } else {
    counts <- LS[, .(feature = unlist(path_features, use.names = FALSE)),
                 by = .(Tree = as.integer(Tree), leaf_id = as.integer(leaf_id))]
    counts <- counts[!is.na(feature) & nzchar(feature)][, .(count = .N), by = .(Tree, leaf_id, feature)]
  }
  
  # Share of top_k features along path
  setorder(counts, Tree, leaf_id, -count)
  perleaf <- counts[, {
    s <- sum(count)
    if (s == 0L) list(topk_share = NA_real_) else list(topk_share = sum(head(count, min(top_k, .N))) / s)
  }, by = .(Tree, leaf_id)]
  
  if (nrow(perleaf) > progress_every) {
    message(sprintf("[per_leaf_topk_share] computed for %s leaves.", format(nrow(perleaf), big.mark=",")))
  }
  setkey(perleaf, Tree, leaf_id)
  perleaf[]
}


# Map per-leaf topk_share → per-SNP by row-mean across trees in chunks.
feature_concentration_per_snp <- function(leaves,
                                          perleaf_share,
                                          chunk_rows = 2000L,
                                          progress_every = 10000L) {
  .info("[feature_concentration_per_snp] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  stopifnot(is.matrix(leaves) || inherits(leaves, "dgCMatrix"))
  P <- data.table::as.data.table(perleaf_share)
  stopifnot(all(c("Tree","leaf_id","topk_share") %in% names(P)))
  data.table::setkey(P, Tree, leaf_id)
  
  n  <- nrow(leaves); Tm <- ncol(leaves)
  res <- numeric(n)
  
  for (start in seq(1L, n, by = chunk_rows)) {
    end <- min(n, start + chunk_rows - 1L)
    idx <- start:end
    
    m <- matrix(NA_real_, nrow = length(idx), ncol = Tm)
    for (tt in seq_len(Tm)) {
      tr <- tt - 1L
      li <- as.integer(leaves[idx, tt])
      dt <- data.table::data.table(Tree = tr, leaf_id = li, ord = seq_along(li))
      data.table::setkey(dt, Tree, leaf_id)
      j <- P[dt]
      m[, tt] <- j$topk_share[order(j$ord)]
    }
    
    res[idx] <- rowMeans(m, na.rm = TRUE)
    if (((end - (start - 1L)) %% progress_every) == 0L || end == n) {
      message(sprintf("[feature_concentration_per_snp] %s / %s rows", format(end, big.mark=","), format(n, big.mark=",")))
    }
    rm(m); gc(FALSE)
  }
  
  data.table::data.table(row = seq_len(n), feature_conc_top2_share = res)
}

#============================#
# Enriched split conditions
#============================#
# Requires: build_keys_from_leaf_steps() output (LS_keys, keys_tbl)
compute_split_enrichment <- function(leaves,
                                     LS_keys,
                                     keys_tbl,
                                     extreme_rows,
                                     min_support_snp = 10L,
                                     rows_per_chunk  = 600L,
                                     progress_every  = 5000L) {
  .info("[compute_split_enrichment] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  stopifnot(is.matrix(leaves) || inherits(leaves, "dgCMatrix"))
  stopifnot(data.table::is.data.table(LS_keys), data.table::is.data.table(keys_tbl))
  
  data.table::setkey(LS_keys, Tree, leaf_id)
  
  n  <- nrow(leaves); Tm <- ncol(leaves)
  all_rows <- seq_len(n)
  extreme_rows <- as.integer(unique(extreme_rows))
  is_extreme <- logical(n); is_extreme[extreme_rows] <- TRUE
  
  N_total   <- n
  N_extreme <- length(extreme_rows)
  N_bg      <- N_total - N_extreme
  if (N_extreme <= 0L || N_bg <= 0L) stop("[drivers] Need both extreme and background SNPs.")
  
  # Compact accumulators (size = number of unique conditions)
  K <- nrow(keys_tbl)
  n_ext     <- integer(K)
  n_bg      <- integer(K)
  sum_d_ext <- numeric(K)
  sum_d_bg  <- numeric(K)
  count_d_ext <- integer(K)
  count_d_bg  <- integer(K)
  dmin_all  <- rep.int(.Machine$integer.max, K)
  dmax_all  <- rep.int(-.Machine$integer.max, K)
  
  # Stream SNPs in chunks
  row_chunks <- split(all_rows, ceiling(seq_along(all_rows) / rows_per_chunk))
  processed <- 0L
  
  for (rows_vec in row_chunks) {
    # (row, Tree, leaf_id) scaffold
    long_list <- vector("list", Tm)
    for (tt in seq_len(Tm)) {
      tr <- tt - 1L
      long_list[[tt]] <- data.table::data.table(
        row = rows_vec, Tree = tr, leaf_id = as.integer(leaves[rows_vec, tt])
      )
    }
    R <- data.table::rbindlist(long_list, use.names = TRUE)
    
    # Join to steps → keys
    S <- LS_keys[R, on = .(Tree, leaf_id), nomatch = 0L, allow.cartesian = TRUE]
    if (!nrow(S)) {
      processed <- processed + length(rows_vec)
      if (processed %% progress_every == 0L || processed >= n) {
        message(sprintf("[compute_split_enrichment] %s / %s SNPs (%.1f%%)",
                        format(processed, big.mark=","), format(n, big.mark=","), 100*processed/n))
      }
      next
    }
    if (!("depth" %in% names(S)) || all(is.na(S$depth))) {
      S[, depth := sequence(.N) - 1L, by = .(row, Tree, leaf_id)]
    }
    
    # Collapse to one (row, key_id) with mean depth
    A <- S[, .(mean_depth = mean(depth)), by = .(row, key_id)]
    if (!nrow(A)) next
    A[, is_ext := is_extreme[row]]
    
    # EXTREME
    A_ext <- A[is_ext == TRUE]
    if (nrow(A_ext)) {
      count <- A_ext[, .(count = .N), by = key_id]
      n_ext[count$key_id] <- n_ext[count$key_id] + count$count
      
      dep <- A_ext[, .(sum_depth = sum(mean_depth),
                       count_depth = .N,
                       dmin = floor(min(mean_depth)),
                       dmax = ceiling(max(mean_depth))), by = key_id]
      sum_d_ext[dep$key_id]  <- sum_d_ext[dep$key_id]  + dep$sum_depth
      count_d_ext[dep$key_id] <- count_d_ext[dep$key_id] + dep$count_depth
      dmin_all[dep$key_id]   <- pmin(dmin_all[dep$key_id], dep$dmin)
      dmax_all[dep$key_id]   <- pmax(dmax_all[dep$key_id], dep$dmax)
    }
    
    # BACKGROUND
    A_bg <- A[is_ext == FALSE]
    if (nrow(A_bg)) {
      count <- A_bg[, .(count = .N), by = key_id]
      n_bg[count$key_id] <- n_bg[count$key_id] + count$count
      
      dep <- A_bg[, .(sum_depth = sum(mean_depth),
                      count_depth = .N,
                      dmin = floor(min(mean_depth)),
                      dmax = ceiling(max(mean_depth))), by = key_id]
      sum_d_bg[dep$key_id]   <- sum_d_bg[dep$key_id]   + dep$sum_depth
      count_d_bg[dep$key_id]  <- count_d_bg[dep$key_id]  + dep$count_depth
      dmin_all[dep$key_id]   <- pmin(dmin_all[dep$key_id], dep$dmin)
      dmax_all[dep$key_id]   <- pmax(dmax_all[dep$key_id], dep$dmax)
    }
    
    processed <- processed + length(rows_vec)
    if (processed %% progress_every == 0L || processed >= n) {
      message(sprintf("[compute_split_enrichment] %s / %s SNPs (%.1f%%)",
                      format(processed, big.mark=","), format(n, big.mark=","), 100*processed/n))
      gc(FALSE)
    }
  }
  
  # Materialize
  dt <- data.table::data.table(
    key_id     = seq_len(K),
    n_extreme  = n_ext,
    n_bg       = n_bg,
    mean_depth_extreme = ifelse(n_ext > 0L, sum_d_ext / pmax(1L, count_d_ext), NA_real_),
    mean_depth_bg      = ifelse(n_bg  > 0L, sum_d_bg  / pmax(1L, count_d_bg),  NA_real_),
    depth_min          = ifelse(dmin_all == .Machine$integer.max, NA_integer_, dmin_all),
    depth_max          = ifelse(dmax_all == - .Machine$integer.max, NA_integer_, dmax_all)
  )
  
  # Filter by support
  dt <- dt[(n_extreme + n_bg) >= as.integer(min_support_snp)]
  if (!nrow(dt)) {
    warning("[compute_split_enrichment] No conditions passed min_support_snp.")
    return(dt)
  }
  
  # One-sided Fisher via hypergeometric
  N_extreme <- as.integer(N_extreme); N_bg <- as.integer(N_bg)
  m_support <- dt$n_extreme + dt$n_bg
  pval <- phyper(q = dt$n_extreme - 1L, m = N_extreme, n = N_bg, k = m_support, lower.tail = FALSE)
  
  setkey(keys_tbl, key_id)
  drivers <- keys_tbl[dt, on = .(key_id)]
  base_rate <- N_extreme / (N_extreme + N_bg)
  
  drivers[, `:=`(
    odds_ratio = ((n_extreme / pmax(1L, N_extreme - n_extreme)) /
                    ( n_bg      / pmax(1L, N_bg      - n_bg)) ),
    pval = pval,
    precision = n_extreme / (n_extreme + n_bg + 1e-12),
    lift      = if (base_rate > 0) (n_extreme / (n_extreme + n_bg + 1e-12)) / base_rate else NA_real_
  )]
  
  drivers[, qval := p.adjust(pval, "BH")]
  setorder(drivers, qval, -odds_ratio)
  drivers[]
}

#============================#
# Feature-level enrichment (fast view)
#============================#
# Build per-leaf feature view from keys, optionally by direction and depth limit.
build_feature_view_from_keys <- function(LS_keys, keys_tbl,
                                         by_direction = TRUE,
                                         depth_limit  = Inf) {
  stopifnot(data.table::is.data.table(LS_keys), data.table::is.data.table(keys_tbl))
  L  <- data.table::copy(LS_keys)[, .(Tree = as.integer(Tree),
                                      leaf_id = as.integer(leaf_id),
                                      depth = as.integer(depth),
                                      key_id = as.integer(key_id))]
  K  <- data.table::copy(keys_tbl)[, .(key_id = as.integer(key_id),
                                       feature = as.character(feature),
                                       direction = as.character(direction))]
  data.table::setkey(L, key_id)
  data.table::setkey(K, key_id)
  
  X <- K[L, nomatch = 0L][, .(Tree, leaf_id, depth, feature, direction)]
  if (is.finite(depth_limit)) X <- X[depth <= as.integer(depth_limit)]
  X <- X[!is.na(feature) & nzchar(feature)]
  if (!by_direction) X[, direction := NA_character_]
  
  feats_tbl <- unique(X[, .(feature, direction)])
  setorder(feats_tbl, feature, direction)
  feats_tbl[, feat_id := .I]
  setkey(feats_tbl, feature, direction)
  
  X  <- feats_tbl[X, on = .(feature, direction), nomatch = 0L][
    , .(mean_depth = mean(as.numeric(depth), na.rm = TRUE)),
    by = .(Tree = as.integer(Tree), leaf_id = as.integer(leaf_id), feat_id = as.integer(feat_id))]
  setkey(X, Tree, leaf_id)
  
  list(LS_feat = X[], feats_tbl = feats_tbl[, .(feat_id, feature, direction)][])
}

# Renamed from: compute_feature_enrichment_from_keys() (but now uses feature view)
# Streams SNPs and tallies feature presence + depth stats with fixed-size accumulators.
compute_feature_enrichment <- function(leaves,
                                       LS_feat,
                                       feats_tbl,
                                       extreme_rows,
                                       min_support_snp = 10L,
                                       rows_per_chunk  = 1000L,
                                       progress_every  = 5000L) {
  .info("[compute_feature_enrichment] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  stopifnot(is.matrix(leaves) || inherits(leaves, "dgCMatrix"))
  stopifnot(data.table::is.data.table(LS_feat), data.table::is.data.table(feats_tbl))
  setkey(LS_feat, Tree, leaf_id)
  
  n  <- nrow(leaves); Tm <- ncol(leaves)
  all_rows <- seq_len(n)
  extreme_rows <- as.integer(unique(extreme_rows))
  is_extreme <- logical(n); is_extreme[extreme_rows] <- TRUE
  N_extreme <- sum(is_extreme); N_bg <- n - N_extreme
  if (N_extreme <= 0L || N_bg <= 0L) stop("[compute_feature_enrichment] Need both extreme and background SNPs.")
  
  Fm <- max(LS_feat$feat_id)
  
  # Accumulators
  n_ext     <- integer(Fm)
  n_bg      <- integer(Fm)
  sum_d_ext <- numeric(Fm)
  sum_d_bg  <- numeric(Fm)
  count_d_ext <- integer(Fm)
  count_d_bg  <- integer(Fm)
  dmin_all  <- rep.int(.Machine$integer.max, Fm)
  dmax_all  <- rep.int(-.Machine$integer.max, Fm)
  
  row_chunks <- split(all_rows, ceiling(seq_along(all_rows) / rows_per_chunk))
  processed <- 0L
  
  for (rows_vec in row_chunks) {
    long_list <- vector("list", Tm)
    for (tt in seq_len(Tm)) {
      tr <- tt - 1L
      long_list[[tt]] <- data.table::data.table(
        row = rows_vec, Tree = tr, leaf_id = as.integer(leaves[rows_vec, tt])
      )
    }
    R <- data.table::rbindlist(long_list, use.names = TRUE)
    
    S <- LS_feat[R, on = .(Tree, leaf_id), nomatch = 0L, allow.cartesian = TRUE]
    if (!nrow(S)) {
      processed <- processed + length(rows_vec)
      if (processed %% progress_every == 0L || processed >= n) {
        message(sprintf("[compute_feature_enrichment] %s / %s SNPs (%.1f%%)",
                        format(processed, big.mark=","), format(n, big.mark=","), 100*processed/n))
      }
      next
    }
    
    A <- S[, .(mean_depth = mean(mean_depth, na.rm = TRUE)), by = .(row, feat_id)]
    if (!nrow(A)) {
      processed <- processed + length(rows_vec)
      if (processed %% progress_every == 0L || processed >= n) {
        message(sprintf("[compute_feature_enrichment] %s / %s SNPs (%.1f%%)",
                        format(processed, big.mark=","), format(n, big.mark=","), 100*processed/n))
      }
      next
    }
    A[, is_ext := is_extreme[row]]
    
    # EXTREME
    A_ext <- A[is_ext == TRUE]
    if (nrow(A_ext)) {
      count <- A_ext[, .(count = .N), by = feat_id]
      n_ext[count$feat_id] <- n_ext[count$feat_id] + count$count
      
      dep <- A_ext[, .(sum_depth = sum(mean_depth),
                       count_depth = .N,
                       dmin = floor(min(mean_depth)),
                       dmax = ceiling(max(mean_depth))), by = feat_id]
      sum_d_ext[dep$feat_id]   <- sum_d_ext[dep$feat_id] + dep$sum_depth
      count_d_ext[dep$feat_id] <- count_d_ext[dep$feat_id] + dep$count_depth
      dmin_all[dep$feat_id]    <- pmin(dmin_all[dep$feat_id], dep$dmin)
      dmax_all[dep$feat_id]    <- pmax(dmax_all[dep$feat_id], dep$dmax)
    }
    
    # BACKGROUND
    A_bg <- A[is_ext == FALSE]
    if (nrow(A_bg)) {
      count <- A_bg[, .(count = .N), by = feat_id]
      n_bg[count$feat_id] <- n_bg[count$feat_id] + count$count
      
      dep <- A_bg[, .(sum_depth = sum(mean_depth),
                      count_depth = .N,
                      dmin = floor(min(mean_depth)),
                      dmax = ceiling(max(mean_depth))), by = feat_id]
      sum_d_bg[dep$feat_id]    <- sum_d_bg[dep$feat_id] + dep$sum_depth
      count_d_bg[dep$feat_id]  <- count_d_bg[dep$feat_id] + dep$count_depth
      dmin_all[dep$feat_id]    <- pmin(dmin_all[dep$feat_id], dep$dmin)
      dmax_all[dep$feat_id]    <- pmax(dmax_all[dep$feat_id], dep$dmax)
    }
    
    processed <- processed + length(rows_vec)
    if (processed %% progress_every == 0L || processed >= n) {
      message(sprintf("[compute_feature_enrichment] %s / %s SNPs (%.1f%%)",
                      format(processed, big.mark=","), format(n, big.mark=","), 100*processed/n))
      gc(FALSE)
    }
  }
  
  out <- data.table::data.table(
    feat_id            = seq_len(Fm),
    n_extreme          = n_ext,
    n_bg               = n_bg,
    mean_depth_extreme = ifelse(count_d_ext > 0L, sum_d_ext / pmax(1L, count_d_ext), NA_real_),
    mean_depth_bg      = ifelse(count_d_bg  > 0L, sum_d_bg  / pmax(1L, count_d_bg),  NA_real_),
    depth_min          = ifelse(dmin_all == .Machine$integer.max, NA_integer_, dmin_all),
    depth_max          = ifelse(dmax_all == - .Machine$integer.max, NA_integer_, dmax_all)
  )
  
  out <- out[(n_extreme + n_bg) >= as.integer(min_support_snp)]
  if (!nrow(out)) {
    warning("[compute_feature_enrichment] No features passed min_support_snp.")
    return(out)
  }
  
  # Stats
  base_rate <- N_extreme / (N_extreme + N_bg)
  m_support <- out$n_extreme + out$n_bg
  pval <- phyper(q = out$n_extreme - 1L, m = N_extreme, n = N_bg, k = m_support, lower.tail = FALSE)
  
  out[, `:=`(
    odds_ratio = ((n_extreme / pmax(1L, N_extreme - n_extreme)) /
                    ( n_bg      / pmax(1L, N_bg      - n_bg)) ),
    pval       = pval,
    precision  = n_extreme / pmax(1e-12, n_extreme + n_bg),
    lift       = if (base_rate > 0) (n_extreme / pmax(1e-12, n_extreme + n_bg)) / base_rate else NA_real_,
    qval       = p.adjust(pval, "BH"),
    N_extreme_total = N_extreme,
    N_bg_total      = N_bg
  )]
  
  setkey(feats_tbl, feat_id)
  out <- feats_tbl[out, on = .(feat_id)]
  setorder(out, qval, -odds_ratio)
  out[]
}

#============================#
# Condition context → feature summaries
#============================#
# How many extreme SNPs hit each (Tree, leaf)
make_extreme_leaf_weights <- function(leaves, extreme_rows) {
  stopifnot(is.matrix(leaves) || inherits(leaves, "dgCMatrix"))
  extreme_rows <- as.integer(unique(extreme_rows))
  Tm <- ncol(leaves)
  out <- vector("list", Tm)
  for (tt in seq_len(Tm)) {
    tr <- tt - 1L
    tab <- table(as.integer(leaves[extreme_rows, tt]))
    if (length(tab)) {
      out[[tt]] <- data.table(Tree = tr,
                              leaf_id = as.integer(names(tab)),
                              w_ext = as.integer(tab))
    }
  }
  ans <- rbindlist(out, use.names = TRUE, fill = TRUE)
  setkey(ans, Tree, leaf_id)
  ans[]
}

# Per-condition context summary 
# leaf_steps must have: Tree, leaf_id, depth, feature, direction, thresh_bin
# Fast context summarizer using non-equi joins + chunked streaming
# Fast neighbor context summarizer: path scans, no joins.
# Counts labels that occur within +/- hops of each target condition along the same path.
summarize_condition_context <- function(drivers,
                                        leaf_steps,
                                        leaves,
                                        extreme_rows,
                                        context_granularity = c("feat_dir","feature","full"),
                                        hops             = 1L,   # window half-width around target (1=immediate)
                                        top_k            = 4L,
                                        rows_per_chunk   = 1500L,
                                        progress_every   = 5000L) {
  context_granularity <- match.arg(context_granularity)
  hops <- as.integer(hops); if (hops < 1L) hops <- 1L
  
  # ---- Normalize inputs ----
  LS <- data.table::as.data.table(leaf_steps)
  need <- c("Tree","leaf_id","depth","feature","direction")
  miss <- setdiff(need, names(LS))
  if (length(miss)) stop(sprintf("[neighbors-fast] leaf_steps missing: %s", paste(miss, collapse=", ")))
  if (!("thresh_bin" %in% names(LS))) {
    if (!("split" %in% names(LS))) stop("[neighbors-fast] need either 'thresh_bin' or 'split' for full granularity")
    LS[, thresh_bin := as.numeric(split)]
  }
  LS[, `:=`(
    Tree      = as.integer(Tree),
    leaf_id   = as.integer(leaf_id),
    depth     = as.integer(depth),
    feature   = as.character(feature),
    direction = as.character(direction)
  )]
  data.table::setorder(LS, Tree, leaf_id, depth)  # ensure path order
  
  D <- data.table::as.data.table(drivers)
  needD <- c("feature","direction","thresh_bin","n_extreme","n_bg","odds_ratio","qval")
  if (!all(needD %in% names(D))) {
    stop(sprintf("[neighbors-fast] drivers missing: %s", paste(setdiff(needD, names(D)), collapse=", ")))
  }
  
  # ---- Collapse drivers to one row per (feature, direction) (best bin) ----
  data.table::setorder(D, feature, direction, qval, -odds_ratio)
  Dbest <- D[, .SD[1L], by = .(feature, direction)]
  
  # ---- Label mapping per granularity ----
  mk_label <- switch(
    context_granularity,
    feature  = function(f, d, b) f,
    feat_dir = function(f, d, b) paste0(f, " ", d),
    full     = function(f, d, b) paste0(f, " ", d, " ", format(b, trim=TRUE))
  )
  
  LS[, ctx_label := mk_label(feature, direction, thresh_bin)]
  Dbest[, tgt_label := mk_label(feature, direction, thresh_bin)]
  
  # Integer encode all labels for speed
  all_labels <- sort(unique(c(LS$ctx_label, Dbest$tgt_label)))
  label_to_id <- structure(seq_along(all_labels), names = all_labels)
  id_to_label <- structure(all_labels, names = seq_along(all_labels))
  
  LS[, `:=`(
    ctx_id  = as.integer(label_to_id[ctx_label]),
    path_id = paste(Tree, leaf_id, sep = "#")
  )]
  Dbest[, tgt_id := as.integer(label_to_id[tgt_label])]
  Dbest[, tgt_key := paste0(feature, "\t", direction)]
  
  # Map target ids -> display keys (feature \t direction)
  id_to_tgt_key <- rep(NA_character_, length(all_labels))
  id_to_tgt_key[Dbest$tgt_id] <- Dbest$tgt_key
  
  # Optional: pre-slim LS to features that appear in drivers (speeds things further)
  keep_feats <- unique(Dbest$feature)
  LS <- LS[feature %in% keep_feats]
  
  # ---- Precompute path sequences per (Tree, leaf_id) for O(1) lookup ----
  # List-of-environments: one per tree, mapping leaf_id -> integer vector of ctx_id
  trees <- sort(unique(LS$Tree))
  seq_maps <- vector("list", length(trees)); names(seq_maps) <- as.character(trees)
  for (tr in trees) {
    LTr <- LS[Tree == tr, .(leaf_id, ctx_id)]
    # pack as list of integer vectors per leaf_id (already depth-sorted)
    split_map <- split(LTr$ctx_id, LTr$leaf_id)
    env <- new.env(parent = emptyenv())
    # store as integer vectors
    for (nm in names(split_map)) env[[nm]] <- as.integer(split_map[[nm]])
    seq_maps[[as.character(tr)]] <- env
  }
  
  # ---- Stream extreme rows, scan paths ----
  n <- nrow(leaves); Tm <- ncol(leaves)
  extreme_rows <- as.integer(unique(extreme_rows))
  extreme_rows <- extreme_rows[extreme_rows >= 1L & extreme_rows <= n]
  if (!length(extreme_rows)) stop("[neighbors-fast] extreme_rows empty")
  
  row_chunks <- split(extreme_rows, ceiling(seq_along(extreme_rows)/rows_per_chunk))
  
  # tallies
  occ_env  <- new.env(parent = emptyenv())  # tgt_key -> occurrences
  pre_env  <- new.env(parent = emptyenv())  # paste(tgt_key, ctx_id, sep="\t") -> count
  post_env <- new.env(parent = emptyenv())
  
  bump1 <- function(env, key, inc = 1L) {
    v <- env[[key]]
    if (is.null(v)) env[[key]] <- inc else env[[key]] <- v + inc
  }
  
  processed <- 0L
  for (rows_vec in row_chunks) {
    for (tt in seq_len(Tm)) {
      tr <- tt - 1L
      env <- seq_maps[[as.character(tr)]]
      if (is.null(env)) next
      li <- as.integer(leaves[rows_vec, tt])
      
      for (i in seq_along(rows_vec)) {
        s <- env[[as.character(li[i])]]
        if (is.null(s) || !length(s)) next
        
        # positions where sequence item is a target
        pos <- which(!is.na(id_to_tgt_key[s]))
        if (!length(pos)) next
        
        # scan neighbors around each target position
        for (p in pos) {
          tgt_key <- id_to_tgt_key[s[p]]
          bump1(occ_env, tgt_key, 1L)
          
          # predecessor window
          if (p > 1L) {
            start <- max(1L, p - hops)
            pre_ids <- unique(s[start:(p-1L)])
            pre_ids <- pre_ids[pre_ids != s[p]]
            if (length(pre_ids)) {
              for (cid in pre_ids) bump1(pre_env, paste0(tgt_key, "\t", cid), 1L)
            }
          }
          # successor window
          if (p < length(s)) {
            end <- min(length(s), p + hops)
            post_ids <- unique(s[(p+1L):end])
            post_ids <- post_ids[post_ids != s[p]]
            if (length(post_ids)) {
              for (cid in post_ids) bump1(post_env, paste0(tgt_key, "\t", cid), 1L)
            }
          }
        }
      }
    }
    processed <- processed + length(rows_vec)
    if (processed %% progress_every == 0L || processed >= length(extreme_rows)) {
      message(sprintf("[neighbors-fast] %s / %s extreme SNPs (%.1f%%)",
                      format(processed, big.mark=","), format(length(extreme_rows), big.mark=","),
                      100*processed/length(extreme_rows)))
      gc(FALSE)
    }
  }
  
  # ---- Materialize tallies ----
  tk <- ls(envir = occ_env)
  if (!length(tk)) {
    warning("[neighbors-fast] no target occurrences observed")
    return(data.table::data.table())
  }
  OCC <- data.table::data.table(
    tgt_key = tk,
    total_occ = vapply(tk, function(k) occ_env[[k]], integer(1))
  )
  
  unpack_side <- function(env) {
    ks <- ls(envir = env)
    if (!length(ks)) return(data.table::data.table(tgt_key=character(), ctx_id=integer(), count=integer()))
    sp <- strsplit(ks, "\t", fixed = TRUE)
    data.table::data.table(
      tgt_key = vapply(sp, `[[`, character(1), 1L),
      ctx_id  = as.integer(vapply(sp, `[[`, character(1), 2L)),
      count   = vapply(ks, function(k) env[[k]], integer(1))
    )
  }
  PRE  <- unpack_side(pre_env)
  POST <- unpack_side(post_env)
  if (nrow(PRE))  PRE[,  ctx_label := id_to_label[as.character(ctx_id)]]
  if (nrow(POST)) POST[, ctx_label := id_to_label[as.character(ctx_id)]]
  
  summarize_top <- function(dt, occ) {
    if (!nrow(dt)) return(data.table::data.table(tgt_key=character(), summary=character()))
    M <- merge(dt[, .(tgt_key, ctx_label, count)], occ, by = "tgt_key", all.x = TRUE, sort = FALSE)
    M[, pct := 100 * count / pmax(1L, total_occ)]
    M <- M[order(tgt_key, -pct)]
    M[, rn := seq_len(.N), by = tgt_key]
    M <- M[rn <= top_k]
    M[, summary := paste0(ctx_label, " (", sprintf("%.0f%%", pct), ")", collapse = " | "), by = tgt_key]
    unique(M[, .(tgt_key, summary)])
  }
  PRE_top  <- summarize_top(PRE,  OCC); data.table::setnames(PRE_top,  "summary", "common_predecessors")
  POST_top <- summarize_top(POST, OCC); data.table::setnames(POST_top, "summary", "common_successors")
  
  OUT <- merge(OCC, PRE_top,  by = "tgt_key", all.x = TRUE, sort = FALSE)
  OUT <- merge(OUT, POST_top, by = "tgt_key", all.x = TRUE, sort = FALSE)
  
  # Split tgt_key back out
  sp <- tstrsplit(OUT$tgt_key, "\t", fixed = TRUE)
  OUT[, `:=`(feature = sp[[1]], direction = sp[[2]])]
  OUT[, tgt_key := NULL]
  
  # Merge best stats from Dbest so you see OR/q/lift next to context
  keep <- c("feature","direction","odds_ratio","qval","lift","precision","n_extreme","n_bg")
  keep <- keep[keep %in% names(Dbest)]
  OUT <- Dbest[, ..keep][OUT, on = .(feature, direction)]
  data.table::setorder(OUT, qval, -odds_ratio, feature, direction)
  OUT[]
}     


# Harvest exact tree-path rules (one rule per (Tree, leaf_id)), then score them.
# Arguments:
#   leaves         : integer matrix n x Tm from predict(..., predleaf=TRUE)
#   leaf_steps     : data.table with columns: Tree, leaf_id, depth, feature, direction, split/thresh_bin
#   extreme_rows   : integer vector of extreme SNP row ids (1..n)
#   max_depth      : keep only paths with <= max_depth splits (depth counted as number of conditions)
#   min_support    : keep only rules seen in at least this many SNPs (ext+bg)
#   trees_subset   : optional integer vector of tree indices (0-based) to limit harvesting (e.g., top-KL trees)
#   pool_identical : if TRUE, merge rules with identical ordered condition strings across trees
#   progress_every : progress heartbeat
harvest_path_rules <- function(leaves,
                               leaf_steps,
                               extreme_rows,
                               max_depth      = 4L,
                               min_support    = 20L,
                               trees_subset   = NULL,   # e.g., select_trees_by_separation(...)$Tree
                               pool_identical = TRUE,
                               progress_every = 500L) {
  stopifnot(is.matrix(leaves) || inherits(leaves, "dgCMatrix"))
  LS <- data.table::as.data.table(leaf_steps)
  need <- c("Tree","leaf_id","depth","feature","direction")
  if (!all(need %in% names(LS))) stop(sprintf("[harvest] leaf_steps missing: %s", paste(setdiff(need, names(LS)), collapse=", ")))
  # Need a bin label; prefer thresh_bin; fall back to split if needed
  if (!("thresh_bin" %in% names(LS))) {
    if (!("split" %in% names(LS))) stop("[harvest] leaf_steps needs either 'thresh_bin' or raw 'split'")
    LS[, thresh_bin := as.numeric(split)]
  } else {
    LS[, thresh_bin := as.numeric(thresh_bin)]
  }
  LS[, `:=`(Tree = as.integer(Tree),
            leaf_id = as.integer(leaf_id),
            depth   = as.integer(depth),
            feature = as.character(feature),
            direction = as.character(direction))]
  data.table::setorder(LS, Tree, leaf_id, depth)
  data.table::setkey(LS, Tree, leaf_id)
  
  n  <- nrow(leaves); Tm <- ncol(leaves)
  extreme_rows <- as.integer(unique(extreme_rows))
  is_extreme <- logical(n); is_extreme[extreme_rows] <- TRUE
  N_extreme <- sum(is_extreme); N_bg <- n - N_extreme
  if (N_extreme <= 0L || N_bg <= 0L) stop("[harvest] Need both extreme and background SNPs.")
  
  # Optional tree restriction (0-based)
  all_trees0 <- 0:(Tm-1L)
  if (!is.null(trees_subset)) {
    trees_subset <- as.integer(sort(unique(trees_subset)))
    trees_subset <- trees_subset[trees_subset %in% all_trees0]
    if (!length(trees_subset)) stop("[harvest] trees_subset empty after filtering")
    use_tt <- (trees_subset + 1L) # convert to 1-based column index into leaves
  } else {
    use_tt <- seq_len(Tm)
  }
  
  # Build per-(Tree, leaf) ordered rule strings once (depth-filtered)
  # Rule string preserves order: "feat dir bin > feat dir bin > ..."
  LS_rule <- LS[depth < max_depth]  # keep <= (max_depth-1) depth indices => <= max_depth splits
  if (!nrow(LS_rule)) stop("[harvest] No steps under the requested max_depth.")
  LS_rule[, cond_lbl := paste0(feature, " ", direction, " ", format(thresh_bin, trim = TRUE))]
  PATHS <- LS_rule[, .(
    rule_len = .N,
    rule_str = paste(cond_lbl, collapse = " > "),
    depth_min = min(depth),
    depth_max = max(depth)
  ), by = .(Tree, leaf_id)]
  data.table::setkey(PATHS, Tree, leaf_id)
  
  # Collect stats per path
  out_list <- vector("list", length(use_tt))
  oi <- 1L
  
  for (tt in use_tt) {
    tr <- tt - 1L  # 0-based tree index
    # SNP-to-leaf column
    li_all <- as.integer(leaves[, tt])
    
    # counts for ALL and EXTREME on this tree
    tab_all <- table(li_all)
    li_ext  <- as.integer(leaves[is_extreme, tt])
    tab_ext <- table(li_ext)
    
    # Assemble counts into a DT keyed by (Tree, leaf_id)
    leaf_ids <- as.integer(names(tab_all))
    cnt_all  <- as.integer(unname(tab_all))
    cnt_ext  <- integer(length(leaf_ids))
    if (length(tab_ext)) {
      match_idx <- match(leaf_ids, as.integer(names(tab_ext)))
      hit <- !is.na(match_idx)
      cnt_ext[hit] <- as.integer(unname(tab_ext)[match_idx[hit]])
    }
    DTc <- data.table::data.table(Tree = tr,
                                  leaf_id = leaf_ids,
                                  n_all = cnt_all,
                                  n_extreme = cnt_ext)
    DTc[, n_bg := pmax(0L, n_all - n_extreme)]
    
    # Join to the rule strings for this tree
    DT <- PATHS[DTc, on = .(Tree, leaf_id), nomatch = 0L]
    if (!nrow(DT)) { oi <- oi + 1L; next }
    
    # Support filter
    DT <- DT[(n_extreme + n_bg) >= as.integer(min_support)]
    if (!nrow(DT)) { oi <- oi + 1L; next }
    
    # Stats per path
    base_rate <- N_extreme / (N_extreme + N_bg)  # scalar
    
    support <- DT$n_extreme + DT$n_bg
    # one-sided hypergeom tail: P[X >= n_extreme]
    pval <- phyper(q = DT$n_extreme - 1L, m = N_extreme, n = N_bg, k = support, lower.tail = FALSE)
    
    # Haldane–Anscombe corrected OR for stability
    a <- DT$n_extreme
    b <- N_extreme - DT$n_extreme
    c <- DT$n_bg
    d <- N_bg - DT$n_bg
    or_ha <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    
    # SAFE sequential assignments (no self-references in same :=)
    DT[, precision := n_extreme / pmax(1e-12, n_extreme + n_bg)]
    DT[, lift      := if (base_rate > 0) precision / base_rate else rep(NA_real_, .N)]
    DT[, `:=`(odds_ratio = or_ha, pval = pval)]
    
    out_list[[oi]] <- DT[, .(Tree, leaf_id, rule_len, rule_str,
                             n_extreme, n_bg, precision, lift, odds_ratio,
                             pval, depth_min, depth_max)]
    oi <- oi + 1L
    
    if (((oi-1L) %% progress_every) == 0L || (tt == tail(use_tt, 1L))) {
      message(sprintf("[harvest] processed tree %d / %d (%.1f%%)",
                      tt, Tm, 100*(oi-1L)/length(use_tt)))
      gc(FALSE)
    }
  }
  
  rules <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (!nrow(rules)) {
    warning("[harvest] No rules after filtering; try lowering min_support or increasing max_depth.")
    return(data.table::data.table())
  }
  
  # Optional: pool identical rule strings across trees (sum across trees), then recount unique SNPs
  if (isTRUE(pool_identical)) {
    # 1) Gather the set of (Tree, leaf_id) per rule_str
    RL <- rules[, .(Tree, leaf_id, rule_str)]
    # 2) For each rule_str, compute the union of rows across all (Tree,leaf_id) that realize it
    pooled_list <- vector("list", length = uniqueN(RL$rule_str))
    ui <- 1L
    for (rs in unique(RL$rule_str)) {
      pairs <- RL[rule_str == rs, .(Tree, leaf_id)]
      if (!nrow(pairs)) next
      
      # Collect row indices that fall into any of these leaves across their trees
      rows_hit <- integer(0)
      for (ii in seq_len(nrow(pairs))) {
        tr0 <- pairs$Tree[ii]                  # 0-based tree id
        tt  <- tr0 + 1L                        # column in 'leaves'
        lid <- pairs$leaf_id[ii]
        rows_hit <- c(rows_hit, which(as.integer(leaves[, tt]) == lid))
      }
      if (length(rows_hit)) rows_hit <- unique(rows_hit)
      
      # Unique support counts
      n_all <- length(rows_hit)
      if (n_all == 0L) {
        n_ext <- 0L; n_bg <- 0L
      } else {
        ext_flags <- is_extreme[rows_hit]
        n_ext <- sum(ext_flags)
        n_bg  <- n_all - n_ext
      }
      
      pooled_list[[ui]] <- data.table::data.table(
        rule_str  = rs,
        rule_len  = unique(rules[rule_str == rs, rule_len])[1L],
        n_extreme = as.integer(n_ext),
        n_bg      = as.integer(n_bg),
        depth_min = min(rules[rule_str == rs, depth_min], na.rm = TRUE),
        depth_max = max(rules[rule_str == rs, depth_max], na.rm = TRUE)
      )
      ui <- ui + 1L
    }
    rules <- data.table::rbindlist(pooled_list, use.names = TRUE, fill = TRUE)
    if (!nrow(rules)) {
      warning("[harvest] No rules after pooling recount; try lowering min_support or increasing max_depth.")
      return(data.table::data.table())
    }
    
    # Recompute stats SAFELY on unique counts
    base_rate <- N_extreme / (N_extreme + N_bg)
    support   <- rules$n_extreme + rules$n_bg
    
    # Guard against any lingering invalids
    ok <- is.finite(support) & support >= 0 &
      is.finite(rules$n_extreme) & rules$n_extreme >= 0 &
      is.finite(N_extreme) & is.finite(N_bg) & N_extreme >= 0 & N_bg >= 0 &
      support <= (N_extreme + N_bg) &
      rules$n_extreme <= support
    
    pval <- rep(NA_real_, nrow(rules))
    if (any(ok)) {
      pval[ok] <- phyper(q = rules$n_extreme[ok] - 1L,
                         m = N_extreme, n = N_bg, k = support[ok],
                         lower.tail = FALSE)
    }
    
    # Haldane–Anscombe OR and precision/lift on the unique counts
    a <- rules$n_extreme
    b <- N_extreme - rules$n_extreme
    c <- rules$n_bg
    d <- N_bg - rules$n_bg
    or_ha <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    
    rules[, precision := n_extreme / pmax(1e-12, n_extreme + n_bg)]
    rules[, lift      := if (base_rate > 0) precision / base_rate else NA_real_]
    rules[, `:=`(odds_ratio = or_ha, pval = pval)]
  }
  
  # FDR and ordering 
  rules[, qval := p.adjust(pval, method = "BH")]
  data.table::setorderv(
    rules,
    cols  = c("qval", "lift", "odds_ratio", "precision"),
    order = c( 1L,   -1L,     -1L,          -1L)
  )
  rules[]
}


#============================#
# Train-reference similarity 
#============================#
# Build a per-tree reference from TRAIN extremes (counts + modal path).
build_extreme_reference_from_leaves <- function(leaves_train,
                                                extreme_rows_train,
                                                helpers,
                                                tdt) {
  stopifnot(is.matrix(leaves_train) || inherits(leaves_train, "dgCMatrix"))
  extreme_rows_train <- as.integer(unique(extreme_rows_train))
  n_ext <- length(extreme_rows_train)
  if (n_ext <= 0L) stop("[ref] No training extremes provided.")
  
  Tm <- ncol(leaves_train)
  counts_by_leaf <- vector("list", Tm)
  modal          <- vector("list", Tm)
  
  for (tt in seq_len(Tm)) {
    tr <- tt - 1L
    lt_ext <- as.integer(leaves_train[extreme_rows_train, tt])
    if (length(lt_ext)) {
      tab <- sort(table(lt_ext), decreasing = TRUE)
      count_vec <- as.integer(tab)
      names(count_vec) <- names(tab)
      counts_by_leaf[[tt]] <- count_vec
      modal_leaf <- as.integer(names(tab)[1L])
      modal[[tt]] <- helpers$get_leaf_path(tr, modal_leaf)
    } else {
      counts_by_leaf[[tt]] <- integer(0)
      modal[[tt]] <- list(nodes = integer(0), gain = numeric(0), gcover = numeric(0))
    }
  }
  
  list(
    n_extreme      = n_ext,
    counts_by_leaf = counts_by_leaf,
    modal          = modal,
    Tm             = Tm
  )
}

# Compute per-SNP similarity/LCP metrics for NEW rows using a TRAIN reference.
compute_similarity_against_reference <- function(leaves_new,
                                                 tdt,
                                                 helpers,
                                                 ref_ctx,
                                                 extreme_rows_eval = NULL,
                                                 show_progress = TRUE,
                                                 progress_every = 100L) {
  .info("[compute_similarity_against_reference] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  stopifnot(is.matrix(leaves_new) || inherits(leaves_new, "dgCMatrix"))
  n  <- nrow(leaves_new)
  Tm <- ncol(leaves_new)
  if (Tm != ref_ctx$Tm) stop("[ref-score] Tree count mismatch vs reference.")
  
  n_ext_train <- as.numeric(ref_ctx$n_extreme)
  if (n_ext_train <= 0) stop("[ref-score] Reference has zero extremes.")
  
  score_raw    <- numeric(n)
  score_depthW <- numeric(n)
  lcp_sum      <- numeric(n)
  lcp_g_sum    <- numeric(n)
  lcp_gc_sum   <- numeric(n)
  lcp_n        <- integer(n)
  
  pb <- if (isTRUE(show_progress)) utils::txtProgressBar(min = 0, max = Tm, style = 3) else NULL
  on.exit(try(close(pb), silent = TRUE), add = TRUE)
  
  for (tt in seq_len(Tm)) {
    tr <- tt - 1L
    lt <- as.integer(leaves_new[, tt])
    
    # similarity from TRAIN extreme counts in this tree
    count_vec_train <- ref_ctx$counts_by_leaf[[tt]]
    if (length(count_vec_train)) {
      leaf_ids_train <- as.integer(names(count_vec_train))
      idx <- match(lt, leaf_ids_train)
      count_vec <- integer(n); hit <- !is.na(idx)
      if (any(hit)) count_vec[hit] <- count_vec_train[idx[hit]]
      prop_ext <- count_vec / n_ext_train
    } else {
      prop_ext <- numeric(n); hit <- rep(FALSE, n)
    }
    
    score_raw <- score_raw + prop_ext
    
    # depth weighting uses test leaf depth
    dep_vals <- helpers$get_depth_map(tr)
    dep_max  <- helpers$get_max_depth(tr)
    depthW   <- if (!is.null(dep_vals) && dep_max > 0L) (dep_vals[lt + 1L] / dep_max) else rep(0, n)
    depthW[!hit] <- 0
    score_depthW <- score_depthW + depthW * prop_ext
    
    # LCP vs TRAIN-modal path
    modal_nodes <- ref_ctx$modal[[tt]]$nodes
    modal_g     <- ref_ctx$modal[[tt]]$gain
    modal_gc    <- ref_ctx$modal[[tt]]$gcover
    sum_modal_g  <- if (length(modal_g))  sum(modal_g,  na.rm = TRUE) else 0
    sum_modal_gc <- if (length(modal_gc)) sum(modal_gc, na.rm = TRUE) else 0
    
    if (any(hit) && length(modal_nodes)) {
      uniq_nids <- unique(lt[hit])
      for (nid in uniq_nids) {
        rows_i <- which(lt == nid)
        lp     <- helpers$get_leaf_path(tr, nid)
        pn     <- lp$nodes
        l <- 0L
        if (length(pn)) {
          k <- min(length(pn), length(modal_nodes))
          while (l < k && pn[l + 1L] == modal_nodes[l + 1L]) l <- l + 1L
        }
        lcp_norm <- if (dep_max > 0L) (l / dep_max) else 0
        if (l > 0L) {
          sum_g_common  <- sum(lp$gain[seq_len(l)],   na.rm = TRUE)
          sum_gc_common <- sum(lp$gcover[seq_len(l)], na.rm = TRUE)
        } else {
          sum_g_common <- 0; sum_gc_common <- 0
        }
        lcp_gain_norm   <- if (sum_modal_g  > 0) sum_g_common  / sum_modal_g  else 0
        lcp_gcover_norm <- if (sum_modal_gc > 0) sum_gc_common / sum_modal_gc else 0
        
        lcp_sum[rows_i]    <- lcp_sum[rows_i]    + lcp_norm
        lcp_g_sum[rows_i]  <- lcp_g_sum[rows_i]  + lcp_gain_norm
        lcp_gc_sum[rows_i] <- lcp_gc_sum[rows_i] + lcp_gcover_norm
        lcp_n[rows_i]      <- lcp_n[rows_i] + 1L
      }
    }
    
    if (!is.null(pb)) utils::setTxtProgressBar(pb, tt)
    if (isTRUE(show_progress) && (tt %% progress_every == 0L || tt == Tm)) {
      message(sprintf("[ref-score] %d / %d trees (%.1f%%)", tt, Tm, 100*tt/Tm))
      if ((tt %% (5L*progress_every)) == 0L) gc(FALSE)
    }
  }
  
  sim_raw        <- score_raw    / Tm
  sim_depthW     <- score_depthW / Tm
  path_lcp_mean  <- ifelse(lcp_n > 0L, lcp_sum    / lcp_n, NA_real_)
  path_lcpG_mean <- ifelse(lcp_n > 0L, lcp_g_sum  / lcp_n, NA_real_)
  path_lcpGC_mean<- ifelse(lcp_n > 0L, lcp_gc_sum / lcp_n, NA_real_)
  
  p_emp <- rep(NA_real_, n)
  if (!is.null(extreme_rows_eval) && length(extreme_rows_eval)) {
    non_ext <- setdiff(seq_len(n), as.integer(unique(extreme_rows_eval)))
    if (length(non_ext)) {
      vals <- sim_raw[non_ext]
      ranks_desc <- rank(-vals, ties.method = "max")
      p_emp[non_ext] <- (ranks_desc + 1) / (length(non_ext) + 1)
    }
  }
  
  data.table::data.table(
    row                         = seq_len(n),
    similarity                  = sim_raw,
    similarity_depthW           = sim_depthW,
    path_overlap_lcp_mean       = path_lcp_mean,
    lcp_gain_weighted_mean      = path_lcpG_mean,
    lcp_gainCover_weighted_mean = path_lcpGC_mean,
    p_emp                       = p_emp
  )
}

# Convenience wrapper to build a TRAIN reference directly from model + train data.
build_reference_from_training <- function(model, X_train, f_train,
                                          lower_tail = FALSE, extreme_k = 1,
                                          tdt = NULL, helpers = NULL) {
  if (is.null(tdt)) {
    tdt <- parse_xgb_tree(model)
  }
  if (is.null(helpers)) {
    helpers <- build_structure_helpers(tdt)
  }
  mu  <- mean(f_train, na.rm = TRUE)
  sdv <- sd(f_train,  na.rm = TRUE)
  thr <- if (!lower_tail) mu + extreme_k * sdv else mu - extreme_k * sdv
  extreme_rows_train <- if (!lower_tail) which(f_train >= thr) else which(f_train <= thr)
  if (!length(extreme_rows_train)) stop("[ref-train] No training extremes at chosen threshold.")
  
  leaves_train <- predict(model, if (inherits(X_train, "dgCMatrix") || is.matrix(X_train)) X_train else as.matrix(X_train),
                          predleaf = TRUE)
  if (is.null(dim(leaves_train))) leaves_train <- matrix(leaves_train, ncol = 1L)
  storage.mode(leaves_train) <- "integer"
  
  ref <- build_extreme_reference_from_leaves(leaves_train, extreme_rows_train, helpers, tdt)
  list(ref = ref, tdt = tdt, helpers = helpers, extreme_rows_train = extreme_rows_train, threshold = thr)
}

# Compute tree-level separation between extreme and background leaf distributions.
# Returns one row per tree (0-based Tree ids).
compute_tree_separation <- function(leaves,
                                    extreme_rows,
                                    smooth = 1e-9,
                                    compute_js = FALSE) {
  stopifnot(is.matrix(leaves) || inherits(leaves, "dgCMatrix"))
  n  <- nrow(leaves)
  Tm <- ncol(leaves)
  extreme_rows <- as.integer(unique(extreme_rows))
  extreme_rows <- extreme_rows[extreme_rows >= 1L & extreme_rows <= n]
  if (!length(extreme_rows)) stop("[tree-sep] extreme_rows is empty or out of range.")
  
  is_ext <- logical(n); is_ext[extreme_rows] <- TRUE
  N_ext  <- sum(is_ext); N_bg <- n - N_ext
  if (N_ext <= 0L || N_bg <= 0L) stop("[tree-sep] Need both extreme and background SNPs.")
  
  # helpers
  kl_div <- function(p, q) {
    # both p and q strictly positive and sum to 1
    sum(p * (log(p) - log(q)))
  }
  js_div <- function(p, q) {
    m <- 0.5 * (p + q)
    0.5 * kl_div(p, m) + 0.5 * kl_div(q, m)
  }
  norm1 <- function(x) { s <- sum(x); if (s == 0) x else (x / s) }
  
  out <- data.table::data.table(
    Tree       = integer(Tm),
    KL_ext_bg  = numeric(Tm),
    KL_bg_ext  = numeric(Tm),
    JS         = if (isTRUE(compute_js)) numeric(Tm) else NULL,
    n_leaves   = integer(Tm)
  )
  
  for (tt in seq_len(Tm)) {
    tr <- tt - 1L
    lt_all <- as.integer(leaves[, tt])
    # counts
    tab_all <- table(lt_all)
    tab_ext <- table(lt_all[is_ext])
    # align support
    leaf_ids <- sort(unique(as.integer(c(names(tab_all), names(tab_ext)))))
    idx_all <- match(leaf_ids, as.integer(names(tab_all))); cnt_all <- integer(length(leaf_ids)); if (length(idx_all)) cnt_all[!is.na(idx_all)] <- as.integer(tab_all[idx_all[!is.na(idx_all)]])
    idx_ext <- match(leaf_ids, as.integer(names(tab_ext))); cnt_ext <- integer(length(leaf_ids)); if (length(idx_ext)) cnt_ext[!is.na(idx_ext)] <- as.integer(tab_ext[idx_ext[!is.na(idx_ext)]])
    cnt_bg  <- pmax(0L, cnt_all - cnt_ext)
    
    # probabilities with tiny smoothing to avoid zeros
    p_ext <- norm1(cnt_ext + smooth)
    p_bg  <- norm1(cnt_bg  + smooth)
    
    kl_e_b <- kl_div(p_ext, p_bg)
    kl_b_e <- kl_div(p_bg,  p_ext)
    
    out[tt, `:=`(
      Tree      = tr,
      KL_ext_bg = kl_e_b,
      KL_bg_ext = kl_b_e,
      n_leaves  = length(leaf_ids)
    )]
    if (isTRUE(compute_js)) out$JS[tt] <- js_div(p_ext, p_bg)
  }
  out[]
}

# Select the top trees by a chosen score column (default = KL_ext_bg).
# You can control how many via top_frac / min_trees / max_trees or a quantile cutoff.
select_trees_by_separation <- function(tree_scores,
                                       score_col   = c("KL_ext_bg", "JS"),
                                       top_frac    = 0.20,     # keep top 20% by default
                                       min_trees   = 50L,      # but at least 50 trees
                                       max_trees   = NULL,     # or cap if you want
                                       q_cutoff    = NULL      # alternatively, keep >= quantile
) {
  stopifnot(is.data.frame(tree_scores))
  score_col <- match.arg(score_col)
  DT <- data.table::as.data.table(tree_scores)
  if (!(score_col %in% names(DT))) stop(sprintf("[select-trees] score_col '%s' not found.", score_col))
  data.table::setorder(DT, -get(score_col))
  
  Tm <- nrow(DT)
  if (!is.null(q_cutoff)) {
    thr <- stats::quantile(DT[[score_col]], probs = q_cutoff, na.rm = TRUE)
    sel <- DT[get(score_col) >= thr, Tree]
  } else {
    k <- max(min_trees, ceiling(top_frac * Tm))
    if (!is.null(max_trees)) k <- min(k, as.integer(max_trees))
    sel <- DT$Tree[seq_len(min(k, Tm))]
  }
  sort(unique(as.integer(sel)))
}


run_snp_level_pipeline <- function(model, # trained XGBoost model
                                   X,     # training feature matrix
                                   f_vec, # training lable vector (F)
                                   save_prefix = "leafsim",
                                   snpcodes    = NULL, # optional; dataframe with columns "snpcode" and "row"; used to annotate results
                                   lower_tail  = FALSE,
                                   extreme_k   = 1,
                                   target_bins = 8L,
                                   min_bin_size = 50L,
                                   winsor_prob  = 0.001,
                                   chunk_size  = 200L,
                                   feat_top_k  = 2L,
                                   progress_every = 500L,
                                   save_checkpoint = TRUE
) {
  # --- Path similarity metrics per‑SNP; keeps leaves/tdt/helpers) ---
  res <- compute_snp_similarity(
    model = model, 
    X = train.dat, 
    f_vec = f.train,
    lower_tail     = lower_tail,
    extreme_k      = extreme_k,
    progress_every = progress_every,
    show_progress  = TRUE
  )
  
  # --- Precompute per‑leaf steps ---
  .info("Precomputing per‑leaf steps...")
  res$leaf_steps <- build_leaf_steps(
    res$leaves,
    res$tdt, 
    res$helpers,
    binning         = "adaptive",
    target_bins     = target_bins,
    min_per_bin     = min_bin_size,
    winsor_prob     = winsor_prob,
    method          = "fd",      # or "quantile"
    trees_per_batch = chunk_size
  )
  
  # --- Feature concentration (per‑leaf → per‑SNP), from leaf_steps) ---
  .info("Computing feature concentration (top-%d) per SNP...", feat_top_k)
  res$perleaf <- per_leaf_topk_share(
    res$leaf_steps, 
    top_k = feat_top_k,
    progress_every = progress_every
  )
  
  fc_per_snp <- feature_concentration_per_snp(
    res$leaves, 
    res$perleaf, 
    chunk_rows = chunk_size
  )
  
  DTs <- list(res$summaries, fc_per_snp, train.codes)
  summaries <- Reduce(function(x, y) merge(x, y, by = "row", all.x = TRUE, sort = FALSE), DTs)
  
  res$summaries <- summaries
  
  # --- Save per‑SNP metrics to disk (small) ---
  out_metrics <- paste0(save_prefix, "_perSNP_metrics_", tail_dir, ".csv")
  fwrite(res$summaries, out_metrics)
  .info("Wrote per‑SNP metrics: %s (rows=%d)", out_metrics, nrow(summaries))
  
  # --- Cache checkpoint  ---
  checkpoint_file <- NA_character_
  if (isTRUE(save_checkpoint)) {
    chk_file <- timestamped_file(paste0(save_prefix, "_persnp_checkpoint_", tail_dir), ".rds")
    .info("Saving checkpoint: %s", chk_file)
    saveRDS(list(
      summaries = res$summaries,
      leaves    = res$leaves,
      tdt       = res$tdt,
      helpers   = res$helpers,
      leaf_steps= res$leaf_steps,
      extreme   = res$extreme,
      tail_dir  = res$tail_dir), file = chk_file)
    .info("Checkpoint saved: %s", chk_file)
  }
  
  # --- Return ---
  list(
    summaries     = res$summaries,        # per‑SNP table (incl. LCP + feature_conc)
    leaves        = res$leaves,           # n x Tm integer matrix
    tdt           = res$tdt,              # compact model tree table
    helpers       = res$helpers,          # cached structure helpers
    leaf_steps    = res$leaf_steps,       # (Tree, leaf_id) → steps (feature,dir,thresh_bin,depth)
    perleaf       = res$perleaf,          # per‑leaf top‑k share (used for FC)
    extreme_idx   = res$extreme_idx,
    tail_dir      = res$tail_dir
  )
}
