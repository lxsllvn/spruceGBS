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
# Path similarity (per-SNP)
#============================#
# Returns per-SNP a log-ratio based similarity metric, mean leaf depth, mean leaf cover, and mean leaf value, plus reusable leaves/tdt/helpers.
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
  
  ## Trees + helpers -----------------------------------------------------------
  tdt <- parse_xgb_tree(model)
  if ((max(tdt$Tree) + 1L) != Tm) {
    .warn_bad("Tree count mismatch: dump=%d vs leaves=%d.", max(tdt$Tree)+1L, Tm)
  }
  .info("Building structure helpers...")
  helpers <- build_structure_helpers(tdt)
  
  ## Accumulators --------------------------------------------------------------
  loglr_sum      <- numeric(n)  # avg(log(prop_ext / p_leaf))
  leafdepth_sum  <- numeric(n)  # avg leaf depth
  leafcover_sum  <- numeric(n)  # avg leaf cover
  leafval_sum    <- numeric(n)  # avg leaf value
  
  m <- length(extreme_idx)
  
  ## Per-tree loop -------------------------------------------------------------
  for (tt in seq_len(Tm)) {
    tr <- tt - 1L
    lt <- leaves[, tt]
    
    # counts of extreme per leaf in this tree
    ft_tab <- table(lt[extreme_idx])
    
    # extreme SNP frequency per leaf (vectorized to rows)
    leaf_ids_ext <- as.integer(names(ft_tab))
    idx_ext <- match(lt, leaf_ids_ext)
    cnt_ext <- integer(n); hit_ext <- !is.na(idx_ext)
    if (any(hit_ext)) cnt_ext[hit_ext] <- as.integer(ft_tab[idx_ext[hit_ext]])
    prop_ext <- cnt_ext / m
    
    # total SNP frequency per leaf
    tab_all <- table(lt)
    leaf_ids_all <- as.integer(names(tab_all))
    idx_all <- match(lt, leaf_ids_all)
    cnt_all <- integer(n); hit_all <- !is.na(idx_all)
    if (any(hit_all)) cnt_all[hit_all] <- as.integer(tab_all[idx_all[hit_all]])
    p_leaf  <- cnt_all / n
    
    # log-ratio similarity; treat zeros as 0 contribution (instead of -Inf)
    llr <- ifelse(prop_ext > 0 & p_leaf > 0, log(prop_ext / p_leaf), 0)
    loglr_sum <- loglr_sum + llr
    
    # depth map (index by node id = leaf id)
    dep_vals <- helpers$get_depth_map(tr)
    if (!is.null(dep_vals)) {
      dvec <- dep_vals[lt + 1L]
      dvec[!is.finite(dvec)] <- 0
      leafdepth_sum <- leafdepth_sum + dvec
    }
    
    # leaf cover & value
    leaf_rows <- tdt[Tree == tr & Leaf == TRUE, .(ID, LeafVal, LeafCover = Cover)]
    if (nrow(leaf_rows)) {
      data.table::setkey(leaf_rows, ID)
      LR <- leaf_rows[.(as.integer(lt))]
      # coalesce NA to 0 before summing
      cov <- suppressWarnings(as.numeric(LR$LeafCover)); cov[!is.finite(cov)] <- 0
      val <- suppressWarnings(as.numeric(LR$LeafVal));   val[!is.finite(val)] <- 0
      leafcover_sum <- leafcover_sum + cov
      leafval_sum   <- leafval_sum   + val
    }
    
    if (show_progress && (tt %% progress_every == 0L || tt == Tm)) {
      .info("processed %d/%d trees (%.1f%%)", tt, Tm, 100*tt/Tm)
      if ((tt %% 50L) == 0L) gc(FALSE)
    }
  }
  
  ## Final per-SNP metrics -----------------------------------------------------
  summaries <- data.table::data.table(
    row             = seq_len(n),
    F               = as.numeric(f_vec),
    sim_logLR       = loglr_sum     / Tm,
    mean_leaf_depth = leafdepth_sum / Tm,
    mean_leaf_cover = leafcover_sum / Tm,
    mean_leaf_value = leafval_sum   / Tm
  )
  
  list(
    summaries   = summaries,
    leaves      = leaves,
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

