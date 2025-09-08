#===============================================================================
# XGBoost Rule Harvesting & Validation 
#===============================================================================


#============================
# Shared validators and helpers
#============================
.require_numeric_vec <- function(x, name) {
  if (!is.numeric(x)) stop(sprintf("[%s] must be numeric.", name))
}

.require_finite_center_scale <- function(y, name) {
  m1 <- mean(y, na.rm = TRUE)
  m2 <- stats::median(y, na.rm = TRUE)
  sdv <- stats::sd(y, na.rm = TRUE)
  md  <- stats::mad(y, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(m1)) stop(sprintf("[%s] mean is NA/Inf.", name))
  if (!is.finite(m2)) stop(sprintf("[%s] median is NA/Inf.", name))
  if (!is.finite(sdv) || sdv <= 0) stop(sprintf("[%s] sd <= 0 or NA/Inf.", name))
  if (!is.finite(md)  || md  <= 0) stop(sprintf("[%s] mad <= 0 or NA/Inf.", name))
}

#============================
# Small shared helpers
#============================
# Derive extremes + background indices from yvar, based on explicit input or statistical banding around the center.
derive_yvar_idx <- function(
    y, n,
    lower_tail  = FALSE,
    extreme_k   = 1,
    extreme_idx = NULL,
    bg_idx      = NULL,
    bg_band_k   = NULL,
    band_center = c("mean","median"),
    band_scale  = c("sd","mad"),
    check_y_stats = NULL
) {
  band_center <- match.arg(band_center)
  band_scale  <- match.arg(band_scale)
  
  # Decide if we need y-stat checks:
  # default: TRUE unless BOTH idx sets are provided and no banding is requested
  if (is.null(check_y_stats)) {
    check_y_stats <- (is.null(extreme_idx) || is.null(bg_idx) || !is.null(bg_band_k))
  }
  
  # Only validate y if needed
  if (isTRUE(check_y_stats)) {
    if (!is.numeric(y)) stop("[derive_yvar_idx] y must be numeric.")
    if (length(y) != n) stop("[derive_yvar_idx] length(y) must equal n.")
    mu  <- mean(y, na.rm = TRUE);  if (!is.finite(mu))  stop("[derive_yvar_idx] mean(y) not finite.")
    md  <- stats::median(y, na.rm = TRUE); if (!is.finite(md))  stop("[derive_yvar_idx] median(y) not finite.")
    sdv <- stats::sd(y, na.rm = TRUE);     if (!is.finite(sdv) || sdv <= 0) stop("[derive_yvar_idx] sd(y) <= 0 or non-finite.")
    madv<- stats::mad(y, constant = 1.4826, na.rm = TRUE); if (!is.finite(madv) || madv <= 0) stop("[derive_yvar_idx] mad(y) <= 0 or non-finite.")
  }
  
  # Build extremes
  if (is.null(extreme_idx)) {
    mu  <- mean(y, na.rm = TRUE); sdv <- sd(y, na.rm = TRUE)
    thr <- if (!lower_tail) mu + extreme_k * sdv else mu - extreme_k * sdv
    extr_idx <- if (!lower_tail) which(y >= thr) else which(y <= thr)
  } else {
    extr_idx <- sort(unique(as.integer(extreme_idx)))
    extr_idx <- extr_idx[extr_idx >= 1L & extr_idx <= n]
  }
  
  # Build background
  if (!is.null(bg_idx)) {
    bg_idx <- sort(unique(as.integer(bg_idx)))
    bg_idx <- bg_idx[bg_idx >= 1L & bg_idx <= n]
    if (length(intersect(extr_idx, bg_idx))) bg_idx <- setdiff(bg_idx, extr_idx)
  } else if (!is.null(bg_band_k)) {
    ctr <- if (band_center == "mean") mean(y, na.rm = TRUE) else stats::median(y, na.rm = TRUE)
    scl <- if (band_scale  == "sd")   stats::sd(y, na.rm = TRUE) else stats::mad(y, constant = 1.4826, na.rm = TRUE)
    lo  <- ctr - bg_band_k * scl; hi <- ctr + bg_band_k * scl
    bg_idx <- which(y >= lo & y <= hi)
    if (length(extr_idx)) bg_idx <- setdiff(bg_idx, extr_idx)
  } else {
    bg_idx <- setdiff(seq_len(n), extr_idx)
  }
  
  # Final checks
  N_extr <- length(extr_idx); N_bg <- length(bg_idx)
  if (N_extr <= 0L) stop("[derive_yvar_idx] No rows in extreme set.")
  if (N_bg   <= 0L) stop("[derive_yvar_idx] No rows in background set.")
  if (length(intersect(extr_idx, bg_idx)) > 0L) stop("[derive_yvar_idx] extreme and background indices overlap.")
  
  is_extreme <- logical(n); is_extreme[extr_idx] <- TRUE
  is_bg      <- logical(n); is_bg[bg_idx]       <- TRUE
  
  list(
    extr_idx  = extr_idx,
    bg_idx    = bg_idx,
    N_extr    = N_extr,
    N_bg      = N_bg,
    base_rate = N_extr / (N_extr + N_bg),
    is_extreme = is_extreme,
    is_bg      = is_bg
  )
}

# Precomputes the number of leaf bins (`nb`) for each tree column in a leaves matrix.
precompute_nb_vec <- function(leaves) {
  Tm <- ncol(leaves)
  # dgCMatrix fast-path via Matrix::colMaxs (not matrixStats)
  if (inherits(leaves, "dgCMatrix")) {
    if (requireNamespace("Matrix", quietly = TRUE)) {
      mx <- Matrix::colMaxs(leaves, na.rm = TRUE)
      mx[!is.finite(mx)] <- -1
      return(as.integer(mx + 1L))
    }
  }
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    mx <- matrixStats::colMaxs(leaves, na.rm = TRUE)
    mx[!is.finite(mx)] <- -1
    return(as.integer(mx + 1L))
  }
  # Fallback loop (no suppressWarnings; handle NA explicitly)
  res <- integer(Tm)
  for (tt in seq_len(Tm)) {
    col <- leaves[, tt]
    mx  <- if (anyNA(col)) max(col, na.rm = TRUE) else max(col)
    if (!is.finite(mx)) res[tt] <- 0L else res[tt] <- as.integer(mx) + 1L
  }
  res
}

# Tabulates counts of a given index set into leaf bins of a single tree.
count_by_leaf <- function(lt, idx, nb) {
  if (!length(idx) || nb <= 0L) return(integer(max(1L, nb)))
  tabulate(lt[idx] + 1L, nbins = nb)
}

# Learn a fixed extreme/background rule on TRAIN and reapply it to any y.
# Canonical names, disjoint sets, and selectable center/scale for BG band.
learn_sigma_rule <- function(yvar_train,
                             extreme_k    = 1,
                             lower_tail   = FALSE,
                             bg_band_k    = NULL,                 # NULL → background = !extreme
                             band_center  = c("mean","median"),   # BG center learned on TRAIN
                             band_scale   = c("sd","mad")) {      # BG scale learned on TRAIN
  band_center <- match.arg(band_center)
  band_scale  <- match.arg(band_scale)
  
  # --- Validate TRAIN y -------------------------------------------------------
  if (!is.numeric(yvar_train)) stop("[learn_sigma_rule] yvar_train must be numeric.")
  if (!length(yvar_train))     stop("[learn_sigma_rule] yvar_train is empty.")
  
  center_val <- if (band_center == "mean") {
    m <- mean(yvar_train, na.rm = TRUE)
    if (!is.finite(m)) stop("[learn_sigma_rule] mean(yvar_train) is not finite.")
    m
  } else {
    med <- stats::median(yvar_train, na.rm = TRUE)
    if (!is.finite(med)) stop("[learn_sigma_rule] median(yvar_train) is not finite.")
    med
  }
  
  # Scale for extreme threshold is always SD (to mirror common sigma rules)
  sd_train <- stats::sd(yvar_train, na.rm = TRUE)
  if (!is.finite(sd_train) || sd_train <= 0)
    stop("[learn_sigma_rule] sd(yvar_train) <= 0 or not finite.")
  
  # Scale for BG band can be SD or MAD (robust), learned on TRAIN
  scale_val <- if (band_scale == "sd") {
    s <- stats::sd(yvar_train, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) stop("[learn_sigma_rule] sd(yvar_train) for BG is not finite/positive.")
    s
  } else {
    md <- stats::mad(yvar_train, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(md) || md <= 0) stop("[learn_sigma_rule] mad(yvar_train) for BG is not finite/positive.")
    md
  }
  
  # --- Fixed thresholds learned on TRAIN --------------------------------------
  thr_ext <- if (!lower_tail) center_val + extreme_k * sd_train
  else             center_val - extreme_k * sd_train
  
  thr_bg  <- if (is.null(bg_band_k)) NULL else {
    lo <- center_val - bg_band_k * scale_val
    hi <- center_val + bg_band_k * scale_val
    c(lo, hi)
  }
  
  # --- Return an object with an apply() that enforces the SAME rule -----------
  list(
    # Apply learned thresholds to any numeric vector f (train/test/heldout)
    apply = function(f) {
      if (!is.numeric(f)) stop("[learn_sigma_rule/apply] input vector must be numeric.")
      n <- length(f); if (!n) stop("[learn_sigma_rule/apply] input vector is empty.")
      
      # extremes by the fixed TRAIN threshold
      extr_idx <- if (!lower_tail) which(f >= thr_ext) else which(f <= thr_ext)
      
      # background either as band (TRAIN-learned) or !extreme
      if (is.null(thr_bg)) {
        bg_idx <- setdiff(seq_len(n), extr_idx)
      } else {
        bg_idx <- which(f >= thr_bg[1] & f <= thr_bg[2])
        # enforce disjointness explicitly
        if (length(extr_idx)) bg_idx <- setdiff(bg_idx, extr_idx)
      }
      
      # Canonical outputs (sorted, unique)
      extr_idx <- sort(unique(as.integer(extr_idx)))
      bg_idx   <- sort(unique(as.integer(bg_idx)))
      
      N_extr <- length(extr_idx)
      N_bg   <- length(bg_idx)
      
      if (N_extr <= 0L) stop("[learn_sigma_rule/apply] No rows in extreme set under learned rule.")
      if (N_bg   <= 0L) stop("[learn_sigma_rule/apply] No rows in background set under learned rule.")
      if (length(intersect(extr_idx, bg_idx)) > 0L)
        stop("[learn_sigma_rule/apply] extreme and background indices overlap (should be disjoint).")
      
      list(
        extr_idx = extr_idx,
        bg_idx   = bg_idx,
        N_extr   = N_extr,
        N_bg     = N_bg
      )
    },
    
    # Expose learned parameters for reproducibility/auditing
    thr_ext      = thr_ext,                 # scalar extreme threshold (TRAIN)
    thr_bg       = thr_bg,                  # length-2 band [lo, hi] or NULL
    center_name  = band_center,
    scale_name   = band_scale,
    center_val   = center_val,              # learned on TRAIN
    sd_train     = sd_train,                # used for extreme threshold
    scale_val    = scale_val,               # used for BG band
    extreme_k    = extreme_k,
    bg_band_k    = bg_band_k,
    lower_tail   = lower_tail
  )
}

#============================
# Parse model → compact tree table (model → tdt)
#============================
# Converts an `xgb.Booster` into a tidy, numeric, per-node `data.table` with consistent schema (Tree, ID, children, split, gain, leaf values, etc.).
parse_xgb_tree <- function(model) {
  stopifnot(inherits(model, "xgb.Booster"))
  dt <- xgboost::xgb.model.dt.tree(model = model)
  data.table::setDT(dt)
  
  to_int_child <- function(x) {
    if (is.null(x)) return(NA_integer_)
    y <- suppressWarnings(as.integer(sub(".*-", "", x)))
    y[is.na(x)] <- NA_integer_
    y
  }
  to_num <- function(x) suppressWarnings(as.numeric(x))
  
  if (!"Node" %in% names(dt)) stop("[parse_xgb_tree] column 'Node' not found in xgb dump")
  
  dt[, `:=`(Tree = as.integer(Tree), ID = as.integer(Node))]
  
  for (col in c("Yes","No","Missing")) {
    if (col %in% names(dt)) dt[, (col) := to_int_child(get(col))] else dt[, (col) := NA_integer_]
  }
  
  if (!"Split" %in% names(dt)) dt[, Split := NA_real_]
  if (!"Cover" %in% names(dt)) dt[, Cover := NA_real_]
  dt[, `:=`(Split = to_num(Split), Cover = to_num(Cover))]
  
  dt[, Leaf := (Feature == "Leaf")]
  
  has_quality <- "Quality" %in% names(dt)
  has_gaincol <- "Gain" %in% names(dt)
  Q <- if (has_quality) to_num(dt$Quality) else rep(NA_real_, nrow(dt))
  G <- if (has_gaincol) to_num(dt$Gain)    else rep(NA_real_, nrow(dt))
  
  dt[, LeafVal := ifelse(Leaf, Q, NA_real_)]
  dt[Leaf & is.na(LeafVal), LeafVal := Split]
  dt[, Gain := ifelse(!Leaf, ifelse(!is.na(G), G, Q), NA_real_)]
  
  drop_cols <- intersect(c("Node","Quality"), names(dt))
  if (length(drop_cols)) dt[, (drop_cols) := NULL]
  
  data.table::setkey(dt, Tree, ID)
  data.table::setorder(dt, Tree, ID)
  dt[, .(Tree, ID, Feature, Split, Yes, No, Missing, Leaf, LeafVal, Gain, Cover)]
}


#============================#
# Structure helpers
#============================#
# Constructs parent/child maps and depth/path lookup functions for tree navigation.
build_structure_helpers <- function(tdt) {
  stopifnot(is.data.frame(tdt),
            all(c("Tree","ID","Yes","No","Missing","Leaf","Gain","Cover") %in% names(tdt)))
  tdt <- data.table::as.data.table(tdt)
  trees <- sort(unique(tdt$Tree)) #
  
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
  data.table::setkey(split_gc, Tree, ID)
  
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
    cur <- leaf_id; k <- 0L
    while (!is.na(cur) && cur != 0L) { par <- pmap[cur + 1L]; if (is.na(par)) break; k <- k + 1L; cur <- par }
    if (k == 0L) {
      out <- list(nodes=integer(0), gain=numeric(0), gcover=numeric(0))
      leaf_cache[[tkey]][[as.character(leaf_id)]] <- out; return(out)
    }
    nodes <- integer(k); cur <- leaf_id
    for (pos in k:1) { par <- pmap[cur + 1L]; nodes[pos] <- par; cur <- par; if (is.na(cur) || cur == 0L) break }
    sgc <- split_gc[.(tr, nodes), .(Gain, SplitCover)]
    out <- list(nodes = nodes, gain = as.numeric(sgc$Gain), gcover = as.numeric(sgc$SplitCover))
    leaf_cache[[tkey]][[as.character(leaf_id)]] <- out
    out
  }
  
  list(get_depth_map = get_depth_map,
       get_max_depth = get_max_depth,
       get_leaf_path = get_leaf_path)
}

#============================
# Build per-leaf steps
#============================
# Generates per-leaf decision paths (feature, direction, threshold bins) from tree structure and leaves
build_leaf_steps <- function(leaves,
                             tdt,
                             helpers,
                             target_bins     = 10L,
                             min_per_bin     = 50L,
                             winsor_prob     = 0.01,
                             method          = c("fd","quantile"),
                             trees_per_batch = 250L,
                             progress_every  = 1000L) {
  message(sprintf("[build_leaf_steps] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  stopifnot(is.matrix(leaves) || inherits(leaves, "dgCMatrix"))
  method <- match.arg(method)
  Tm <- ncol(leaves)
  
  tdt <- data.table::as.data.table(tdt)
  data.table::setkey(tdt, Tree, ID)
  
  build_raw_steps_batch <- function(batch_start, batch_end) {
    total_steps <- 0L
    meta <- vector("list", batch_end - batch_start + 1L); mi <- 1L
    for (tt in batch_start:batch_end) {
      tr <- tt - 1L
      lt <- leaves[, tt]
      uleaf <- sort(unique(as.integer(lt)))
      meta[[mi]] <- list(tt = tt, tr = tr, uleaf = uleaf)
      if (length(uleaf)) {
        for (leaf_id in uleaf) {
          pn <- helpers$get_leaf_path(tr, leaf_id)$nodes
          if (length(pn)) total_steps <- total_steps + length(pn)
        }
      }
      mi <- mi + 1L
    }
    if (total_steps == 0L) {
      return(data.table::data.table(Tree = integer(0), leaf_id = integer(0),
                                    depth = integer(0), feature = character(0),
                                    direction = character(0), split_val = numeric(0)))
    }
    
    Tree_v  <- integer(total_steps)
    Leaf_v  <- integer(total_steps)
    Depth_v <- integer(total_steps)
    Feat_v  <- character(total_steps)
    Dir_v   <- character(total_steps)
    Split_v <- numeric(total_steps)
    cursor  <- 0L
    
    for (m in meta) {
      tt <- m$tt; tr <- m$tr; uleaf <- m$uleaf
      if (!length(uleaf)) next
      tt_dt <- tdt[.(tr)]
      if (!nrow(tt_dt)) next
      
      maxid <- max(tt_dt$ID, na.rm = TRUE)
      yesA  <- rep.int(NA_integer_,  maxid + 1L);  yesA [tt_dt$ID + 1L] <- tt_dt$Yes
      noA   <- rep.int(NA_integer_,  maxid + 1L);  noA  [tt_dt$ID + 1L] <- tt_dt$No
      featA <- rep.int(NA_character_,maxid + 1L);  featA[tt_dt$ID + 1L] <- as.character(tt_dt$Feature)
      spltA <- rep.int(NA_real_,     maxid + 1L);  spltA[tt_dt$ID + 1L] <- suppressWarnings(as.numeric(tt_dt$Split))
      
      for (leaf_id in uleaf) {
        pn <- helpers$get_leaf_path(tr, leaf_id)$nodes
        k  <- length(pn)
        if (!k) next
        child_along <- integer(k)
        if (k > 1L) child_along[1:(k-1L)] <- pn[2:k]
        child_along[k] <- leaf_id
        parents <- pn + 1L
        yes  <- yesA [parents]
        no   <- noA  [parents]
        feat <- featA[parents]
        splt <- spltA[parents]
        dir <- ifelse(yes == child_along, "<", ifelse(no == child_along, ">=", "missing"))
        
        rng <- (cursor + 1L):(cursor + k)
        Tree_v [rng] <- tr
        Leaf_v [rng] <- leaf_id
        Depth_v[rng] <- 0L:(k - 1L)
        Feat_v [rng] <- feat
        Dir_v  [rng] <- dir
        Split_v[rng] <- splt
        cursor <- cursor + k
      }
    }
    
    data.table::data.table(
      Tree      = Tree_v,
      leaf_id   = Leaf_v,
      depth     = Depth_v,
      feature   = Feat_v,
      direction = Dir_v,
      split_val = Split_v
    )
  }
  
  out_list <- vector("list", ceiling(Tm / trees_per_batch)); oi <- 1L
  for (batch_start in seq(1L, Tm, by = trees_per_batch)) {
    batch_end <- min(batch_start + trees_per_batch - 1L, Tm)
    out_list[[oi]] <- build_raw_steps_batch(batch_start, batch_end); oi <- oi + 1L
    if (!is.null(progress_every) && progress_every > 0L &&
        ((batch_end %% progress_every) == 0L || batch_end == Tm)) {
      message(sprintf("[build_leaf_steps] built steps up to tree %d / %d (%.1f%%)",
                      batch_end, Tm, 100*batch_end/Tm))
    }
  }
  
  LS <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (!nrow(LS)) {
    data.table::setkey(LS, Tree, leaf_id)
    return(LS[, .(Tree, leaf_id, depth, feature, direction, thresh_bin = numeric())])
  }
  
  LS <- LS[!is.na(feature) & nzchar(feature)]
  LS[, `:=`(
    Tree      = as.integer(Tree),
    leaf_id   = as.integer(leaf_id),
    depth     = as.integer(depth),
    feature   = as.character(feature),
    direction = as.character(direction),
    split_val = as.numeric(split_val)
  )]
  
  message(sprintf("[build_leaf_steps] adaptive bins per feature: target=%d, min_per_bin=%d, method=%s, winsor=%.3g",
        as.integer(target_bins), as.integer(min_per_bin), method, as.numeric(winsor_prob)))
  
  LS_split <- LS[is.finite(split_val)]
  feats <- sort(unique(LS_split$feature))
  breaks_env <- new.env(parent = emptyenv()); mids_env <- new.env(parent = emptyenv())
  
  propose_breaks <- function(v, max_bins, method, winsor) {
    v <- sort(v[is.finite(v)]); if (length(v) <= 1L) { r <- range(v, na.rm=TRUE); if(!all(is.finite(r))) r <- c(0,0); return(unique(c(r[1], r[2]))) }
    lo <- suppressWarnings(stats::quantile(v, probs = winsor, names = FALSE))
    hi <- suppressWarnings(stats::quantile(v, probs = 1 - winsor, names = FALSE))
    v  <- v[v >= lo & v <= hi]; if (length(v) <= 1L) { r <- range(v, na.rm=TRUE); if(!all(is.finite(r))) r <- c(0,0); return(unique(c(r[1], r[2]))) }
    if (method == "fd") {
      h <- 2 * stats::IQR(v) / (length(v)^(1/3)); if (!is.finite(h) || h <= 0) h <- (max(v) - min(v)) / max(2L, max_bins)
      nb <- ceiling((max(v) - min(v)) / max(h, .Machine$double.eps)); nb <- min(max(nb, 2L), max_bins)
      br <- pretty(v, n = nb)
    } else {
      nb <- max(2L, as.integer(max_bins))
      probs <- seq(0, 1, length.out = nb + 1L)
      br <- unique(stats::quantile(v, probs = probs, names = FALSE, type = 7))
      if (length(br) < 3L) br <- pretty(v, n = min(max_bins, 3L))
    }
    br <- sort(unique(as.numeric(br))); if (length(br) < 2L) br <- unique(c(min(v), max(v))); br
  }
  enforce_min_bin <- function(v, br, target_bins, min_per_bin) {
    if (length(br) < 2L) return(br)
    repeat {
      idx <- findInterval(v, br, all.inside = TRUE)
      tab <- tabulate(idx, nbins = length(br) - 1L)
      if ((length(tab) <= target_bins) && all(tab >= min_per_bin | tab == 0L)) break
      k <- if (length(tab)) which.min(tab) else 1L
      if (length(tab) <= 1L) break
      rm_pos <- if (k == 1L) 2L else if (k == length(tab)) length(tab) else if (tab[k - 1L] <= tab[k + 1L]) k else k + 1L
      br <- br[-rm_pos]; if (length(br) < 2L) { br <- br[1:2]; break }
    }
    br
  }
  
  for (f in feats) {
    v  <- LS_split[feature == f, split_val]; if (!length(v)) next
    br <- propose_breaks(v, max_bins = as.integer(target_bins), method = method, winsor = as.numeric(winsor_prob))
    br <- enforce_min_bin(v, br, target_bins = as.integer(target_bins), min_per_bin = as.integer(min_per_bin))
    if (length(br) < 2L) { r <- range(v, na.rm=TRUE); if(!all(is.finite(r))) r <- c(0,0); br <- unique(c(r[1], r[2])) }
    mids <- (br[-1L] + br[-length(br)]) / 2
    breaks_env[[f]] <- br; mids_env[[f]] <- mids
  }
  
  LS[, thresh_bin := {
    br <- breaks_env[[feature[1L]]]; md <- mids_env[[feature[1L]]]
    if (is.null(br) || is.null(md) || !is.finite(split_val[1L])) rep(NA_real_, .N) else {
      idx <- findInterval(split_val, br, all.inside = TRUE); as.numeric(md[idx])
    }
  }, by = feature]
  
  LS[, split_val := NULL]
  data.table::setkey(LS, Tree, leaf_id)
  LS[]
}

#============================
# Build (Tree,leaf) rule strings
#============================
.format_num <- function(x) formatC(as.numeric(x), format = "e", digits = 2)

# Converts per-leaf decision paths into canonical rule strings, optionally tightening redundant splits.
build_path_rule_strings <- function(leaf_steps, max_depth = 4L, tighten_monotone = TRUE) {
  LS <- data.table::as.data.table(leaf_steps)
  need <- c("Tree","leaf_id","depth","feature","direction","thresh_bin")
  miss <- setdiff(need, names(LS))
  if (length(miss)) stop(sprintf("[paths] leaf_steps missing: %s", paste(miss, collapse=", ")))
  LS[, `:=`(
    Tree       = as.integer(Tree),
    leaf_id    = as.integer(leaf_id),
    depth      = as.integer(depth),
    feature    = as.character(feature),
    direction  = as.character(direction),
    thresh_bin = as.numeric(thresh_bin)
  )]
  LS <- LS[depth < as.integer(max_depth)]
  data.table::setorder(LS, Tree, leaf_id, depth)
  
  if (isTRUE(tighten_monotone)) {
    eps <- 1e-12
    LS_tight <- LS[, {
      out <- list(cond = character(0), depth = integer(0))
      best_ge <- new.env(parent = emptyenv()); best_lt <- new.env(parent = emptyenv())
      for (i in seq_len(.N)) {
        f <- feature[i]; d <- direction[i]; b <- thresh_bin[i]; keep <- TRUE
        if (d == ">=") { cur <- best_ge[[f]]; if (is.null(cur) || b > cur + eps) best_ge[[f]] <- b else keep <- FALSE }
        else if (d == "<") { cur <- best_lt[[f]]; if (is.null(cur) || b < cur - eps) best_lt[[f]] <- b else keep <- FALSE }
        if (keep) { out$cond <- c(out$cond, paste0(f, " ", d, " ", .format_num(b))); out$depth <- c(out$depth, depth[i]) }
      }
      if (length(out$cond)) .(cond = out$cond, depth = out$depth) else NULL
    }, by = .(Tree, leaf_id)]
  } else {
    LS_tight <- LS[, .(cond = paste0(feature, " ", direction, " ", .format_num(thresh_bin)),
                       depth = depth),
                   by = .(Tree, leaf_id)]
  }
  
  PATHS <- LS_tight[, .(
    rule_len  = .N,
    rule_str  = paste(cond, collapse = " | "),
    depth_max = max(depth)
  ), by = .(Tree, leaf_id)]
  data.table::setkey(PATHS, Tree, leaf_id)
  PATHS[]
}

#============================
# Path similarity (per-SNP)
#============================
# Computes per-row similarity scores (log-likelihood ratios, mean depth, cover, leaf value) between extreme and background sets.
compute_snp_similarity <- function(model, X, yvar_train,
                                   lower_tail     = FALSE,
                                   extreme_k      = 1,
                                   extreme_idx    = NULL,
                                   bg_idx         = NULL,
                                   bg_band_k      = NULL,
                                   band_center    = c("mean","median"),
                                   band_scale     = c("sd","mad"),
                                   progress_every = NULL,
                                   leaves_override    = NULL,
                                   compute_depth      = TRUE,
                                   compute_leaf_stats = TRUE) {
  message(sprintf("[compute_snp_similarity] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  stopifnot(inherits(model, "xgb.Booster"))
  X <- if (inherits(X, "dgCMatrix") || is.matrix(X)) X else as.matrix(X)
  .require_numeric_vec(yvar_train, "compute_snp_similarity:yvar_train")
  if (length(yvar_train) != nrow(X)) stop("[compute_snp_similarity] length(yvar_train) != nrow(X).")
  .require_finite_center_scale(yvar_train, "compute_snp_similarity:yvar_train")
  band_center <- match.arg(band_center); band_scale <- match.arg(band_scale)
  
  n <- nrow(X)
  idx <- derive_yvar_idx(
    y             = yvar_train,
    n             = n,
    lower_tail    = lower_tail,         
    extreme_k     = extreme_k,             
    extreme_idx   = extreme_idx,
    bg_idx        = bg_idx,          
    bg_band_k     = bg_band_k,
    band_center   = band_center,
    band_scale    = band_scale
  )
  
  extr_idx   <- idx$extr_idx
  bg_idx     <- idx$bg_idx
  N_extr     <- idx$N_extr
  N_bg       <- idx$N_bg

  inv_N_extr <- if (N_extr > 0L) 1 / N_extr else 0
  inv_N_bg   <- if (N_bg > 0L) 1 / N_bg else 0
  eps <- 1e-12
  
  leaves <- if (is.null(leaves_override)) predict(model, X, predleaf = TRUE) else leaves_override
  if (is.null(dim(leaves))) leaves <- matrix(leaves, ncol = 1L)
  if (nrow(leaves) != n) stop("[compute_snp_similarity] nrow(leaves) must equal nrow(X).")
  storage.mode(leaves) <- "integer"
  Tm <- ncol(leaves)
  tdt <- parse_xgb_tree(model)
  if ((max(tdt$Tree) + 1L) != Tm) stop(sprintf("[compute_snp_similarity] Tree count mismatch: dump=%d vs leaves=%d.", max(tdt$Tree)+1L, Tm))
  
  helpers <- build_structure_helpers(tdt)
  nb_vec <- precompute_nb_vec(leaves)
  
  depth_map <- if (isTRUE(compute_depth)) {
    out <- vector("list", Tm); for (tt in seq_len(Tm)) out[[tt]] <- helpers$get_depth_map(tt - 1L); out
  } else NULL
  
  leaf_cover <- leaf_value <- NULL
  if (isTRUE(compute_leaf_stats)) {
    leaf_cover <- vector("list", Tm); leaf_value <- vector("list", Tm)
    for (tr in sort(unique(tdt$Tree))) {
      tt_dt <- tdt[Tree == tr & Leaf == TRUE, .(ID, Cover, LeafVal)]
      if (!nrow(tt_dt)) next
      nb <- max(tt_dt$ID) + 1L
      cov <- numeric(nb); val <- numeric(nb)
      cov[tt_dt$ID + 1L] <- as.numeric(tt_dt$Cover)
      val[tt_dt$ID + 1L] <- as.numeric(tt_dt$LeafVal)
      leaf_cover[[tr + 1L]] <- cov; leaf_value[[tr + 1L]] <- val
    }
  }
  
  loglr_sum     <- numeric(n)
  leafdepth_sum <- if (isTRUE(compute_depth))      numeric(n) else NULL
  leafcover_sum <- if (isTRUE(compute_leaf_stats)) numeric(n) else NULL
  leafval_sum   <- if (isTRUE(compute_leaf_stats)) numeric(n) else NULL
  
  for (tt in seq_len(Tm)) {
    tr <- tt - 1L
    lt <- leaves[, tt]
    nb <- nb_vec[tt]; if (nb <= 0L) next
    count_extr_vec <- count_by_leaf(lt, extr_idx, nb)
    count_bg_vec   <- count_by_leaf(lt, bg_idx,   nb)
    pE_leaf  <- count_extr_vec * inv_N_extr
    pB_leaf  <- count_bg_vec   * inv_N_bg
    llr_leaf <- log((pE_leaf + eps) / (pB_leaf + eps))
    loglr_sum <- loglr_sum + llr_leaf[lt + 1L]
    
    if (isTRUE(compute_depth)) {
      dep_vals <- depth_map[[tt]]
      if (!is.null(dep_vals)) {
        dv <- dep_vals[lt + 1L]; dv[!is.finite(dv)] <- 0
        leafdepth_sum <- leafdepth_sum + dv
      }
    }
    if (isTRUE(compute_leaf_stats)) {
      if (!is.null(leaf_cover[[tt]])) leafcover_sum <- leafcover_sum + leaf_cover[[tt]][lt + 1L]
      if (!is.null(leaf_value[[tt]])) leafval_sum   <- leafval_sum   + leaf_value[[tt]][lt + 1L]
    }
    
    if (!is.null(progress_every) && progress_every > 0L &&
        (tt %% progress_every == 0L || tt == Tm)) {
      message(sprintf("[compute_snp_similarity] processed %d/%d trees (%.1f%%)", 
                       tt, Tm, 100 * tt / Tm))
    }
  }
  
  summaries <- data.table::data.table(
    row             = seq_len(n),
    yvar            = as.numeric(yvar_train),
    sim_logLR       = loglr_sum     / Tm,
    mean_leaf_depth = if (isTRUE(compute_depth))      leafdepth_sum / Tm else NA_real_,
    mean_leaf_cover = if (isTRUE(compute_leaf_stats)) leafcover_sum / Tm else NA_real_,
    mean_leaf_value = if (isTRUE(compute_leaf_stats)) leafval_sum   / Tm else NA_real_
  )
  
  list(
    summaries   = summaries,
    leaves      = leaves,
    tdt         = tdt,
    helpers     = helpers,
    extreme_idx = extr_idx,
    bg_idx      = bg_idx,
    tail_dir    = if (lower_tail) "low" else "high"
  )
}

#============================
# Feature concentration 
#============================
# Calculates how concentrated each leaf’s decision path is in its top-k features.
per_leaf_topk_share <- function(leaf_steps, top_k = 2L, progress_every = 5000L) {
  stopifnot(data.table::is.data.table(leaf_steps) || is.data.frame(leaf_steps))
  LS <- data.table::as.data.table(leaf_steps)
  has_long <- "feature" %in% names(LS); has_list <- "path_features" %in% names(LS)
  if (!has_long && !has_list) stop("[per_leaf_topk_share] Need 'feature' (long) or 'path_features' (list).")
  
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
  
  data.table::setorder(counts, Tree, leaf_id, -count)
  perleaf <- counts[, {
    s <- sum(count)
    if (s == 0L) list(topk_share = NA_real_) else list(topk_share = sum(head(count, min(top_k, .N))) / s)
  }, by = .(Tree, leaf_id)]
  
  if (nrow(perleaf) > progress_every) {
    message(sprintf("[per_leaf_topk_share] computed for %s leaves.", format(nrow(perleaf), big.mark=",")))
  }
  data.table::setkey(perleaf, Tree, leaf_id)
  perleaf[]
}

# Aggregates per-leaf feature concentration into per-row (SNP) averages across trees.
feature_concentration_per_snp <- function(leaves, perleaf_share, chunk_rows = 2000L, progress_every = 10000L) {
  message(sprintf("[feature_concentration_per_snp] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
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

#============================
# Harvest path rules (TRAIN)
#============================
# Extracts all candidate rules from training leaves, scoring them on extreme vs. background rows with statistics and medians.
harvest_path_rules <- function(
    leaves,
    leaf_steps,
    yvar_train        = NULL,  # optional if extreme_idx & bg_idx provided
    extreme_idx       = NULL,
    lower_tail        = FALSE,
    extreme_k         = 1,
    bg_idx            = NULL,
    bg_band_k         = NULL,
    band_center       = c("mean","median"),
    band_scale        = c("sd","mad"),
    max_depth         = 4L,
    min_support       = 20L,
    trees_subset      = NULL,   # 0-based tree ids
    pool_identical    = TRUE,
    trees_per_batch   = 250L,
    progress_every    = 1000L,
    tighten_monotone  = TRUE
) {
  message(sprintf("[harvest_path_rules] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  stopifnot(is.matrix(leaves) || inherits(leaves, "dgCMatrix"))
  
  band_center <- match.arg(band_center); band_scale <- match.arg(band_scale)
  
  # yvar is:
  # - optional if both extr/bg idx are supplied (used only for medians if provided)
  # - mandatory if extr/bg idx are not supplied
  if (is.null(extreme_idx) || is.null(bg_idx)) {
    if (is.null(yvar_train)) stop("[harvest_path_rules] yvar_train is required when extreme_idx/bg_idx are not provided.")
    .require_numeric_vec(yvar_train, "harvest_path_rules:yvar_train")
    if (length(yvar_train) != nrow(leaves)) stop("[harvest_path_rules] length(yvar_train) != nrow(leaves).")
    .require_finite_center_scale(yvar_train, "harvest_path_rules:yvar_train")
  } else if (!is.null(yvar_train)) {
    # yvar only used for medians; basic numeric/length check if present
    .require_numeric_vec(yvar_train, "harvest_path_rules:yvar_train")
    if (length(yvar_train) != nrow(leaves)) stop("[harvest_path_rules] length(yvar_train) != nrow(leaves).")
  }
  
  # Globals
  n  <- nrow(leaves); Tm <- ncol(leaves)
  
  # Canonical PATHS (drop unused cols early)
  PATHS_full <- build_path_rule_strings(leaf_steps, max_depth = max_depth, tighten_monotone = tighten_monotone)
  PATHS <- PATHS_full[, .(Tree, leaf_id, rule_str)]
  data.table::setkey(PATHS, Tree, leaf_id)
  
  # Derive/validate indices
  idx <- derive_yvar_idx(
    y           = if (is.null(extreme_idx) || is.null(bg_idx)) yvar_train else numeric(n), # not used if both idx given
    n           = n,
    lower_tail  = lower_tail,
    extreme_k   = extreme_k,
    extreme_idx = extreme_idx,
    bg_idx      = bg_idx,
    bg_band_k   = bg_band_k,
    band_center = band_center,
    band_scale  = band_scale,
    check_y_stats = (is.null(extreme_idx) || is.null(bg_idx) || !is.null(bg_band_k))
  )
  extr_idx  <- idx$extr_idx
  bg_idx    <- idx$bg_idx
  N_extr    <- idx$N_extr
  N_bg      <- idx$N_bg
  base_rate <- idx$base_rate
  
  # Tree subset → column indices
  all_trees0 <- 0:(Tm - 1L)
  if (!is.null(trees_subset)) {
    trees_subset <- as.integer(sort(unique(trees_subset)))
    trees_subset <- trees_subset[trees_subset %in% all_trees0]
    if (!length(trees_subset)) stop("[harvest] trees_subset empty after filtering")
    use_tt <- trees_subset + 1L
  } else {
    use_tt <- seq_len(Tm)
  }
  
  # Precompute nb per tree & rule_len by (Tree,leaf)
  nb_vec <- precompute_nb_vec(leaves)
  
  rule_len_by_leaf <- vector("list", Tm)
  split_paths <- split(PATHS_full[, .(leaf_id, rule_len)], PATHS_full$Tree)
  for (tr_chr in names(split_paths)) {
    tr <- as.integer(tr_chr)
    dt <- split_paths[[tr_chr]]
    if (!nrow(dt)) next
    nb <- max(dt$leaf_id, na.rm = TRUE) + 1L
    rlen <- rep.int(NA_integer_, nb)
    rlen[dt$leaf_id + 1L] <- as.integer(dt$rule_len)
    rule_len_by_leaf[[tr + 1L]] <- rlen
  }
  rm(PATHS_full)
  
  pool_beat <- if (!is.null(progress_every) && progress_every > 0L) progress_every * 5L else NULL
  y_num <- if (!is.null(yvar_train)) as.numeric(yvar_train) else NULL
  
  if (isTRUE(pool_identical)) {
    # ---- POOL: produce R ------------------------------------------------------
    rules_min_list <- vector("list", ceiling(length(use_tt) / trees_per_batch)); rmi <- 1L
    pairs_list     <- vector("list", ceiling(length(use_tt) / trees_per_batch)); pli <- 1L
    
    for (batch_start in seq(1L, length(use_tt), by = trees_per_batch)) {
      batch_end <- min(batch_start + trees_per_batch - 1L, length(use_tt))
      batch_rules_min <- vector("list", batch_end - batch_start + 1L); bi <- 1L
      batch_pairs     <- vector("list", batch_end - batch_start + 1L); bj <- 1L
      
      for (tt in use_tt[batch_start:batch_end]) {
        tr <- tt - 1L
        lt <- as.integer(leaves[, tt])
        nb <- nb_vec[tt]; if (!is.finite(nb) || nb <= 0L) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        count_extr_vec <- count_by_leaf(lt, extr_idx, nb)
        count_bg_vec   <- count_by_leaf(lt, bg_idx,   nb)
        
        leaf_ids <- which((count_extr_vec + count_bg_vec) > 0L) - 1L
        if (!length(leaf_ids)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        rlen <- rule_len_by_leaf[[tt]]
        if (is.null(rlen)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        DT0 <- data.table::data.table(
          Tree      = tr,
          leaf_id   = leaf_ids,
          rule_len  = rlen[leaf_ids + 1L],
          n_extreme = count_extr_vec[leaf_ids + 1L],
          n_bg      = count_bg_vec [leaf_ids + 1L]
        )
        
        DT0 <- DT0[(n_extreme + n_bg) >= as.integer(min_support)]
        if (!nrow(DT0)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        J <- PATHS[DT0, on = .(Tree, leaf_id), nomatch = 0L]
        if (!nrow(J)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        batch_rules_min[[bi]] <- J[, .(rule_str, rule_len)]
        batch_pairs    [[bj]] <- J[, .(Tree, leaf_id, rule_str)]
        bi <- bi + 1L; bj <- bj + 1L
      }
      
      rules_min_list[[rmi]] <- data.table::rbindlist(batch_rules_min, use.names = TRUE, fill = TRUE); rmi <- rmi + 1L
      pairs_list    [[pli]] <- data.table::rbindlist(batch_pairs,     use.names = TRUE, fill = TRUE); pli <- pli + 1L
      
      if (!is.null(progress_every) && progress_every > 0L &&
          ((batch_end %% progress_every) == 0L || batch_end == length(use_tt))) {
        message(sprintf("[harvest] processed tree %d / %d (%.1f%%)",
                        use_tt[batch_end], Tm, 100*batch_end/length(use_tt)))
      }
    }
    
    rules_min <- data.table::rbindlist(rules_min_list, use.names = TRUE, fill = TRUE)
    if (!nrow(rules_min)) {
      warning("[harvest] No rules after filtering; try lowering min_support or increasing max_depth.")
      return(data.table::data.table())
    }
    rule_len_map <- rules_min[, .(rule_len = max(rule_len, na.rm = TRUE)), by = rule_str]
    data.table::setkey(rule_len_map, rule_str)
    
    pairs_all <- data.table::rbindlist(pairs_list, use.names = TRUE, fill = TRUE)
    if (!nrow(pairs_all)) {
      warning("[harvest] No rule pairs to pool after filtering.")
      return(data.table::data.table())
    }
    
    # Row-index caches for pooling
    leaf_rows_index_ext <- vector("list", length = Tm)
    leaf_rows_index_bg  <- vector("list", length = Tm)
    for (tt in use_tt) {
      lt <- as.integer(leaves[, tt])
      leaf_rows_index_ext[[tt]] <- split(extr_idx, lt[extr_idx], drop = TRUE)
      leaf_rows_index_bg [[tt]] <- split(bg_idx,   lt[bg_idx],   drop = TRUE)
    }
    
    uniq_rules <- sort(unique(pairs_all$rule_str))
    pooled <- vector("list", length(uniq_rules))
    
    last_beat <- 0L
    for (i in seq_along(uniq_rules)) {
      rs <- uniq_rules[i]
      pairs <- unique(pairs_all[rule_str == rs, .(Tree, leaf_id)])
      if (!nrow(pairs)) next
      
      rows_ext <- integer(0); rows_bg <- integer(0)
      for (j in seq_len(nrow(pairs))) {
        tt  <- pairs$Tree[j] + 1L
        lid <- as.character(pairs$leaf_id[j])
        re  <- leaf_rows_index_ext[[tt]][[lid]]
        rb  <- leaf_rows_index_bg [[tt]][[lid]]
        if (!is.null(re)) rows_ext <- c(rows_ext, re)
        if (!is.null(rb)) rows_bg  <- c(rows_bg,  rb)
      }
      rows_ext <- unique(rows_ext)
      rows_bg  <- unique(rows_bg)
      
      n_e <- length(rows_ext)
      n_b <- length(rows_bg)
      support <- n_e + n_b
      
      precision <- n_e / pmax(1e-12, support)
      recall    <- if (N_extr > 0) n_e / N_extr else NA_real_
      lift      <- if (base_rate  > 0) precision / base_rate else NA_real_
      pval      <- stats::phyper(q = max(0L, n_e - 1L), m = N_extr, n = N_bg, k = support, lower.tail = FALSE)
      or_ha     <- ((n_e + 0.5) * (N_bg - n_b + 0.5)) / ((N_extr - n_e + 0.5) * (n_b + 0.5))
      
      med_e <- if (!is.null(y_num) && n_e) median(y_num[rows_ext], na.rm = TRUE) else NA_real_
      med_b <- if (!is.null(y_num) && n_b) median(y_num[rows_bg ], na.rm = TRUE) else NA_real_
      med_o <- if (!is.null(y_num) && support > 0L) median(c(y_num[rows_ext], y_num[rows_bg]), na.rm = TRUE) else NA_real_
      
      rl <- rule_len_map[rs, rule_len]
      
      pooled[[i]] <- data.table::data.table(
        rule_str      = rs,
        rule_len      = as.integer(rl),
        n_extreme     = n_e,
        n_bg          = n_b,
        support       = support,
        precision     = precision,
        recall        = recall,
        lift          = lift,
        odds_ratio    = or_ha,
        pval          = pval,
        med_y_extreme = med_e,
        med_y_bg      = med_b,
        med_y_overall = med_o
      )
      
      if (!is.null(pool_beat) && (i - last_beat) >= pool_beat) {
        message(sprintf("[harvest/pool] %d of %d rule strings", i, length(uniq_rules)))
        last_beat <- i
      }
    }
    
    R <- data.table::rbindlist(pooled, use.names = TRUE, fill = TRUE)
    
  } else {
    # ---- UNPOOLED: produce rules, then R <- rules -----------------------------
    out_list <- vector("list", ceiling(length(use_tt) / trees_per_batch)); oi <- 1L
    
    for (batch_start in seq(1L, length(use_tt), by = trees_per_batch)) {
      batch_end <- min(batch_start + trees_per_batch - 1L, length(use_tt))
      batch_rules <- vector("list", batch_end - batch_start + 1L); bi <- 1L
      
      for (tt in use_tt[batch_start:batch_end]) {
        tr <- tt - 1L
        lt <- as.integer(leaves[, tt])
        nb <- nb_vec[tt]; if (!is.finite(nb) || nb <= 0L) { bi <- bi + 1L; next }
        
        count_extr_vec <- count_by_leaf(lt, extr_idx, nb)
        count_bg_vec   <- count_by_leaf(lt, bg_idx,   nb)
        
        leaf_ids <- which((count_extr_vec + count_bg_vec) > 0L) - 1L
        if (!length(leaf_ids)) { bi <- bi + 1L; next }
        
        rlen <- rule_len_by_leaf[[tt]]
        if (is.null(rlen)) { bi <- bi + 1L; next }
        
        DT0 <- data.table::data.table(
          Tree      = tr,
          leaf_id   = leaf_ids,
          rule_len  = rlen[leaf_ids + 1L],
          n_extreme = count_extr_vec[leaf_ids + 1L],
          n_bg      = count_bg_vec [leaf_ids + 1L]
        )
        
        DT0 <- DT0[(n_extreme + n_bg) >= as.integer(min_support)]
        if (!nrow(DT0)) { bi <- bi + 1L; next }
        
        J <- PATHS[DT0, on = .(Tree, leaf_id), nomatch = 0L]
        if (!nrow(J)) { bi <- bi + 1L; next }
        
        support <- J$n_extreme + J$n_bg
        pval    <- stats::phyper(q = pmax(0L, J$n_extreme - 1L), m = N_extr, n = N_bg, k = support, lower.tail = FALSE)
        a <- J$n_extreme; b <- N_extr - J$n_extreme; c <- J$n_bg; d <- N_bg - J$n_bg
        or_ha <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        
        DT <- J[, .(
          rule_str,
          rule_len,
          n_extreme,
          n_bg,
          support,
          precision  = n_extreme / pmax(1e-12, n_extreme + n_bg),
          recall     = if (N_extr > 0) n_extreme / N_extr else NA_real_,
          lift       = if (base_rate  > 0) (n_extreme / pmax(1e-12, n_extreme + n_bg)) / base_rate else NA_real_,
          odds_ratio = or_ha,
          pval       = pval
        )]
        
        batch_rules[[bi]] <- DT; bi <- bi + 1L
      }
      
      out_list[[oi]] <- data.table::rbindlist(batch_rules, use.names = TRUE, fill = TRUE); oi <- oi + 1L
      if (!is.null(progress_every) && progress_every > 0L &&
          ((batch_end %% progress_every) == 0L || batch_end == length(use_tt))) {
        message(sprintf("[harvest] processed tree %d / %d (%.1f%%)",
                        use_tt[batch_end], Tm, 100*batch_end/length(use_tt)))
      }
    }
    
    rules <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
    if (!nrow(rules)) {
      warning("[harvest] No rules after filtering; try lowering min_support or increasing max_depth.")
      return(data.table::data.table())
    }
    R <- rules
  }
  
  # q-values & ordering on R (always)
  R[, qval := p.adjust(pval, method = "BH")]
  data.table::setorderv(R,
                        cols  = c("qval", "lift", "odds_ratio", "precision", "recall"),
                        order = c( 1L,   -1L,     -1L,          -1L,         -1L))
  R[]
}

#===============================================================================
# XGBoost Rule Harvesting & Validation — Canonical Schema & Conventions
#===============================================================================
# All functions in this module MUST conform to these names, shapes, formulas,
# and join/counting styles. If any contract is violated, the function should
# stop() with an informative, prefixed error message.
#-------------------------------------------------------------------------------

#---------------------------
# Global dimensions
#---------------------------
# n      : number of rows (SNPs / observations)
# Tm     : number of trees in the ensemble (ncol(leaves))

#---------------------------
# Indices
#---------------------------
# tt     : 1-based tree index used in R loops (tt ∈ 1..Tm)
# tr     : 0-based tree index (tr = tt - 1), matches xgboost dump column 'Tree'
# extr_idx : integer vector, 1-based row indices for the extreme set (sorted, unique)
# bg_idx   : integer vector, 1-based row indices for the background set (sorted, unique, disjoint from extr_idx)

#---------------------------
# Counts
#---------------------------
# N_extr : length(extr_idx)
# N_bg   : length(bg_idx)

#---------------------------
# Leaves
#---------------------------
# leaves : integer matrix (n × Tm), 0-based leaf IDs; from predict(model, X, predleaf=TRUE)
# lt     : integer vector, leaf IDs for rows in tree tt (lt <- as.integer(leaves[, tt]))
# nb     : integer, number of leaf bins in tree tt (nb = max(lt) + 1L)
# nb_vec : integer vector of nb per tree; precompute once with precompute_nb_vec(leaves)

#---------------------------
# Per-leaf counts (single counting style)
#---------------------------
# count_extr_vec : tabulate(lt[extr_idx] + 1L, nbins = nb)     # via count_by_leaf()
# count_bg_vec   : tabulate(lt[bg_idx]   + 1L, nbins = nb)

#---------------------------
# Canonical per-tree tables
#---------------------------
# DT0 : base table per (Tree, leaf) BEFORE adding rule_str:
#       columns: Tree, leaf_id, rule_len, n_extreme, n_bg
# PATHS : canonical mapping (Tree, leaf_id) → rule_str, built from leaf_steps
# J   : DT0 joined with PATHS on (Tree, leaf_id), adds rule_str
#       columns after join: rule_str + all DT0 columns
# pairs_all : ledger of (Tree, leaf_id, rule_str) AFTER min_support filter (used for pooling)

#---------------------------
# Rule sets (outputs)
#---------------------------
# rules : unpooled per-(Tree,leaf) rules with stats (if pool_identical = FALSE)
# R     : pooled rules (union over rows across trees) with full stats (final return)

#---------------------------
# Stats (single formulas — identical everywhere)
#---------------------------
# support    = n_extreme + n_bg
# precision  = n_extreme / pmax(eps, support)        # eps = 1e-12; ONLY for denominators/log guards
# recall     = n_extreme / N_extr
# base_rate  = N_extr / (N_extr + N_bg)
# lift       = precision / base_rate
# odds_ratio = ((a+0.5)*(d+0.5)) / ((b+0.5)*(c+0.5)), where:
#              a = n_extreme; b = N_extr - n_extreme; c = n_bg; d = N_bg - n_bg
# pval       = phyper(q = max(0L, n_extreme - 1L), m = N_extr, n = N_bg, k = support, lower.tail = FALSE)
# qval       = p.adjust(pval, method = "BH")

#---------------------------
# Column/shape contracts (must be identical everywhere)
#---------------------------
# leaf_steps : Tree, leaf_id, depth, feature, direction, thresh_bin
# PATHS      : Tree, leaf_id, rule_len, rule_str, depth_max
# DT0        : Tree, leaf_id, rule_len, n_extreme, n_bg
# J          : rule_str + (Tree, leaf_id, rule_len, n_extreme, n_bg)  # keep Tree/leaf_id only for ledgers
# pairs_all  : Tree, leaf_id, rule_str
# Final rules / R :
#   rule_str, rule_len, n_extreme, n_bg, support, precision, recall,
#   lift, odds_ratio, pval, qval[, med_y_extreme, med_y_bg, med_y_overall]
# Final sort order: qval (asc), then -lift, -odds_ratio, -precision, -recall

#---------------------------
# Data/type contracts
#---------------------------
# leaves        : integer matrix (n × Tm), 0-based leaf IDs
# extr_idx/bg_idx : integer vectors (1-based indices), sorted, unique, disjoint
# yvar_*        : numeric vectors (length = n) with finite mean/median; sd(y) > 0; mad(y, 1.4826) > 0
# rule_len      : integer ≥ 1
# leaf_id       : integer ≥ 0
# Tree          : integer ≥ 0
# Booleans      : prefix with is_ (e.g., is_extreme, is_bg), length n

#---------------------------
# Join & counting policy
#---------------------------
# - ONE join style: data.table keyed joins only:
#     setkey(PATHS, Tree, leaf_id); J <- PATHS[DT0, on = .(Tree, leaf_id), nomatch = 0L]
# - ONE counting style: ALWAYS tabulate(lt[idx] + 1L, nbins = nb) via count_by_leaf()
# - Ledger: ALWAYS record (Tree, leaf_id, rule_str) AFTER min_support filter; pool using that ledger
# - When pool_identical = TRUE:
#     • In the batch loop: build DT0 → filter by support → join for rule_str
#       Append to a minimal rules_min: (rule_str, rule_len) only (no stats)
#       Append to pairs_all: (Tree, leaf_id, rule_str)
#     • In pooling pass: reconstruct n_e, n_b via union of rows across trees,
#       then compute support/precision/recall/lift/pval/odds_ratio/medians;
#       take rule_len = max(rule_len) from rules_min[rule_str]

#---------------------------
# Options / knobs (consistent names across functions)
#---------------------------
# tighten_monotone, bg_band_k, band_center ∈ {"mean","median"}, band_scale ∈ {"sd","mad"}
# extreme_k, lower_tail, trees_subset
# trees_per_batch, progress_every (pool_progress_every = 5 × progress_every)
# max_depth, min_support, target_bins, min_per_bin, winsor_prob
# eps = 1e-12 (ONLY for denominators/log guards — NEVER to “fix” pvals)

#---------------------------
# Error/validation conventions
#---------------------------
# Prefix all messages with function tag in brackets, e.g. "[validate] ..."
# stop() on:
#   - Non-numeric yvar; length mismatch to n; non-finite mean/median
#   - sd(y) <= 0 or mad(y, 1.4826) <= 0
#   - No rows in extreme or background set
#   - Overlap between extreme and background indices
#   - Tree count mismatch (xgb dump vs leaves)
#   - trees_subset empty after filtering
#   - candidate_rules missing 'rule_str'
#   - leaf_steps missing required columns
#   - Any NA from phyper
#   - Non-numeric medians
# warn() on:
#   - No rules after filters; no pooled support

#---------------------------
# Determinism & indices
#---------------------------
# - Tree IDs in dumps are 0-based (tr); loop index is 1-based (tt); always convert with tr = tt - 1
# - Row indices extr_idx/bg_idx are 1-based R indices
# - Avoid silent recycling/clamping; fail fast with informative errors

#---------------------------
# Performance notes
#---------------------------
# - precompute_nb_vec(leaves) is a single pass to obtain nb per tree
# - Use batched tree loops; precompute rule_len_by_leaf
# - For pooling, rely on pairs_all ledger instead of scanning full rules tables
# - Before loops, drop unused PATHS columns (keep only Tree, leaf_id, rule_str)

#---------------------------
# Logging
#---------------------------
# - Use message() with a consistent template:
#     "[stage] processed %d/%d (%.1f%%)"
#     "[harvest/pool] %d of %d rule strings"
#===============================================================================

#============================
# Shared validators and helpers
#============================
.require_numeric_vec <- function(x, name) {
  if (!is.numeric(x)) stop(sprintf("[%s] must be numeric.", name))
}

.require_finite_center_scale <- function(y, name) {
  m1 <- mean(y, na.rm = TRUE)
  m2 <- stats::median(y, na.rm = TRUE)
  sdv <- stats::sd(y, na.rm = TRUE)
  md  <- stats::mad(y, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(m1)) stop(sprintf("[%s] mean is NA/Inf.", name))
  if (!is.finite(m2)) stop(sprintf("[%s] median is NA/Inf.", name))
  if (!is.finite(sdv) || sdv <= 0) stop(sprintf("[%s] sd <= 0 or NA/Inf.", name))
  if (!is.finite(md)  || md  <= 0) stop(sprintf("[%s] mad <= 0 or NA/Inf.", name))
}

#============================
# Small shared helpers
#============================
# Derive extremes + background indices from yvar, based on explicit input or statistical banding around the center.
derive_yvar_idx <- function(
    y, n,
    lower_tail  = FALSE,
    extreme_k   = 1,
    extreme_idx = NULL,
    bg_idx      = NULL,
    bg_band_k   = NULL,
    band_center = c("mean","median"),
    band_scale  = c("sd","mad"),
    check_y_stats = NULL
) {
  band_center <- match.arg(band_center)
  band_scale  <- match.arg(band_scale)
  
  # Decide if we need y-stat checks:
  # default: TRUE unless BOTH idx sets are provided and no banding is requested
  if (is.null(check_y_stats)) {
    check_y_stats <- (is.null(extreme_idx) || is.null(bg_idx) || !is.null(bg_band_k))
  }
  
  # Only validate y if needed
  if (isTRUE(check_y_stats)) {
    if (!is.numeric(y)) stop("[derive_yvar_idx] y must be numeric.")
    if (length(y) != n) stop("[derive_yvar_idx] length(y) must equal n.")
    mu  <- mean(y, na.rm = TRUE);  if (!is.finite(mu))  stop("[derive_yvar_idx] mean(y) not finite.")
    md  <- stats::median(y, na.rm = TRUE); if (!is.finite(md))  stop("[derive_yvar_idx] median(y) not finite.")
    sdv <- stats::sd(y, na.rm = TRUE);     if (!is.finite(sdv) || sdv <= 0) stop("[derive_yvar_idx] sd(y) <= 0 or non-finite.")
    madv<- stats::mad(y, constant = 1.4826, na.rm = TRUE); if (!is.finite(madv) || madv <= 0) stop("[derive_yvar_idx] mad(y) <= 0 or non-finite.")
  }
  
  # Build extremes
  if (is.null(extreme_idx)) {
    mu  <- mean(y, na.rm = TRUE); sdv <- sd(y, na.rm = TRUE)
    thr <- if (!lower_tail) mu + extreme_k * sdv else mu - extreme_k * sdv
    extr_idx <- if (!lower_tail) which(y >= thr) else which(y <= thr)
  } else {
    extr_idx <- sort(unique(as.integer(extreme_idx)))
    extr_idx <- extr_idx[extr_idx >= 1L & extr_idx <= n]
  }
  
  # Build background
  if (!is.null(bg_idx)) {
    bg_idx <- sort(unique(as.integer(bg_idx)))
    bg_idx <- bg_idx[bg_idx >= 1L & bg_idx <= n]
    if (length(intersect(extr_idx, bg_idx))) bg_idx <- setdiff(bg_idx, extr_idx)
  } else if (!is.null(bg_band_k)) {
    ctr <- if (band_center == "mean") mean(y, na.rm = TRUE) else stats::median(y, na.rm = TRUE)
    scl <- if (band_scale  == "sd")   stats::sd(y, na.rm = TRUE) else stats::mad(y, constant = 1.4826, na.rm = TRUE)
    lo  <- ctr - bg_band_k * scl; hi <- ctr + bg_band_k * scl
    bg_idx <- which(y >= lo & y <= hi)
    if (length(extr_idx)) bg_idx <- setdiff(bg_idx, extr_idx)
  } else {
    bg_idx <- setdiff(seq_len(n), extr_idx)
  }
  
  # Final checks
  N_extr <- length(extr_idx); N_bg <- length(bg_idx)
  if (N_extr <= 0L) stop("[derive_yvar_idx] No rows in extreme set.")
  if (N_bg   <= 0L) stop("[derive_yvar_idx] No rows in background set.")
  if (length(intersect(extr_idx, bg_idx)) > 0L) stop("[derive_yvar_idx] extreme and background indices overlap.")
  
  is_extreme <- logical(n); is_extreme[extr_idx] <- TRUE
  is_bg      <- logical(n); is_bg[bg_idx]       <- TRUE
  
  list(
    extr_idx  = extr_idx,
    bg_idx    = bg_idx,
    N_extr    = N_extr,
    N_bg      = N_bg,
    base_rate = N_extr / (N_extr + N_bg),
    is_extreme = is_extreme,
    is_bg      = is_bg
  )
}

# Precomputes the number of leaf bins (`nb`) for each tree column in a leaves matrix.
precompute_nb_vec <- function(leaves) {
  Tm <- ncol(leaves)
  # dgCMatrix fast-path via Matrix::colMaxs (not matrixStats)
  if (inherits(leaves, "dgCMatrix")) {
    if (requireNamespace("Matrix", quietly = TRUE)) {
      mx <- Matrix::colMaxs(leaves, na.rm = TRUE)
      mx[!is.finite(mx)] <- -1
      return(as.integer(mx + 1L))
    }
  }
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    mx <- matrixStats::colMaxs(leaves, na.rm = TRUE)
    mx[!is.finite(mx)] <- -1
    return(as.integer(mx + 1L))
  }
  # Fallback loop (no suppressWarnings; handle NA explicitly)
  res <- integer(Tm)
  for (tt in seq_len(Tm)) {
    col <- leaves[, tt]
    mx  <- if (anyNA(col)) max(col, na.rm = TRUE) else max(col)
    if (!is.finite(mx)) res[tt] <- 0L else res[tt] <- as.integer(mx) + 1L
  }
  res
}

# Tabulates counts of a given index set into leaf bins of a single tree.
count_by_leaf <- function(lt, idx, nb) {
  if (!length(idx) || nb <= 0L) return(integer(max(1L, nb)))
  tabulate(lt[idx] + 1L, nbins = nb)
}

# Learn a fixed extreme/background rule on TRAIN and reapply it to any y.
# Canonical names, disjoint sets, and selectable center/scale for BG band.
learn_sigma_rule <- function(yvar_train,
                             extreme_k    = 1,
                             lower_tail   = FALSE,
                             bg_band_k    = NULL,                 # NULL → background = !extreme
                             band_center  = c("mean","median"),   # BG center learned on TRAIN
                             band_scale   = c("sd","mad")) {      # BG scale learned on TRAIN
  band_center <- match.arg(band_center)
  band_scale  <- match.arg(band_scale)
  
  # --- Validate TRAIN y -------------------------------------------------------
  if (!is.numeric(yvar_train)) stop("[learn_sigma_rule] yvar_train must be numeric.")
  if (!length(yvar_train))     stop("[learn_sigma_rule] yvar_train is empty.")
  
  center_val <- if (band_center == "mean") {
    m <- mean(yvar_train, na.rm = TRUE)
    if (!is.finite(m)) stop("[learn_sigma_rule] mean(yvar_train) is not finite.")
    m
  } else {
    med <- stats::median(yvar_train, na.rm = TRUE)
    if (!is.finite(med)) stop("[learn_sigma_rule] median(yvar_train) is not finite.")
    med
  }
  
  # Scale for extreme threshold is always SD (to mirror common sigma rules)
  sd_train <- stats::sd(yvar_train, na.rm = TRUE)
  if (!is.finite(sd_train) || sd_train <= 0)
    stop("[learn_sigma_rule] sd(yvar_train) <= 0 or not finite.")
  
  # Scale for BG band can be SD or MAD (robust), learned on TRAIN
  scale_val <- if (band_scale == "sd") {
    s <- stats::sd(yvar_train, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) stop("[learn_sigma_rule] sd(yvar_train) for BG is not finite/positive.")
    s
  } else {
    md <- stats::mad(yvar_train, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(md) || md <= 0) stop("[learn_sigma_rule] mad(yvar_train) for BG is not finite/positive.")
    md
  }
  
  # --- Fixed thresholds learned on TRAIN --------------------------------------
  thr_ext <- if (!lower_tail) center_val + extreme_k * sd_train
  else             center_val - extreme_k * sd_train
  
  thr_bg  <- if (is.null(bg_band_k)) NULL else {
    lo <- center_val - bg_band_k * scale_val
    hi <- center_val + bg_band_k * scale_val
    c(lo, hi)
  }
  
  # --- Return an object with an apply() that enforces the SAME rule -----------
  list(
    # Apply learned thresholds to any numeric vector f (train/test/heldout)
    apply = function(f) {
      if (!is.numeric(f)) stop("[learn_sigma_rule/apply] input vector must be numeric.")
      n <- length(f); if (!n) stop("[learn_sigma_rule/apply] input vector is empty.")
      
      # extremes by the fixed TRAIN threshold
      extr_idx <- if (!lower_tail) which(f >= thr_ext) else which(f <= thr_ext)
      
      # background either as band (TRAIN-learned) or !extreme
      if (is.null(thr_bg)) {
        bg_idx <- setdiff(seq_len(n), extr_idx)
      } else {
        bg_idx <- which(f >= thr_bg[1] & f <= thr_bg[2])
        # enforce disjointness explicitly
        if (length(extr_idx)) bg_idx <- setdiff(bg_idx, extr_idx)
      }
      
      # Canonical outputs (sorted, unique)
      extr_idx <- sort(unique(as.integer(extr_idx)))
      bg_idx   <- sort(unique(as.integer(bg_idx)))
      
      N_extr <- length(extr_idx)
      N_bg   <- length(bg_idx)
      
      if (N_extr <= 0L) stop("[learn_sigma_rule/apply] No rows in extreme set under learned rule.")
      if (N_bg   <= 0L) stop("[learn_sigma_rule/apply] No rows in background set under learned rule.")
      if (length(intersect(extr_idx, bg_idx)) > 0L)
        stop("[learn_sigma_rule/apply] extreme and background indices overlap (should be disjoint).")
      
      list(
        extr_idx = extr_idx,
        bg_idx   = bg_idx,
        N_extr   = N_extr,
        N_bg     = N_bg
      )
    },
    
    # Expose learned parameters for reproducibility/auditing
    thr_ext      = thr_ext,                 # scalar extreme threshold (TRAIN)
    thr_bg       = thr_bg,                  # length-2 band [lo, hi] or NULL
    center_name  = band_center,
    scale_name   = band_scale,
    center_val   = center_val,              # learned on TRAIN
    sd_train     = sd_train,                # used for extreme threshold
    scale_val    = scale_val,               # used for BG band
    extreme_k    = extreme_k,
    bg_band_k    = bg_band_k,
    lower_tail   = lower_tail
  )
}

#============================
# Parse model → compact tree table (model → tdt)
#============================
# Converts an `xgb.Booster` into a tidy, numeric, per-node `data.table` with consistent schema (Tree, ID, children, split, gain, leaf values, etc.).
parse_xgb_tree <- function(model) {
  stopifnot(inherits(model, "xgb.Booster"))
  dt <- xgboost::xgb.model.dt.tree(model = model)
  data.table::setDT(dt)
  
  to_int_child <- function(x) {
    if (is.null(x)) return(NA_integer_)
    y <- suppressWarnings(as.integer(sub(".*-", "", x)))
    y[is.na(x)] <- NA_integer_
    y
  }
  to_num <- function(x) suppressWarnings(as.numeric(x))
  
  if (!"Node" %in% names(dt)) stop("[parse_xgb_tree] column 'Node' not found in xgb dump")
  
  dt[, `:=`(Tree = as.integer(Tree), ID = as.integer(Node))]
  
  for (col in c("Yes","No","Missing")) {
    if (col %in% names(dt)) dt[, (col) := to_int_child(get(col))] else dt[, (col) := NA_integer_]
  }
  
  if (!"Split" %in% names(dt)) dt[, Split := NA_real_]
  if (!"Cover" %in% names(dt)) dt[, Cover := NA_real_]
  dt[, `:=`(Split = to_num(Split), Cover = to_num(Cover))]
  
  dt[, Leaf := (Feature == "Leaf")]
  
  has_quality <- "Quality" %in% names(dt)
  has_gaincol <- "Gain" %in% names(dt)
  Q <- if (has_quality) to_num(dt$Quality) else rep(NA_real_, nrow(dt))
  G <- if (has_gaincol) to_num(dt$Gain)    else rep(NA_real_, nrow(dt))
  
  dt[, LeafVal := ifelse(Leaf, Q, NA_real_)]
  dt[Leaf & is.na(LeafVal), LeafVal := Split]
  dt[, Gain := ifelse(!Leaf, ifelse(!is.na(G), G, Q), NA_real_)]
  
  drop_cols <- intersect(c("Node","Quality"), names(dt))
  if (length(drop_cols)) dt[, (drop_cols) := NULL]
  
  data.table::setkey(dt, Tree, ID)
  data.table::setorder(dt, Tree, ID)
  dt[, .(Tree, ID, Feature, Split, Yes, No, Missing, Leaf, LeafVal, Gain, Cover)]
}


#============================#
# Structure helpers
#============================#
# Constructs parent/child maps and depth/path lookup functions for tree navigation.
build_structure_helpers <- function(tdt) {
  stopifnot(is.data.frame(tdt),
            all(c("Tree","ID","Yes","No","Missing","Leaf","Gain","Cover") %in% names(tdt)))
  tdt <- data.table::as.data.table(tdt)
  trees <- sort(unique(tdt$Tree)) #
  
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
  data.table::setkey(split_gc, Tree, ID)
  
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
    cur <- leaf_id; k <- 0L
    while (!is.na(cur) && cur != 0L) { par <- pmap[cur + 1L]; if (is.na(par)) break; k <- k + 1L; cur <- par }
    if (k == 0L) {
      out <- list(nodes=integer(0), gain=numeric(0), gcover=numeric(0))
      leaf_cache[[tkey]][[as.character(leaf_id)]] <- out; return(out)
    }
    nodes <- integer(k); cur <- leaf_id
    for (pos in k:1) { par <- pmap[cur + 1L]; nodes[pos] <- par; cur <- par; if (is.na(cur) || cur == 0L) break }
    sgc <- split_gc[.(tr, nodes), .(Gain, SplitCover)]
    out <- list(nodes = nodes, gain = as.numeric(sgc$Gain), gcover = as.numeric(sgc$SplitCover))
    leaf_cache[[tkey]][[as.character(leaf_id)]] <- out
    out
  }
  
  list(get_depth_map = get_depth_map,
       get_max_depth = get_max_depth,
       get_leaf_path = get_leaf_path)
}

#============================
# Build per-leaf steps
#============================
# Generates per-leaf decision paths (feature, direction, threshold bins) from tree structure and leaves
build_leaf_steps <- function(leaves,
                             tdt,
                             helpers,
                             target_bins     = 10L,
                             min_per_bin     = 50L,
                             winsor_prob     = 0.01,
                             method          = c("fd","quantile"),
                             trees_per_batch = 250L,
                             progress_every  = 1000L) {
  message(sprintf("[build_leaf_steps] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  stopifnot(is.matrix(leaves) || inherits(leaves, "dgCMatrix"))
  method <- match.arg(method)
  Tm <- ncol(leaves)
  
  tdt <- data.table::as.data.table(tdt)
  data.table::setkey(tdt, Tree, ID)
  
  build_raw_steps_batch <- function(batch_start, batch_end) {
    total_steps <- 0L
    meta <- vector("list", batch_end - batch_start + 1L); mi <- 1L
    for (tt in batch_start:batch_end) {
      tr <- tt - 1L
      lt <- leaves[, tt]
      uleaf <- sort(unique(as.integer(lt)))
      meta[[mi]] <- list(tt = tt, tr = tr, uleaf = uleaf)
      if (length(uleaf)) {
        for (leaf_id in uleaf) {
          pn <- helpers$get_leaf_path(tr, leaf_id)$nodes
          if (length(pn)) total_steps <- total_steps + length(pn)
        }
      }
      mi <- mi + 1L
    }
    if (total_steps == 0L) {
      return(data.table::data.table(Tree = integer(0), leaf_id = integer(0),
                                    depth = integer(0), feature = character(0),
                                    direction = character(0), split_val = numeric(0)))
    }
    
    Tree_v  <- integer(total_steps)
    Leaf_v  <- integer(total_steps)
    Depth_v <- integer(total_steps)
    Feat_v  <- character(total_steps)
    Dir_v   <- character(total_steps)
    Split_v <- numeric(total_steps)
    cursor  <- 0L
    
    for (m in meta) {
      tt <- m$tt; tr <- m$tr; uleaf <- m$uleaf
      if (!length(uleaf)) next
      tt_dt <- tdt[.(tr)]
      if (!nrow(tt_dt)) next
      
      maxid <- max(tt_dt$ID, na.rm = TRUE)
      yesA  <- rep.int(NA_integer_,  maxid + 1L);  yesA [tt_dt$ID + 1L] <- tt_dt$Yes
      noA   <- rep.int(NA_integer_,  maxid + 1L);  noA  [tt_dt$ID + 1L] <- tt_dt$No
      featA <- rep.int(NA_character_,maxid + 1L);  featA[tt_dt$ID + 1L] <- as.character(tt_dt$Feature)
      spltA <- rep.int(NA_real_,     maxid + 1L);  spltA[tt_dt$ID + 1L] <- suppressWarnings(as.numeric(tt_dt$Split))
      
      for (leaf_id in uleaf) {
        pn <- helpers$get_leaf_path(tr, leaf_id)$nodes
        k  <- length(pn)
        if (!k) next
        child_along <- integer(k)
        if (k > 1L) child_along[1:(k-1L)] <- pn[2:k]
        child_along[k] <- leaf_id
        parents <- pn + 1L
        yes  <- yesA [parents]
        no   <- noA  [parents]
        feat <- featA[parents]
        splt <- spltA[parents]
        dir <- ifelse(yes == child_along, "<", ifelse(no == child_along, ">=", "missing"))
        
        rng <- (cursor + 1L):(cursor + k)
        Tree_v [rng] <- tr
        Leaf_v [rng] <- leaf_id
        Depth_v[rng] <- 0L:(k - 1L)
        Feat_v [rng] <- feat
        Dir_v  [rng] <- dir
        Split_v[rng] <- splt
        cursor <- cursor + k
      }
    }
    
    data.table::data.table(
      Tree      = Tree_v,
      leaf_id   = Leaf_v,
      depth     = Depth_v,
      feature   = Feat_v,
      direction = Dir_v,
      split_val = Split_v
    )
  }
  
  out_list <- vector("list", ceiling(Tm / trees_per_batch)); oi <- 1L
  for (batch_start in seq(1L, Tm, by = trees_per_batch)) {
    batch_end <- min(batch_start + trees_per_batch - 1L, Tm)
    out_list[[oi]] <- build_raw_steps_batch(batch_start, batch_end); oi <- oi + 1L
    if (!is.null(progress_every) && progress_every > 0L &&
        ((batch_end %% progress_every) == 0L || batch_end == Tm)) {
      message(sprintf("[build_leaf_steps] built steps up to tree %d / %d (%.1f%%)",
                      batch_end, Tm, 100*batch_end/Tm))
    }
  }
  
  LS <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (!nrow(LS)) {
    data.table::setkey(LS, Tree, leaf_id)
    return(LS[, .(Tree, leaf_id, depth, feature, direction, thresh_bin = numeric())])
  }
  
  LS <- LS[!is.na(feature) & nzchar(feature)]
  LS[, `:=`(
    Tree      = as.integer(Tree),
    leaf_id   = as.integer(leaf_id),
    depth     = as.integer(depth),
    feature   = as.character(feature),
    direction = as.character(direction),
    split_val = as.numeric(split_val)
  )]
  
  message(sprintf("[build_leaf_steps] adaptive bins per feature: target=%d, min_per_bin=%d, method=%s, winsor=%.3g",
        as.integer(target_bins), as.integer(min_per_bin), method, as.numeric(winsor_prob)))
  
  LS_split <- LS[is.finite(split_val)]
  feats <- sort(unique(LS_split$feature))
  breaks_env <- new.env(parent = emptyenv()); mids_env <- new.env(parent = emptyenv())
  
  propose_breaks <- function(v, max_bins, method, winsor) {
    v <- sort(v[is.finite(v)]); if (length(v) <= 1L) { r <- range(v, na.rm=TRUE); if(!all(is.finite(r))) r <- c(0,0); return(unique(c(r[1], r[2]))) }
    lo <- suppressWarnings(stats::quantile(v, probs = winsor, names = FALSE))
    hi <- suppressWarnings(stats::quantile(v, probs = 1 - winsor, names = FALSE))
    v  <- v[v >= lo & v <= hi]; if (length(v) <= 1L) { r <- range(v, na.rm=TRUE); if(!all(is.finite(r))) r <- c(0,0); return(unique(c(r[1], r[2]))) }
    if (method == "fd") {
      h <- 2 * stats::IQR(v) / (length(v)^(1/3)); if (!is.finite(h) || h <= 0) h <- (max(v) - min(v)) / max(2L, max_bins)
      nb <- ceiling((max(v) - min(v)) / max(h, .Machine$double.eps)); nb <- min(max(nb, 2L), max_bins)
      br <- pretty(v, n = nb)
    } else {
      nb <- max(2L, as.integer(max_bins))
      probs <- seq(0, 1, length.out = nb + 1L)
      br <- unique(stats::quantile(v, probs = probs, names = FALSE, type = 7))
      if (length(br) < 3L) br <- pretty(v, n = min(max_bins, 3L))
    }
    br <- sort(unique(as.numeric(br))); if (length(br) < 2L) br <- unique(c(min(v), max(v))); br
  }
  enforce_min_bin <- function(v, br, target_bins, min_per_bin) {
    if (length(br) < 2L) return(br)
    repeat {
      idx <- findInterval(v, br, all.inside = TRUE)
      tab <- tabulate(idx, nbins = length(br) - 1L)
      if ((length(tab) <= target_bins) && all(tab >= min_per_bin | tab == 0L)) break
      k <- if (length(tab)) which.min(tab) else 1L
      if (length(tab) <= 1L) break
      rm_pos <- if (k == 1L) 2L else if (k == length(tab)) length(tab) else if (tab[k - 1L] <= tab[k + 1L]) k else k + 1L
      br <- br[-rm_pos]; if (length(br) < 2L) { br <- br[1:2]; break }
    }
    br
  }
  
  for (f in feats) {
    v  <- LS_split[feature == f, split_val]; if (!length(v)) next
    br <- propose_breaks(v, max_bins = as.integer(target_bins), method = method, winsor = as.numeric(winsor_prob))
    br <- enforce_min_bin(v, br, target_bins = as.integer(target_bins), min_per_bin = as.integer(min_per_bin))
    if (length(br) < 2L) { r <- range(v, na.rm=TRUE); if(!all(is.finite(r))) r <- c(0,0); br <- unique(c(r[1], r[2])) }
    mids <- (br[-1L] + br[-length(br)]) / 2
    breaks_env[[f]] <- br; mids_env[[f]] <- mids
  }
  
  LS[, thresh_bin := {
    br <- breaks_env[[feature[1L]]]; md <- mids_env[[feature[1L]]]
    if (is.null(br) || is.null(md) || !is.finite(split_val[1L])) rep(NA_real_, .N) else {
      idx <- findInterval(split_val, br, all.inside = TRUE); as.numeric(md[idx])
    }
  }, by = feature]
  
  LS[, split_val := NULL]
  data.table::setkey(LS, Tree, leaf_id)
  LS[]
}

#============================
# Build (Tree,leaf) rule strings
#============================
.format_num <- function(x) formatC(as.numeric(x), format = "e", digits = 2)

# Converts per-leaf decision paths into canonical rule strings, optionally tightening redundant splits.
build_path_rule_strings <- function(leaf_steps, max_depth = 4L, tighten_monotone = TRUE) {
  LS <- data.table::as.data.table(leaf_steps)
  need <- c("Tree","leaf_id","depth","feature","direction","thresh_bin")
  miss <- setdiff(need, names(LS))
  if (length(miss)) stop(sprintf("[paths] leaf_steps missing: %s", paste(miss, collapse=", ")))
  LS[, `:=`(
    Tree       = as.integer(Tree),
    leaf_id    = as.integer(leaf_id),
    depth      = as.integer(depth),
    feature    = as.character(feature),
    direction  = as.character(direction),
    thresh_bin = as.numeric(thresh_bin)
  )]
  LS <- LS[depth < as.integer(max_depth)]
  data.table::setorder(LS, Tree, leaf_id, depth)
  
  if (isTRUE(tighten_monotone)) {
    eps <- 1e-12
    LS_tight <- LS[, {
      out <- list(cond = character(0), depth = integer(0))
      best_ge <- new.env(parent = emptyenv()); best_lt <- new.env(parent = emptyenv())
      for (i in seq_len(.N)) {
        f <- feature[i]; d <- direction[i]; b <- thresh_bin[i]; keep <- TRUE
        if (d == ">=") { cur <- best_ge[[f]]; if (is.null(cur) || b > cur + eps) best_ge[[f]] <- b else keep <- FALSE }
        else if (d == "<") { cur <- best_lt[[f]]; if (is.null(cur) || b < cur - eps) best_lt[[f]] <- b else keep <- FALSE }
        if (keep) { out$cond <- c(out$cond, paste0(f, " ", d, " ", .format_num(b))); out$depth <- c(out$depth, depth[i]) }
      }
      if (length(out$cond)) .(cond = out$cond, depth = out$depth) else NULL
    }, by = .(Tree, leaf_id)]
  } else {
    LS_tight <- LS[, .(cond = paste0(feature, " ", direction, " ", .format_num(thresh_bin)),
                       depth = depth),
                   by = .(Tree, leaf_id)]
  }
  
  PATHS <- LS_tight[, .(
    rule_len  = .N,
    rule_str  = paste(cond, collapse = " | "),
    depth_max = max(depth)
  ), by = .(Tree, leaf_id)]
  data.table::setkey(PATHS, Tree, leaf_id)
  PATHS[]
}

#============================
# Path similarity (per-SNP)
#============================
# Computes per-row similarity scores (log-likelihood ratios, mean depth, cover, leaf value) between extreme and background sets.
compute_snp_similarity <- function(model, X, yvar_train,
                                   lower_tail     = FALSE,
                                   extreme_k      = 1,
                                   extreme_idx    = NULL,
                                   bg_idx         = NULL,
                                   bg_band_k      = NULL,
                                   band_center    = c("mean","median"),
                                   band_scale     = c("sd","mad"),
                                   progress_every = NULL,
                                   leaves_override    = NULL,
                                   compute_depth      = TRUE,
                                   compute_leaf_stats = TRUE) {
  message(sprintf("[compute_snp_similarity] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  stopifnot(inherits(model, "xgb.Booster"))
  X <- if (inherits(X, "dgCMatrix") || is.matrix(X)) X else as.matrix(X)
  .require_numeric_vec(yvar_train, "compute_snp_similarity:yvar_train")
  if (length(yvar_train) != nrow(X)) stop("[compute_snp_similarity] length(yvar_train) != nrow(X).")
  .require_finite_center_scale(yvar_train, "compute_snp_similarity:yvar_train")
  band_center <- match.arg(band_center); band_scale <- match.arg(band_scale)
  
  n <- nrow(X)
  idx <- derive_yvar_idx(
    y             = yvar_train,
    n             = n,
    lower_tail    = lower_tail,         
    extreme_k     = extreme_k,             
    extreme_idx   = extreme_idx,
    bg_idx        = bg_idx,          
    bg_band_k     = bg_band_k,
    band_center   = band_center,
    band_scale    = band_scale
  )
  
  extr_idx   <- idx$extr_idx
  bg_idx     <- idx$bg_idx
  N_extr     <- idx$N_extr
  N_bg       <- idx$N_bg

  inv_N_extr <- if (N_extr > 0L) 1 / N_extr else 0
  inv_N_bg   <- if (N_bg > 0L) 1 / N_bg else 0
  eps <- 1e-12
  
  leaves <- if (is.null(leaves_override)) predict(model, X, predleaf = TRUE) else leaves_override
  if (is.null(dim(leaves))) leaves <- matrix(leaves, ncol = 1L)
  if (nrow(leaves) != n) stop("[compute_snp_similarity] nrow(leaves) must equal nrow(X).")
  storage.mode(leaves) <- "integer"
  Tm <- ncol(leaves)
  tdt <- parse_xgb_tree(model)
  if ((max(tdt$Tree) + 1L) != Tm) stop(sprintf("[compute_snp_similarity] Tree count mismatch: dump=%d vs leaves=%d.", max(tdt$Tree)+1L, Tm))
  
  helpers <- build_structure_helpers(tdt)
  nb_vec <- precompute_nb_vec(leaves)
  
  depth_map <- if (isTRUE(compute_depth)) {
    out <- vector("list", Tm); for (tt in seq_len(Tm)) out[[tt]] <- helpers$get_depth_map(tt - 1L); out
  } else NULL
  
  leaf_cover <- leaf_value <- NULL
  if (isTRUE(compute_leaf_stats)) {
    leaf_cover <- vector("list", Tm); leaf_value <- vector("list", Tm)
    for (tr in sort(unique(tdt$Tree))) {
      tt_dt <- tdt[Tree == tr & Leaf == TRUE, .(ID, Cover, LeafVal)]
      if (!nrow(tt_dt)) next
      nb <- max(tt_dt$ID) + 1L
      cov <- numeric(nb); val <- numeric(nb)
      cov[tt_dt$ID + 1L] <- as.numeric(tt_dt$Cover)
      val[tt_dt$ID + 1L] <- as.numeric(tt_dt$LeafVal)
      leaf_cover[[tr + 1L]] <- cov; leaf_value[[tr + 1L]] <- val
    }
  }
  
  loglr_sum     <- numeric(n)
  leafdepth_sum <- if (isTRUE(compute_depth))      numeric(n) else NULL
  leafcover_sum <- if (isTRUE(compute_leaf_stats)) numeric(n) else NULL
  leafval_sum   <- if (isTRUE(compute_leaf_stats)) numeric(n) else NULL
  
  for (tt in seq_len(Tm)) {
    tr <- tt - 1L
    lt <- leaves[, tt]
    nb <- nb_vec[tt]; if (nb <= 0L) next
    count_extr_vec <- count_by_leaf(lt, extr_idx, nb)
    count_bg_vec   <- count_by_leaf(lt, bg_idx,   nb)
    pE_leaf  <- count_extr_vec * inv_N_extr
    pB_leaf  <- count_bg_vec   * inv_N_bg
    llr_leaf <- log((pE_leaf + eps) / (pB_leaf + eps))
    loglr_sum <- loglr_sum + llr_leaf[lt + 1L]
    
    if (isTRUE(compute_depth)) {
      dep_vals <- depth_map[[tt]]
      if (!is.null(dep_vals)) {
        dv <- dep_vals[lt + 1L]; dv[!is.finite(dv)] <- 0
        leafdepth_sum <- leafdepth_sum + dv
      }
    }
    if (isTRUE(compute_leaf_stats)) {
      if (!is.null(leaf_cover[[tt]])) leafcover_sum <- leafcover_sum + leaf_cover[[tt]][lt + 1L]
      if (!is.null(leaf_value[[tt]])) leafval_sum   <- leafval_sum   + leaf_value[[tt]][lt + 1L]
    }
    
    if (!is.null(progress_every) && progress_every > 0L &&
        (tt %% progress_every == 0L || tt == Tm)) {
      message(sprintf("[compute_snp_similarity] processed %d/%d trees (%.1f%%)", 
                       tt, Tm, 100 * tt / Tm))
    }
  }
  
  summaries <- data.table::data.table(
    row             = seq_len(n),
    yvar            = as.numeric(yvar_train),
    sim_logLR       = loglr_sum     / Tm,
    mean_leaf_depth = if (isTRUE(compute_depth))      leafdepth_sum / Tm else NA_real_,
    mean_leaf_cover = if (isTRUE(compute_leaf_stats)) leafcover_sum / Tm else NA_real_,
    mean_leaf_value = if (isTRUE(compute_leaf_stats)) leafval_sum   / Tm else NA_real_
  )
  
  list(
    summaries   = summaries,
    leaves      = leaves,
    tdt         = tdt,
    helpers     = helpers,
    extreme_idx = extr_idx,
    bg_idx      = bg_idx,
    tail_dir    = if (lower_tail) "low" else "high"
  )
}

#============================
# Feature concentration 
#============================
# Calculates how concentrated each leaf’s decision path is in its top-k features.
per_leaf_topk_share <- function(leaf_steps, top_k = 2L, progress_every = 5000L) {
  stopifnot(data.table::is.data.table(leaf_steps) || is.data.frame(leaf_steps))
  LS <- data.table::as.data.table(leaf_steps)
  has_long <- "feature" %in% names(LS); has_list <- "path_features" %in% names(LS)
  if (!has_long && !has_list) stop("[per_leaf_topk_share] Need 'feature' (long) or 'path_features' (list).")
  
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
  
  data.table::setorder(counts, Tree, leaf_id, -count)
  perleaf <- counts[, {
    s <- sum(count)
    if (s == 0L) list(topk_share = NA_real_) else list(topk_share = sum(head(count, min(top_k, .N))) / s)
  }, by = .(Tree, leaf_id)]
  
  if (nrow(perleaf) > progress_every) {
    message(sprintf("[per_leaf_topk_share] computed for %s leaves.", format(nrow(perleaf), big.mark=",")))
  }
  data.table::setkey(perleaf, Tree, leaf_id)
  perleaf[]
}

# Aggregates per-leaf feature concentration into per-row (SNP) averages across trees.
feature_concentration_per_snp <- function(leaves, perleaf_share, chunk_rows = 2000L, progress_every = 10000L) {
  message(sprintf("[feature_concentration_per_snp] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
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

#============================
# Harvest path rules (TRAIN)
#============================
# Extracts all candidate rules from training leaves, scoring them on extreme vs. background rows with statistics and medians.
harvest_path_rules <- function(
    leaves,
    leaf_steps,
    yvar_train        = NULL,  # optional if extreme_idx & bg_idx provided
    extreme_idx       = NULL,
    lower_tail        = FALSE,
    extreme_k         = 1,
    bg_idx            = NULL,
    bg_band_k         = NULL,
    band_center       = c("mean","median"),
    band_scale        = c("sd","mad"),
    max_depth         = 4L,
    min_support       = 20L,
    trees_subset      = NULL,   # 0-based tree ids
    pool_identical    = TRUE,
    trees_per_batch   = 250L,
    progress_every    = 1000L,
    tighten_monotone  = TRUE,
    return_ledger     = TRUE
) {
  message(sprintf("[harvest_path_rules] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  stopifnot(is.matrix(leaves) || inherits(leaves, "dgCMatrix"))
  
  band_center <- match.arg(band_center); band_scale <- match.arg(band_scale)
  
  # yvar is:
  # - optional if both extr/bg idx are supplied (used only for medians if provided)
  # - mandatory if extr/bg idx are not supplied
  if (is.null(extreme_idx) || is.null(bg_idx)) {
    if (is.null(yvar_train)) stop("[harvest_path_rules] yvar_train is required when extreme_idx/bg_idx are not provided.")
    .require_numeric_vec(yvar_train, "harvest_path_rules:yvar_train")
    if (length(yvar_train) != nrow(leaves)) stop("[harvest_path_rules] length(yvar_train) != nrow(leaves).")
    .require_finite_center_scale(yvar_train, "harvest_path_rules:yvar_train")
  } else if (!is.null(yvar_train)) {
    # yvar only used for medians; basic numeric/length check if present
    .require_numeric_vec(yvar_train, "harvest_path_rules:yvar_train")
    if (length(yvar_train) != nrow(leaves)) stop("[harvest_path_rules] length(yvar_train) != nrow(leaves).")
  }
  
  # Globals
  n  <- nrow(leaves); Tm <- ncol(leaves)
  
  # Canonical PATHS (drop unused cols early)
  PATHS_full <- build_path_rule_strings(leaf_steps, max_depth = max_depth, tighten_monotone = tighten_monotone)
  PATHS <- PATHS_full[, .(Tree, leaf_id, rule_str)]
  data.table::setkey(PATHS, Tree, leaf_id)
  
  # Derive/validate indices
  idx <- derive_yvar_idx(
    y           = if (is.null(extreme_idx) || is.null(bg_idx)) yvar_train else numeric(n), # not used if both idx given
    n           = n,
    lower_tail  = lower_tail,
    extreme_k   = extreme_k,
    extreme_idx = extreme_idx,
    bg_idx      = bg_idx,
    bg_band_k   = bg_band_k,
    band_center = band_center,
    band_scale  = band_scale,
    check_y_stats = (is.null(extreme_idx) || is.null(bg_idx) || !is.null(bg_band_k))
  )
  extr_idx  <- idx$extr_idx
  bg_idx    <- idx$bg_idx
  N_extr    <- idx$N_extr
  N_bg      <- idx$N_bg
  base_rate <- idx$base_rate
  
  # Tree subset → column indices
  all_trees0 <- 0:(Tm - 1L)
  if (!is.null(trees_subset)) {
    trees_subset <- as.integer(sort(unique(trees_subset)))
    trees_subset <- trees_subset[trees_subset %in% all_trees0]
    if (!length(trees_subset)) stop("[harvest] trees_subset empty after filtering")
    use_tt <- trees_subset + 1L
  } else {
    use_tt <- seq_len(Tm)
  }
  
  # Precompute nb per tree & rule_len by (Tree,leaf)
  nb_vec <- precompute_nb_vec(leaves)
  
  rule_len_by_leaf <- vector("list", Tm)
  split_paths <- split(PATHS_full[, .(leaf_id, rule_len)], PATHS_full$Tree)
  for (tr_chr in names(split_paths)) {
    tr <- as.integer(tr_chr)
    dt <- split_paths[[tr_chr]]
    if (!nrow(dt)) next
    nb <- max(dt$leaf_id, na.rm = TRUE) + 1L
    rlen <- rep.int(NA_integer_, nb)
    rlen[dt$leaf_id + 1L] <- as.integer(dt$rule_len)
    rule_len_by_leaf[[tr + 1L]] <- rlen
  }
  rm(PATHS_full)
  
  pool_beat <- if (!is.null(progress_every) && progress_every > 0L) progress_every * 5L else NULL
  y_num <- if (!is.null(yvar_train)) as.numeric(yvar_train) else NULL
  
  pairs_all <- NULL
  
  if (isTRUE(pool_identical)) {
    # ---- POOL: produce R ------------------------------------------------------
    rules_min_list <- vector("list", ceiling(length(use_tt) / trees_per_batch)); rmi <- 1L
    pairs_list     <- vector("list", ceiling(length(use_tt) / trees_per_batch)); pli <- 1L
    
    for (batch_start in seq(1L, length(use_tt), by = trees_per_batch)) {
      batch_end <- min(batch_start + trees_per_batch - 1L, length(use_tt))
      batch_rules_min <- vector("list", batch_end - batch_start + 1L); bi <- 1L
      batch_pairs     <- vector("list", batch_end - batch_start + 1L); bj <- 1L
      
      for (tt in use_tt[batch_start:batch_end]) {
        tr <- tt - 1L
        lt <- as.integer(leaves[, tt])
        nb <- nb_vec[tt]; if (!is.finite(nb) || nb <= 0L) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        count_extr_vec <- count_by_leaf(lt, extr_idx, nb)
        count_bg_vec   <- count_by_leaf(lt, bg_idx,   nb)
        
        leaf_ids <- which((count_extr_vec + count_bg_vec) > 0L) - 1L
        if (!length(leaf_ids)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        rlen <- rule_len_by_leaf[[tt]]
        if (is.null(rlen)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        DT0 <- data.table::data.table(
          Tree      = tr,
          leaf_id   = leaf_ids,
          rule_len  = rlen[leaf_ids + 1L],
          n_extreme = count_extr_vec[leaf_ids + 1L],
          n_bg      = count_bg_vec [leaf_ids + 1L]
        )
        
        DT0 <- DT0[(n_extreme + n_bg) >= as.integer(min_support)]
        if (!nrow(DT0)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        J <- PATHS[DT0, on = .(Tree, leaf_id), nomatch = 0L]
        if (!nrow(J)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        batch_rules_min[[bi]] <- J[, .(rule_str, rule_len)]
        batch_pairs    [[bj]] <- J[, .(Tree, leaf_id, rule_str)]
        bi <- bi + 1L; bj <- bj + 1L
      }
      
      rules_min_list[[rmi]] <- data.table::rbindlist(batch_rules_min, use.names = TRUE, fill = TRUE); rmi <- rmi + 1L
      pairs_list    [[pli]] <- data.table::rbindlist(batch_pairs,     use.names = TRUE, fill = TRUE); pli <- pli + 1L
      
      if (!is.null(progress_every) && progress_every > 0L &&
          ((batch_end %% progress_every) == 0L || batch_end == length(use_tt))) {
        message(sprintf("[harvest] processed tree %d / %d (%.1f%%)",
                        use_tt[batch_end], Tm, 100*batch_end/length(use_tt)))
      }
    }
    
    rules_min <- data.table::rbindlist(rules_min_list, use.names = TRUE, fill = TRUE)
    if (!nrow(rules_min)) {
      warning("[harvest] No rules after filtering; try lowering min_support or increasing max_depth.")
      return(data.table::data.table())
    }
    rule_len_map <- rules_min[, .(rule_len = max(rule_len, na.rm = TRUE)), by = rule_str]
    data.table::setkey(rule_len_map, rule_str)
    
    pairs_all <- data.table::rbindlist(pairs_list, use.names = TRUE, fill = TRUE)
    if (!nrow(pairs_all)) {
      warning("[harvest] No rule pairs to pool after filtering.")
      return(data.table::data.table())
    }
    
    # Row-index caches for pooling
    leaf_rows_index_ext <- vector("list", length = Tm)
    leaf_rows_index_bg  <- vector("list", length = Tm)
    for (tt in use_tt) {
      lt <- as.integer(leaves[, tt])
      leaf_rows_index_ext[[tt]] <- split(extr_idx, lt[extr_idx], drop = TRUE)
      leaf_rows_index_bg [[tt]] <- split(bg_idx,   lt[bg_idx],   drop = TRUE)
    }
    
    uniq_rules <- sort(unique(pairs_all$rule_str))
    pooled <- vector("list", length(uniq_rules))
    
    last_beat <- 0L
    for (i in seq_along(uniq_rules)) {
      rs <- uniq_rules[i]
      pairs <- unique(pairs_all[rule_str == rs, .(Tree, leaf_id)])
      if (!nrow(pairs)) next
      
      rows_ext <- integer(0); rows_bg <- integer(0)
      for (j in seq_len(nrow(pairs))) {
        tt  <- pairs$Tree[j] + 1L
        lid <- as.character(pairs$leaf_id[j])
        re  <- leaf_rows_index_ext[[tt]][[lid]]
        rb  <- leaf_rows_index_bg [[tt]][[lid]]
        if (!is.null(re)) rows_ext <- c(rows_ext, re)
        if (!is.null(rb)) rows_bg  <- c(rows_bg,  rb)
      }
      rows_ext <- unique(rows_ext)
      rows_bg  <- unique(rows_bg)
      
      n_e <- length(rows_ext)
      n_b <- length(rows_bg)
      support <- n_e + n_b
      
      precision <- n_e / pmax(1e-12, support)
      recall    <- if (N_extr > 0) n_e / N_extr else NA_real_
      lift      <- if (base_rate  > 0) precision / base_rate else NA_real_
      pval      <- stats::phyper(q = max(0L, n_e - 1L), m = N_extr, n = N_bg, k = support, lower.tail = FALSE)
      or_ha     <- ((n_e + 0.5) * (N_bg - n_b + 0.5)) / ((N_extr - n_e + 0.5) * (n_b + 0.5))
      
      med_e <- if (!is.null(y_num) && n_e) median(y_num[rows_ext], na.rm = TRUE) else NA_real_
      med_b <- if (!is.null(y_num) && n_b) median(y_num[rows_bg ], na.rm = TRUE) else NA_real_
      med_o <- if (!is.null(y_num) && support > 0L) median(c(y_num[rows_ext], y_num[rows_bg]), na.rm = TRUE) else NA_real_
      
      rl <- rule_len_map[rs, rule_len]
      
      pooled[[i]] <- data.table::data.table(
        rule_str      = rs,
        rule_len      = as.integer(rl),
        n_extreme     = n_e,
        n_bg          = n_b,
        support       = support,
        precision     = precision,
        recall        = recall,
        lift          = lift,
        odds_ratio    = or_ha,
        pval          = pval,
        med_y_extreme = med_e,
        med_y_bg      = med_b,
        med_y_overall = med_o
      )
      
      if (!is.null(pool_beat) && (i - last_beat) >= pool_beat) {
        message(sprintf("[harvest/pool] %d of %d rule strings", i, length(uniq_rules)))
        last_beat <- i
      }
    }
    
    R <- data.table::rbindlist(pooled, use.names = TRUE, fill = TRUE)
    
  } else {
    # ---- UNPOOLED: produce rules, then R <- rules -----------------------------
    out_list <- vector("list", ceiling(length(use_tt) / trees_per_batch)); oi <- 1L
    
    for (batch_start in seq(1L, length(use_tt), by = trees_per_batch)) {
      batch_end <- min(batch_start + trees_per_batch - 1L, length(use_tt))
      batch_rules <- vector("list", batch_end - batch_start + 1L); bi <- 1L
      
      for (tt in use_tt[batch_start:batch_end]) {
        tr <- tt - 1L
        lt <- as.integer(leaves[, tt])
        nb <- nb_vec[tt]; if (!is.finite(nb) || nb <= 0L) { bi <- bi + 1L; next }
        
        count_extr_vec <- count_by_leaf(lt, extr_idx, nb)
        count_bg_vec   <- count_by_leaf(lt, bg_idx,   nb)
        
        leaf_ids <- which((count_extr_vec + count_bg_vec) > 0L) - 1L
        if (!length(leaf_ids)) { bi <- bi + 1L; next }
        
        rlen <- rule_len_by_leaf[[tt]]
        if (is.null(rlen)) { bi <- bi + 1L; next }
        
        DT0 <- data.table::data.table(
          Tree      = tr,
          leaf_id   = leaf_ids,
          rule_len  = rlen[leaf_ids + 1L],
          n_extreme = count_extr_vec[leaf_ids + 1L],
          n_bg      = count_bg_vec [leaf_ids + 1L]
        )
        
        DT0 <- DT0[(n_extreme + n_bg) >= as.integer(min_support)]
        if (!nrow(DT0)) { bi <- bi + 1L; next }
        
        J <- PATHS[DT0, on = .(Tree, leaf_id), nomatch = 0L]
        if (!nrow(J)) { bi <- bi + 1L; next }
        
        support <- J$n_extreme + J$n_bg
        pval    <- stats::phyper(q = pmax(0L, J$n_extreme - 1L), m = N_extr, n = N_bg, k = support, lower.tail = FALSE)
        a <- J$n_extreme; b <- N_extr - J$n_extreme; c <- J$n_bg; d <- N_bg - J$n_bg
        or_ha <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        
        DT <- J[, .(
          rule_str,
          rule_len,
          n_extreme,
          n_bg,
          support,
          precision  = n_extreme / pmax(1e-12, n_extreme + n_bg),
          recall     = if (N_extr > 0) n_extreme / N_extr else NA_real_,
          lift       = if (base_rate  > 0) (n_extreme / pmax(1e-12, n_extreme + n_bg)) / base_rate else NA_real_,
          odds_ratio = or_ha,
          pval       = pval
        )]
        
        batch_rules[[bi]] <- DT; bi <- bi + 1L
      }
      
      out_list[[oi]] <- data.table::rbindlist(batch_rules, use.names = TRUE, fill = TRUE); oi <- oi + 1L
      if (!is.null(progress_every) && progress_every > 0L &&
          ((batch_end %% progress_every) == 0L || batch_end == length(use_tt))) {
        message(sprintf("[harvest] processed tree %d / %d (%.1f%%)",
                        use_tt[batch_end], Tm, 100*batch_end/length(use_tt)))
      }
    }
    
    rules <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
    if (!nrow(rules)) {
      warning("[harvest] No rules after filtering; try lowering min_support or increasing max_depth.")
      return(data.table::data.table())
    }
    R <- rules
  }
  
  # q-values & ordering on R (always)
  R[, qval := p.adjust(pval, method = "BH")]
  data.table::setorderv(R,
                        cols  = c("qval", "lift", "odds_ratio", "precision", "recall"),
                        order = c( 1L,   -1L,     -1L,          -1L,         -1L))
  if (isTRUE(pool_identical) && isTRUE(return_ledger)) {
    return(list(R = R[], pairs_all = pairs_all))
  } else {
    return(R[])
  }
}

#============================
# Validate rules on test set
#============================
# Applies candidate rules to a test set, recomputing support and metrics (identical schema to harvest).
validate_rules_on_test <- function(
    leaves_test,
    leaf_steps,
    yvar_test          = NULL,   # optional if extreme_idx & bg_idx provided
    lower_tail         = FALSE,
    extreme_k          = 1,
    extreme_idx        = NULL,
    bg_idx             = NULL,
    bg_band_k          = NULL,
    band_center        = c("mean","median"),
    band_scale         = c("sd","mad"),
    candidate_rules,
    max_depth          = 4L,
    min_support        = 1L,
    pool_identical     = TRUE,
    trees_per_batch    = 250L,
    progress_every     = 1000L,
    tighten_monotone   = TRUE,
    return_ledger     = TRUE
) {
  message(sprintf("[validate_rules_on_test] start: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  stopifnot(is.matrix(leaves_test) || inherits(leaves_test, "dgCMatrix"))
  stopifnot(is.data.frame(candidate_rules) || data.table::is.data.table(candidate_rules))
  if (!("rule_str" %in% names(candidate_rules))) stop("[validate_rules_on_test] candidate_rules must contain 'rule_str'")
  
  band_center <- match.arg(band_center); band_scale <- match.arg(band_scale)
  
  # yvar is:
  # - optional if both extr/bg idx are supplied (used only for medians if provided)
  # - mandatory if extr/bg idx are not supplied
  if (is.null(extreme_idx) || is.null(bg_idx)) {
    if (is.null(yvar_test)) stop("[validate_rules_on_test] yvar_test is required when extreme_idx/bg_idx are not provided.")
    .require_numeric_vec(yvar_test, "validate_rules_on_test:yvar_test")
    if (length(yvar_test) != nrow(leaves_test)) stop("[validate_rules_on_test] length(yvar_test) != nrow(leaves_test).")
    .require_finite_center_scale(yvar_test, "validate_rules_on_test:yvar_test")
  } else if (!is.null(yvar_test)) {
    .require_numeric_vec(yvar_test, "validate_rules_on_test:yvar_test")
    if (length(yvar_test) != nrow(leaves_test)) stop("[validate_rules_on_test] length(yvar_test) != nrow(leaves_test).")
  }
  
  # Globals / Sets
  n  <- nrow(leaves_test); Tm <- ncol(leaves_test)
  
  idx <- derive_yvar_idx(
    y           = if (is.null(extreme_idx) || is.null(bg_idx)) yvar_test else numeric(n),
    n           = n,
    lower_tail  = lower_tail,
    extreme_k   = extreme_k,
    extreme_idx = extreme_idx,
    bg_idx      = bg_idx,
    bg_band_k   = bg_band_k,
    band_center = band_center,
    band_scale  = band_scale,
    check_y_stats = (is.null(extreme_idx) || is.null(bg_idx) || !is.null(bg_band_k))
  )
  extr_idx  <- idx$extr_idx
  bg_idx    <- idx$bg_idx
  N_extr    <- idx$N_extr
  N_bg      <- idx$N_bg
  base_rate <- idx$base_rate
  
  # Canonical rule strings; filter to candidates; drop unused cols early
  PATHS_full <- build_path_rule_strings(leaf_steps, max_depth = max_depth, tighten_monotone = tighten_monotone)
  data.table::setDT(PATHS_full); data.table::setkey(PATHS_full, rule_str)
  CR <- data.table::as.data.table(unique(candidate_rules[, "rule_str", drop = FALSE]))
  data.table::setkey(CR, rule_str)
  PATHS_full <- PATHS_full[CR, nomatch = 0L]
  PATHS <- PATHS_full[, .(Tree, leaf_id, rule_str)]
  data.table::setkey(PATHS, Tree, leaf_id)
  
  # Precompute nb per tree & rule_len by (Tree,leaf) (only for pooling’s max)
  nb_vec <- precompute_nb_vec(leaves_test)
  
  rule_len_by_leaf <- vector("list", Tm)
  split_paths <- split(PATHS_full[, .(leaf_id, rule_len)], PATHS_full$Tree)
  for (tr_chr in names(split_paths)) {
    tr <- as.integer(tr_chr)
    dt <- split_paths[[tr_chr]]
    if (!nrow(dt)) next
    nb <- max(dt$leaf_id, na.rm = TRUE) + 1L
    rlen <- rep.int(NA_integer_, nb)
    rlen[dt$leaf_id + 1L] <- as.integer(dt$rule_len)
    rule_len_by_leaf[[tr + 1L]] <- rlen
  }
  rm(PATHS_full)
  
  y_num <- if (!is.null(yvar_test)) as.numeric(yvar_test) else NULL
  pool_beat <- if (!is.null(progress_every) && progress_every > 0L) progress_every * 5L else NULL
  
  pairs_all <- NULL
  
  if (isTRUE(pool_identical)) {
    # ---- POOL: produce R ------------------------------------------------------
    rules_min_list <- vector("list", ceiling(Tm / trees_per_batch)); rmi <- 1L
    pairs_list     <- vector("list", ceiling(Tm / trees_per_batch)); pli <- 1L
    
    for (batch_start in seq(1L, Tm, by = trees_per_batch)) {
      batch_end <- min(batch_start + trees_per_batch - 1L, Tm)
      batch_rules_min <- vector("list", batch_end - batch_start + 1L); bi <- 1L
      batch_pairs     <- vector("list", batch_end - batch_start + 1L); bj <- 1L
      
      for (tt in batch_start:batch_end) {
        tr <- tt - 1L
        lt <- as.integer(leaves_test[, tt])
        nb <- nb_vec[tt]; if (!is.finite(nb) || nb <= 0L) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        count_extr_vec <- count_by_leaf(lt, extr_idx, nb)
        count_bg_vec   <- count_by_leaf(lt, bg_idx,   nb)
        
        leaf_ids <- which((count_extr_vec + count_bg_vec) > 0L) - 1L
        if (!length(leaf_ids)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        rlen <- rule_len_by_leaf[[tt]]
        if (is.null(rlen)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        DT0 <- data.table::data.table(
          Tree      = tr,
          leaf_id   = leaf_ids,
          rule_len  = rlen[leaf_ids + 1L],
          n_extreme = count_extr_vec[leaf_ids + 1L],
          n_bg      = count_bg_vec [leaf_ids + 1L]
        )
        
        DT0 <- DT0[(n_extreme + n_bg) >= as.integer(min_support)]
        if (!nrow(DT0)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        J <- PATHS[DT0, on = .(Tree, leaf_id), nomatch = 0L]
        if (!nrow(J)) { bi <- bi + 1L; bj <- bj + 1L; next }
        
        batch_rules_min[[bi]] <- J[, .(rule_str, rule_len)]
        batch_pairs    [[bj]] <- J[, .(Tree, leaf_id, rule_str)]
        bi <- bi + 1L; bj <- bj + 1L
      }
      
      rules_min_list[[rmi]] <- data.table::rbindlist(batch_rules_min, use.names = TRUE, fill = TRUE); rmi <- rmi + 1L
      pairs_list    [[pli]] <- data.table::rbindlist(batch_pairs,     use.names = TRUE, fill = TRUE); pli <- pli + 1L
      
      if (!is.null(progress_every) && progress_every > 0L &&
          ((batch_end %% progress_every) == 0L || batch_end == Tm)) {
        message(sprintf("[validate] processed tree %d / %d (%.1f%%)", batch_end, Tm, 100*batch_end/Tm))
      }
    }
    
    rules_min <- data.table::rbindlist(rules_min_list, use.names = TRUE, fill = TRUE)
    if (!nrow(rules_min)) {
      warning("[validate] No rule support on the test set under the given settings.")
      return(data.table::data.table())
    }
    rule_len_map <- rules_min[, .(rule_len = max(rule_len, na.rm = TRUE)), by = rule_str]
    data.table::setkey(rule_len_map, rule_str)
    
    pairs_all <- data.table::rbindlist(pairs_list, use.names = TRUE, fill = TRUE)
    if (!nrow(pairs_all)) {
      warning("[validate] No rule pairs to pool on the test set.")
      return(data.table::data.table())
    }
    
    # Row-index caches on TEST for pooling
    leaf_rows_index_ext <- vector("list", length = Tm)
    leaf_rows_index_bg  <- vector("list", length = Tm)
    for (tt in seq_len(Tm)) {
      lt <- as.integer(leaves_test[, tt])
      leaf_rows_index_ext[[tt]] <- split(extr_idx, lt[extr_idx], drop = TRUE)
      leaf_rows_index_bg [[tt]] <- split(bg_idx,   lt[bg_idx],   drop = TRUE)
    }
    
    uniq_rules <- sort(unique(pairs_all$rule_str))
    pooled <- vector("list", length(uniq_rules))
    
    last_beat <- 0L
    for (i in seq_along(uniq_rules)) {
      rs <- uniq_rules[i]
      pairs <- unique(pairs_all[rule_str == rs, .(Tree, leaf_id)])
      if (!nrow(pairs)) next
      
      rows_ext <- integer(0); rows_bg <- integer(0)
      for (j in seq_len(nrow(pairs))) {
        tt  <- pairs$Tree[j] + 1L
        lid <- as.character(pairs$leaf_id[j])
        re  <- leaf_rows_index_ext[[tt]][[lid]]
        rb  <- leaf_rows_index_bg [[tt]][[lid]]
        if (!is.null(re)) rows_ext <- c(rows_ext, re)
        if (!is.null(rb)) rows_bg  <- c(rows_bg,  rb)
      }
      rows_ext <- unique(rows_ext)
      rows_bg  <- unique(rows_bg)
      
      n_e <- length(rows_ext)
      n_b <- length(rows_bg)
      support <- n_e + n_b
      
      precision <- n_e / pmax(1e-12, support)
      recall    <- if (N_extr > 0) n_e / N_extr else NA_real_
      lift      <- if (base_rate  > 0) precision / base_rate else NA_real_
      pval      <- stats::phyper(q = max(0L, n_e - 1L), m = N_extr, n = N_bg, k = support, lower.tail = FALSE)
      or_ha     <- ((n_e + 0.5) * (N_bg - n_b + 0.5)) / ((N_extr - n_e + 0.5) * (n_b + 0.5))
      
      rl <- rule_len_map[rs, rule_len]
      
      med_e <- if (!is.null(y_num <- y_num <- if (!is.null(yvar_test)) as.numeric(yvar_test) else NULL) && n_e) median(yvar_test[rows_ext], na.rm = TRUE) else NA_real_
      med_b <- if (!is.null(yvar_test) && n_b) median(yvar_test[rows_bg ], na.rm = TRUE) else NA_real_
      med_o <- if (!is.null(yvar_test) && support > 0L) median(c(yvar_test[rows_ext], yvar_test[rows_bg]), na.rm = TRUE) else NA_real_
      
      pooled[[i]] <- data.table::data.table(
        rule_str      = rs,
        rule_len      = as.integer(rl),
        n_extreme     = n_e,
        n_bg          = n_b,
        support       = support,
        precision     = precision,
        recall        = recall,
        lift          = lift,
        odds_ratio    = or_ha,
        pval          = pval,
        med_y_extreme = med_e,
        med_y_bg      = med_b,
        med_y_overall = med_o
      )
      
      if (!is.null(pool_beat) && (i - last_beat) >= pool_beat) {
        message(sprintf("[validate/pool] %d of %d rule strings", i, length(uniq_rules)))
        last_beat <- i
      }
    }
    
    R <- data.table::rbindlist(pooled, use.names = TRUE, fill = TRUE)
    
  } else {
    # ---- UNPOOLED: produce rules, then R <- rules -----------------------------
    out_list <- vector("list", ceiling(Tm / trees_per_batch)); oi <- 1L
    
    for (batch_start in seq(1L, Tm, by = trees_per_batch)) {
      batch_end <- min(batch_start + trees_per_batch - 1L, Tm)
      batch_rules <- vector("list", batch_end - batch_start + 1L); bi <- 1L
      
      for (tt in batch_start:batch_end) {
        tr <- tt - 1L
        lt <- as.integer(leaves_test[, tt])
        nb <- nb_vec[tt]; if (!is.finite(nb) || nb <= 0L) { bi <- bi + 1L; next }
        
        count_extr_vec <- count_by_leaf(lt, extr_idx, nb)
        count_bg_vec   <- count_by_leaf(lt, bg_idx,   nb)
        
        leaf_ids <- which((count_extr_vec + count_bg_vec) > 0L) - 1L
        if (!length(leaf_ids)) { bi <- bi + 1L; next }
        
        rlen <- rule_len_by_leaf[[tt]]
        if (is.null(rlen)) { bi <- bi + 1L; next }
        
        DT0 <- data.table::data.table(
          Tree      = tr,
          leaf_id   = leaf_ids,
          rule_len  = rlen[leaf_ids + 1L],
          n_extreme = count_extr_vec[leaf_ids + 1L],
          n_bg      = count_bg_vec [leaf_ids + 1L]
        )
        
        DT0 <- DT0[(n_extreme + n_bg) >= as.integer(min_support)]
        if (!nrow(DT0)) { bi <- bi + 1L; next }
        
        J <- PATHS[DT0, on = .(Tree, leaf_id), nomatch = 0L]
        if (!nrow(J)) { bi <- bi + 1L; next }
        
        support <- J$n_extreme + J$n_bg
        pval    <- stats::phyper(q = pmax(0L, J$n_extreme - 1L), m = N_extr, n = N_bg, k = support, lower.tail = FALSE)
        a <- J$n_extreme; b <- N_extr - J$n_extreme; c <- J$n_bg; d <- N_bg - J$n_bg
        or_ha <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        
        DT <- J[, .(
          rule_str,
          rule_len,
          n_extreme,
          n_bg,
          support,
          precision  = n_extreme / pmax(1e-12, n_extreme + n_bg),
          recall     = if (N_extr > 0) n_extreme / N_extr else NA_real_,
          lift       = if (base_rate  > 0) (n_extreme / pmax(1e-12, n_extreme + n_bg)) / base_rate else NA_real_,
          odds_ratio = or_ha,
          pval       = pval
        )]
        
        batch_rules[[bi]] <- DT; bi <- bi + 1L
      }
      
      out_list[[oi]] <- data.table::rbindlist(batch_rules, use.names = TRUE, fill = TRUE); oi <- oi + 1L
      if (!is.null(progress_every) && progress_every > 0L &&
          ((batch_end %% progress_every) == 0L || batch_end == Tm)) {
        message(sprintf("[validate] processed tree %d / %d (%.1f%%)", batch_end, Tm, 100*batch_end/Tm))
      }
    }
    
    rules <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
    if (!nrow(rules)) {
      warning("[validate] No rule support on the test set under the given settings.")
      return(data.table::data.table())
    }
    R <- rules
  }
  
  # q-values & ordering on R (always)
  R[, qval := p.adjust(pval, method = "BH")]
  data.table::setorderv(R,
                        cols  = c("qval", "lift", "odds_ratio", "precision", "recall"),
                        order = c( 1L,   -1L,     -1L,          -1L,         -1L))
  if (isTRUE(pool_identical) && isTRUE(return_ledger)) {
    return(list(R = R[], pairs_all = pairs_all))
  } else {
    return(R[])
  }
}
