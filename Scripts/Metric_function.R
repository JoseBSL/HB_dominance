#Function to compute metrics

compute_metrics = function(df,
                           focal_species = "Apis mellifera",
                           return_observed = TRUE,
                           eps = 1e-12) {
  
  ## ---------- helpers (mirror bipartite::dfun for d′ finite d_min) ----------
  .d_row_kl <- function(x, q) {
    s <- sum(x); if (s <= 0) return(0)
    p <- x / s
    nz <- p > 0
    sum(p[nz] * log(p[nz] / q[nz]))
  }
  
  .dmin_bip <- function(x, q, cs = NULL) {
    n <- sum(x)
    if (n <= 0) return(0)
    expec <- floor(q * n)
    rest  <- n - sum(expec)
    xnew  <- expec
    for (j in seq_len(rest)) {
      i.vec <- if (is.null(cs)) seq_along(xnew) else which(xnew < cs)
      xsum <- sum(xnew)
      p1 <- xnew / (xsum + 1)
      dstep1 <- ifelse(xnew != 0, p1 * log(p1 / q), 0)
      dcheck <- rep(Inf, length(xnew))
      for (i in i.vec) {
        pi1 <- (xnew[i] + 1) / (xsum + 1)
        dcheck[i] <- pi1 * log(pi1 / q[i]) + sum(dstep1[-i])
      }
      i.best <- which.min(dcheck)[1]
      xnew[i.best] <- xnew[i.best] + 1L
    }
    .d_row_kl(xnew, q)
  }
  
  ## ---------- basic checks & matrix ----------
  req <- c("pollinator", "plant", "interaction")
  stopifnot(all(req %in% names(df)))
  df <- as.data.frame(df)
  df <- df[stats::complete.cases(df[req]), req, drop = FALSE]
  
  if (nrow(df) == 0) return(NULL)
  
  W <- xtabs(interaction ~ pollinator + plant, df)  # fills unobserved pairs with 0
  if (!(focal_species %in% rownames(W))) return(NULL)
  
  keep_r <- rowSums(W) > 0
  keep_c <- colSums(W) > 0
  if (!keep_r[focal_species]) return(NULL)
  W <- W[keep_r, keep_c, drop = FALSE]
  
  others <- setdiff(rownames(W), focal_species)
  if (!length(others)) return(NULL)
  
  out <- data.frame(pollinator = others, stringsAsFactors = FALSE)
  
  ## ------------ normalized degree (observed adjacency) ------------
  A <- W > 0
  deg <- rowSums(A)
  n_plants <- ncol(W)
  out$degree         <- as.numeric(deg[others])
  out$norm_degree    <- as.numeric(deg[others] / max(n_plants, 1L))
  out$n_plants_total <- n_plants
  
  ## ------------ PDI (Poisot) standardized to 0–1 ------------
  if (n_plants <= 1) {
    out$pdi <- NA_real_
  } else {
    pdi_all <- apply(W, 1, function(v) {
      if (all(v <= 0)) return(NA_real_)
      vmax <- max(v)
      if (vmax <= 0) return(NA_real_)
      vnorm <- v / vmax                      # scale by strongest link
      imax  <- which.max(v)[1]               # one instance of the max
      rest  <- vnorm[-imax]
      (sum(1 - rest)) / (length(v) - 1)      # = average (1 - v_rest)
    })
    out$pdi <- as.numeric(pdi_all[others])
  }
  
  ## ------------ observed overlaps (resource use) ------------
  if (return_observed) {
    f      <- as.numeric(W[focal_species, ])
    f_sum  <- sum(f)
    f_sq   <- sum(f^2)
    f_d    <- f_sq / (f_sum^2)
    
    O      <- W[others, , drop = FALSE]
    O_sum  <- rowSums(O)
    O_sq   <- rowSums(O^2)
    cross  <- as.numeric(O %*% f)
    
    out$pianka <- cross / sqrt(pmax(f_sq * O_sq, eps))
    out$morisita_horn <- (2 * cross) /
      ((f_d + (O_sq / pmax(O_sum^2, eps))) * f_sum * pmax(O_sum, eps))
  }
  
  ## ------------ preference overlaps (closed-form E = r*c/N) ------------
  r  <- rowSums(W); c <- colSums(W); N <- sum(W)
  E      <- (r %o% c) / max(N, eps)
  Pref   <- W / pmax(E, eps)
  Pref[E == 0] <- 0
  rs     <- rowSums(Pref)
  Pref   <- Pref / pmax(rs, eps)   # row-normalize
  
  Fp     <- as.numeric(Pref[focal_species, ])
  Op     <- Pref[others, , drop = FALSE]
  numer  <- as.numeric(Op %*% Fp)
  Sf     <- sum(Fp^2)
  So     <- rowSums(Op^2)
  
  out$pref_pianka        <- numer / sqrt(pmax(Sf * So, eps))
  out$pref_morisita_horn <- (2 * numer) / pmax(Sf + So, eps)
  
  ## ------------ Blüthgen's d' (standard, all species in q; optional keep) ------------
  P_use <- W / pmax(rowSums(W), eps)
  Q_all <- colSums(W); Q_all <- Q_all / pmax(sum(Q_all), eps); Q_all[Q_all < eps] <- eps
  log_ratio_all <- log(P_use / matrix(Q_all, nrow(P_use), ncol(P_use), byrow = TRUE))
  log_ratio_all[!is.finite(log_ratio_all)] <- 0
  D_all <- rowSums(P_use * log_ratio_all)
  dmax_all <- if (any(Q_all > 0)) log(1 / min(Q_all[Q_all > 0])) else NA_real_
  Dprime_all <- if (is.finite(dmax_all) && dmax_all > 0) pmin(pmax(D_all / dmax_all, 0), 1) else NA_real_
  out$d_prime <- as.numeric(Dprime_all[others])
  
  ## ------------ d' with Apis excluded from availability (bipartite-consistent) ------------
  nonapis_rows <- setdiff(rownames(W), focal_species)
  Wq <- if (length(nonapis_rows)) W[nonapis_rows, , drop = FALSE] else W * 0
  q  <- colSums(Wq)
  keep_cols <- q > 0
  if (!any(keep_cols)) {
    out$d_prime_no_apis_finite <- NA_real_
  } else {
    q  <- q[keep_cols]; q <- q / sum(q)
    cs <- colSums(Wq[, keep_cols, drop = FALSE])   # column caps
    m  <- sum(Wq)                                  # total (non-Apis) interactions
    X  <- W[, keep_cols, drop = FALSE]
    
    Ai    <- rowSums(X)
    draw  <- apply(X, 1, .d_row_kl, q = q)
    dmin  <- vapply(seq_len(nrow(X)), function(i) .dmin_bip(X[i, ], q, cs = cs), numeric(1))
    dmax  <- log(pmax(m / pmax(Ai, 1), 1))        # bipartite::dfun style
    
    dprime_excl <- (draw - dmin) / pmax(dmax - dmin, eps)
    dprime_excl <- pmin(pmax(dprime_excl, 0), 1)
    out$d_prime_no_apis_finite <- as.numeric(dprime_excl[others])
  }
  
  ## ------------ Müller's PAC: exposure of i to Apis (directional) ------------
  # s_apis(k) = alpha_{Apis,k} / sum_m alpha_{mk}
  s_apis <- as.numeric(W[focal_species, ]) / pmax(colSums(W), eps)
  
  # p_i(k) = alpha_{i,k} / sum_l alpha_{i,l}; PAC(i<-Apis) = sum_k p_i(k) * s_apis(k)
  O      <- W[others, , drop = FALSE]
  O_rs   <- pmax(rowSums(O), eps)
  P_i    <- O / O_rs
  out$pac_from_apis <- as.numeric(P_i %*% s_apis)
  
  rownames(out) <- NULL
  out
}

###
