#Function to compute metrics

compute_metrics = function(df,
                           focal_species = "Apis mellifera",
                           return_observed = TRUE,
                           eps = 1e-12) {

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
  
  out$pref_morisita_horn <- (2 * numer) / pmax(Sf + So, eps)
  
  rownames(out) <- NULL
  out
}

###

