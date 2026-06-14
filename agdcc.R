# agdcc.R
#
# Three-stage QMLE for the DCC/ADCC/G-DCC/AG-DCC family
# (Cappiello, Engle & Sheppard 2006; Engle & Sheppard 2001 NBER WP 8554).
#
# Two VCV options via vcv_method:
#   "none"   — no standard errors
#   "3stage" — Full 3-block sandwich A0^{-1} B0 A0^{-1}' / T (Engle & Sheppard 2001)

.build_backend <- function(cpp_file = "agdcc_core.cpp") {
  if (!exists("qt_recursion_cpp", envir = .GlobalEnv)) {
    if (!file.exists(cpp_file)) stop("C++ backend '", cpp_file, "' not found.")
    Rcpp::sourceCpp(cpp_file)
    message("agdcc_core.cpp compiled.")
  }
}

.load_deps <- function() {
  for (p in c("Rcpp","RcppArmadillo","rugarch","numDeriv","parallel","MASS")) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org")
      library(p, character.only = TRUE)
    }
  }
}

DCC_SCALAR  <- 1L; ADCC_SCALAR <- 2L
GDCC_DIAG   <- 3L; AGDCC_DIAG  <- 4L

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

cov2cor_safe <- function(S) {
  d <- sqrt(pmax(diag(S), .Machine$double.eps))
  R <- S / outer(d, d); diag(R) <- 1
  list(R = R, sds = d)
}

is_pd <- function(M) tryCatch({ chol(M); TRUE }, error = function(e) FALSE)

.sym <- function(M) (M + t(M)) / 2   # force numerical symmetry after Jacobian accumulation

# Ridge and pinv fallbacks: A22 can be near-singular near the constraint boundary.
.safe_solve <- function(A, ridge = 1e-8) {
  A <- .sym(A)
  tryCatch(solve(A), error = function(e)
    tryCatch(solve(A + diag(max(abs(diag(A))) * ridge, nrow(A))),
             error = function(e2) MASS::ginv(A)))
}

# Newey-West HAC covariance for a raw T x p score matrix.
covnw <- function(X, L = NULL) {
  T  <- nrow(X)
  if (is.null(L)) L <- min(floor(1.2 * T^(1/3)), T)
  Xc <- scale(X, scale = FALSE)
  B  <- crossprod(Xc) / T
  for (h in seq_len(L)) {
    w  <- 1 - h / (L + 1)
    Xh <- crossprod(Xc[(h+1):T, , drop=FALSE],
                    Xc[1:(T-h),  , drop=FALSE]) / T
    B  <- B + w * (Xh + t(Xh))
  }
  B
}

# Pair index helpers for intercept parameter ordering (column-major, lower triangle).
.pairs_lower <- function(k) {
  do.call(rbind, lapply(seq_len(k-1), function(j) cbind(row=(j+1):k, col=j)))
}
.pairs_lower_diag <- function(k) {
  do.call(rbind, lapply(seq_len(k), function(j) cbind(row=j:k, col=j)))
}

# ---------------------------------------------------------------------------
# Stage 1
# ---------------------------------------------------------------------------

fit_stage1 <- function(data, stage1_fit = NULL, n_cores = 1L, return_fits = FALSE) {
  k <- ncol(data); T <- nrow(data)

  extract_one <- function(fit_obj, ret_i, spec_obj = NULL) {
    ht  <- as.numeric(rugarch::sigma(fit_obj))^2
    res <- list(h=ht, e=ret_i/sqrt(ht))
    if (return_fits) {
      res$spec  <- if (!is.null(spec_obj)) spec_obj else rugarch::getspec(fit_obj)
      res$theta <- coef(fit_obj)
    }
    res
  }

  if (inherits(stage1_fit, "uGARCHmultifit")) {
    # rugarch >= 1.4: individual fits stored in @fit slot as named list
    ind <- tryCatch(stage1_fit@fit,
                    error = function(e) lapply(seq_len(k), function(i) stage1_fit[[i]]))
    if (is.null(ind) || length(ind) < k)
      stop("Cannot extract individual fits from ugarchMultifit.")
    univar <- lapply(seq_len(k), function(i) extract_one(ind[[i]], data[,i]))

  } else if (is.list(stage1_fit) &&
             all(vapply(stage1_fit, inherits, logical(1), "uGARCHfit"))) {
    if (length(stage1_fit) != k)
      stop("stage1_fit length (", length(stage1_fit), ") != k (", k, ")")
    univar <- lapply(seq_len(k), function(i) extract_one(stage1_fit[[i]], data[,i]))

  } else {
    spec <- rugarch::ugarchspec(
      variance.model     = list(model="sGARCH", garchOrder=c(1,1)),
      mean.model         = list(armaOrder=c(0,0), include.mean=FALSE),
      distribution.model = "norm")
    fit_one <- function(i) {
      fo <- rugarch::ugarchfit(spec, data[,i], solver="solnp",
                               solver.control=list(trace=0))
      extract_one(fo, data[,i], spec_obj = spec)   # pass spec directly (avoids @model accessor)
    }
    if (n_cores > 1L && k > 1L) {
      cl <- parallel::makeCluster(min(n_cores, k))
      on.exit(parallel::stopCluster(cl), add=TRUE)
      parallel::clusterExport(cl, c("spec","data","extract_one","return_fits"),
                              envir=environment())
      parallel::clusterEvalQ(cl, { library(rugarch); Rcpp::sourceCpp("agdcc_core.cpp") })
      univar <- parallel::parLapply(cl, seq_len(k), fit_one)
    } else {
      univar <- lapply(seq_len(k), fit_one)
    }
  }

  H     <- do.call(cbind, lapply(univar, `[[`, "h"))
  e_std <- do.call(cbind, lapply(univar, `[[`, "e"))
  out   <- list(H=H, e_std=e_std)
  if (return_fits) {
    out$specs  <- lapply(univar, `[[`, "spec")
    out$thetas <- lapply(univar, `[[`, "theta")
  }
  out
}

# ---------------------------------------------------------------------------
# Stage 2
# ---------------------------------------------------------------------------

compute_intercepts <- function(e_std) {
  T <- nrow(e_std)
  n_std <- e_std * (e_std < 0)
  list(P=cov2cor_safe(crossprod(e_std)/T)$R,
       N=crossprod(n_std)/T, n_std=n_std)
}

build_outer_cubes <- function(e_std, n_std) {
  T <- nrow(e_std); k <- ncol(e_std)
  ee <- array(0, c(k,k,T)); nn <- array(0, c(k,k,T))
  for (t in seq_len(T)) { ee[,,t] <- tcrossprod(e_std[t,]); nn[,,t] <- tcrossprod(n_std[t,]) }
  list(ee=ee, nn=nn)
}

compute_delta <- function(P, N) {
  Psq <- tryCatch({ ch <- chol(P); backsolve(ch, diag(nrow(P))) },
                  error=function(e) diag(nrow(P)))
  max(eigen(Psq %*% N %*% Psq, symmetric=TRUE, only.values=TRUE)$values)
}

# ---------------------------------------------------------------------------
# Stage 3
# ---------------------------------------------------------------------------

.stat_constr <- function(psi, model, k, delta) {
  if (model==DCC_SCALAR)  return(psi[1]^2 + psi[2]^2)
  if (model==ADCC_SCALAR) return(psi[1]^2 + delta*psi[2]^2 + psi[3]^2)
  if (model==GDCC_DIAG)   return(max(psi[1:k]^2 + psi[(k+1):(2*k)]^2))
  max(psi[1:k]^2 + delta*psi[(k+1):(2*k)]^2 + psi[(2*k+1):(3*k)]^2)
}

# Gradient of stat_constr w.r.t. psi (for penalty gradient in analytical gr)
.stat_grad_raw <- function(psi, model, k, delta) {
  np <- length(psi); g <- numeric(np)
  if (model==DCC_SCALAR) {
    g[1] <- 2*psi[1]; g[2] <- 2*psi[2]
  } else if (model==ADCC_SCALAR) {
    g[1] <- 2*psi[1]; g[2] <- 2*delta*psi[2]; g[3] <- 2*psi[3]
  } else if (model==GDCC_DIAG) {
    vals <- psi[1:k]^2 + psi[(k+1):(2*k)]^2; idx <- which.max(vals)
    g[idx] <- 2*psi[idx]; g[k+idx] <- 2*psi[k+idx]
  } else {
    vals <- psi[1:k]^2 + delta*psi[(k+1):(2*k)]^2 + psi[(2*k+1):(3*k)]^2
    idx  <- which.max(vals)
    g[idx] <- 2*psi[idx]; g[k+idx] <- 2*delta*psi[k+idx]; g[2*k+idx] <- 2*psi[2*k+idx]
  }
  g
}

stage3_obj <- function(psi, P, N, e_std, ee, nn, bc, bca, model, delta, k) {
  ll    <- stage3_qll_fast_cpp(psi, P, N, e_std, ee, nn, bc, bca, model)
  cv    <- .stat_constr(psi, model, k, delta)
  pen   <- if (cv >= 0.9999) 1e6*(cv-0.9999)^2 else 0
  ll + pen
}

stage3_grid_search <- function(P, N, e_std, ee, nn, bc, bca,
                                model, delta, k, n_cores=1L) {
  if (model==DCC_SCALAR) {
    grid <- as.matrix(expand.grid(a=c(0.01,0.05,0.10), b=c(0.85,0.90,0.95)))
  } else if (model==ADCC_SCALAR) {
    grid <- as.matrix(expand.grid(a=c(0.01,0.05), g=c(0.01,0.03), b=c(0.85,0.92)))
  } else if (model==GDCC_DIAG) {
    # Multiple scalar candidates replicated to k-dim: better coverage than a single point
    sv <- expand.grid(a=c(0.02,0.05,0.10), b=c(0.85,0.90,0.95))
    grid <- do.call(rbind, lapply(seq_len(nrow(sv)), function(i)
      c(rep(sv$a[i],k), rep(sv$b[i],k))))
  } else {
    sv <- expand.grid(a=c(0.02,0.05,0.10), g=c(0.01,0.03), b=c(0.85,0.90,0.95))
    grid <- do.call(rbind, lapply(seq_len(nrow(sv)), function(i)
      c(rep(sv$a[i],k), rep(sv$g[i],k), rep(sv$b[i],k))))
  }

  ok   <- apply(grid, 1, function(r) .stat_constr(r,model,k,delta) < 0.9999)
  grid <- grid[ok,,drop=FALSE]
  if (nrow(grid)==0)
    grid <- matrix(switch(as.character(model),
      "1"=c(0.05,0.90),"2"=c(0.05,0.02,0.90),
      "3"=c(rep(0.05,k),rep(0.90,k)),"4"=c(rep(0.05,k),rep(0.02,k),rep(0.90,k))), nrow=1)

  lls <- stage3_qll_grid_cpp(grid,P,N,e_std,ee,nn,bc,bca,model,as.integer(n_cores))
  as.numeric(grid[which.min(lls),])
}

optimise_stage3 <- function(P, N, e_std, ee, nn, bc, bca,
                             model, delta, k, n_cores=1L,
                             tol=1e-8, maxit=500L) {
  sv <- stage3_grid_search(P,N,e_std,ee,nn,bc,bca,model,delta,k,n_cores)
  np <- length(sv)
  T  <- nrow(e_std)
  LB <- rep(1e-6, np); UB <- rep(0.9999, np)
  it1 <- min(max(maxit, np * 20L), 2000L)
  it2 <- min(max(maxit * 2L, np * 40L), 4000L)
  obj <- function(p) stage3_obj(p,P,N,e_std,ee,nn,bc,bca,model,delta,k)
  # Analytical gradient: qt_and_mean_scores_cpp returns (1/T)*Σ dℓ_t/dψ
  # (mean log-lik score). The objective is Σ(-ℓ_t)+penalty, so:
  #   d(obj)/dψ = -T * mean_score + d(penalty)/dψ
  # n_cores parallelises the inner j-loop over parameters (score per param);
  # for large np (diagonal models) this gives a ~n_cores-fold speedup per
  # gradient call. The L-BFGS-B loop itself is single-threaded.
  gr_fn <- function(p) {
    ms <- qt_and_mean_scores_cpp(p, P, N, e_std, ee, nn, bc, bca, model, as.integer(n_cores))
    cv <- .stat_constr(p, model, k, delta)
    pg <- if (cv >= 0.9999) 2e6*(cv-0.9999)*(.stat_grad_raw(p,model,k,delta)) else rep(0,np)
    -T * ms + pg
  }
  bfgs_ctrl <- list(maxit=it1, factr=tol/.Machine$double.eps)
  f1 <- optim(sv, obj, gr=gr_fn, method="L-BFGS-B", lower=LB, upper=UB, control=bfgs_ctrl)

  # Nelder-Mead is useful for small-np scalar models to escape local optima
  # but degrades badly for diagonal models (np grows with k, simplex needs
  # k+1 vertices in np-dimensional space). Skip it for np > 5.
  if (np <= 5L) {
    f2  <- optim(f1$par, obj, method="Nelder-Mead",
                 control=list(maxit=it2, reltol=tol))
    sv3 <- if (f2$value < f1$value) f2$par else f1$par
  } else {
    f2  <- f1; sv3 <- f1$par
  }

  f3 <- optim(sv3, obj, gr=gr_fn, method="L-BFGS-B", lower=LB, upper=UB, control=bfgs_ctrl)
  f  <- list(f1, f2, f3)[[which.min(c(f1$value, f2$value, f3$value))]]
  list(psi_D=f$par, ll=-f$value, convergence=f$convergence)
}

# Analytic per-obs scores for the Stage-2 intercept pseudo-objective.
# Returns T x n2 matrix: R̄ block (K(K-1)/2 cols) + optional N̄ block (K(K+1)/2 cols).
intercept_perobs_scores <- function(e_std, n_std, P, N, model) {
  T <- nrow(e_std); k <- ncol(e_std)
  scales <- colMeans(e_std^2)
  pairs  <- .pairs_lower(k)
  n2_sym <- nrow(pairs)
  S2_R   <- matrix(0, T, n2_sym)
  for (p in seq_len(n2_sym)) {
    i <- pairs[p, 1]; j <- pairs[p, 2]
    S2_R[, p] <- e_std[, i] * e_std[, j] / sqrt(scales[i] * scales[j]) - P[i, j]
  }
  if (model %in% c(ADCC_SCALAR, AGDCC_DIAG)) {
    pairs_n <- .pairs_lower_diag(k)
    S2_N    <- matrix(0, T, nrow(pairs_n))
    for (p in seq_len(nrow(pairs_n))) {
      i <- pairs_n[p, 1]; j <- pairs_n[p, 2]
      S2_N[, p] <- n_std[, i] * n_std[, j] - N[i, j]
    }
    cbind(S2_R, S2_N)
  } else {
    S2_R
  }
}

# ---------------------------------------------------------------------------
# VCV: Full 3-block cross-stage sandwich (Engle & Sheppard 2001 Theorem 2 in full)
#
# Three parameter blocks:
#   Block 1 (φ): GARCH parameters  — n1 = Σ nᵢ
#   Block 2 (z): Intercepts P̄,N̄   — n2 = K(K-1)/2 [+ K(K+1)/2 for ADCC/AGDCC]
#   Block 3 (ψ): DCC dynamics       — n3 = length(psi_D)
#
# A₀ = [ A11    0     0  ]     lower-triangular (n_tot × n_tot)
#      [ A21   A22    0  ]     A22 = I_{n2} (intercept pseudo-Hessian is identity)
#      [ A31   A32   A33 ]
#
# B₀ = covnw([S1 | S2 | S3])  — NW-HAC, matching Sheppard dcc.m 3-stage convention
# VCV(ψ̂) = bottom-right n3×n3 block of A₀⁻¹ B₀ (A₀⁻¹)' / T
# ---------------------------------------------------------------------------

# Fused computation of S1 (GARCH per-obs scores), A11 (GARCH info matrix),
# A21 (intercept x GARCH cross-Jacobian), A31 (DCC x GARCH cross-Jacobian).
# Shares each ugarchfilter call across all three outputs, cutting filter calls
# from 3*n1 to n1 + k (k baseline calls for S1, n1 perturbed calls shared).
# Parallel over assets via parLapply when n_cores > 1.
compute_garch_blocks <- function(e_std, P, N, psi_D, model, n_cores,
                                 specs, thetas, data, s3_base, h = 1e-5) {
  k       <- ncol(data); T <- nrow(data); n3 <- length(psi_D)
  n1      <- sum(vapply(thetas, length, integer(1)))
  is_asym <- model %in% c(ADCC_SCALAR, AGDCC_DIAG)
  pairs   <- .pairs_lower(k)
  n2_sym  <- nrow(pairs)
  n2      <- n2_sym + if (is_asym) k*(k+1L)/2L else 0L
  pairs_n <- if (is_asym) .pairs_lower_diag(k) else NULL
  scales  <- colMeans(e_std^2)
  col_end   <- cumsum(vapply(thetas, length, integer(1)))
  col_start <- c(1L, col_end[-k] + 1L)

  # Per-asset worker: returns S1i (T x ni), A21i (n2 x ni), A31i (n3 x ni).
  # Free variables accessed via lexical scoping (serial) or clusterExport (parallel).
  do_asset <- function(task) {
    i  <- task$i
    ni <- length(task$theta_i)
    S1i  <- matrix(0, T, ni)
    A21i <- matrix(0, n2, ni)
    A31i <- matrix(0, n3, ni)

    # Baseline LL for S1 forward difference
    spec_f <- task$spec_i; spec_f@model$fixed.pars <- as.list(task$theta_i)
    filt0  <- rugarch::ugarchfilter(spec_f, task$data_i)
    ll0    <- filt0@filter$log.likelihoods; ll0[is.na(ll0)] <- 0

    r_pairs_i <- which(pairs[, 1] == i | pairs[, 2] == i)
    n_pairs_i <- if (is_asym) which(pairs_n[, 1] == i | pairs_n[, 2] == i) else integer(0)

    for (l in seq_len(ni)) {
      th_p   <- task$theta_i; th_p[l] <- th_p[l] + h
      sp_p   <- task$spec_i; sp_p@model$fixed.pars <- as.list(th_p)
      filt_p <- rugarch::ugarchfilter(sp_p, task$data_i)
      e_p    <- task$data_i / sqrt(pmax(as.numeric(rugarch::sigma(filt_p))^2,
                                        .Machine$double.eps))
      ll_p   <- filt_p@filter$log.likelihoods; ll_p[is.na(ll_p)] <- 0

      # S1 column
      S1i[, l] <- (ll_p - ll0) / h

      # A21: R̄ block (pairs involving asset i; m_base = 0 at variance-targeting solution)
      scale_ip <- mean(e_p^2)
      for (p in r_pairs_i) {
        ii <- pairs[p, 1]; jj <- pairs[p, 2]
        e_ii <- if (ii == i) e_p else e_std[, ii]
        e_jj <- if (jj == i) e_p else e_std[, jj]
        sc_i <- if (ii == i) scale_ip else scales[ii]
        sc_j <- if (jj == i) scale_ip else scales[jj]
        A21i[p, l] <- -mean(e_ii * e_jj) / sqrt(sc_i * sc_j) / h
      }
      # A21: N̄ block (ADCC/AGDCC)
      if (is_asym) {
        n_p <- e_p * (e_p < 0)
        for (p in n_pairs_i) {
          ii   <- pairs_n[p, 1]; jj <- pairs_n[p, 2]
          n_ii <- if (ii == i) n_p else e_std[, ii] * (e_std[, ii] < 0)
          n_jj <- if (jj == i) n_p else e_std[, jj] * (e_std[, jj] < 0)
          A21i[n2_sym + p, l] <- -mean(n_ii * n_jj) / h
        }
      }

      # A31: DCC mean-score Jacobian (Qt recursion on perturbed e_std[:,i])
      ep   <- e_std; ep[, i] <- e_p
      np_  <- ep * (ep < 0)
      cb_p <- build_outer_cubes_cpp(ep, np_)
      bc_p <- ewma_backcast_cpp(cb_p$ee); bca_p <- ewma_backcast_cpp(cb_p$nn)
      s3_p <- qt_and_mean_scores_cpp(psi_D, P, N, ep, cb_p$ee, cb_p$nn,
                                     bc_p, bca_p, model, 1L)
      A31i[, l] <- -(s3_p - s3_base) / h
    }
    list(S1i = S1i, A21i = A21i, A31i = A31i)
  }

  tasks <- lapply(seq_len(k), function(i)
    list(i = i, spec_i = specs[[i]], theta_i = thetas[[i]], data_i = data[, i]))

  if (n_cores > 1L && k > 1L) {
    cl <- parallel::makeCluster(min(n_cores, k))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    env_vars <- c("do_asset", "T", "n3", "n2", "n2_sym", "is_asym",
                  "pairs", "scales", "h", "psi_D", "P", "N", "model",
                  "e_std", "s3_base", "ADCC_SCALAR", "AGDCC_DIAG")
    if (is_asym) env_vars <- c(env_vars, "pairs_n")
    parallel::clusterExport(cl, env_vars, envir = environment())
    parallel::clusterEvalQ(cl, { library(rugarch); Rcpp::sourceCpp("agdcc_core.cpp") })
    res_list <- parallel::parLapply(cl, tasks, do_asset)
  } else {
    res_list <- lapply(tasks, do_asset)
  }

  S1  <- matrix(0, T, n1)
  A21 <- matrix(0, n2, n1)
  A31 <- matrix(0, n3, n1)
  for (i in seq_len(k)) {
    cols        <- col_start[i]:col_end[i]
    S1[, cols]  <- res_list[[i]]$S1i
    A21[, cols] <- res_list[[i]]$A21i
    A31[, cols] <- res_list[[i]]$A31i
  }
  S1[is.na(S1)] <- 0

  A11 <- matrix(0, n1, n1)
  for (i in seq_len(k)) {
    idx <- col_start[i]:col_end[i]
    A11[idx, idx] <- crossprod(S1[, idx]) / T
  }
  list(S1 = S1, A11 = A11, A21 = A21, A31 = A31)
}

vcv_3stage <- function(psi_D, P, N, e_std, ee, nn, bc, bca, model, n_cores,
                       specs, thetas, data) {
  T   <- nrow(e_std); k <- ncol(e_std)
  n3  <- length(psi_D)
  n1  <- sum(vapply(thetas, length, integer(1)))
  n_std   <- e_std * (e_std < 0)
  is_asym <- model %in% c(ADCC_SCALAR, AGDCC_DIAG)
  n2      <- nrow(.pairs_lower(k)) + if (is_asym) k*(k+1L)/2L else 0L
  n_tot   <- n1 + n2 + n3

  if (n_tot >= T)
    warning(sprintf(
      "3stage VCV: n_tot=%d >= T=%d; B0 rank-deficient. DCC (psi) SEs remain valid.",
      n_tot, T))

  # A33 + S3 (DCC dynamics block)
  message("  3stage: A33 + S3 (DCC scores)...")
  rec <- qt_recursion_cpp(psi_D, P, N, ee, nn, bc, bca, model)
  sc  <- analytic_scores_parallel_cpp(psi_D, P, N, e_std, ee, nn,
                                      bc, bca, model, rec$Qt, rec$Rt,
                                      as.integer(n_cores))
  S3  <- sc$scores
  sfn <- function(p)
    qt_and_mean_scores_cpp(p, P, N, e_std, ee, nn, bc, bca, model, as.integer(n_cores))
  jac_method <- if (n3 <= 5L) "Richardson" else "simple"
  A33 <- .sym(-numDeriv::jacobian(sfn, psi_D, method = jac_method))

  # Baseline DCC mean-scores for A31 finite differences (reuse cubes already in memory)
  s3_base <- qt_and_mean_scores_cpp(psi_D, P, N, e_std, ee, nn,
                                    bc, bca, model, as.integer(n_cores))

  # S1 / A11 / A21 / A31 — fused, parallel over assets
  message("  3stage: S1/A11/A21/A31 (fused, parallel over k assets)...")
  gb  <- compute_garch_blocks(e_std, P, N, psi_D, model, n_cores,
                              specs, thetas, data, s3_base)
  S1  <- gb$S1; A11 <- gb$A11; A21 <- gb$A21; A31 <- gb$A31

  # S2 (analytic intercept scores) + A22 = I_{n2}
  message("  3stage: S2 (analytic)...")
  S2      <- intercept_perobs_scores(e_std, n_std, P, N, model)
  A22_int <- diag(n2)

  # A32 — C++ OpenMP over n2 pairs (cubes shared read-only; avoids parLapply overhead)
  message("  3stage: A32 (C++ OpenMP over n2 pairs)...")
  A32 <- compute_A32_cpp(psi_D, P, N, e_std, ee, nn, bc, bca,
                          as.integer(model), as.integer(n_cores))

  # Assemble lower-triangular A0
  A0 <- matrix(0, n_tot, n_tot)
  A0[1:n1,               1:n1]            <- A11
  A0[(n1+1):(n1+n2),     1:n1]            <- A21
  A0[(n1+1):(n1+n2),     (n1+1):(n1+n2)] <- A22_int
  A0[(n1+n2+1):n_tot,    1:n1]            <- A31
  A0[(n1+n2+1):n_tot,    (n1+1):(n1+n2)] <- A32
  A0[(n1+n2+1):n_tot,    (n1+n2+1):n_tot]<- A33

  # Full sandwich: VCV(psi) = [A0_inv B0 A0_inv']_{psi,psi} / T
  message("  3stage: B0 = covnw(S_full), full sandwich...")
  S_full   <- cbind(S1, S2, S3)               # T x n_tot
  B0       <- covnw(S_full)                   # n_tot x n_tot
  A0inv    <- .safe_solve(A0)                 # n_tot x n_tot
  VCV_full <- A0inv %*% B0 %*% t(A0inv) / T  # n_tot x n_tot
  idx_psi  <- (n1+n2+1):n_tot
  VCV      <- .sym(VCV_full[idx_psi, idx_psi])
  se       <- sqrt(pmax(diag(VCV), 0))

  list(VCV = VCV, se = se,
       A11 = A11, A21 = A21, A31 = A31, A32 = A32, A33 = A33,
       B0 = B0,
       n1 = n1, n2 = n2, n3 = n3,
       scores = S_full, method = "3stage")
}

# ---------------------------------------------------------------------------
# Post-hoc VCV computation
# ---------------------------------------------------------------------------

#' Add or replace VCV on an existing agdcc_fit object.
#'
#' Useful when the model was fitted with vcv_method="none" and inference is
#' added separately. Requires fit$data (stored automatically by fit_agdcc).
#'
#' @param fit        agdcc_fit object (any vcv_method)
#' @param vcv_method "3stage"
#' @param n_cores    integer; NULL → detectCores()-1
#' @return updated agdcc_fit with VCV, se, tstat, pval, vcv_detail filled in
compute_vcv <- function(fit, vcv_method = "3stage", n_cores = NULL) {
  stopifnot(inherits(fit, "agdcc_fit"))
  vcv_method <- match.arg(vcv_method, "3stage")
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
  n_cores <- max(1L, as.integer(n_cores))

  psi_D   <- fit$psi_D
  P       <- fit$P
  N_mat   <- fit$N
  e_std   <- fit$e_std
  model   <- fit$model_int
  n_std   <- e_std * (e_std < 0)
  cb      <- build_outer_cubes_cpp(e_std, n_std)
  bc      <- ewma_backcast_cpp(cb$ee)
  bca     <- ewma_backcast_cpp(cb$nn)

  message(sprintf("VCV (3stage) for %s...", fit$model))
  if (is.null(fit$data))
    stop("fit$data is NULL — re-fit with fit_agdcc() to populate fit$data")
  s1 <- fit$stage1
  if (is.null(s1$specs))
    stop("fit$stage1$specs is NULL — re-fit with fit_agdcc() (specs always stored now)")
  vcv_res <- vcv_3stage(psi_D, P, N_mat, e_std, cb$ee, cb$nn,
                        bc, bca, model, n_cores,
                        s1$specs, s1$thetas, fit$data)

  se    <- vcv_res$se; names(se) <- names(psi_D)
  tstat <- psi_D / se
  pval  <- 2 * pnorm(-abs(tstat))

  fit$VCV        <- vcv_res$VCV
  fit$se         <- se
  fit$tstat      <- tstat
  fit$pval       <- pval
  fit$vcv_method <- vcv_method
  fit$vcv_detail <- vcv_res
  fit
}

#' Recompute VCV from stored vcv_detail without re-running estimation.
#'
#' Rebuilds A0 from cached blocks (A11, A21, A31, A32, A33), recomputes B0
#' from stored scores using the current covnw bandwidth, and re-extracts
#' the (psi,psi) VCV block. All A-matrix blocks remain unchanged.
#'
#' @param fit  agdcc_fit with vcv_detail populated (vcv_method="3stage")
#' @return updated agdcc_fit with VCV, se, tstat, pval, vcv_detail$B0 refreshed
recompute_vcv <- function(fit) {
  stopifnot(inherits(fit, "agdcc_fit"),
            !is.null(fit$vcv_detail),
            !is.null(fit$vcv_detail$scores))
  d     <- fit$vcv_detail
  n1    <- d$n1; n2 <- d$n2; n3 <- d$n3
  n_tot <- n1 + n2 + n3
  T     <- nrow(fit$e_std)

  A0 <- matrix(0, n_tot, n_tot)
  A0[1:n1,              1:n1]             <- d$A11
  A0[(n1+1):(n1+n2),    1:n1]             <- d$A21
  A0[(n1+1):(n1+n2),    (n1+1):(n1+n2)]  <- diag(n2)
  A0[(n1+n2+1):n_tot,   1:n1]             <- d$A31
  A0[(n1+n2+1):n_tot,   (n1+1):(n1+n2)]  <- d$A32
  A0[(n1+n2+1):n_tot,   (n1+n2+1):n_tot] <- d$A33

  B0       <- covnw(d$scores)
  A0inv    <- .safe_solve(A0)
  VCV_full <- A0inv %*% B0 %*% t(A0inv) / T
  idx_psi  <- (n1+n2+1):n_tot
  VCV      <- .sym(VCV_full[idx_psi, idx_psi])
  se       <- sqrt(pmax(diag(VCV), 0))
  names(se) <- names(fit$psi_D)
  tstat     <- fit$psi_D / se
  pval      <- 2 * pnorm(-abs(tstat))

  fit$VCV              <- VCV
  fit$se               <- se
  fit$tstat            <- tstat
  fit$pval             <- pval
  fit$vcv_detail$B0    <- B0
  fit$vcv_detail$VCV   <- VCV
  fit$vcv_detail$se    <- se
  fit
}

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Estimate DCC / ADCC / G-DCC / AG-DCC.
#'
#' @param data        T x k zero-mean returns matrix
#' @param model       "DCC" | "ADCC" | "GDCC" | "AGDCC"
#' @param stage1_fit  NULL | ugarchMultifit | list<ugarchfit>
#' @param vcv_method  "none" | "3stage"
#'   "none"   — no standard errors
#'   "3stage" — Full 3-block sandwich (E&S 2001 Theorem 2 complete): lower-triangular
#'              A0=[A11,0,0;A21,A22,0;A31,A32,A33], B0=covnw([S1|S2|S3]);
#'              VCV(psi) = bottom-right n3x n3 block of A0^{-1} B0 (A0^{-1})'/T.
#'              Requires stage1_fit to provide specs/thetas (or re-fits sGARCH(1,1)).
#' @param n_cores     integer; cores for parallel operations
#' @param cpp_file    path to agdcc_core.cpp
#' @return object of class "agdcc_fit"
fit_agdcc <- function(data,
                      model       = "AGDCC",
                      stage1_fit  = NULL,
                      vcv_method  = c("none", "3stage"),
                      n_cores     = 1L,
                      cpp_file    = "agdcc_core.cpp") {

  .load_deps(); .build_backend(cpp_file)
  vcv_method <- match.arg(vcv_method)

  model_int <- switch(toupper(model),
    "DCC"=DCC_SCALAR,"ADCC"=ADCC_SCALAR,"GDCC"=GDCC_DIAG,"AGDCC"=AGDCC_DIAG,
    stop("model must be DCC | ADCC | GDCC | AGDCC"))

  if (!is.matrix(data)) data <- as.matrix(data)
  T <- nrow(data); k <- ncol(data)
  n_cores <- max(1L, as.integer(n_cores))
  message(sprintf("=== %s | T=%d k=%d cores=%d vcv=%s ===",
                  toupper(model), T, k, n_cores, vcv_method))

  message("Stage 1...")
  s1 <- fit_stage1(data, stage1_fit, n_cores, return_fits = TRUE)

  message("Stage 2...")
  ic  <- compute_intercepts(s1$e_std)
  P   <- ic$P; N_mat <- ic$N
  cb  <- build_outer_cubes_cpp(s1$e_std, ic$n_std)
  bc  <- ewma_backcast_cpp(cb$ee); bca <- ewma_backcast_cpp(cb$nn)
  delta <- tryCatch(compute_delta(P, N_mat), error=function(e) 0.5)

  message("Stage 3...")
  opt   <- optimise_stage3(P, N_mat, s1$e_std, cb$ee, cb$nn,
                           bc, bca, model_int, delta, k, n_cores)
  psi_D <- opt$psi_D
  message(sprintf("  OK | conv=%d | LL=%.4f", opt$convergence, opt$ll))

  rec    <- qt_recursion_cpp(psi_D, P, N_mat, cb$ee, cb$nn, bc, bca, model_int)
  Qt_arr <- rec$Qt; Rt_arr <- rec$Rt

  Ht_arr <- array(0, c(k,k,T))
  for (t in seq_len(T)) {
    Dt <- diag(sqrt(s1$H[t,]), k)
    Ht_arr[,,t] <- Dt %*% Rt_arr[,,t] %*% Dt
  }

  asset_names <- if (!is.null(colnames(data))) colnames(data) else as.character(seq_len(k))
  pnames <- switch(as.character(model_int),
    "1"=c("a","b"), "2"=c("a","g","b"),
    "3"=c(sprintf("[%s].a1", asset_names), sprintf("[%s].b1", asset_names)),
    "4"=c(sprintf("[%s].a1", asset_names), sprintf("[%s].g1", asset_names),
          sprintf("[%s].b1", asset_names)))
  names(psi_D) <- pnames

  vcv_res <- NULL
  if (vcv_method == "3stage") {
    message("Inference (E&S 2001 full 3-stage sandwich)...")
    vcv_res <- vcv_3stage(psi_D, P, N_mat, s1$e_std, cb$ee, cb$nn,
                          bc, bca, model_int, n_cores,
                          s1$specs, s1$thetas, data)
  }

  se    <- if (!is.null(vcv_res)) vcv_res$se    else rep(NA_real_, length(psi_D))
  names(se) <- names(psi_D)
  tstat <- psi_D / se
  pval  <- 2 * pnorm(-abs(tstat))

  structure(
    list(psi_D=psi_D, ll=opt$ll, Ht=Ht_arr, Rt=Rt_arr, Qt=Qt_arr,
         VCV=if(!is.null(vcv_res)) vcv_res$VCV else NULL,
         se=se, tstat=tstat, pval=pval,
         H_univ=s1$H, e_std=s1$e_std, P=P, N=N_mat, delta=delta,
         model=toupper(model), model_int=model_int,
         vcv_method=vcv_method, T=T, k=k, stage1=s1, data=data,
         vcv_detail=vcv_res, convergence=opt$convergence),
    class="agdcc_fit")
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' @export
print.agdcc_fit <- function(x, ...) {
  ic <- AIC(x)
  cat(sprintf("\n=== %s  [k=%d, T=%d, vcv=%s] ===\n",
              x$model, x$k, x$T, x$vcv_method))
  cat(sprintf("Log-likelihood : %.4f\nAIC / BIC      : %.2f / %.2f\nConverged      : %s\n\n",
              x$ll, ic["AIC"], ic["BIC"], x$convergence==0))
  np          <- length(x$psi_D)
  pname_width <- max(12L, max(nchar(names(x$psi_D))) + 2L)
  ws <- c(pname_width, 12L, 12L, 10L, 10L)
  hs <- c("Parameter","Estimate","Std.Err.","t-stat","p-val")
  cat(paste(mapply(formatC, hs, width=ws, MoreArgs=list(flag="-")), collapse=""), "\n")
  cat(strrep("-", sum(ws)), "\n")
  for (j in seq_len(np)) {
    st <- if (!is.na(x$pval[j]))
            ifelse(x$pval[j]<0.01,"***",ifelse(x$pval[j]<0.05,"**",
              ifelse(x$pval[j]<0.10,"*",""))) else ""
    cat(formatC(names(x$psi_D)[j], width=pname_width, flag="-"),
        formatC(x$psi_D[j], format="f", digits=6, width=12L),
        formatC(x$se[j],    format="f", digits=6, width=12L),
        formatC(x$tstat[j], format="f", digits=4, width=10L),
        formatC(x$pval[j],  format="f", digits=4, width=10L),
        st, "\n")
  }
  cat(strrep("-", sum(ws)), "\n*** p<.01, ** p<.05, * p<.10\n\n")
  invisible(x)
}

#' @export
AIC.agdcc_fit <- function(object, ..., k=2) {
  np <- length(object$psi_D)
  c(AIC=-2*object$ll+k*np, BIC=-2*object$ll+np*log(object$T))
}

#' @export
coef.agdcc_fit <- function(object, ...) object$psi_D

#' @export
vcov.agdcc_fit <- function(object, ...) object$VCV

#' @export
conditional_correlations <- function(fit, i, j) {
  stopifnot(inherits(fit, "agdcc_fit"))
  as.numeric(fit$Rt[i, j, ])
}

