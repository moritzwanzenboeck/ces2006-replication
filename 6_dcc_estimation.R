# 6_dcc_estimation.R
# Fits all four DCC variants (DCC, ADCC, GDCC, AGDCC) to the 39-asset USD
# returns matrix using pre-fitted Stage 1 rugarch objects.
#
# Two-phase workflow:
#   Phase 1 — optimisation only (vcv_method="none"); cached per-model.
#   Phase 2 — 3-stage VCV via compute_vcv(); cached per-model.
#
# Feasibility (k=39, T=1355):
#   DCC  : n_tot ~938   < T — full rank
#   GDCC : n_tot ~1014  < T — full rank
#   ADCC : n_tot ~1719  > T — B0 rank-deficient; vcv_3stage() warns, proceeds
#   AGDCC: n_tot ~1833  > T — B0 rank-deficient; vcv_3stage() warns, proceeds
#   (DCC psi-block SEs remain valid; S3 has rank n3 << T in all cases)
#
# Outputs:
#   fits       : named list of 4 agdcc_fit objects (with VCV)

rm(list = ls())
tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)),
  error = function(e) {
    setwd("C:\\Users\\morit\\OneDrive\\Uni\\Universität Zürich\\FS 2026\\Topics in Time Series Analysis\\Seminar Paper\\implementation")
  }
)

# -------------------------------------------------------------------------
# 0. TIMESTAMPED LOG FILE
# -------------------------------------------------------------------------

dir.create("logs", showWarnings = FALSE)
log_file <- file.path("logs", paste0("dcc_estimation_",
                                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
con_log <- file(log_file, open = "wt")
sink(con_log, split = TRUE)
sink(con_log, type = "message")
on.exit({ sink(type = "message"); sink(); close(con_log) }, add = TRUE)

ms <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

ms("=== DCC Estimation Log — started ===")

packages <- c("Rcpp", "RcppArmadillo", "rugarch", "parallel", "dplyr", "rmgarch")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# -------------------------------------------------------------------------
# 1. LOAD AG-DCC ESTIMATOR AND DATA
# -------------------------------------------------------------------------

ms("Compiling C++ backend...")
Rcpp::sourceCpp("agdcc_core.cpp")
source("agdcc.R")
ms("Backend ready.")

if (!file.exists("data/weekly_returns_usd.RData"))
  stop("Run 2_returns_compute.R first.")
if (!file.exists("data/stage1_fits.RData"))
  stop("Run 3_stage1_garch_fit.R first.")

load("data/weekly_returns_usd.RData")   # returns_mat, returns_df, asset_meta
load("data/stage1_fits.RData")          # stage1_best

T_obs   <- nrow(returns_mat)
k_obs   <- ncol(returns_mat)
n_cores <- max(1L, detectCores() - 1L)
ms(sprintf("Data: T=%d  k=%d  cores=%d", T_obs, k_obs, n_cores))

model_names <- c("DCC", "ADCC", "GDCC", "AGDCC")
fits        <- setNames(vector("list", length(model_names)), model_names)

# -------------------------------------------------------------------------
# 2. PHASE 1 — OPTIMISATION (vcv_method="none")
# -------------------------------------------------------------------------

ms(sprintf("\n=== PHASE 1: Optimisation (%d models, no VCV) ===", length(model_names)))

for (m in model_names) {
  cache_none <- sprintf("data/dcc_fit_%s_none_tmp.rds", m)

  if (file.exists(cache_none)) {
    fits[[m]] <- readRDS(cache_none)
    ms(sprintf("%-6s: loaded from cache | conv=%d | LL=%.4f",
               m, fits[[m]]$convergence, fits[[m]]$ll))
    next
  }

  ms(sprintf("%-6s: starting optimisation...", m))
  t0 <- proc.time()
  fits[[m]] <- tryCatch(
    fit_agdcc(returns_mat,
              model      = m,
              stage1_fit = stage1_best,
              vcv_method = "none",
              n_cores    = n_cores),
    error = function(e) { ms(sprintf("ERROR: %s", e$message)); NULL }
  )
  elapsed <- (proc.time() - t0)["elapsed"]

  if (!is.null(fits[[m]])) {
    ms(sprintf("%-6s: done [%.0fs] | conv=%d | LL=%.4f",
               m, elapsed, fits[[m]]$convergence, fits[[m]]$ll))
    cat(sprintf("  Parameters: %s\n",
                paste(sprintf("%s=%.5f", names(fits[[m]]$psi_D), fits[[m]]$psi_D),
                      collapse="  ")))
    saveRDS(fits[[m]], cache_none)
    ms(sprintf("%-6s: cached → %s", m, cache_none))
  }
}

ms("\nPhase 1 complete — all models optimised.")

# Quick LL summary
cat("\n--- Phase 1 LL summary ---\n")
for (m in model_names) {
  f <- fits[[m]]
  if (!is.null(f)) cat(sprintf("  %-6s  LL=%10.4f  conv=%d\n", m, f$ll, f$convergence))
}

# -------------------------------------------------------------------------
# 3. PHASE 2 — 3-STAGE VCV
# -------------------------------------------------------------------------

ms(sprintf("\n=== PHASE 2: 3-stage VCV (%d models) ===", length(model_names)))
ms("(compute_A32_cpp uses OpenMP; GARCH perturbations via parLapply)")

for (m in model_names) {
  cache_vcv <- sprintf("data/dcc_fit_%s_3stage_tmp.rds", m)

  if (file.exists(cache_vcv)) {
    fits[[m]] <- readRDS(cache_vcv)
    ms(sprintf("%-6s: VCV loaded from cache | SE range [%.6f, %.6f]",
               m, min(fits[[m]]$se, na.rm=TRUE), max(fits[[m]]$se, na.rm=TRUE)))
    print(fits[[m]])
    next
  }

  if (is.null(fits[[m]])) {
    ms(sprintf("%-6s: skipped (Phase 1 fit missing)", m)); next
  }

  ms(sprintf("%-6s: starting 3-stage VCV...", m))
  t0 <- proc.time()
  fits[[m]] <- tryCatch(
    compute_vcv(fits[[m]], vcv_method = "3stage", n_cores = n_cores),
    error = function(e) { ms(sprintf("ERROR in VCV: %s", e$message)); fits[[m]] }
  )
  elapsed <- (proc.time() - t0)["elapsed"]

  ms(sprintf("%-6s: VCV done [%.0fs] | SE range [%.6f, %.6f]",
             m, elapsed,
             min(fits[[m]]$se, na.rm=TRUE), max(fits[[m]]$se, na.rm=TRUE)))
  saveRDS(fits[[m]], cache_vcv)
  ms(sprintf("%-6s: cached → %s", m, cache_vcv))
  print(fits[[m]])
}

ms("\nPhase 2 complete — all VCVs computed.")

# -------------------------------------------------------------------------
# 4. POST-ESTIMATION CHECKS
# -------------------------------------------------------------------------

cat("\n--- Post-estimation checks ---\n")
for (m in model_names) {
  f <- fits[[m]]
  if (is.null(f)) { cat(m, ": fit failed\n"); next }
  psi   <- f$psi_D
  conv  <- f$convergence
  ll    <- f$ll
  se_ok <- if (is.null(f$se) || all(is.na(f$se))) "NA" else
             as.character(all(is.finite(f$se) & f$se > 0))
  cat(sprintf("%-6s: conv=%d | LL=%10.4f | psi in [%.5f, %.5f] | SE ok=%s\n",
              m, conv, ll, min(psi), max(psi), se_ok))
}

# -------------------------------------------------------------------------
# 5. LIKELIHOOD RATIO TESTS
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# 6. SAVE
# -------------------------------------------------------------------------

dir.create("data", showWarnings = FALSE)
save(fits, file = "data/dcc_fits.RData")
ms("Saved data/dcc_fits.RData")

# -------------------------------------------------------------------------
# 7. RMGARCH COMPARISON (DCC and ADCC only)
# -------------------------------------------------------------------------
# Parameterisation note:
#   Our estimator stores  a  and uses a^2 in the Qt recursion.
#   rmgarch stores dcca1 = a^2 directly.  So our_a^2 approx rmgarch_dcca1.
#   Similarly our_g^2 approx rmgarch_dccg1, and our_b approx rmgarch_dccb1
#   (rmgarch does NOT square b; it uses dccb1 directly).
# Stage-1 note: rmgarch uses sGARCH(1,1) with normal errors below;
#   our Stage 1 uses BIC-selected heterogeneous specs, so LL values differ.

if (!requireNamespace("rmgarch", quietly = TRUE)) {
  cat("\nrmgarch not installed — skipping comparison.\n")
} else {
  library(rmgarch)
  cat("\n--- rmgarch DCC / ADCC comparison ---\n")
  cat("(rmgarch Stage 1: sGARCH(1,1), normal; ours: BIC-selected heterogeneous)\n")
  cat("Parameterisation: our a^2 ~ dcca1, our g^2 ~ dccg1, our b ~ dccb1\n\n")

  uspec_uni <- ugarchspec(
    variance.model    = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model        = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  k_rm <- ncol(returns_mat)

  spec_dcc  <- dccspec(
    uspec = multispec(replicate(k_rm, uspec_uni)),
    dccOrder = c(1, 1), distribution = "mvnorm"
  )
  spec_adcc <- dccspec(
    uspec = multispec(replicate(k_rm, uspec_uni)),
    dccOrder = c(1, 1), distribution = "mvnorm", model = "aDCC"
  )

  rfit_dcc <- tryCatch(
    dccfit(spec_dcc,  data = returns_mat, solver = "solnp",
           fit.control = list(eval.se = FALSE)),
    error = function(e) { cat("rmgarch DCC error:", e$message, "\n"); NULL }
  )
  rfit_adcc <- tryCatch(
    dccfit(spec_adcc, data = returns_mat, solver = "solnp",
           fit.control = list(eval.se = FALSE)),
    error = function(e) { cat("rmgarch ADCC error:", e$message, "\n"); NULL }
  )

  if (!is.null(rfit_dcc)) {
    cp  <- coef(rfit_dcc)
    ra1 <- cp[grep("dcca1", names(cp))]
    rb1 <- cp[grep("dccb1", names(cp))]
    rll <- likelihood(rfit_dcc)
    cat(sprintf("rmgarch DCC :  dcca1=%.6f  dccb1=%.6f  LL=%.2f\n", ra1, rb1, rll))
    if (!is.null(fits[["DCC"]])) {
      our <- fits[["DCC"]]$psi_D
      cat(sprintf("Our DCC     :  a=%.6f  a^2=%.6f  b=%.6f  LL=%.2f\n",
                  our["a"], our["a"]^2, our["b"], fits[["DCC"]]$ll))
    }
  }

  if (!is.null(rfit_adcc)) {
    cp  <- coef(rfit_adcc)
    ra1 <- cp[grep("dcca1", names(cp))]
    rg1 <- cp[grep("dccg1", names(cp))]
    rb1 <- cp[grep("dccb1", names(cp))]
    rll <- likelihood(rfit_adcc)
    cat(sprintf("\nrmgarch ADCC:  dcca1=%.6f  dccg1=%.6f  dccb1=%.6f  LL=%.2f\n",
                ra1, rg1, rb1, rll))
    if (!is.null(fits[["ADCC"]])) {
      our <- fits[["ADCC"]]$psi_D
      cat(sprintf("Our ADCC    :  a=%.6f  a^2=%.6f  g=%.6f  g^2=%.6f  b=%.6f  LL=%.2f\n",
                  our["a"], our["a"]^2, our["g"], our["g"]^2, our["b"],
                  fits[["ADCC"]]$ll))
    }
  }
}

ms("=== Log ended ===")


