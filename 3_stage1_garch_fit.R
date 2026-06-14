# 3_stage1_garch_fit.R
# Fits all 8 GARCH specs to each of the 39 return series (21 equity + 18 bond)
# using Windows-safe parallel::parLapply, selects the best model per series by
# BIC, and saves:
#   - stage1_best    : named list of 39 uGARCHfit objects (BIC-selected model)
#   - stage1_all     : named list of 39 lists of 8 uGARCHfit objects (all fits)
#   - stage1_selected: character vector of selected model names (length 39)
#   - stage1_bic     : 39×8 matrix of BIC values
#   - e_std          : T×39 matrix of standardised residuals ε̂_{it} = r_{it}/√ĥ_{it}
#   - h_mat          : T×39 matrix of conditional variances ĥ_{it}

rm(list = ls())
tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)),
  error = function(e) {
    setwd("C:\\Users\\morit\\OneDrive\\Uni\\Universität Zürich\\FS 2026\\Topics in Time Series Analysis\\Seminar Paper\\implementation")
  }
)

packages <- c("rugarch", "parallel", "dplyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# -------------------------------------------------------------------------
# 1. LOAD DATA
# -------------------------------------------------------------------------

if (!file.exists("data/weekly_returns_usd.RData"))
  stop("Run 2_returns_compute.R first.")

load("data/weekly_returns_usd.RData")   # returns_mat, returns_df, asset_meta

series_names <- colnames(returns_mat)
n_series     <- length(series_names)
T_obs        <- nrow(returns_mat)
n_specs      <- 8L

cat("Fitting", n_specs, "GARCH specs to", n_series, "series (T =", T_obs, ")...\n")

# -------------------------------------------------------------------------
# 2. FIT ONE SPEC TO ONE SERIES
# -------------------------------------------------------------------------

fit_one <- function(spec, ret_vec) {
  tryCatch({
    fit <- ugarchfit(spec = spec, data = ret_vec, solver = "solnp",
                     solver.control = list(trace = 0))
    conv <- fit@fit$convergence
    if (conv != 0 || !all(is.finite(coef(fit)))) return(NULL)
    fit
  }, error = function(e) NULL)
}

get_bic <- function(fit) {
  if (is.null(fit)) return(Inf)
  ic <- tryCatch(infocriteria(fit), error = function(e) NULL)
  if (is.null(ic)) return(Inf)
  # infocriteria returns matrix; BIC is row 2 (Bayes)
  as.numeric(ic["Bayes",])
}

# -------------------------------------------------------------------------
# 3. PARALLEL FITTING (Windows-safe parLapply)
# -------------------------------------------------------------------------

source("3_stage1_garch_specs.R")   # loads garch_specs

n_cores <- max(1L, detectCores() - 1L)
cat("Using", n_cores, "cores\n")
cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(rugarch))
})
clusterExport(cl, c("garch_specs", "returns_mat", "fit_one", "get_bic"), envir = environment())

# Build task list as a list of (si, spec_name) pairs — passed directly to workers
tasks <- do.call(c, lapply(names(garch_specs), function(sn)
  lapply(seq_len(n_series), function(si) list(si = si, spec_name = sn))))

cat("Running", length(tasks), "fits...\n")
results_flat <- parLapply(cl, tasks, function(task) {
  si        <- task$si
  spec_name <- task$spec_name
  ret_vec   <- returns_mat[, si]
  spec      <- garch_specs[[spec_name]]
  fit       <- fit_one(spec, ret_vec)
  bic       <- get_bic(fit)
  list(si = si, spec_name = spec_name, fit = fit, bic = bic)
})

stopCluster(cl)

# -------------------------------------------------------------------------
# 4. REORGANISE RESULTS
# -------------------------------------------------------------------------

# stage1_all: list of 34, each a named list of 8 fits
stage1_all <- lapply(seq_len(n_series), function(si) {
  fits_i <- results_flat[sapply(results_flat, `[[`, "si") == si]
  out <- setNames(
    lapply(fits_i, `[[`, "fit"),
    sapply(fits_i, `[[`, "spec_name")
  )
  out[names(garch_specs)]   # enforce canonical order
})
names(stage1_all) <- series_names

# BIC matrix: rows = series, cols = specs
stage1_bic <- matrix(Inf, nrow = n_series, ncol = n_specs,
                     dimnames = list(series_names, names(garch_specs)))
for (res in results_flat) {
  stage1_bic[res$si, res$spec_name] <- res$bic
}

# Best model per series
stage1_selected <- apply(stage1_bic, 1, function(row) {
  if (all(is.infinite(row))) return(NA_character_)
  names(which.min(row))
})

stage1_best <- setNames(
  lapply(seq_len(n_series), function(si) stage1_all[[si]][[stage1_selected[si]]]),
  series_names
)

# -------------------------------------------------------------------------
# 5. EXTRACT RESIDUALS AND CONDITIONAL VARIANCES
# -------------------------------------------------------------------------

e_std <- matrix(NA_real_, nrow = T_obs, ncol = n_series,
                dimnames = list(NULL, series_names))
h_mat <- matrix(NA_real_, nrow = T_obs, ncol = n_series,
                dimnames = list(NULL, series_names))

for (si in seq_len(n_series)) {
  fit <- stage1_best[[si]]
  if (is.null(fit)) next
  h_mat[, si] <- as.numeric(sigma(fit))^2
  e_std[, si] <- as.numeric(residuals(fit, standardize = TRUE))
}

# -------------------------------------------------------------------------
# 6. DIAGNOSTICS
# -------------------------------------------------------------------------

cat("\n--- Stage 1 model selection ---\n")
sel_tbl <- table(stage1_selected)
print(sel_tbl)

n_failed <- sum(is.na(stage1_selected))
if (n_failed > 0)
  cat("WARNING:", n_failed, "series had all specs fail to converge\n")

n_asym <- sum(stage1_selected %in% c("AVGARCH","eGARCH","TGARCH","gjrGARCH","apARCH","NAGARCH","NGARCH"),
              na.rm = TRUE)
cat(n_asym, "of", n_series - n_failed, "series selected an asymmetric model\n")

cat("\nAnnualized SD check (should be ~15-25% equity, ~3-8% bonds):\n")
ann_sd <- apply(returns_mat, 2, sd, na.rm = TRUE) * sqrt(52)
print(round(sort(ann_sd) * 100, 1))

# -------------------------------------------------------------------------
# 7. SAVE
# -------------------------------------------------------------------------

save(stage1_best, stage1_all, stage1_selected, stage1_bic, e_std, h_mat,
     file = "data/stage1_fits.RData")

cat("\nSaved data/stage1_fits.RData\n")
