# 5_bootstrap.R
# Stationary bootstrap significance test for median pairwise correlations
# between asset classes (CES 2006, footnote 9).
#
# The paper reports that median bond-bond, equity-equity, and equity-bond
# correlations are each significantly different from the other two at the 1%
# level. Bootstrap distribution tabulated via Politis & Romano (1994)
# stationary bootstrap with mean block length l = 13 weeks (as stated in
# footnote 9), B = 1000 replications.
#
# Test: two-sided percentile p-value — p = 2*min(P(d_b* <= 0), P(d_b* > 0)),
# where d_b* = group1 median - group2 median. Matches CES 2006's "statistically
# different from" language, which is non-directional.
#
# Outputs: data/bootstrap_corr.RData (obs_medians, boot_medians, p_values)
#          and prints a summary table.

rm(list = ls())
tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)),
  error = function(e) {
    setwd("C:\\Users\\morit\\OneDrive\\Uni\\Universität Zürich\\FS 2026\\Topics in Time Series Analysis\\Seminar Paper\\implementation")
  }
)

packages <- c("modernBoot")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

if (!file.exists("data/weekly_returns_usd.RData"))
  stop("Run 2_returns_compute.R first.")

load("data/weekly_returns_usd.RData")   # returns_mat, asset_meta

T_obs <- nrow(returns_mat)
nms   <- colnames(returns_mat)
eq_idx <- which(startsWith(nms, "eq_"))
bd_idx <- which(startsWith(nms, "bd_"))

B      <- 1000L
L_MEAN <- 13L        # mean block length (CES footnote 8)
p      <- 1 / L_MEAN # geometric success probability (modernBoot convention)

cat("Stationary bootstrap: B =", B, ", l =", L_MEAN, "weeks\n")
cat("T =", T_obs, ", k =", ncol(returns_mat),
    "(", length(eq_idx), "equity,", length(bd_idx), "bond)\n\n")

# -------------------------------------------------------------------------
# 1. GROUP MEDIAN STATISTIC
# -------------------------------------------------------------------------

get_medians <- function(C) {
  ee <- C[eq_idx, eq_idx]; bb <- C[bd_idx, bd_idx]
  c(
    med_ee = median(ee[upper.tri(ee)]),
    med_bb = median(bb[upper.tri(bb)]),
    med_eb = median(C[eq_idx, bd_idx])   # all cross pairs (no diagonal)
  )
}

C_obs    <- cor(returns_mat)
obs      <- get_medians(C_obs)

cat("Observed median correlations:\n")
cat(sprintf("  Equity-Equity : %.4f  (%d pairs)\n",
            obs["med_ee"], sum(upper.tri(C_obs[eq_idx, eq_idx]))))
cat(sprintf("  Bond-Bond     : %.4f  (%d pairs)\n",
            obs["med_bb"], sum(upper.tri(C_obs[bd_idx, bd_idx]))))
cat(sprintf("  Equity-Bond   : %.4f  (%d pairs)\n",
            obs["med_eb"], length(eq_idx) * length(bd_idx)))
cat("\n")

# -------------------------------------------------------------------------
# 2. STATIONARY BOOTSTRAP INDICES (Politis & Romano 1994)
#    Use modernBoot::stationary_boot on a dummy 1:T series to obtain B
#    sets of integer indices; apply each to all columns simultaneously so
#    the cross-sectional correlation structure is preserved.
# -------------------------------------------------------------------------

set.seed(42)
cat("Generating bootstrap index sets via modernBoot::stationary_boot...\n")
idx_list <- modernBoot::stationary_boot(seq_len(T_obs), p = p, R = B)

# -------------------------------------------------------------------------
# 3. BOOTSTRAP
# -------------------------------------------------------------------------

boot_mat <- matrix(NA_real_, B, 3L,
                   dimnames = list(NULL, c("med_ee", "med_bb", "med_eb")))

pb_step <- B %/% 10L
for (b in seq_len(B)) {
  if (b %% pb_step == 0L)
    cat(sprintf("  bootstrap %4d / %d\n", b, B))
  idx           <- as.integer(round(idx_list[[b]]))
  boot_mat[b, ] <- get_medians(cor(returns_mat[idx, ]))
}

# -------------------------------------------------------------------------
# 4. PAIRWISE TESTS (two-sided percentile p-value)
# -------------------------------------------------------------------------

# CES 2006: "statistically different from" — no direction assumed.
# p = 2 * min(P(d_b* <= 0), P(d_b* > 0)), where d_b* = group1 - group2.

two_sided_p <- function(d_boot) {
  p_lo <- mean(d_boot <= 0)
  2 * min(p_lo, 1 - p_lo)
}

d_bb_ee <- boot_mat[, "med_bb"] - boot_mat[, "med_ee"]
d_bb_eb <- boot_mat[, "med_bb"] - boot_mat[, "med_eb"]
d_ee_eb <- boot_mat[, "med_ee"] - boot_mat[, "med_eb"]

p_values <- c(
  bb_vs_ee = two_sided_p(d_bb_ee),
  bb_vs_eb = two_sided_p(d_bb_eb),
  ee_vs_eb = two_sided_p(d_ee_eb)
)

# -------------------------------------------------------------------------
# 5. SUMMARY
# -------------------------------------------------------------------------

sig_label <- function(p) {
  if (is.na(p)) "NA"
  else if (p < 0.01) "***"
  else if (p < 0.05) "**"
  else if (p < 0.10) "*"
  else ""
}

cat("\n--- Bootstrap significance (two-sided, l = 13) ---\n")
cat(sprintf("  Bond-Bond vs Equity-Equity  : d_obs = %+.4f  p = %.4f %s\n",
            obs["med_bb"] - obs["med_ee"],
            p_values["bb_vs_ee"], sig_label(p_values["bb_vs_ee"])))
cat(sprintf("  Bond-Bond vs Equity-Bond    : d_obs = %+.4f  p = %.4f %s\n",
            obs["med_bb"] - obs["med_eb"],
            p_values["bb_vs_eb"], sig_label(p_values["bb_vs_eb"])))
cat(sprintf("  Equity-Equity vs Equity-Bond: d_obs = %+.4f  p = %.4f %s\n",
            obs["med_ee"] - obs["med_eb"],
            p_values["ee_vs_eb"], sig_label(p_values["ee_vs_eb"])))
cat("(*** p<0.01, ** p<0.05, * p<0.10)\n\n")

# Bootstrap 95% CIs for each group median
boot_ci <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)
cat("--- Bootstrap 95% CIs for group medians ---\n")
for (nm in c("med_ee", "med_bb", "med_eb")) {
  ci <- boot_ci(boot_mat[, nm])
  cat(sprintf("  %-12s : obs = %.4f  95%% CI = [%.4f, %.4f]\n",
              nm, obs[nm], ci[1], ci[2]))
}
cat("\n")

# -------------------------------------------------------------------------
# 6. SAVE
# -------------------------------------------------------------------------

bootstrap_corr <- list(
  obs_medians  = obs,
  boot_medians = boot_mat,
  p_values     = p_values,
  params       = list(B = B, l = L_MEAN, p = p, seed = 42)
)

save(bootstrap_corr, file = "data/bootstrap_corr.RData")
cat("Saved: data/bootstrap_corr.RData\n")
