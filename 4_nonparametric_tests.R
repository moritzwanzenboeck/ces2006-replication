# 4_nonparametric_tests.R
# Two nonparametric asymmetry tests following CES 2006 (Section 2 / p.546-548).
# All regressions use Newey-West HAC standard errors (sandwich package).
#
# TEST 1 — VARIANCE ASYMMETRY (per series, CES p.546-548):
#   r²_{it} = alpha + beta*I(r_{i,t-1}<0) + eps_t
#   H0: beta = 0  vs  H1: beta > 0  (large-sample normal, one-sided)
#
# TEST 2 — COVARIANCE ASYMMETRY (pairwise, CES Table 4):
#   q_{ij,t} = alpha + beta*I(e_i<0,e_j<0) + gamma*I(e_i>0,e_j>0) + eps_t
#   H0: beta = gamma  vs  H1: beta ≠ gamma  (two-sided)
#   t-stat = (beta-hat - gamma-hat) / NW-SE(beta-gamma)
#
# Outputs:
#   var_results    : data.frame — series-level variance asymmetry stats
#   test_results   : data.frame — pair-level covariance asymmetry stats
#   table4_summary : proportion of pairs significant at 10% and 20%

rm(list = ls())
tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)),
  error = function(e) {
    setwd("C:\\Users\\morit\\OneDrive\\Uni\\Universität Zürich\\FS 2026\\Topics in Time Series Analysis\\Seminar Paper\\implementation")
  }
)

packages <- c("dplyr", "purrr", "tidyr", "sandwich")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# -------------------------------------------------------------------------
# 1. LOAD
# -------------------------------------------------------------------------

if (!file.exists("data/stage1_fits.RData"))
  stop("Run 3_stage1_garch_fit.R first.")
if (!file.exists("data/weekly_returns_usd.RData"))
  stop("Run 2_returns_compute.R first.")

load("data/stage1_fits.RData")        # e_std, h_mat, stage1_best, stage1_selected
load("data/weekly_returns_usd.RData") # returns_mat, returns_df, asset_meta

series_names <- colnames(e_std)
n            <- ncol(e_std)
T_obs        <- nrow(e_std)

# -------------------------------------------------------------------------
# 2. TEST 1 — VARIANCE ASYMMETRY (CES 2006, p.546-548)
# -------------------------------------------------------------------------
# Regress r²_it on I(r_{i,t-1} < 0) with NW-HAC SE.
# r_{i,t-1} < 0 → larger r²_it means asymmetric volatility (leverage / vol-feedback).

test_var_asym <- function(nm) {
  r     <- returns_mat[, nm]
  r_sq  <- r^2
  I_lag <- c(NA_real_, as.integer(r[-length(r)] < 0))
  valid <- !is.na(I_lag)
  fit   <- lm(r_sq[valid] ~ I_lag[valid])
  vcv   <- sandwich::NeweyWest(fit)
  cf    <- coef(fit)
  b     <- cf["I_lag[valid]"]
  se    <- sqrt(vcv["I_lag[valid]", "I_lag[valid]"])
  tstat <- b / se
  data.frame(series = nm, beta = unname(b),
             tstat_var = unname(tstat),
             p_var = pnorm(tstat, lower.tail = FALSE),
             stringsAsFactors = FALSE)
}

cat("Running variance asymmetry test for", ncol(returns_mat), "series...\n")
var_results <- do.call(rbind, lapply(colnames(returns_mat), test_var_asym))

is_eq  <- startsWith(var_results$series, "eq_")
is_bd  <- startsWith(var_results$series, "bd_")
cat(sprintf("\n--- Test 1: Variance asymmetry (NW-HAC, one-sided H1: beta>0) ---\n"))
report_var <- function(label, mask) {
  cat(sprintf("%-7s %d/%d beta>0; sig: %d at 10%%, %d at 5%%, %d at 1%% (of %d)\n",
              label,
              sum(var_results$beta[mask] > 0), sum(mask),
              sum(var_results$p_var[mask] < 0.10, na.rm = TRUE),
              sum(var_results$p_var[mask] < 0.05, na.rm = TRUE),
              sum(var_results$p_var[mask] < 0.01, na.rm = TRUE),
              sum(mask)))
}
report_var("All:",    rep(TRUE, nrow(var_results)))
report_var("Equity:", is_eq)
report_var("Bond:",   is_bd)

# -------------------------------------------------------------------------
# 3. REGIONAL GROUPINGS (matches CES 2006 Table 4)
# -------------------------------------------------------------------------

region_map <- list(
  # Equity
  eq_australasia = c("eq_australia", "eq_hong_kong", "eq_japan",
                     "eq_new_zealand", "eq_singapore"),
  eq_europe      = c("eq_austria", "eq_belgium", "eq_denmark", "eq_france",
                     "eq_germany", "eq_ireland", "eq_italy", "eq_netherlands",
                     "eq_norway", "eq_spain", "eq_sweden", "eq_switzerland",
                     "eq_united_kingdom"),
  eq_namerica    = c("eq_canada", "eq_mexico", "eq_usa"),
  # Bond
  bd_australasia = c("bd_australia", "bd_japan", "bd_new_zealand"),
  bd_europe      = c("bd_austria", "bd_belgium", "bd_denmark", "bd_france",
                     "bd_germany", "bd_ireland", "bd_italy", "bd_netherlands",
                     "bd_norway", "bd_spain", "bd_sweden", "bd_switzerland",
                     "bd_united_kingdom"),
  bd_namerica    = c("bd_canada", "bd_united_states")
)

# Filter to only series present in e_std
region_map <- lapply(region_map, function(s) intersect(s, series_names))

eq_series <- series_names[startsWith(series_names, "eq_")]
bd_series <- series_names[startsWith(series_names, "bd_")]

# Build all pairs and classify them
all_pairs <- combn(series_names, 2, simplify = FALSE)

classify_pair <- function(s1, s2) {
  is_eq1 <- s1 %in% eq_series; is_eq2 <- s2 %in% eq_series

  # Intra-equity / intra-bond / inter
  type <- if (is_eq1 && is_eq2) "intra_equity" else
          if (!is_eq1 && !is_eq2) "intra_bond" else "inter"

  # Regional sub-group for same-asset pairs
  region <- NA_character_
  if (type != "inter") {
    for (rn in names(region_map)) {
      if (s1 %in% region_map[[rn]] && s2 %in% region_map[[rn]]) {
        region <- rn; break
      }
    }
  }
  list(type = type, region = region)
}

# -------------------------------------------------------------------------
# 4. TEST 2 — COVARIANCE ASYMMETRY (CES 2006, p.548, Table 4)
# -------------------------------------------------------------------------

run_pair_test <- function(i, j) {
  ei   <- e_std[, i]; ej <- e_std[, j]
  q    <- ei * ej
  I_mm <- as.integer(ei < 0 & ej < 0)
  I_pp <- as.integer(ei > 0 & ej > 0)

  if (sum(I_mm) < 2L || sum(I_pp) < 2L)
    return(list(tstat = NA_real_, pval = NA_real_))

  fit  <- lm(q ~ I_mm + I_pp)
  vcv  <- sandwich::NeweyWest(fit)   # HAC SE: autocorrelation-robust
  cf   <- coef(fit)
  diff  <- cf["I_mm"] - cf["I_pp"]
  se    <- sqrt(vcv["I_mm","I_mm"] + vcv["I_pp","I_pp"] - 2 * vcv["I_mm","I_pp"])
  tstat <- diff / se
  pval  <- 2 * pnorm(abs(tstat), lower.tail = FALSE)  # two-sided, large-sample normal
  list(tstat = as.numeric(tstat), pval = as.numeric(pval))
}

cat("Running", length(all_pairs), "pairwise tests...\n")

test_results <- map_dfr(all_pairs, function(pair) {
  s1 <- pair[1]; s2 <- pair[2]
  i  <- which(series_names == s1)
  j  <- which(series_names == s2)
  cl <- classify_pair(s1, s2)
  tt <- run_pair_test(i, j)

  data.frame(
    series_i = s1, series_j = s2,
    type     = cl$type, region = cl$region,
    tstat    = tt$tstat, p_reg = tt$pval,
    stringsAsFactors = FALSE
  )
})

# -------------------------------------------------------------------------
# 5. TABLE 4 SUMMARY
# -------------------------------------------------------------------------

summarise_group <- function(df, label) {
  n_total <- nrow(df)
  n10 <- sum(df$p_reg < 0.10, na.rm = TRUE)
  n20 <- sum(df$p_reg < 0.20, na.rm = TRUE)
  data.frame(group = label, n_pairs = n_total,
             pct_sig10 = n10 / n_total, pct_sig20 = n20 / n_total)
}

all_eq_pairs <- test_results %>% filter(type == "intra_equity")
all_bd_pairs <- test_results %>% filter(type == "intra_bond")
all_cr_pairs <- test_results %>% filter(type == "inter")

table4_summary <- bind_rows(
  summarise_group(test_results,                                "All pairs"),
  summarise_group(all_eq_pairs,                               "Intra-equity"),
  summarise_group(all_bd_pairs,                               "Intra-bond"),
  summarise_group(all_cr_pairs,                               "Inter (eq-bd)"),
  # Regional equity
  map_dfr(c("eq_australasia","eq_europe","eq_namerica"), function(r) {
    sub <- test_results %>% filter(type == "intra_equity",
                                   !is.na(region) & region == r)
    if (nrow(sub) == 0) return(NULL)
    summarise_group(sub, paste0("Eq: ", r))
  }),
  # Regional bond
  map_dfr(c("bd_australasia","bd_europe","bd_namerica"), function(r) {
    sub <- test_results %>% filter(type == "intra_bond",
                                   !is.na(region) & region == r)
    if (nrow(sub) == 0) return(NULL)
    summarise_group(sub, paste0("Bd: ", r))
  })
)

cat("\n--- Test 2: Covariance asymmetry (NW-HAC regression, CES 2006 Table 4) ---\n")
print(table4_summary, digits = 3)

# -------------------------------------------------------------------------
# 6. SAVE
# -------------------------------------------------------------------------

save(var_results, test_results, table4_summary, file = "data/nonparam_tests.RData")
cat("\nSaved data/nonparam_tests.RData\n")
