# count_garch_params.R
# Counts univariate GARCH parameters (n_1) across the 39 BIC-selected Stage-1
# fits, excluding any mean term ("mu") used to demean the data.

packages <- c("rugarch")  # registers the S4 coef() method for uGARCHfit
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
load("data/stage1_fits.RData")  # provides stage1_best, stage1_selected

# Per-asset parameter count, dropping the mean parameter if present
n_par <- vapply(stage1_best, function(fit) {
  nm <- names(coef(fit))
  sum(nm != "mu")
}, integer(1))

n1 <- sum(n_par)

# Per-asset breakdown
res <- data.frame(
  asset    = names(stage1_best),
  model    = stage1_selected,
  n_params = n_par,
  has_mu   = vapply(stage1_best, function(f) "mu" %in% names(coef(f)), logical(1)),
  row.names = NULL
)
print(res)

cat(sprintf("\nTotal univariate GARCH parameters (n_1, excl. mu): %d\n", n1))
cat(sprintf("Assets carrying a 'mu' term: %d\n", sum(res$has_mu)))
cat("\nParameters by model type:\n")
print(tapply(res$n_params, res$model, sum))

# ----------------------------------------------------------------------------
# Rank of the three-stage score covariance B0 for each DCC specification.
# B0 = T^-1 S'S is (n_tot x n_tot) but has rank <= min(T, n_tot); it is
# rank-deficient when n_tot = n1 + n2 + n3 exceeds T, with at least
# n_tot - T exact zero eigenvalues.
# ----------------------------------------------------------------------------
T_obs <- nrow(e_std)   # 1,355 observations
K     <- ncol(e_std)   # 39 assets

models <- c("DCC", "ADCC", "GDCC", "AGDCC")
asym   <- c(FALSE, TRUE, FALSE, TRUE)            # estimates N-bar target as well
n2     <- ifelse(asym, K^2, K * (K - 1) / 2)     # Stage-2 variance-targeting intercepts
n3     <- c(2, 3, 2 * K, 3 * K)                  # Stage-3 DCC dynamics
n_tot  <- n1 + n2 + n3
zero_eig  <- pmax(0L, n_tot - T_obs)
deficient <- n_tot > T_obs

rank_tbl <- data.frame(
  model            = models,
  n1               = n1,
  n2               = n2,
  n3               = n3,
  n_tot            = n_tot,
  zero_eigenvalues = zero_eig,
  rank_deficient   = ifelse(deficient, "Yes", "No"),
  row.names        = NULL
)
cat(sprintf("\nThree-stage VCV rank (T = %d):\n", T_obs))
print(rank_tbl)

# Emit a booktabs LaTeX table matching the style of table6b.
ch <- function(x) paste0("\\multicolumn{1}{c}{", x, "}")
tex <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Rank of the three-stage score covariance $\\hat{\\mathbf{B}}_0$}",
  "\\label{tab:vcv_rank}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lSSSSc}",
  "\\toprule",
  paste("Model &", ch("$n_1$"), "&", ch("$n_2$"), "&", ch("$n_3$"), "&",
        ch("$n_{\\mathrm{tot}}$"), "& Rank-deficient \\\\"),
  "\\midrule"
)
for (i in seq_along(models)) {
  tex <- c(tex, sprintf("%s & %d & %d & %d & %d & %s \\\\",
                        models[i], n1, n2[i], n3[i], n_tot[i],
                        rank_tbl$rank_deficient[i]))
}
tex <- c(tex,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\\small",
  sprintf(paste0("\\item $n_{\\mathrm{tot}} = n_1 + n_2 + n_3$, with $T = %s$. ",
                 "The mean term $\\mu$ is excluded from $n_1$. ",
                 "$\\hat{\\mathbf{B}}_0 = T^{-1}S'S$ has rank $\\leq \\min(T, n_{\\mathrm{tot}})$, ",
                 "so it is rank-deficient when $n_{\\mathrm{tot}} > T$, with at least ",
                 "$n_{\\mathrm{tot}} - T$ zero eigenvalues (here %d for ADCC, %d for AGDCC)."),
          format(T_obs, big.mark = "{,}"), zero_eig[2], zero_eig[4]),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
out_tex <- "output/tables/table_rank_deficiency.tex"
writeLines(tex, out_tex)
cat(sprintf("\nWrote %s\n", out_tex))
