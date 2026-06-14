# 2_returns_inspect.R  (run manually — not part of the main pipeline)
# Quick sanity checks on the weekly USD returns matrix.
# Annualized moments, correlation heatmap, rolling volatility.

rm(list = ls())
tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)),
  error = function(e) {
    setwd("C:\\Users\\morit\\OneDrive\\Uni\\Universität Zürich\\FS 2026\\Topics in Time Series Analysis\\Seminar Paper\\implementation")
  }
)

packages <- c("dplyr", "tidyr", "ggplot2", "moments")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

load("data/weekly_returns_usd.RData")   # returns_mat, returns_df, asset_meta

T_obs  <- nrow(returns_mat)
series <- colnames(returns_mat)
eq_idx <- startsWith(series, "eq_")

# -------------------------------------------------------------------------
# 1. ANNUALIZED MOMENTS TABLE
# -------------------------------------------------------------------------

moments_tbl <- data.frame(
  series   = series,
  class    = ifelse(eq_idx, "equity", "bond"),
  mean_ann = apply(returns_mat, 2, mean) * 52 * 100,
  sd_ann   = apply(returns_mat, 2, sd)   * sqrt(52) * 100,
  skew     = apply(returns_mat, 2, skewness),
  kurt     = apply(returns_mat, 2, kurtosis)
)

cat("=== Annualized moments (%) ===\n")
print(moments_tbl %>% arrange(class, series), digits = 3)

cat("\nEquity SD range: [",
    round(min(moments_tbl$sd_ann[eq_idx]), 1), ",",
    round(max(moments_tbl$sd_ann[eq_idx]), 1), "]\n")
cat("Bond   SD range: [",
    round(min(moments_tbl$sd_ann[!eq_idx]), 1), ",",
    round(max(moments_tbl$sd_ann[!eq_idx]), 1), "]\n")

# -------------------------------------------------------------------------
# 2. CORRELATION HEATMAP (equity only for readability)
# -------------------------------------------------------------------------

eq_mat  <- returns_mat[, eq_idx]
corr_eq <- cor(eq_mat)

corr_long <- as.data.frame(as.table(corr_eq)) %>%
  rename(s1 = Var1, s2 = Var2, corr = Freq)

p_heat <- ggplot(corr_long, aes(x = s1, y = s2, fill = corr)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title  = element_blank()) +
  labs(title = "Equity return correlations (full sample)", fill = "ρ")

print(p_heat)

# -------------------------------------------------------------------------
# 3. ROLLING 52-WEEK VOLATILITY (aggregate equity)
# -------------------------------------------------------------------------

eq_avg_ret <- rowMeans(returns_mat[, eq_idx])
dates      <- returns_df$date

roll_vol <- sapply(53:T_obs, function(t) sd(eq_avg_ret[(t - 52):t]) * sqrt(52) * 100)

roll_df <- data.frame(
  date    = dates[53:T_obs],
  vol_ann = roll_vol
)

p_vol <- ggplot(roll_df, aes(x = date, y = vol_ann)) +
  geom_line(colour = "steelblue") +
  labs(title  = "Rolling 52-week equity volatility (average across 21 markets)",
       x = NULL, y = "Annualized SD (%)") +
  theme_minimal()

print(p_vol)
