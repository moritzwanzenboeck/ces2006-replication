# 2_returns_compute.R
# Aligns all USD price series to Thursday-to-Thursday weekly frequency,
# computes log returns, and outputs a complete-case returns matrix.
# Equity: 21 MSCI IMI USD series. Bond: 18 TR Government Benchmark USD series.

rm(list = ls())
tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)),
  error = function(e) {
    setwd("C:\\Users\\morit\\OneDrive\\Uni\\Universität Zürich\\FS 2026\\Topics in Time Series Analysis\\Seminar Paper\\implementation")
  }
)

packages <- c("dplyr", "tidyr", "readr", "lubridate", "purrr", "fs")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# -------------------------------------------------------------------------
# 1. LOAD USD PRICE SERIES
# -------------------------------------------------------------------------

msci_path <- "msci_country_imi_indeces/msci_imi_usd_calculated.csv"
tr_path   <- "tr_10yr_gb_indices/tr_gb_usd_calculated.csv"

if (!file.exists(msci_path) || !file.exists(tr_path)) {
  stop("USD data not found. Run 0_data_msci_imi.R, 0_data_tr_gb.R, and 1_data_manipulation.R first.")
}

msci_raw <- read_csv(msci_path, show_col_types = FALSE)
tr_raw   <- read_csv(tr_path,   show_col_types = FALSE)

# Keep USD series only; prefer original USD (is_calculated == FALSE) when both exist
keep_usd <- function(df) {
  df %>%
    filter(currency == "USD") %>%
    group_by(index_name, date) %>%
    arrange(is_calculated) %>%          # FALSE < TRUE, so original USD first
    slice(1) %>%
    ungroup() %>%
    select(date, index_name, value)
}

msci_usd <- keep_usd(msci_raw)
tr_usd   <- keep_usd(tr_raw)

all_usd <- bind_rows(msci_usd, tr_usd) %>%
  filter(date >= as.Date("2000-01-01"), date <= as.Date("2025-12-31"))

cat("Series available:", n_distinct(all_usd$index_name), "\n")

# -------------------------------------------------------------------------
# 2. THURSDAY ALIGNMENT
# -------------------------------------------------------------------------
# For each calendar Thursday in 2000-2025, take the most recent available
# observation on or before that Thursday (prior-or-equal business day fill).

thursdays <- seq(
  from = as.Date("2000-01-06"),   # first Thursday on or after 2000-01-01
  to   = as.Date("2025-12-31"),
  by   = "week"
)

align_to_thursdays <- function(df_series) {
  df_series <- df_series %>% arrange(date)
  # findInterval: for each Thursday, the index of the last date <= that Thursday.
  # O(T·log N) vs the previous O(T·N) loop.
  idx    <- findInterval(as.numeric(thursdays), as.numeric(df_series$date))
  values <- ifelse(idx > 0L, df_series$value[idx], NA_real_)
  data.frame(date = thursdays, value = values)
}

cat("Aligning to Thursday dates...\n")
prices_thu <- all_usd %>%
  group_by(index_name) %>%
  group_modify(~ align_to_thursdays(.x)) %>%
  ungroup()

# -------------------------------------------------------------------------
# 3. FORWARD-FILL GAPS ≤ 2 WEEKS, THEN COMPUTE LOG RETURNS
# -------------------------------------------------------------------------

prices_wide <- prices_thu %>%
  pivot_wider(names_from = index_name, values_from = value) %>%
  arrange(date)

# Forward-fill gaps of ≤ 2 consecutive weeks per column
fill_max2 <- function(x) {
  n   <- length(x)
  out <- x
  run <- 0L
  for (i in seq_len(n)) {
    if (is.na(out[i])) {
      run <- run + 1L
      if (run <= 2L && i > 1L && !is.na(out[i - 1L])) out[i] <- out[i - 1L]
    } else {
      run <- 0L
    }
  }
  out
}

price_cols <- setdiff(names(prices_wide), "date")
prices_wide[price_cols] <- lapply(prices_wide[price_cols], fill_max2)

# Log returns: r_t = log(P_t / P_{t-1}), drops first row (no lag)
returns_wide <- prices_wide %>%
  mutate(across(all_of(price_cols), ~ log(.x / lag(.x)))) %>%
  slice(-1)

# -------------------------------------------------------------------------
# 4. COMPLETE-CASE MATRIX
# -------------------------------------------------------------------------

# Report series with high NA rates before dropping
na_rates <- returns_wide %>%
  summarise(across(all_of(price_cols), ~ mean(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "series", values_to = "na_rate") %>%
  filter(na_rate > 0.05)

if (nrow(na_rates) > 0) {
  cat("\nSeries with >5% NA in returns (before complete-case drop):\n")
  print(na_rates, n = Inf)
}

returns_complete <- returns_wide %>% filter(if_all(all_of(price_cols), ~ !is.na(.x)))

cat("\n--- Returns matrix summary ---\n")
cat("T (weeks)  :", nrow(returns_complete), "\n")
cat("k (series) :", length(price_cols), "\n")
cat("Date range :", format(min(returns_complete$date)), "to", format(max(returns_complete$date)), "\n")

# -------------------------------------------------------------------------
# 5. ASSET METADATA
# -------------------------------------------------------------------------

asset_meta <- data.frame(
  index_name  = price_cols,
  asset_class = ifelse(startsWith(price_cols, "eq_"), "equity", "bond"),
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------------
# 6. SAVE OUTPUTS
# -------------------------------------------------------------------------

dir_create("data")

returns_mat <- as.matrix(returns_complete[, price_cols])
rownames(returns_mat) <- format(returns_complete$date)
returns_df  <- returns_complete

write_csv(returns_df, "data/weekly_returns_usd.csv")
save(returns_mat, returns_df, asset_meta, file = "data/weekly_returns_usd.RData")

cat("Saved data/weekly_returns_usd.csv and data/weekly_returns_usd.RData\n")
