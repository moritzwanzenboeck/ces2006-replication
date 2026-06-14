# 1_data_manipulation.R
# This script loads the cleaned datasets from the source scripts (0_data_msci_imi.R and 0_data_tr_gb.R),
# standardizes missing currency values and computes USD conversions for local currency representations.

rm(list = ls())
tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), 
  error = function(e) {
    setwd("C:\\Users\\morit\\OneDrive\\Uni\\Universität Zürich\\FS 2026\\Topics in Time Series Analysis\\Seminar Paper\\implementation")
  }
)

# Load required packages
packages <- c("dplyr", "readr", "purrr", "stringr", "quantmod", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# -------------------------------------------------------------------------
# 1. LOAD CLEANED DATA
# -------------------------------------------------------------------------

msci_path <- "msci_country_imi_indeces/msci_imi_clean_data.csv"
tr_path <- "tr_10yr_gb_indices/tr_gb_clean_data.csv"

# Check if data exists
if (!file.exists(msci_path) | !file.exists(tr_path)) {
  stop("Cleaned data not found. Please run '0_data_msci_imi.R' and '0_data_tr_gb.R' first.")
}

msci_df <- read_csv(msci_path, show_col_types = FALSE)
tr_df <- read_csv(tr_path, show_col_types = FALSE)

# -------------------------------------------------------------------------
# 2. CALCULATE CURRENCY CONVERSIONS
# -------------------------------------------------------------------------

# Generic function to convert local currencies to USD utilizing FRED exchange rates
convert_to_usd <- function(df) {
  unique_currencies <- unique(df$currency[df$currency != "USD" & !is.na(df$currency)])
  
  if (length(unique_currencies) == 0) return(df)
  
  # Define FRED mapping for currencies
  # invert = TRUE implies the rate is expressed as Foreign Currency per 1 USD
  fred_mapping <- list(
    "AUD" = list(sym = "DEXUSAL", invert = FALSE), # 1 AUD in USD
    "EUR" = list(sym = "DEXUSEU", invert = FALSE), # 1 EUR in USD
    "GBP" = list(sym = "DEXUSUK", invert = FALSE), # 1 GBP in USD
    "NZD" = list(sym = "DEXUSNZ", invert = FALSE), # 1 NZD in USD
    "CAD" = list(sym = "DEXCAUS", invert = TRUE),  # CAD per 1 USD
    "CHF" = list(sym = "DEXSZUS", invert = TRUE),
    "JPY" = list(sym = "DEXJPUS", invert = TRUE),
    "CNY" = list(sym = "DEXCHUS", invert = TRUE),
    "SEK" = list(sym = "DEXSDUS", invert = TRUE),
    "DKK" = list(sym = "DEXDNUS", invert = TRUE),
    "NOK" = list(sym = "DEXNOUS", invert = TRUE),
    "HKD" = list(sym = "DEXHKUS", invert = TRUE),
    "SGD" = list(sym = "DEXSIUS", invert = TRUE),
    "MXN" = list(sym = "DEXMXUS", invert = TRUE),
    "KRW" = list(sym = "DEXKOUS", invert = TRUE),
    "INR" = list(sym = "DEXINUS", invert = TRUE),
    "BRL" = list(sym = "DEXBZUS", invert = TRUE),
    "ZAR" = list(sym = "DEXSFUS", invert = TRUE)
  )
  
  # Fetch exchange rates
  rates_list <- map(unique_currencies, function(curr) {
    if (!(curr %in% names(fred_mapping))) {
      warning(paste("No FRED symbol mapping configured for currency:", curr))
      return(NULL)
    }
    
    symbol <- fred_mapping[[curr]]$sym
    tryCatch({
      getSymbols(symbol, src = "FRED", auto.assign = FALSE)
    }, error = function(e) {
      warning(paste("Could not fetch rates for", curr))
      return(NULL)
    })
  })
  names(rates_list) <- unique_currencies
  
  # Process conversions across the unique currency list
  usd_versions <- map_dfr(unique_currencies, function(curr) {
    if (is.null(rates_list[[curr]])) return(NULL)
    
    curr_data <- df %>% filter(currency == curr)
    rates <- rates_list[[curr]]
    
    rates_df <- data.frame(date = as.Date(index(rates)), rate = as.numeric(rates[, 1]))

    # Expand to every calendar day and forward-fill so weekly bond dates
    # (which may fall on weekends) always find a rate.
    rates_filled <- data.frame(
      date = seq(min(rates_df$date), max(rates_df$date), by = "day")
    ) %>%
      left_join(rates_df, by = "date") %>%
      arrange(date) %>%
      fill(rate, .direction = "down") %>%
      fill(rate, .direction = "up")   # backward fill for the very first gap

    curr_data %>%
      left_join(rates_filled, by = "date") %>%
      fill(rate, .direction = "downup") %>%
      mutate(
        rate = if (fred_mapping[[curr]]$invert) 1 / rate else rate,
        value = value * rate,
        currency = "USD",
        is_calculated = TRUE
      ) %>%
      select(-rate)
  })
  
  # Return original combined with USD variants
  bind_rows(df, usd_versions)
}

# Custom function dedicated to constructing Local Currency versions for MSCI assets reported strictly in USD
calculate_missing_lcy <- function(df) {
  # Map index sub-tokens to their intrinsic home currencies
  country_currency_map <- list(
    "Australia" = "AUD",
    "EMU" = "EUR", "Europe" = "EUR", "Germany" = "EUR", "France" = "EUR", 
    "Italy" = "EUR", "Spain" = "EUR", "Netherlands" = "EUR",
    "United Kingdom" = "GBP", "UK" = "GBP",
    "New Zealand" = "NZD",
    "Canada" = "CAD",
    "Switzerland" = "CHF",
    "Japan" = "JPY",
    "China" = "CNY",
    "Sweden" = "SEK",
    "Denmark" = "DKK",
    "Norway" = "NOK",
    "Hong Kong" = "HKD",
    "Singapore" = "SGD",
    "Mexico" = "MXN",
    "Korea" = "KRW",
    "India" = "INR",
    "Brazil" = "BRL",
    "South Africa" = "ZAR"
  )
  
  usd_only_indices <- df %>%
    group_by(index_name) %>%
    summarise(has_only_usd = all(currency == "USD"), .groups = "drop") %>%
    filter(has_only_usd) %>%
    pull(index_name)
  
  if (length(usd_only_indices) == 0) return(df)
  
  usd_only_data <- df %>% filter(index_name %in% usd_only_indices & !is_calculated)
  
  fred_mapping <- list(
    "AUD" = list(sym = "DEXUSAL", invert = FALSE),
    "EUR" = list(sym = "DEXUSEU", invert = FALSE),
    "GBP" = list(sym = "DEXUSUK", invert = FALSE),
    "NZD" = list(sym = "DEXUSNZ", invert = FALSE),
    "CAD" = list(sym = "DEXCAUS", invert = TRUE),
    "CHF" = list(sym = "DEXSZUS", invert = TRUE),
    "JPY" = list(sym = "DEXJPUS", invert = TRUE),
    "CNY" = list(sym = "DEXCHUS", invert = TRUE),
    "SEK" = list(sym = "DEXSDUS", invert = TRUE),
    "DKK" = list(sym = "DEXDNUS", invert = TRUE),
    "NOK" = list(sym = "DEXNOUS", invert = TRUE),
    "HKD" = list(sym = "DEXHKUS", invert = TRUE),
    "SGD" = list(sym = "DEXSIUS", invert = TRUE),
    "MXN" = list(sym = "DEXMXUS", invert = TRUE),
    "KRW" = list(sym = "DEXKOUS", invert = TRUE),
    "INR" = list(sym = "DEXINUS", invert = TRUE),
    "BRL" = list(sym = "DEXBZUS", invert = TRUE),
    "ZAR" = list(sym = "DEXSFUS", invert = TRUE)
  )
  
  lcy_versions <- map_dfr(usd_only_indices, function(idx_name) {
    implied_curr <- NA
    for (country in names(country_currency_map)) {
      if (str_detect(idx_name, regex(country, ignore_case = TRUE))) {
        implied_curr <- country_currency_map[[country]]
        break
      }
    }
    
    if (is.na(implied_curr) || implied_curr == "USD") return(NULL)
    
    symbol <- fred_mapping[[implied_curr]]$sym
    if (is.null(symbol)) return(NULL)
    
    idx_data <- usd_only_data %>% filter(index_name == idx_name)
    
    tryCatch({
      rates <- getSymbols(symbol, src = "FRED", auto.assign = FALSE)
      rates_df <- data.frame(date = as.Date(index(rates)), rate = as.numeric(rates[, 1]))

      rates_filled <- data.frame(
        date = seq(min(rates_df$date), max(rates_df$date), by = "day")
      ) %>%
        left_join(rates_df, by = "date") %>%
        arrange(date) %>%
        fill(rate, .direction = "down") %>%
        fill(rate, .direction = "up")

      idx_data %>%
        left_join(rates_filled, by = "date") %>%
        fill(rate, .direction = "downup") %>%
        mutate(
          exchange_rate = if (fred_mapping[[implied_curr]]$invert) 1 / rate else rate,
          value = value / exchange_rate,
          currency = implied_curr,
          is_calculated = TRUE
        ) %>%
        select(-rate, -exchange_rate)
    }, error = function(e) {
      warning(paste("Could not fetch rates to calculate LCY for", idx_name))
      return(NULL)
    })
  })
  
  bind_rows(df, lcy_versions)
}

# -------------------------------------------------------------------------
# 3. APPLY MANIPULATIONS
# -------------------------------------------------------------------------

cat("Converting TR Gov Benchmark missing currencies to USD...\n")
tr_df_usd <- convert_to_usd(tr_df)

cat("Generating missing Local Currency variants for MSCI assets...\n")
msci_df_extended <- calculate_missing_lcy(msci_df)

cat("Converting available non-USD MSCI assets to USD...\n")
msci_df_usd <- convert_to_usd(msci_df_extended)

# -------------------------------------------------------------------------
# 4. SAVE MANIPULATED OUTPUTS
# -------------------------------------------------------------------------

write.csv(tr_df_usd, "tr_10yr_gb_indices/tr_gb_usd_calculated.csv", row.names = FALSE)
write.csv(msci_df_usd, "msci_country_imi_indeces/msci_imi_usd_calculated.csv", row.names = FALSE)

save(tr_df_usd, file = "tr_10yr_gb_indices/tr_gb_manipulated.RData")
save(msci_df_usd, file = "msci_country_imi_indeces/msci_imi_manipulated.RData")

cat("Manipulations successful. Process completed.\n")
