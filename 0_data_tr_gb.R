# 0_data_tr_gb.R
# Reads Thomson Reuters 10-Year Government Benchmark TR index data.
# All CSV files (including (1)/(2) variants) are ingested — the extra files
# typically cover different time periods, not duplicate observations.
# Deduplicates by (date, index_name), filters to Mon–Fri, forward-fills
# gaps ≤ 7 working days, and saves outputs in the pipeline schema.

rm(list = ls())
tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)),
  error = function(e) {
    setwd("C:\\Users\\morit\\OneDrive\\Uni\\Universität Zürich\\FS 2026\\Topics in Time Series Analysis\\Seminar Paper\\implementation")
  }
)

packages <- c("readr", "dplyr", "purrr", "stringr", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# -------------------------------------------------------------------------
# 1. CONFIGURATION
# -------------------------------------------------------------------------

folder_path <- "tr_10yr_gb_indices"

# Map extracted country name → standardized bd_* index name
country_name_map <- c(
  "Australia"   = "bd_australia",
  "Austria"     = "bd_austria",
  "Belgium"     = "bd_belgium",
  "Canada"      = "bd_canada",
  "Denmark"     = "bd_denmark",
  "France"      = "bd_france",
  "Germany"     = "bd_germany",
  "Ireland"     = "bd_ireland",
  "Italy"       = "bd_italy",
  "Japan"       = "bd_japan",
  "Netherlands" = "bd_netherlands",
  "New Zealand" = "bd_new_zealand",
  "Norway"      = "bd_norway",
  "Spain"       = "bd_spain",
  "Sweden"      = "bd_sweden",
  "Switzerland" = "bd_switzerland",
  "UK"          = "bd_united_kingdom",
  "US"          = "bd_united_states"
)

currency_dict <- c(
  "bd_australia"      = "AUD",
  "bd_austria"        = "EUR",
  "bd_belgium"        = "EUR",
  "bd_canada"         = "CAD",
  "bd_denmark"        = "DKK",
  "bd_france"         = "EUR",
  "bd_germany"        = "EUR",
  "bd_ireland"        = "EUR",
  "bd_italy"          = "EUR",
  "bd_japan"          = "JPY",
  "bd_netherlands"    = "EUR",
  "bd_new_zealand"    = "NZD",
  "bd_norway"         = "NOK",
  "bd_spain"          = "EUR",
  "bd_sweden"         = "SEK",
  "bd_switzerland"    = "CHF",
  "bd_united_kingdom" = "GBP",
  "bd_united_states"  = "USD"
)

# -------------------------------------------------------------------------
# 2. READ ALL CSV FILES
# -------------------------------------------------------------------------

# Match original TR benchmark files only (not generated output files).
# The (1)/(2) variants are different date ranges and are all needed.
files <- list.files(folder_path, pattern = "^Thomson Reuters.*\\.csv$", full.names = TRUE)
cat("Found", length(files), "TR benchmark CSV files.\n")

read_benchmark_data <- function(file_path) {
  file_name    <- basename(file_path)
  country_name <- str_extract(file_name, "(?<=Reuters\\s).*(?=\\s10\\sYear)")

  if (is.na(country_name)) {
    warning(paste("Could not extract country name from:", file_name))
    return(NULL)
  }
  if (!(country_name %in% names(country_name_map))) {
    warning(paste("Unknown country name:", country_name, "in", file_name))
    return(NULL)
  }

  idx  <- country_name_map[[country_name]]
  curr <- currency_dict[[idx]]

  df <- read_csv(file_path, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
    transmute(
      date          = as.Date(Date, format = "%m/%d/%Y"),
      value         = parse_number(Price),       # handles "1,283.295"
      index_name    = idx,
      currency      = curr,
      index_level   = "Total Return",
      is_calculated = FALSE
    ) %>%
    filter(!is.na(date), !is.na(value))

  df
}

joint_raw <- files %>%
  map(read_benchmark_data) %>%
  compact() %>%
  bind_rows()

cat("Raw rows across all files:", nrow(joint_raw), "\n")

# -------------------------------------------------------------------------
# 3. DEDUPLICATE AND FILTER TO WORKING DAYS (Mon–Fri)
# -------------------------------------------------------------------------

# %u gives ISO weekday: 1=Mon … 5=Fri, 6=Sat, 7=Sun
is_workday <- function(d) as.integer(format(d, "%u")) <= 5L

joint_dedup <- joint_raw %>%
  distinct(date, index_name, .keep_all = TRUE) %>%
  filter(is_workday(date))

cat("After dedup + working-day filter:", nrow(joint_dedup), "\n")
cat("Indices found:", paste(sort(unique(joint_dedup$index_name)), collapse = ", "), "\n")

# -------------------------------------------------------------------------
# 4. EXPAND TO FULL WORKING-DAY GRID AND FORWARD-FILL ≤ 7 CONSECUTIVE DAYS
# -------------------------------------------------------------------------

forward_fill_max7 <- function(x) {
  run <- 0L
  for (i in seq_along(x)) {
    if (is.na(x[i])) {
      run <- run + 1L
      if (run <= 7L && i > 1L && !is.na(x[i - 1L])) x[i] <- x[i - 1L]
    } else {
      run <- 0L
    }
  }
  x
}

meta_cols <- joint_dedup %>%
  select(index_name, currency, index_level, is_calculated) %>%
  distinct()

joint_tr_df <- joint_dedup %>%
  group_by(index_name) %>%
  group_modify(function(.x, .y) {
    idx_min  <- min(.x$date)
    idx_max  <- max(.x$date)
    all_days <- seq(idx_min, idx_max, by = "day")
    wd       <- data.frame(date = all_days[is_workday(all_days)])

    expanded       <- left_join(wd, select(.x, date, value), by = "date")
    expanded$value <- forward_fill_max7(expanded$value)
    expanded
  }) %>%
  ungroup() %>%
  filter(!is.na(value)) %>%
  left_join(meta_cols, by = "index_name") %>%
  select(date, value, index_name, currency, index_level, is_calculated) %>%
  arrange(index_name, date)

cat("Final rows after working-day expansion + forward-fill:", nrow(joint_tr_df), "\n")

# -------------------------------------------------------------------------
# 5. DATA DICTIONARY
# -------------------------------------------------------------------------

joint_tr_df_lag <- joint_tr_df %>%
  group_by(index_name) %>%
  arrange(date) %>%
  mutate(gap_days = as.numeric(difftime(date, lag(date), units = "days"))) %>%
  ungroup()

data_dictionary_tr <- joint_tr_df_lag %>%
  group_by(index_name, currency, index_level, is_calculated) %>%
  summarise(
    start_date   = min(date, na.rm = TRUE),
    end_date     = max(date, na.rm = TRUE),
    observations = n(),
    nas          = sum(is.na(value)),
    max_gap      = max(gap_days, na.rm = TRUE),
    .groups = "drop"
  )

print("TR Data Dictionary:")
print(data_dictionary_tr, n = Inf, width = Inf)

# -------------------------------------------------------------------------
# 6. SAVE OUTPUTS
# -------------------------------------------------------------------------

write.csv(joint_tr_df,        file.path(folder_path, "tr_gb_clean_data.csv"),      row.names = FALSE)
write.csv(data_dictionary_tr, file.path(folder_path, "tr_gb_data_dictionary.csv"), row.names = FALSE)
save(joint_tr_df, data_dictionary_tr, file = file.path(folder_path, "tr_gb_data.RData"))

cat("Successfully parsed Thomson Reuters data and saved outputs.\n")
