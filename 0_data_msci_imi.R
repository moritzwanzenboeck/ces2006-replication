# 0_data_msci_imi.R
# This script reads and parses MSCI IMI index data from Excel files.
# It formats the data, extracts metadata, and saves the cleaned dataset 
# as well as its data dictionary into CSV and RData format for further manipulation.

rm(list = ls())
tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), 
  error = function(e) {
    # Fallback if not running in RStudio
    setwd("C:\\Users\\morit\\OneDrive\\Uni\\Universität Zürich\\FS 2026\\Topics in Time Series Analysis\\Seminar Paper\\implementation")
  }
)

# Load required packages
packages <- c("readxl", "dplyr", "purrr", "stringr", "fs", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# -------------------------------------------------------------------------
# 1. READ MSCI DATA
# -------------------------------------------------------------------------

folder_path <- "msci_country_imi_indeces"
file_list <- dir_ls(folder_path, glob = "*.xls*")

# Function to parse a single excel file
parse_msci_excel <- function(file) {
  # Read the first few rows to extract metadata
  raw_meta <- read_excel(file, col_names = FALSE, n_max = 7, .name_repair = "minimal")
  
  # Extract metadata based on known cell locations
  currency <- as.character(raw_meta[2, 2])
  index_name <- as.character(raw_meta[4, 2])
  index_level <- "Gross"
  
  # Read the actual datasets skipping the first 6 metadata rows
  data <- read_excel(
    file, 
    skip = 6, 
    col_names = c("date", "value"), 
    .name_repair = "unique", 
    range = cell_cols("A:B")
  )
  
  # Clean and format data
  data <- data %>%
    slice(5:n()) %>% # Skip non-data rows
    filter(!is.na(date) & !is.na(value)) %>% # rename from "MSCI Australia IMI Index" to "eq_australia"
    mutate(
      index_name = paste0("eq_", str_replace_all(tolower(index_name), "msci\\s|imi\\s|index", "") %>% str_trim() %>% str_replace_all("\\s+", "_")),
      date = as.Date(as.character(date), format = "%Y-%m-%d"),
      # Ensure value is numeric, replacing comma with dot if necessary
      value = as.double(str_replace_all(as.character(value), ",", ".")),
      currency = currency,
      index_level = index_level,
      is_calculated = FALSE # Identifier for values that were inherently reported, not estimated
    ) %>%
    # Select columns matching the uniform naming scheme
    select(date, value, index_name, currency, index_level, is_calculated)
  
  return(data)
}

# Apply function to all files and bind together into one long dataframe
joint_msci_df <- map_dfr(file_list, parse_msci_excel)

# -------------------------------------------------------------------------
# 2. CREATE DATA DICTIONARY
# -------------------------------------------------------------------------

data_dictionary_msci <- joint_msci_df %>%
  group_by(index_name, currency, index_level, is_calculated) %>%
  summarise(
    start_date = min(date, na.rm = TRUE),
    end_date = max(date, na.rm = TRUE),
    observations = n(),
    nas = sum(is.na(value)),
    max_gap = max(difftime(date, lag(date), units = "days"), na.rm = TRUE),
    .groups = "drop"
  )

print("Joint MSCI Dataframe (Head):")
head(joint_msci_df)
print("MSCI Data Dictionary:")
print(data_dictionary_msci, n = Inf, width = Inf)

# -------------------------------------------------------------------------
# 3. SAVE OUTPUTS
# -------------------------------------------------------------------------

# Save as CSV
write.csv(joint_msci_df, paste0(folder_path, "/msci_imi_clean_data.csv"), row.names = FALSE)
write.csv(data_dictionary_msci, paste0(folder_path, "/msci_imi_data_dictionary.csv"), row.names = FALSE)

# Save as RData
save(joint_msci_df, data_dictionary_msci, file = paste0(folder_path, "/msci_imi_data.RData"))

cat("Successfully parsed MSCI data and saved outputs.\n")
