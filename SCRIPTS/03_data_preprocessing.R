# ==============================================================================
# Script: 03_data_preprocessing.R
# Description: Data cleaning, ICD-10 verification, municipal aggregation, 
#              and Crude Maternal Mortality Ratio (MMR) calculation (2019–2024).
# Inputs: DOMAT_2019_2024_combined.parquet, SINASC_2019_2024_combined.parquet
# Target Output: Municipality-year analytical database for spatial modeling.
# ==============================================================================

# 1. Package Loading -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(arrow, dplyr, tidyr, stringr, readr, lubridate)

# 2. Path Setup ----------------------------------------------------------------
data_raw_dir <- "data_raw"
data_proc_dir <- "data_processed"
if (!dir.exists(data_proc_dir)) dir.create(data_proc_dir, recursive = TRUE)

# 3. Load Raw Datasets ---------------------------------------------------------
message("[INFO] Loading raw datasets...")

sim_file <- file.path(data_raw_dir, "DOMAT_2019_2024_combined.parquet")
sinasc_file <- file.path(data_raw_dir, "sinasc_clean_2019_2024_corrigido.parquet")

if (!file.exists(sim_file)) {
  stop("[ERROR] SIM master file not found. Run 01_fetch_sim_domat.R first.")
}

domat_raw <- read_parquet(sim_file)

# Check if consolidated SINASC exists; if not, read available annual parquet files
if (file.exists(sinasc_file)) {
  sinasc_raw <- read_parquet(sinasc_file)
} else {
  sinasc_files <- list.files(data_raw_dir, pattern = "^SINASC_20[0-9]{2}\\.parquet$", full.names = TRUE)
  if (length(sinasc_files) == 0) {
    stop("[ERROR] No SINASC parquet files found in data_raw/.")
  }
  message(sprintf("[INFO] Found %d annual SINASC files. Combining...", length(sinasc_files)))
  sinasc_raw <- bind_rows(lapply(sinasc_files, read_parquet))
}

# 4. Process SIM Data (Maternal Deaths) ---------------------------------------
message("[PROCESSING] Standardizing SIM (Maternal Mortality) records...")

domat_clean <- domat_raw %>%
  filter(TPMORTEOCO %in% c("1", "2", "4")) %>%
  filter(TPOBITOCOR != "06") %>%
  mutate(
    # Standardize DTOBITO to extract event year (supports DDMMYYYY or YYYY-MM-DD)
    DTOBITO = parse_date_time(as.character(DTOBITO), orders = c("dmy", "ymd", "Y-m-d")),
    year = year(DTOBITO)             # extract year
    )%>% 
  select(CODMUNRES, RACACOR, CAUSABAS, year)

# 5. Process SINASC Data (Live Births) -----------------------------------------
message("[PROCESSING] Standardizing SINASC (Live Births) records...")

sinasc_clean <- sinasc_raw %>%
  select(CODMUNRES, RACACORMAE, year)

################ OVERALL MATERNAL MORTALITY

sinasc_all <- sinasc_clean %>% 
  summarise(all = n(), .groups = "drop")

domat_all <- domat_clean %>% 
  summarise(all = n(), .groups = "drop")

mmr <- domat_all / sinasc_all * 100000

################ MATERNAL MORTALITY AMONG WHITE WOMEN
White_mothers <- domat_clean %>% 
  filter(RACACOR == "1" | RACACOR == 1 | RACACOR == "Branca") %>% 
  summarise(W_mat_death = n(), .groups = "drop")

White_babies <- sinasc_clean %>% 
  filter(RACACORMAE == "1" | RACACORMAE == 1 | RACACORMAE == "Branca")%>% 
  summarise(W_babies = n(), .groups = "drop")

# Calculate Crude Maternal Mortality Ratio (MMR per 100,000 live births)
White_mmr_crude = White_mothers / White_babies * 100000 #52.3

################ MATERNAL MORTALITY AMONG COMBINED BLACK WOMEN

Black_mothers <- domat_clean %>% 
  filter(RACACOR == "2" | RACACOR == "4") %>% 
  summarise(B_mat_death = n(), .groups = "drop")

Black_babies <- sinasc_clean %>% 
  filter(RACACORMAE == "Parda" | RACACORMAE == "Preta")%>% 
  summarise(B_babies = n(), .groups = "drop")

# Calculate Crude Maternal Mortality Ratio (MMR per 100,000 live births)
Black_mmr_crude = Black_mothers / Black_babies * 100000 #55.2

################ Export domat_clean
write_parquet(domat_clean, file.path(data_proc_dir, "domat_clean_2019_2024.parquet"))

################ Export sinasc_clean
write_parquet(sinasc_clean, file.path(data_proc_dir, "sinasc_clean_2019_2024.parquet"))
