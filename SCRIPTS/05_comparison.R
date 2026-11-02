# ==============================================================================
# Script: 05_comparison.R
# Description: Statistics between Black (Preta/Parda) and White (Branca) maternal 
#             mortality and causes of deaths
# Inputs: domat_clean_2019_2024.parquet, sinasc_clean_2019_2024.parquet
# ==============================================================================

# ==============================================================================
# 1. Package Management & Environment Setup
# ==============================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  arrow,
  dplyr,
  tidyr,
  ggplot2,
  stringr,
  SpatialEpi,
  sf,
  geobr,
  tibble
)

# ==============================================================================
# 2. Load Processed Clean Datasets
# ==============================================================================
data_proc_dir <- "data_processed"
sim_path      <- file.path(data_proc_dir, "domat_clean_2019_2024.parquet")
sinasc_path   <- file.path(data_proc_dir, "sinasc_clean_2019_2024.parquet")

if (!file.exists(sim_path) || !file.exists(sinasc_path)) {
  stop("[ERROR] Processed files not found. Run 03_data_preprocessing.R first.")
}

message("[INFO] Loading cleaned datasets...")
sim_clean    <- read_parquet(sim_path)
sinasc_clean <- read_parquet(sinasc_path)

# Vector for Black/Parda race/color codes
black <- c("Preta", "Parda", "2", "4")
white <- c("Branca", "1")

# ==============================================================================
# 3. Causes by Race
# ==============================================================================

# 1. Helper function for cause classification (prevents code duplication)
classify_cause <- function(causabas) {
  cid3 <- substr(causabas, 1, 3)
  case_when(
    cid3 %in% c("O00", "O01", "O02", "O03", "O04", "O05", "O06", "O07", "O08") ~ "Pregnancy with abortive outcome (O00-O08)",
    cid3 %in% paste0("O", 10:16) ~ "Hypertensive disorders (O10-O16)",
    cid3 %in% c("O20", "O44", "O45", "O46", "O67", "O72") | causabas %in% c("O432", "O43.2", "O622", "O62.2") ~ "Obstetric hemorrhage (O20, O43.2, O44-O46, O62.2, O67, O72)",
    cid3 %in% c("O85", "O86", "O87", "O91") ~ "Puerperal infection/sepsis (O85-O87, O91)",
    cid3 %in% c("O88", "O90", "O71", "O74", "O75") ~ "Other direct obstetric complications",
    str_detect(cid3, "^O98|^O99|^B2|^F5|^D3") ~ "Indirect obstetric causes (O98-O99 and preexisting diseases)",
    cid3 == "O95" ~ "Obstetric cause of unspecified origin (O95)",
    TRUE ~ "Other causes"
  )
}

# 2. Overall counts
sim_clean %>% 
  count(CAUSABAS, sort = TRUE) %>% 
  mutate(percentage = (n / sum(n)) * 100)

sim_clean_causes <- sim_clean %>%
  mutate(group_cause = classify_cause(CAUSABAS)) %>%
  group_by(group_cause) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percentage = (n / sum(n)) * 100) %>%
  arrange(desc(n))

# 3. Processing by racial group
black_deaths_causes <- sim_clean %>% 
  filter(RACACOR %in% black) %>%
  mutate(group_cause = classify_cause(CAUSABAS)) %>%
  group_by(group_cause) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percentage = (n / sum(n)) * 100) %>%
  arrange(desc(n))

white_deaths_causes <- sim_clean %>% 
  filter(RACACOR %in% white) %>%
  mutate(group_cause = classify_cause(CAUSABAS)) %>%
  group_by(group_cause) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percentage = (n / sum(n)) * 100) %>%
  arrange(desc(n))

# ==============================================================================
# 4. TEST FOR DIFFERENCES (Dynamic Matrix Construction)
# ==============================================================================
# Merge computed counts and convert to contingency matrix
cause_data_df <- full_join(
  black_deaths_causes %>% select(group_cause, Black = n),
  white_deaths_causes %>% select(group_cause, White = n),
  by = "group_cause"
) %>%
  mutate(
    Black = coalesce(Black, 0L),
    White = coalesce(White, 0L)
  ) %>%
  arrange(desc(Black + White))

# Convert to numeric matrix with dynamic row names
cause_data <- cause_data_df %>%
  column_to_rownames("group_cause") %>%
  as.matrix()

# Display dynamically constructed matrix
print(cause_data)

# ==============================================================================
# 5. Statistical Testing for Cause-Specific Differences
# ==============================================================================

# 5.1 Pearson's Chi-Square Test of Independence
chi_results <- chisq.test(cause_data)
print(chi_results)

# 5.2 Post-hoc Analysis: Adjusted Standardized Residuals (z-scores)
adj_residuals <- chi_results$stdres
round(adj_residuals, 2)

# ==============================================================================
# 6. Overall Rate Comparison & Risk Estimation
# ==============================================================================

# Dynamic extraction of total maternal deaths from sim_clean [1 = Black, 2 = White]
n_deaths_black <- sim_clean %>% filter(RACACOR %in% black) %>% nrow()
n_deaths_white <- sim_clean %>% filter(RACACOR %in% white) %>% nrow()

deaths <- c(n_deaths_black, n_deaths_white)

# Dynamic extraction of live births from sinasc_clean (or hardcoded vector if external)
n_births_black <- sinasc_clean %>% filter(RACACORMAE %in% black) %>% nrow()
n_births_white <- sinasc_clean %>% filter(RACACORMAE %in% white) %>% nrow()

live_births <- c(n_births_black, n_births_white)

# Option 1: Two-sample test for equality of proportions (p-value and CI for difference)
prop.test(x = deaths, n = live_births)

# Option 2: Exact Poisson Test (Exact Rate Ratio and 95% CI)
poisson.test(x = deaths, T = live_births)
