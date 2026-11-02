# ==============================================================================
# Script: 04_mmr_trend.R
# Description: Analyzes Black (Preta/Parda) and White (Branca) maternal 
#             mortality trends (2019–2024), computes annualized MMRs, 
#             and exports high-resolution figures.
# Inputs: domat_clean_2019_2024.parquet, sinasc_clean_2019_2024.parquet
# Outputs: High-resolution Figure 1 (TIFF 600 DPI, PNG)
# ==============================================================================

# 1. Package Management --------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(arrow, dplyr, tidyr, ggplot2, broom)

# 2. Load Processed Clean Datasets ---------------------------------------------
data_proc_dir <- "data_processed"
sim_path      <- file.path(data_proc_dir, "domat_clean_2019_2024.parquet")
sinasc_path   <- file.path(data_proc_dir, "sinasc_clean_2019_2024.parquet")

if (!file.exists(sim_path) || !file.exists(sinasc_path)) {
  stop("[ERROR] Processed files not found. Run 03_data_preprocessing.R first.")
}

message("[INFO] Loading cleaned datasets...")
sim_clean    <- read_parquet(sim_path)
sinasc_clean <- read_parquet(sinasc_path)

# 3. Filter Black Population & Aggregate Trends --------------------------------
message("[PROCESSING] Computing Overall maternal mortality trend ...")

deaths <- sim_clean %>%
  group_by(year) %>%
  summarise(deaths = n(), .groups = "drop")

live_births <- sinasc_clean %>%
  group_by(year) %>%
  summarise(live_births = n(), .groups = "drop")

mmr_trend <- live_births %>%
  left_join(deaths, by = "year") %>%
  mutate(
    deaths = replace_na(deaths, 0),
    mmr = (deaths / live_births) * 100000
  ) %>%
  arrange(year)

total_deaths <- sum(mmr_trend$deaths)
total_births <- sum(mmr_trend$live_births)
period_mmr   <- (total_deaths / total_births) * 100000

message("==========================================================")
message(sprintf("Maternal Deaths (2019–2024): %d", total_deaths))
message(sprintf("Live Births (2019–2024): %d", total_births))
message(sprintf("Overall MMR (2019–2024): %.2f per 100,000", period_mmr))
message("==========================================================")

message("[PROCESSING] Computing combined Black maternal mortality trend ...")

# Filter Black maternal deaths (2 = Preta, 4 = Parda)
deaths_black <- sim_clean %>%
  filter(RACACOR %in% c("2", "4", "Preta", "Parda")) %>%
  group_by(year) %>%
  summarise(deaths = n(), .groups = "drop")

# Filter Black live births denominator
births_black <- sinasc_clean %>%
  filter(RACACORMAE %in% c("2", "4", "Preta", "Parda")) %>%
  group_by(year) %>%
  summarise(live_births = n(), .groups = "drop")

# Merge and calculate annual MMR per 100,000 live births
black_mmr_trend <- births_black %>%
  left_join(deaths_black, by = "year") %>%
  mutate(
    deaths = replace_na(deaths, 0),
    mmr_black = (deaths / live_births) * 100000
  ) %>%
  arrange(year)

# Print Summary Statistics
b_total_deaths <- sum(black_mmr_trend$deaths)
b_total_births <- sum(black_mmr_trend$live_births)
b_period_mmr   <- (b_total_deaths / b_total_births) * 100000

message("==========================================================")
message(sprintf("Black Maternal Deaths (2019–2024): %d", b_total_deaths))
message(sprintf("Black Live Births (2019–2024): %d", b_total_births))
message(sprintf("Overall Black MMR (2019–2024): %.2f per 100,000", b_period_mmr))
message("==========================================================")

message("[PROCESSING] Computing White maternal mortality trend ...")

# Filter White maternal deaths (1 = Branca)
deaths_white <- sim_clean %>%
  filter(RACACOR %in% c("1", "Branca")) %>%
  group_by(year) %>%
  summarise(deaths = n(), .groups = "drop")

# Filter Black live births denominator
births_white <- sinasc_clean %>%
  filter(RACACORMAE %in% c("1", "Branca")) %>%
  group_by(year) %>%
  summarise(live_births = n(), .groups = "drop")

# Merge and calculate annual MMR per 100,000 live births
white_mmr_trend <- births_white %>%
  left_join(deaths_white, by = "year") %>%
  mutate(
    deaths = replace_na(deaths, 0),
    mmr_white = (deaths / live_births) * 100000
  ) %>%
  arrange(year)

# Print Summary Statistics
w_total_deaths <- sum(white_mmr_trend$deaths)
w_total_births <- sum(white_mmr_trend$live_births)
w_period_mmr   <- (w_total_deaths / w_total_births) * 100000

message("==========================================================")
message(sprintf("White Maternal Deaths (2019–2024): %d", w_total_deaths))
message(sprintf("White Live Births (2019–2024): %d", w_total_births))
message(sprintf("Overall White MMR (2019–2024): %.2f per 100,000", w_period_mmr))
message("==========================================================")

# 4. Publication-Ready Figure --------------------------------------------------
message("[VISUALIZATION] Generating Figure 1...")

mmr_data <- white_mmr_trend %>% 
  left_join(black_mmr_trend, by = "year") %>% 
  select(year, mmr_black, mmr_white)

# Convert to long format (required for ggplot2 to map two series)
data_long <- mmr_data %>%
  pivot_longer(
    cols = c(mmr_black, mmr_white),
    names_to = "group",
    values_to = "mmr"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("mmr_black", "mmr_white"),
      labels = c("Black women", "White women (reference)")
    )
  )

# Fine-tuning the positioning of text labels to avoid overlap
data_long <- data_long %>%
  mutate(
    vjust_pos = case_when(
      year %in% c(2019, 2020, 2021, 2022, 2023, 2024) & group == "White women (reference)" ~ 1.8,
      year %in% c(2021) & group == "Black women" ~ 1.8,
      # Abaixo do ponto
      TRUE ~ -0.8 # Acima do ponto (padrão para todos os outros)
    )
  )

# Figure: maternal mortality trends among black and white women
p <- ggplot(data_long, aes(x = year, y = mmr, color = group, linetype = group, shape = group)) +
  # Shaded region for the critical COVID-19 period
  annotate("rect", xmin = 2019.8, xmax = 2021.2, ymin = -Inf, ymax = Inf,
           alpha = 0.12, fill = "gray20") +
  
  # Lines and points with aesthetic mapping by group
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.5) +
  
  # Numeric labels above/below the points
  geom_text(
    aes(label = sprintf("%.1f", mmr), vjust = vjust_pos),
    size = 3.5, 
    fontface = "bold", 
    show.legend = FALSE
  ) +
  
  # Custom scales for visual distinction
  scale_color_manual(values = c("Black women" = "#1B3B6F", "White women (reference)" = "#708090")) +
  scale_linetype_manual(values = c("Black women" = "solid", "White women (reference)" = "dashed")) +
  scale_shape_manual(values = c("Black women" = "circle", "White women (reference)" = "square")) +
  
  # Axis scales and limits
  scale_x_continuous(breaks = 2019:2024) +
  scale_y_continuous(limits = c(0, 135), expand = c(0, 0)) +
  
  # Axis labels
  labs(
    x = "Year",
    y = "Maternal Mortality Ratio\n(per 100,000 live births)",
    color = NULL, linetype = NULL, shape = NULL
  ) +
  
  # Classic theme
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    axis.text.x  = element_text(size = 10, color = "black", face = "bold"),
    axis.text.y  = element_text(size = 10, color = "black"),
    axis.line    = element_line(linewidth = 0.7, color = "black"),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    legend.position = "bottom",
    legend.text = element_text(size = 11, face = "bold"),
    plot.margin  = margin(t = 15, r = 15, b = 15, l = 15)
  )

# Visualizing
print(p)

# Export in high resolution (600 DPI)
# TIFF
ggsave("Figure1_MMR_Trend_new.tiff", plot = p, width = 7, height = 4.5, units = "in", dpi = 600, compression = "lzw")

# PNG
ggsave("Figure1_MMR_Trend_new.png", plot = p, width = 7, height = 4.5, units = "in", dpi = 600)

############ CONTRAFACTUAL 

# Manual Data Frame construction with annual absolute values

# Unifies data for Black women
black_data <- inner_join(deaths_black, births_black, by = "year") %>%
  mutate(race = "Black women")

# Unifies data for White women
white_data <- inner_join(deaths_white, births_white, by = "year") %>%
  mutate(race = "White women")

# Combines both groups and calculates the maternal mortality ratio (MMR)
observed <- bind_rows(black_data, white_data) %>%
    mutate(
    race = factor(race, levels = c("Black women", "White women")),
    mmr_observed = round((deaths / live_births) * 100000, 2)
  )

print(observed)

# Calculation of Counterfactual MMR and Decomposition of Variations
contrafactual <- observed %>%
  group_by(race) %>%
  mutate(
    # Reference values for 2019
    lb_2019 = live_births[year == 2019],
    deaths_2019 = deaths[year == 2019],
    
    # RMM Contrafactual (fixing live births in 2019)
    mmr_contrafactual = round((deaths / lb_2019) * 100000, 2),
    
    # Percentage changes relative to 2019
    var_pct_deaths = round(((deaths - deaths_2019) / deaths_2019) * 100, 2),
    var_pct_live_births = round(((live_births - lb_2019) / lb_2019) * 100, 2),
    
    # Magnitude of the denominator effect (Observed MMR - Counterfactual MMR)
    denominator_effect = round(mmr_observed - mmr_contrafactual, 2)
  ) %>%
  ungroup()

# Comparative table focusing on the peak of the pandemic (2021)
result_2021 <- contrafactual %>%
  filter(year == 2021) %>%
  select(
    race, 
    var_pct_deaths, 
    var_pct_live_births, 
    mmr_observed, 
    mmr_contrafactual, 
    denominator_effect
  )

print(result_2021)

# Poisson test for differences

test_poisson <- poisson.test(
  x = c(b_total_deaths, w_total_deaths),
  T = c(b_total_births, w_total_births) / 100000
)


tidy(test_poisson) %>%
  select(estimate, conf.low, conf.high, p.value) %>%
  rename(
    Rate_Ratio = estimate,
    CI_95_Low  = conf.low,
    CI_95_High = conf.high
  )