# ==============================================================================
# Script: 07_bayesian_spatial_model.R
# Description: Bayesian spatial smoothing model (CARBYM) and integrated 
#              epidemiological risk typology for Black maternal mortality 
#              in Brazil (2019–2024).
# Inputs: domat_clean_2019_2024.parquet, sinasc_clean_2019_2024.parquet, 
#         tab_cluster.parquet
# Output: High-resolution integrated risk scenario map (PNG and PDF)
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
  spdep,
  CARBayes
)

# ==============================================================================
# 2. Load Processed Clean Datasets
# ==============================================================================
data_proc_dir <- "data_processed"
figures_dir   <- "figures"

sim_path     <- file.path(data_proc_dir, "domat_clean_2019_2024.parquet")
sinasc_path  <- file.path(data_proc_dir, "sinasc_clean_2019_2024.parquet")
cluster_path <- file.path(data_proc_dir, "tab_cluster.parquet")

if (!file.exists(sim_path) || !file.exists(sinasc_path) || !file.exists(cluster_path)) {
  stop("[ERROR] Required processed files not found. Ensure previous scripts have executed.")
}

message("[INFO] Loading clean datasets...")
sim_clean    <- read_parquet(sim_path)
sinasc_clean <- read_parquet(sinasc_path)
cluster_df   <- read_parquet(cluster_path)

# Vector for Black/Parda race/color codes
black_races <- c("Preta", "Parda", "2", "4")

# ==============================================================================
# 3. Process Cases & Spatial Population Denominators
# ==============================================================================
# 3.1 Observed Maternal Deaths (Black/Parda)
deaths_grouped <- sim_clean %>%
  filter(RACACOR %in% black_races) %>%
  group_by(CODMUNRES) %>%
  summarise(Y = n(), .groups = "drop")

# 3.2 Live Birth Denominators (Black/Parda Mothers) with Spatial Filtering
pop_spatial <- sinasc_clean %>%
  filter(RACACORMAE %in% black_races) %>%
  mutate(CODMUNRES_str = str_trim(as.character(CODMUNRES))) %>%
  filter(
    !str_ends(CODMUNRES_str, "0000"),
    !is.na(CODMUNRES_str),
    CODMUNRES_str != ""
  ) %>%
  group_by(CODMUNRES) %>%
  summarise(pop = n(), .groups = "drop")

# 3.3 Merge Cases, Population, and Compute Standardized Expected Cases (E)
d <- full_join(deaths_grouped, pop_spatial, by = "CODMUNRES") %>%
  filter(!is.na(CODMUNRES)) %>%
  mutate(
    Y   = coalesce(Y, 0L),
    pop = coalesce(pop, 0L)
  ) %>%
  arrange(CODMUNRES) %>%
  mutate(
    E   = SpatialEpi::expected(population = pop, cases = Y, n.strata = 1),
    SMR = ifelse(E > 0, Y / E, 0)
  )

# ==============================================================================
# 4. Spatial Join with IBGE Boundaries
# ==============================================================================
message("[INFO] Downloading IBGE municipality boundaries...")
muni_sf <- read_municipality(year = 2020, showProgress = FALSE) %>%
  mutate(CODMUNRES = as.character(substr(code_muni, 1, 6)))

d_munic <- muni_sf %>%
  left_join(d, by = "CODMUNRES") %>%
  mutate(
    Y   = coalesce(Y, 0L),
    pop = coalesce(pop, 0L),
    E   = ifelse(is.na(E) | E == 0, 0.001, E) # Small constant offset for Poisson GLM
  )

which(d_munic$pop == 0)

# Zero out isolated island municipalities missing spatial contiguity ("4300001", "4300002") and NA births ("4120309")
d_munic[d_munic$code_muni %in% c("4300001", "4300002", "4120309"), c("Y", "pop", "E")] <- 0

# ==============================================================================
# 5. Spatial Adjacency Matrix Construction (W Matrix)
# ==============================================================================
message("[INFO] Constructing spatial weight matrix (Queen contiguity)...")
nb <- poly2nb(muni_sf, queen = TRUE)

# Connect island municipalities to nearest spatial neighbor
deg0 <- which(card(nb) == 0)
if (length(deg0) > 0) {
  message(sprintf("[INFO] Connecting %d isolated island municipalities...", length(deg0)))
  for (i in deg0) {
    dist_vec <- as.numeric(st_distance(muni_sf[i, ], muni_sf))
    dist_vec[i] <- Inf
    nearest_idx <- which.min(dist_vec)
    
    nb[[i]]           <- as.integer(nearest_idx)
    nb[[nearest_idx]] <- sort(unique(as.integer(c(nb[[nearest_idx]], i))))
  }
}

class(nb) <- "nb"
attr(nb, "region.id") <- as.character(1:nrow(muni_sf))

# Generate symmetric binary adjacency matrix
W <- nb2mat(nb, style = "B", zero.policy = TRUE)
W_sym <- (W + t(W)) / 2
W_final <- matrix(as.numeric(W_sym > 0), nrow = nrow(W), ncol = ncol(W))

if (!isSymmetric(W_final)) {
  stop("[ERROR] Adjacency matrix W is non-symmetric.")
}

# ==============================================================================
# 6. Fit Bayesian CARBYM Model (Besag-York-Mollié)
# ==============================================================================
message("[INFO] Running Bayesian CARBYM model (S.CARbym)...")
set.seed(12345)

d_munic <- d_munic %>%
  mutate(
    E = pmax(E, 1e-6))

bayes_model <- S.CARbym(
  formula  = Y ~ offset(log(E)),
  family   = "poisson",
  data     = d_munic,
  W        = W_final,
  burnin   = 20000,
  n.sample = 120000,
  thin     = 10
)

# Extract MCMC posterior samples
beta_samples <- as.matrix(bayes_model$samples$beta)
psi_samples  <- as.matrix(bayes_model$samples$psi)

# Relative Risk (RR) posterior calculation
rr_samples <- sweep(psi_samples, 1, beta_samples[, 1], "+")

# Posterior metrics calculation
d_munic$SMR_smoothed  <- colMeans(exp(rr_samples))
d_munic$fitted_counts <- colMeans(as.matrix(bayes_model$samples$fitted))
d_munic$prob_post     <- colMeans(exp(rr_samples) > 1)

# Categorical posterior risk probability
d_munic$prob_cat <- cut(
  d_munic$prob_post,
  breaks         = c(0, 0.8, 0.95, 0.99, 1),
  labels         = c("Low", "Moderate", "High (>95%)", "Very High (>99%)"),
  include.lowest = TRUE
)

# ==============================================================================
# 7. Integrated Epidemiological Risk Typology
# ==============================================================================
# Join flexible spatial scan cluster results
cluster_sub <- cluster_df %>%
  dplyr::select(code_muni, cluster)

d_munic <- d_munic %>%
  left_join(cluster_sub, by = "code_muni") %>%
  mutate(
    risk_scenario = case_when(
      (Y == 0) & (prob_post >= 0.95 | cluster %in% c("primary", "secondary")) ~ "Epidemiological Silence (Investigate)",
      (Y > 0) & (prob_post >= 0.95) & (cluster %in% c("primary", "secondary")) ~ "Confirmed Critical Risk",
      (Y > 0) & (prob_post >= 0.95) ~ "Persistent Regional Risk",
      (Y > 0) & (cluster %in% c("primary", "secondary")) ~ "Local Excess Alert",
      TRUE ~ "Expected Risk Level"
    ),
    risk_scenario = factor(
      risk_scenario,
      levels = c(
        "Confirmed Critical Risk",
        "Persistent Regional Risk",
        "Local Excess Alert",
        "Epidemiological Silence (Investigate)",
        "Expected Risk Level"
      )
    )
  )

# Audit scenario frequencies
message("[INFO] Integrated risk scenario summary:")
print(table(d_munic$risk_scenario))

# ==============================================================================
# 8. Plot & Export Integrated Surveillance Map
# ==============================================================================
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

map <- ggplot(d_munic) +
  geom_sf(aes(fill = risk_scenario), color = "grey85", size = 0.05) +
  scale_fill_manual(
    values = c(
      "Confirmed Critical Risk"              = "#b2182b", # Dark Red
      "Persistent Regional Risk"            = "#ef8a62", # Orange
      "Local Excess Alert"                  = "#fdbb2d", # Gold Yellow
      "Epidemiological Silence (Investigate)" = "#9b59b6", # Purple (Underreporting)
      "Expected Risk Level"                 = "#f7f7f7"  # Light Grey
    ),
    name = "Surveillance Priority Level"
  ) +
  labs(
    title    = "Integrated Spatial Risk Assessment for Black Maternal Mortality",
    subtitle = "Synthesis of Bayesian CARBYM Smoothing & Flexible Spatial Scan Statistics (2019–2024)",
    caption  = "Data sources: SIM & SINASC (MS/SVSA/CGIAE). Integrated modeling via CARBayes and rflexscan."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 10),
    legend.text     = element_text(size = 9),
    plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "grey30"),
    plot.caption    = element_text(size = 8, color = "grey40", hjust = 0.5)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# Render map in graphics device
print(map)

# Export high-resolution publication maps
message("[INFO] Exporting map files...")
ggsave(
  filename = file.path(figures_dir, "risk_scenarios.png"),
  plot     = map,
  width    = 10,
  height   = 10,
  dpi      = 600
)

# Export shapefile
st_write(
  obj = d_munic,
  dsn = file.path("data_processed", "risk_black_maternal_mortality.shp"),
  delete_layer = TRUE
)


message("[SUCCESS] Script 07 executed and figures successfully saved.")

# ==============================================================================
# 9. MCMC Convergence Diagnostics (Supplementary Material)
# ==============================================================================
pacman::p_load(coda, lattice)

# Extract MCMC chain samples for hyperparameters
beta_mcmc  <- model$samples$beta
tau2_mcmc  <- model$samples$tau2
sigma2_mcmc <-model$samples$nu2 # or sigma2 depending on CARBayes version

# Export Trace & Density Plots for Main Parameters
figures_dir <- "figures"
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
supp_fig_path <- file.path(figures_dir, "supp_figure_mcmc_diagnostics.png")

png(filename = supp_fig_path, width = 2400, height = 1800, res = 300)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

# 1. Intercept (beta_0)
plot(beta_mcmc, auto.layout = FALSE, main = "Trace: Intercept (beta_0)")
# 2. Spatial Variance (tau2)
plot(tau2_mcmc, auto.layout = FALSE, main = "Trace: Spatial Variance (tau2)")
# 3. Unstructured Variance (sigma2)
if (!is.null(sigma2_mcmc)) {
  plot(sigma2_mcmc, auto.layout = FALSE, main = "Trace: Unstructured Variance (sigma2)")
}

dev.off()

# Numerical Convergence Diagnostic: Geweke's Z-score
# Values between -1.96 and +1.96 indicate convergence (p > 0.05)
geweke_diag <- geweke.diag(beta_mcmc)
message(sprintf("[DIAGNOSTIC] Geweke Z-score for Beta Intercept: %.4f", geweke_diag$z))