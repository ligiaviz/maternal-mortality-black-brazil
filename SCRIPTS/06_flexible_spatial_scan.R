# ==============================================================================
# Script: 06_flexible_spatial_scan.R
# Description: Flexible scan statistics for Black Maternal deaths
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
  geobr
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
target_races <- c("Preta", "Parda", "2", "4")

# ==============================================================================
# 3. Process Cases (SIM) & Population (SINASC)
# ==============================================================================
# 3.1 Cases (Observed Deaths) by Municipality
cases_grouped <- sim_clean %>%
  filter(RACACOR %in% target_races) %>%
  group_by(CODMUNRES) %>%
  summarise(Y = n(), .groups = "drop")

# 3.2 Population (Live Births) Spatial Cleaning & Audit
pop_raw <- sinasc_clean %>%
  filter(RACACORMAE %in% target_races)

pop_spatial <- pop_raw %>%
  mutate(CODMUNRES_str = str_trim(as.character(CODMUNRES))) %>%
  filter(
    !str_ends(CODMUNRES_str, "0000"),
    !is.na(CODMUNRES_str),
    CODMUNRES_str != ""
  )

# Data Loss Audit Report
n_original    <- nrow(pop_raw)
n_filtrado    <- nrow(pop_spatial)
n_removidos   <- n_original - n_filtrado
pct_removidos <- round((n_removidos / n_original) * 100, 2)

cat(sprintf("Original records (SINASC): %d\n", n_original))
cat(sprintf("Records ready for spatial analysis: %d\n", n_filtrado))
cat(sprintf("Unassigned records discarded (ending in 0000): %d (%.2f%%)\n", n_removidos, pct_removidos))

pop_grouped <- pop_spatial %>%
  group_by(CODMUNRES) %>%
  summarise(pop = n(), .groups = "drop")

# ==============================================================================
# 4. Standardized Incidence Ratio (TIP) & Expected Cases
# ==============================================================================
d <- full_join(cases_grouped, pop_grouped, by = "CODMUNRES") %>%
  filter(!is.na(CODMUNRES)) %>%
  mutate(
    Y   = coalesce(Y, 0L),
    pop = coalesce(pop, 0L)
  ) %>%
  arrange(CODMUNRES) %>%
  mutate(
    # Compute expected cases using SpatialEpi
    E   = SpatialEpi::expected(
      population = pop,
      cases      = Y,
      n.strata   = 1
    ),
    # Calculate TIP (Standardized Rate) and MM (Rate per 1,000 live births)
    TIP = ifelse(E > 0, Y / E, 0),
    MM  = ifelse(pop > 0, (Y / pop) * 1000, 0)
  )

# Validation Checks
message(sprintf("[VALIDATION] Total Expected Cases (sum E): %.0f", sum(d$E)))
message(sprintf("[VALIDATION] Total Population (sum pop): %d", sum(d$pop)))
message(sprintf("[VALIDATION] Range TIP: [%.4f, %.4f]", min(d$TIP), max(d$TIP)))
message(sprintf("[VALIDATION] Range MM: [%.4f, %.4f]", min(d$MM), max(d$MM)))

# ==============================================================================
# 5. Spatial Join with IBGE Boundaries
# ==============================================================================
message("[INFO] Downloading IBGE municipality boundaries...")
muni_sf <- read_municipality(year = 2020, showProgress = FALSE) %>%
  mutate(CODMUNRES = as.character(substr(code_muni, 1, 6)))

d_munic <- muni_sf %>%
  left_join(d, by = "CODMUNRES")

# ==============================================================================
# 6. Spatial Scan Statistics (rflexscan) - Maternal Mortality
# ==============================================================================
pacman::p_load(spdep, rflexscan, purrr, stringr)

# Transform coordinates to Cartesian CRS (SIRGAS 2000 / Brazil Polyconic - EPSG:5880)
munic_cartesian <- st_transform(d_munic, crs = 5880)

# Create spatial adjacency weight matrix (Queen contiguity)
munic_cart.nb <- spdep::poly2nb(munic_cartesian)

# Extract centroids and coordinates (Longitude / Latitude)
centroids <- st_centroid(munic_cartesian)
coords    <- st_coordinates(centroids)

centroids <- centroids %>%
  mutate(
    long = coords[, 1],
    lat  = coords[, 2]
  )

# Zero out island / unconnected municipalities missing spatial contiguity
centroids[centroids$code_muni %in% c("4300001", "4300002", "4120309"), c("Y", "pop", "E")] <- 0

# Diagnostic Checks
message(sprintf("[DIAGNOSTIC] Total Population in Centroids: %d", sum(centroids$pop, na.rm = TRUE)))
message(sprintf("[DIAGNOSTIC] Missing Values in Centroid Population: %d", sum(is.na(centroids$pop))))

# ==============================================================================
# 7. Run Flexible Spatial Scan Test
# ==============================================================================
fls <- rflexscan(
  x           = centroids$long,
  y           = centroids$lat,
  observed    = centroids$Y,
  expected    = centroids$E,
  name        = centroids$code_muni,
  clustersize = 14,            # Maximum cluster size (10 to 15 areas)
  nb          = munic_cart.nb  # Adjacency matrix
)

# Display scan test summary
print(fls)
summary(fls)

# ==============================================================================
# 8. Extract Clusters and Spatial Join
# ==============================================================================
# Extract primary and secondary clusters with p-values into a single clean data frame
cluster_data <- map_dfr(seq_along(fls$cluster), function(i) {
  tibble(
    nome_x  = str_replace_all(as.character(fls$cluster[[i]]$name), "_", " "),
    pvalue  = fls$cluster[[i]]$pval,
    cluster = if (i == 1) "primary" else "secondary"
  )
})

# Prepare IBGE spatial boundaries and join cluster statistics
muni_sf_prep <- muni_sf %>%
  mutate(
    nome_x = as.character(code_muni),
    codigo = str_sub(code_muni, 1, 6)
  )

cluster_flexscan_sf <- muni_sf_prep %>%
  left_join(cluster_data, by = "nome_x") %>%
  mutate(
    cluster = case_when(
      is.na(cluster) | pvalue > 0.05 ~ "non-significant",
      TRUE ~ cluster
    )
  ) %>%
  st_as_sf()

# Validation Checks
print(unique(cluster_flexscan_sf$cluster))
class(cluster_flexscan_sf)

if (!dir.exists("data_processed")) dir.create("data_processed", recursive = TRUE)

# Create the table without geometry
tab_cluster <- cluster_flexscan_sf %>% 
  dplyr::select(code_muni, name_muni, name_state, CODMUNRES, cluster) %>%
  st_drop_geometry()

# Export
write_parquet(tab_cluster, file.path("data_processed", "tab_cluster.parquet"))
