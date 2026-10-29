# ==============================================================================
# Script: 02_fetch_sinasc.R
# Description: High-performance parallel extraction pipeline for DATASUS 
#              SINASC (Live Births) data across all 27 Brazilian states (2019–2024).
# Dataset Source: Ministry of Health, Brazil (DATASUS / SINASC - DNRES & DNPRELIM)
# ==============================================================================

# 1. Package Installation and Loading ------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(read.dbc, arrow, furrr, dplyr, purrr)

# 2. Pipeline Configuration ----------------------------------------------------
# Set parallel processing across available CPU cores
plan(multisession, workers = 4)

years <- 2019:2024
ufs <- c(
  "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", 
  "MG", "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", 
  "RO", "RR", "RS", "SC", "SE", "SP", "TO"
)

output_dir <- "data_raw"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 3. Individual State-Year Downloader ------------------------------------------
download_single_sinasc <- function(uf, year) {
  # Set timeout inside worker session
  options(timeout = 300)
  
  yy <- sprintf("%02d", year %% 100)
  yyyy <- as.character(year)
  
  # Determine base folder (Recent years 2023-2024 often reside in DNPRELIM)
  folders <- if (year >= 2023) c("DNPRELIM", "DNRES") else c("DNRES", "DNPRELIM")
  
  file_bases <- c(
    sprintf("DN%s%s", uf, yyyy),
    sprintf("DN%s%s", uf, yy)
  )
  
  # Generate priority URL list
  urls <- c()
  for (folder in folders) {
    for (base in file_bases) {
      urls <- c(urls,
                sprintf("http://ftp.datasus.gov.br/dissemin/publicos/SINASC/NOV/%s/%s.dbc", folder, base),
                sprintf("http://ftp.datasus.gov.br/dissemin/publicos/SINASC/NOV/%s/%s.DBC", folder, base)
      )
    }
  }
  
  temp_dbc <- tempfile(fileext = ".dbc")
  download_success <- FALSE
  
  # Attempt downloading from constructed candidate endpoints
  for (url in urls) {
    status <- tryCatch({
      download.file(url, destfile = temp_dbc, mode = "wb", method = "libcurl", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE)
    
    if (status && file.exists(temp_dbc) && file.size(temp_dbc) > 0) {
      download_success <- TRUE
      break
    }
  }
  
  if (!download_success) {
    warning(sprintf("[WARNING] Could not retrieve SINASC for %s - %d", uf, year))
    return(NULL)
  }
  
  # Parse DBC file
  df <- tryCatch({
    read.dbc(temp_dbc)
  }, error = function(e) {
    warning(sprintf("[ERROR] Failed to parse DBC for %s - %d: %s", uf, year, e$message))
    NULL
  })
  
  if (file.exists(temp_dbc)) unlink(temp_dbc)
  return(df)
}

# 4. Annual Processing Function with Local Caching -----------------------------
process_sinasc_year <- function(year) {
  parquet_filename <- file.path(output_dir, sprintf("SINASC_%d.parquet", year))
  
  # Local Cache Check: Skip download if Parquet file already exists
  if (file.exists(parquet_filename)) {
    message(sprintf("[CACHE] SINASC Year %d already processed. Loading local Parquet...", year))
    return(read_parquet(parquet_filename))
  }
  
  message(sprintf("[DOWNLOADING] Extracting SINASC data for all UFs (Year %d)...", year))
  
  # Retrieve all 27 UFs for the target year
  uf_list <- map(ufs, ~download_single_sinasc(.x, year))
  
  # Combine non-null state data frames
  year_df <- bind_rows(compact(uf_list))
  
  if (nrow(year_df) > 0) {
    write_parquet(year_df, parquet_filename)
    message(sprintf("[SUCCESS] SINASC Year %d processed (%d records).", year, nrow(year_df)))
    return(year_df)
  } else {
    warning(sprintf("[ERROR] No records retrieved for SINASC Year %d.", year))
    return(NULL)
  }
}

# 5. Pipeline Execution --------------------------------------------------------
message("Starting parallel SINASC data extraction pipeline (2019-2024)...")
start_time <- Sys.time()

# Execute annual downloads in parallel
sinasc_list <- future_map(years, process_sinasc_year, .options = furrr_options(seed = TRUE))

# Combine all years into master national dataset
combined_sinasc <- bind_rows(compact(sinasc_list))

# Save consolidated master Parquet file
consolidated_path <- file.path(output_dir, "SINASC_2019_2024_combined.parquet")
write_parquet(combined_sinasc, consolidated_path)

end_time <- Sys.time()
elapsed <- round(as.numeric(end_time - start_time, units = "secs"), 2)

# Summary Report
message("==========================================================")
message(sprintf("Execution finished in %s seconds.", elapsed))
message(sprintf("Consolidated SINASC dataset: %d total live births.", nrow(combined_sinasc)))
message(sprintf("Master file saved to: %s", consolidated_path))
message("==========================================================")
