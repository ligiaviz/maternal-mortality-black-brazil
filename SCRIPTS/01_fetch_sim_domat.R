# ==============================================================================
# Script: 01_fetch_sim_domat.R
# Description: Robust parallel download for DATASUS SIM (DOMAT) 2019-2024
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(read.dbc, arrow, furrr, dplyr, purrr)

# Set parallel execution
plan(multisession, workers = 4)

years <- 2019:2024
output_dir <- "data_raw"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fetch_and_process_domat <- function(year) {
  # FIX: Set timeout inside worker session so background processes do not reset to 60s
  options(timeout = 600)
  
  yy <- sprintf("%02d", year %% 100)
  parquet_filename <- file.path(output_dir, sprintf("DOMAT_%d.parquet", year))
  
  # Local Cache Check
  if (file.exists(parquet_filename)) {
    message(sprintf("[CACHE] Year %d already processed. Loading local Parquet...", year))
    return(read_parquet(parquet_filename))
  }
  
  file_base <- sprintf("DOMAT%s", yy)
  
  # List of DATASUS mirrors and file extension variants (HTTP/FTP, .dbc/.DBC)
  urls <- c(
    sprintf("http://ftp.datasus.gov.br/dissemin/publicos/SIM/CID10/DOFET/%s.dbc", file_base),
    sprintf("http://ftp.datasus.gov.br/dissemin/publicos/SIM/CID10/DOFET/%s.DBC", file_base),
    sprintf("ftp://ftp.datasus.gov.br/dissemin/publicos/SIM/CID10/DOFET/%s.dbc", file_base),
    sprintf("ftp://ftp.datasus.gov.br/dissemin/publicos/SIM/CID10/DOFET/%s.DBC", file_base)
  )
  
  temp_dbc <- tempfile(fileext = ".dbc")
  download_success <- FALSE
  
  message(sprintf("[DOWNLOADING] Fetching data for year %d...", year))
  
  # Loop through mirrors until download succeeds
  for (url in urls) {
    status <- tryCatch({
      download.file(url, destfile = temp_dbc, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE)
    
    if (status && file.exists(temp_dbc) && file.size(temp_dbc) > 0) {
      download_success <- TRUE
      break
    }
  }
  
  if (!download_success) {
    warning(sprintf("[ERROR] Failed to download year %d from all DATASUS mirrors.", year))
    return(NULL)
  }
  
  # Process DBC and save as Parquet
  tryCatch({
    df <- read.dbc(temp_dbc)
    write_parquet(df, parquet_filename)
    if (file.exists(temp_dbc)) unlink(temp_dbc)
    
    message(sprintf("[SUCCESS] Year %d processed successfully.", year))
    return(df)
  }, error = function(e) {
    warning(sprintf("[ERROR] Failed to read DBC for year %d: %s", year, e$message))
    if (file.exists(temp_dbc)) unlink(temp_dbc)
    return(NULL)
  })
}

# Run parallel extraction
message("Starting parallel data extraction pipeline (2019-2024)...")
start_time <- Sys.time()

domat_list <- future_map(years, fetch_and_process_domat, .options = furrr_options(seed = TRUE))
combined_domat <- bind_rows(domat_list)

consolidated_path <- file.path(output_dir, "DOMAT_2019_2024_combined.parquet")
write_parquet(combined_domat, consolidated_path)

end_time <- Sys.time()
elapsed <- round(as.numeric(end_time - start_time, units = "secs"), 2)

message("==========================================================")
message(sprintf("Execution finished in %s seconds.", elapsed))
message(sprintf("Consolidated dataset: %d total records.", nrow(combined_domat)))
message(sprintf("Master file saved to: %s", consolidated_path))
message("==========================================================")
