library(lubridate)
library(stringr)
library(readr)

relabel_mst <- function(path) {
  yr    <- as.integer(str_extract(basename(path), "\\d{4}"))   # basename: robust year extraction
  start <- as.POSIXct(sprintf("%d-06-01 17:00:00", yr), tz = "Etc/GMT+7")  # explicit tz, no DST
  end   <- as.POSIXct(sprintf("%d-08-31 16:00:00", yr), tz = "Etc/GMT+7")
  df <- readRDS(path)
  df$dtuse   <- with_tz(df$dtuse, "Etc/GMT+7")
  df$hour    <- hour(df$dtuse); df$doy <- yday(df$dtuse)
  df$month   <- month(df$dtuse); df$day <- day(df$dtuse)
  df$dt_noyr <- df$dtuse; year(df$dt_noyr) <- 2024
  df %>% filter(dtuse >= start, dtuse <= end)
}


process_weather_files <- function(root_dir = "C:/Users/jmsmi/era5_workflow_new/micro2") {
  
  rds_files <- list.files(
    root_dir,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  for (rds_file in rds_files) {
    
    cat("Processing:", rds_file, "\n")
    
    csv_file <- sub("\\.rds$", ".csv", rds_file)
    
    # Backup original RDS
    rds_backup <- sub("\\.rds$", "_original.rds", rds_file)
    if (!file.exists(rds_backup)) {
      file.rename(rds_file, rds_backup)
    }
    
    # Backup original CSV
    csv_backup <- sub("\\.csv$", "_original.csv", csv_file)
    if (file.exists(csv_file) && !file.exists(csv_backup)) {
      file.rename(csv_file, csv_backup)
    }
    
    # Create updated dataset
    #
    # Option 1: relabel_mst() returns the modified object
    dat <- relabel_mst(rds_backup)
    
    saveRDS(dat, rds_file)
    
    # Recreate CSV from updated RDS
    write_csv(dat, csv_file)
    
    # Option 2 (if relabel_mst writes the file itself):
    # relabel_mst(rds_backup)
    # dat <- readRDS(rds_file)
    # write_csv(dat, csv_file)
    
  }
  
  invisible(rds_files)
}




restore_originals <- function(root_dir = "C:/Users/jmsmi/era5_workflow_new/micro2") {
  
  original_rds <- list.files(
    root_dir,
    pattern = "_original\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  for (orig_rds in original_rds) {
    
    current_rds <- sub("_original\\.rds$", ".rds", orig_rds)
    
    if (file.exists(current_rds)) {
      file.remove(current_rds)
    }
    
    file.rename(orig_rds, current_rds)
    
    orig_csv <- sub("_original\\.rds$", "_original.csv", orig_rds)
    current_csv <- sub("_original\\.rds$", ".csv", orig_rds)
    
    if (file.exists(orig_csv)) {
      
      if (file.exists(current_csv)) {
        file.remove(current_csv)
      }
      
      file.rename(orig_csv, current_csv)
    }
    
    cat("Restored:", basename(current_rds), "\n")
  }
  
  invisible(NULL)
}
