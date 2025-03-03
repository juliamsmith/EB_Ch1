# R/data_functions.R

# Load climate data from RDS files
load_climate_data <- function(site, year, base_path = "../multi_microclimate") {
  # Construct file path
  file_path <- file.path(base_path, site, sprintf("climateuse_combined_%s_%d_6_7_8.rds", site, year))
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop(paste("Climate data file not found:", file_path))
  }
  
  # Load the data
  climate_data <- readRDS(file_path)
  
  return(climate_data)
}

# Load population data
load_pops_data <- function(file_path = "data/pops.rds") {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop(paste("Population data file not found:", file_path))
  }
  
  # Load the data
  pops <- readRDS(file_path)
  
  return(pops)
}

# Get surface roughness data
get_surface_roughness <- function(file_path = "data/surface_roughness.rds") {
  if (file.exists(file_path)) {
    return(readRDS(file_path))
  } else {
    # Default value based on run_era_Eb.R
    warning("Surface roughness file not found, using default value")
    return(0.001) # Default value
  }
}

# Save results to RDS file
save_results <- function(results, output_path, job_id, metadata) {
  # Create directory if it doesn't exist
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Create filename with relevant information
  filename <- file.path(output_path, 
                        sprintf("eb_results_%s_%s_%s_%s_%s_%s_%s.rds", 
                                metadata$species, metadata$year,
                                metadata$site_orig, metadata$site_clim,
                                metadata$height, metadata$shade,
                                job_id))
  
  # Add metadata to results
  results <- results %>%
    mutate(job_id = job_id)
  
  # Save the results
  saveRDS(results, filename)
  
  return(filename)
}