#!/usr/bin/env Rscript

# Process multiple years in one job
args <- commandArgs(trailingOnly = TRUE)
species <- args[1]              # e.g., "MS"
years <- as.numeric(strsplit(args[2], ",")[[1]])  # e.g., "1945,1946,1947"
site_orig <- args[3]            # Original site, e.g., "A1"
site_clim <- args[4]            # Climate site, e.g., "B1" 
sex <- args[5]                  # Sex, e.g., "F" or "M"
output_dir <- args[6]           # Output directory
microclim_path <- if(length(args) > 6) args[7] else "/mmfs1/gscratch/biology/jmsmith/targeted_microclimate3"

# Print arguments for debugging
cat("Arguments received:\n")
cat(sprintf("Species: %s\n", species))
cat(sprintf("Years: %s\n", paste(years, collapse=", ")))
cat(sprintf("Origin site: %s\n", site_orig))
cat(sprintf("Climate site: %s\n", site_clim))
cat(sprintf("Sex: %s\n", sex))
cat(sprintf("Output directory: %s\n", output_dir))

# Load functions
source("R/setup.R")
setup_packages()
source("R/biophysical_functions.R")
source("R/energy_functions.R") 
source("R/data_functions.R")

# Load population data
pops <- load_pops_data()

# Process all years
for (year in years) {
  # Create year/sex output directory
  year_output_dir <- file.path(output_dir, as.character(year), sex)
  dir.create(year_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("\nProcessing: Species=%s, Year=%d, Origin=%s, Climate=%s, Sex=%s\n",
              species, year, site_orig, site_clim, sex))
  
  # Try to load climate data
  tryCatch({
    # Load climate data
    climate_data <- load_climate_data(site_clim, year, microclim_path)
    
    # Calculate energy budget
    result <- calculate_energy_budget(
      climate_data, pops,
      species, sex, site_orig, site_clim, year
    )
    
    # Save results
    job_id <- Sys.getenv("SLURM_JOB_ID")
    if (job_id == "") job_id <- format(Sys.time(), "%Y%m%d%H%M%S")
    
    metadata <- list(
      species = species,
      year = year,
      site_orig = site_orig,
      site_clim = site_clim, 
      sex = sex
    )
    
    output_file <- save_results(result, year_output_dir, job_id, metadata)
    cat(sprintf("Saved results to %s\n", output_file))
  }, error = function(e) {
    cat(sprintf("Error processing year %d: %s\n", year, e$message))
  })
}

cat("\nAll years processed.\n")
