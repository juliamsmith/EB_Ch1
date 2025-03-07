#!/usr/bin/env Rscript

# Calculate Energy Budget
# This script calculates body temperatures and energy budgets
# for grasshoppers based on climate data and population parameters

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
species <- args[1]          # e.g., "MB"
year <- as.numeric(args[2]) # e.g., 2022
site_orig <- args[3]        # Original site, e.g., "A1"
site_clim <- args[4]        # Climate site, e.g., "B1" 
sex <- args[5]              # Sex, e.g., "F" or "M"
output_dir <- args[6]       # Output directory
microclim_path <- if(length(args) > 6) args[7] else "/mmfs1/gscratch/biology/jmsmith/targeted_microclimate"

# Print command line args for debugging
cat("Command line arguments:\n")
print(args)
cat(sprintf("microclim_path: %s\n", microclim_path))

# Load functions
source("R/setup.R")
setup_packages()
source("R/biophysical_functions.R")
source("R/energy_functions.R") 
source("R/data_functions.R")

# Load population data
pops <- load_pops_data()

# Load climate data for the entire period
climate_data <- load_climate_data(site_clim, year, microclim_path)

# No filtering - use all data in the file
cat(sprintf("Processing: Species=%s, Year=%d, Origin=%s, Climate=%s, Sex=%s\n",
            species, year, site_orig, site_clim, sex))
cat(sprintf("Time period: %s to %s\n",
            format(min(climate_data$dtuse), "%Y-%m-%d"),
            format(max(climate_data$dtuse), "%Y-%m-%d")))

# Calculate energy budget
result <- calculate_energy_budget(
  climate_data, pops,
  species, sex, site_orig, site_clim, year
)

# Save results with metadata
job_id <- Sys.getenv("SLURM_JOB_ID")
if (job_id == "") job_id <- format(Sys.time(), "%Y%m%d%H%M%S")

metadata <- list(
  species = species,
  year = year,
  site_orig = site_orig,
  site_clim = site_clim, 
  sex = sex
)

output_file <- save_results(result, output_dir, job_id, metadata)

cat(sprintf("Processing complete. Results saved to %s\n", output_file))
