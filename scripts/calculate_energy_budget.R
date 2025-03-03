#!/usr/bin/env Rscript

# Calculate Energy Budget
# This script calculates body temperatures and energy budgets
# for grasshoppers based on climate data and population parameters

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
species <- args[1]         # e.g., "MB"
year <- as.numeric(args[2]) # e.g., 2011
site_orig <- args[3]       # Original site, e.g., "A1"
site_clim <- args[4]       # Climate site, e.g., "B1" 
height <- as.numeric(args[5])    # Height above ground, e.g., 0.03
shade <- as.numeric(args[6])     # Shade level, e.g., 0.1
output_dir <- args[7]      # Output directory

# Set working directory to the script location
script_dir <- dirname(sys.frame(1)$ofile)
setwd(file.path(script_dir, ".."))

# Load functions
source("R/setup.R")
setup_packages()
source("R/biophysical_functions.R")
source("R/energy_functions.R") 
source("R/data_functions.R")

# Load population data
pops <- load_pops_data()

# Get surface roughness
surface_roughness <- get_surface_roughness()

# Load climate data
climate_data <- load_climate_data(site_clim, year)

cat(sprintf("Processing: Species=%s, Year=%d, Origin=%s, Climate=%s, Height=%.2f, Shade=%.2f\n",
            species, year, site_orig, site_clim, height, shade))

# Calculate energy budget
result <- calculate_energy_budget(
  climate_data, height, shade, surface_roughness, pops,
  species, "F", site_orig, site_clim, year
)

# Save results with metadata
job_id <- Sys.getenv("SLURM_JOB_ID")
if (job_id == "") job_id <- format(Sys.time(), "%Y%m%d%H%M%S")

metadata <- list(
  species = species,
  year = year,
  site_orig = site_orig,
  site_clim = site_clim,
  height = height,
  shade = shade
)

output_file <- save_results(result, output_dir, job_id, metadata)

cat(sprintf("Processing complete. Results saved to %s\n", output_file))