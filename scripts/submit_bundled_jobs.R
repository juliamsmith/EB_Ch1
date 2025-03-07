#!/usr/bin/env Rscript

# Submit Bundled Energy Budget Jobs
# This script generates and submits Slurm jobs for energy budget calculations
# with years bundled together

# Load required packages
library(tidyverse)

# Parse command line arguments for species
args <- commandArgs(trailingOnly = TRUE)

# Define species to process
if (length(args) > 0 && args[1] %in% c("MS", "MB", "both")) {
  if (args[1] == "both") {
    species_list <- c("MS", "MB")
  } else {
    species_list <- args[1]
  }
} else {
  # Default to both species
  species_list <- c("MS", "MB")
}

# Get remaining arguments
args <- args[-1]

sexes <- c("F", "M")

# Define year ranges - put them in bundles
historical_years <- "1945,1946,1947,1948,1949,1950,1951,1952,1953,1954,1955,1956,1957,1958,1959,1960,1961,1962,1963,1964"
recent_years <- "2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024"

# Define valid sites for each species
ms_sites <- c("Eldo", "A1", "B1")
mb_sites <- c("A1", "B1", "C1", "D1")

# Base output directory
output_base_dir <- "output/results"

# Read the Slurm template
template <- readLines("slurm/energy_bundle.sh")

# Define function to create job script
create_job_script <- function(species, years, site_orig, site_clim, sex) {
  # Create job name
  period <- ifelse(startsWith(years, "19"), "hist", "recent")
  job_name <- sprintf("%s_%s_%s_%s_%s", 
                     species, period, site_orig, site_clim, sex)
  script_name <- sprintf("bundle_job_%s.sh", job_name)
  
  # Create base output directory
  output_dir <- file.path(output_base_dir, site_orig)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create job script by substituting parameters
  job_script <- template
  job_script <- gsub("\\{JOB_NAME\\}", job_name, job_script)
  job_script <- gsub("\\{SPECIES\\}", species, job_script)
  job_script <- gsub("\\{YEARS\\}", years, job_script)
  job_script <- gsub("\\{SITE_ORIG\\}", site_orig, job_script)
  job_script <- gsub("\\{SITE_CLIM\\}", site_clim, job_script)
  job_script <- gsub("\\{SEX\\}", sex, job_script)
  job_script <- gsub("\\{OUTPUT_DIR\\}", output_dir, job_script)
  
  return(list(script = job_script, filename = script_name, job_name = job_name))
}

# Generate and submit jobs
submit_job <- function(script_info, submit = TRUE) {
  # Write the script to a file
  writeLines(script_info$script, script_info$filename)
  
  # Make it executable
  system(sprintf("chmod +x %s", script_info$filename))
  
  # Submit the job if requested
  if (submit) {
    cat(sprintf("Submitting job: %s\n", script_info$filename))
    system(sprintf("sbatch %s", script_info$filename))
  } else {
    cat(sprintf("Created job script: %s\n", script_info$filename))
  }
  
  return(script_info$filename)
}

# Check for dry run option
dry_run <- "--dry-run" %in% args

# Handle additional options
run_historical <- TRUE
run_recent <- TRUE

if ("--historical-only" %in% args) {
  run_recent <- FALSE
}
if ("--recent-only" %in% args) {
  run_historical <- FALSE
}

# Determine which year bundles to use
year_bundles <- c()
if (run_historical) year_bundles <- c(year_bundles, historical_years)
if (run_recent) year_bundles <- c(year_bundles, recent_years)

# Generate job scripts and optionally submit them
job_count <- 0
job_files <- list()

for (species in species_list) {
  # Get valid sites for this species
  valid_sites <- if (species == "MS") ms_sites else mb_sites
  
  for (years in year_bundles) {
    for (site_orig in valid_sites) {
      for (site_clim in valid_sites) { # Only consider valid climate sites
        for (sex in sexes) {
          script_info <- create_job_script(species, years, site_orig, site_clim, sex)
          job_files[[length(job_files) + 1]] <- submit_job(script_info, !dry_run)
          job_count <- job_count + 1
        }
      }
    }
  }
}

cat(sprintf("%s %d jobs.\n", 
           ifelse(dry_run, "Created", "Submitted"), 
           job_count))

# If it's a dry run, create a submission script
if (dry_run) {
  submit_script <- "#!/bin/bash\n\n# Submit all job scripts\n"
  for (job_file in job_files) {
    submit_script <- paste0(submit_script, sprintf("sbatch %s\n", job_file))
  }
  
  writeLines(submit_script, sprintf("submit_all_%s_bundles.sh", 
                                   paste(species_list, collapse="_")))
  system(sprintf("chmod +x submit_all_%s_bundles.sh", 
                paste(species_list, collapse="_")))
  cat(sprintf("Created submission script: submit_all_%s_bundles.sh\n", 
             paste(species_list, collapse="_")))
}
