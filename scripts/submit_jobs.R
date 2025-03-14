#!/usr/bin/env Rscript

# Submit Energy Budget Jobs
# This script generates and submits Slurm jobs for energy budget calculations

# Load required packages
library(tidyverse)

# Define job parameters - focus on MS only
species_list <- c("MS")  # Only MS for now
sexes <- c("F", "M")

# Define year ranges
historical_years <- c(1945:1964)
recent_years <- c(2005:2024)

# Sites for MS (they only occur at these sites)
ms_sites <- c("Eldo", "A1", "B1")

# Base output directory
output_base_dir <- "output/results"

# Read the Slurm template
template <- readLines("slurm/energy_budget.sh")

# Define function to create job script
create_job_script <- function(species, year, site_orig, site_clim, sex) {
  # Create job name
  job_name <- sprintf("%s_%d_%s_%s_%s", 
                     species, year, site_orig, site_clim, sex)
  script_name <- sprintf("temp_job_%s.sh", job_name)
  
  # Create output directory structure
  output_dir <- file.path(output_base_dir, sex, site_orig, format(year, digits=4))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create job script by substituting parameters
  job_script <- template
  job_script <- gsub("\\{JOB_NAME\\}", job_name, job_script)
  job_script <- gsub("\\{SPECIES\\}", species, job_script)
  job_script <- gsub("\\{YEAR\\}", year, job_script)
  job_script <- gsub("\\{SITE_ORIG\\}", site_orig, job_script)
  job_script <- gsub("\\{SITE_CLIM\\}", site_clim, job_script)
  job_script <- gsub("\\{SEX\\}", sex, job_script)
  job_script <- gsub("\\{OUTPUT_DIR\\}", output_dir, job_script)
  
  return(list(script = job_script, filename = script_name, job_name = job_name))
}

# Create output directory if it doesn't exist
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("output/logs", recursive = TRUE, showWarnings = FALSE)

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

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
dry_run <- length(args) > 0 && args[1] == "--dry-run"

# Handle additional options
run_historical <- TRUE
run_recent <- TRUE
sample_years <- FALSE

if (length(args) > 1) {
  if ("--historical-only" %in% args) {
    run_recent <- FALSE
  } else if ("--recent-only" %in% args) {
    run_historical <- FALSE
  }
  
  if ("--sample-years" %in% args) {
    sample_years <- TRUE
    # Use a subset of years for testing
    historical_years <- c(1945, 1950, 1955, 1960, 1964)
    recent_years <- c(2005, 2010, 2015, 2020, 2024) 
  }
  
  # Allow specification of specific years
  year_arg <- args[grep("^--years=", args)]
  if (length(year_arg) > 0) {
    years_str <- sub("^--years=", "", year_arg)
    year_list <- as.numeric(strsplit(years_str, ",")[[1]])
    
    # Filter historical and recent years
    historical_years <- historical_years[historical_years %in% year_list]
    recent_years <- recent_years[recent_years %in% year_list]
  }
}

# Determine which years to use
years_to_run <- c()
if (run_historical) years_to_run <- c(years_to_run, historical_years)
if (run_recent) years_to_run <- c(years_to_run, recent_years)

# Generate job scripts and optionally submit them
job_count <- 0
job_files <- list()

for (species in species_list) {
  # Get valid sites for this species
  valid_sites <- ms_sites  # For MS, use the MS sites list
  
  for (year in years_to_run) {
    for (site_orig in valid_sites) {
      for (site_clim in valid_sites) { # Only consider valid climate sites
        for (sex in sexes) {
          script_info <- create_job_script(species, year, site_orig, site_clim, sex)
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
  
  writeLines(submit_script, "submit_all_jobs.sh")
  system("chmod +x submit_all_jobs.sh")
  cat("Created submission script: submit_all_jobs.sh\n")
}
