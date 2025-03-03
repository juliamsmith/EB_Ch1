#!/usr/bin/env Rscript

# Submit Energy Budget Jobs
# This script generates and submits Slurm jobs for energy budget calculations

# Load required packages
library(tidyverse)

# Define job parameters
species_list <- c("MB", "MS")
years <- c(2011, 2022, 2023)  # Adjust based on available data
site_list <- c("Eldo", "A1", "B1", "C1")
heights <- c(0.03, 0.1, 0.3)  # Different heights above ground
shade_levels <- c(0, 0.25, 0.5, 0.75, 0.9)  # Different shade levels
sexes <- c("F", "M")  # Add both sexes
output_dir <- "output/results"

# Read the Slurm template
template <- readLines("slurm/energy_budget.sh")

# Define function to create job script
create_job_script <- function(species, year, site_orig, site_clim, height, shade) {
  # Create job name and filename
  job_name <- sprintf("%s_%d_%s_%s_%.2f_%.2f", 
                      species, year, site_orig, site_clim, height, shade)
  script_name <- sprintf("temp_job_%s.sh", gsub("\\.", "", job_name))
  
  # Create job script by substituting parameters
  job_script <- template
  job_script <- gsub("\\{JOB_NAME\\}", job_name, job_script)
  job_script <- gsub("\\{SPECIES\\}", species, job_script)
  job_script <- gsub("\\{YEAR\\}", year, job_script)
  job_script <- gsub("\\{SITE_ORIG\\}", site_orig, job_script)
  job_script <- gsub("\\{SITE_CLIM\\}", site_clim, job_script) 
  job_script <- gsub("\\{HEIGHT\\}", height, job_script)
  job_script <- gsub("\\{SHADE\\}", shade, job_script)
  job_script <- gsub("\\{OUTPUT_DIR\\}", output_dir, job_script)
  
  return(list(script = job_script, filename = script_name))
}

# Create output directory if it doesn't exist
dir.create("output/logs", recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

# Parse command line arguments for dry run
args <- commandArgs(trailingOnly = TRUE)
dry_run <- length(args) > 0 && args[1] == "--dry-run"

# Generate job scripts and optionally submit them
job_count <- 0
job_files <- list()

for (species in species_list) {
  for (year in years) {
    for (site_orig in site_list) {
      for (site_clim in site_list) {
        for (height in heights) {
          for (shade in shade_levels) {
            script_info <- create_job_script(species, year, site_orig, site_clim, height, shade)
            job_files[[length(job_files) + 1]] <- submit_job(script_info, !dry_run)
            job_count <- job_count + 1
          }
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