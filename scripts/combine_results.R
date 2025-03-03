#!/usr/bin/env Rscript

# Combine Energy Budget Results
# This script combines the results from individual energy budget jobs

# Load required packages
library(tidyverse)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
results_dir <- ifelse(length(args) > 0, args[1], "output/results")
output_file <- ifelse(length(args) > 1, args[2], "output/combined_results.rds")

# Check if results directory exists
if (!dir.exists(results_dir)) {
  stop(paste("Results directory not found:", results_dir))
}

# Get all result files
result_files <- list.files(results_dir, pattern = "eb_results_.*\\.rds", full.names = TRUE)

if (length(result_files) == 0) {
  stop(paste("No result files found in", results_dir))
}

cat(sprintf("Found %d result files to combine\n", length(result_files)))

# Initialize empty dataframe to store combined results
combined_results <- NULL

# Process each file
for (i in seq_along(result_files)) {
  file_path <- result_files[i]
  
  # Load result
  tryCatch({
    result <- readRDS(file_path)
    
    # Add to combined results
    if (is.null(combined_results)) {
      combined_results <- result
    } else {
      combined_results <- bind_rows(combined_results, result)
    }
    
    # Print progress
    if (i %% 100 == 0 || i == length(result_files)) {
      cat(sprintf("Processed %d/%d files\n", i, length(result_files)))
    }
  }, error = function(e) {
    warning(sprintf("Error processing file %s: %s", file_path, e$message))
  })
}

# Create output directory if it doesn't exist
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Save combined results
saveRDS(combined_results, output_file)
cat(sprintf("Combined results saved to %s\n", output_file))

# Optional: Create a summary of the results
summary_file <- file.path(dirname(output_file), "results_summary.csv")
summary_data <- combined_results %>%
  group_by(species, year, site_orig, site_clim, height, shade) %>%
  summarize(
    mean_tb = mean(Tb, na.rm = TRUE),
    max_tb = max(Tb, na.rm = TRUE),
    min_tb = min(Tb, na.rm = TRUE),
    total_net_gain = sum(net_gains, na.rm = TRUE),
    total_gains = sum(gains, na.rm = TRUE),
    total_losses = sum(losses, na.rm = TRUE),
    n_records = n(),
    .groups = "drop"
  )

# Save summary
write_csv(summary_data, summary_file)
cat(sprintf("Summary saved to %s\n", summary_file))