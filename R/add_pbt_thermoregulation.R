#!/usr/bin/env Rscript

# Add PBT Thermoregulation to Energy Budget Results
# This script reads existing energy budget results and adds a PBT (Preferred Body Temperature)
# thermoregulation condition

library(tidyverse)

# Define PBT values for each species
PBT_VALUES <- list(
  MB = 32.82804878,
  MS = 30.62596154 
)

#' Add PBT thermoregulation to energy budget results
#' 
#' @param results_file Path to the results RDS file
#' @param species Species code (MB or MS)
#' @param pops Population data frame
#' @param output_file Optional output file path. If NULL, overwrites input file
#' @return Modified results with PBT thermoregulation added
add_pbt_thermoregulation <- function(results_file, species, pops, output_file = NULL) {
  
  # Load the results
  cat(sprintf("Loading results from %s\n", results_file))
  results <- readRDS(results_file)
  
  # Get PBT for this species
  pbt <- PBT_VALUES[[species]]
  if (is.null(pbt)) {
    stop(sprintf("PBT not defined for species %s", species))
  }
  
  cat(sprintf("Using PBT = %.2f for species %s\n", pbt, species))
  
  # Get unique metadata from the results
  metadata <- results %>%
    select(species, sex, site_orig, site_clim) %>%
    distinct()
  
  if (nrow(metadata) != 1) {
    stop("Results file contains multiple species/sex/site combinations")
  }
  
  # Calculate gains and losses at PBT
  cat("Calculating energy at PBT...\n")
  
  # Get population parameters
  pop_dat <- pops %>%
    filter(spp == species & site == metadata$site_orig & sex == metadata$sex) %>%
    slice(1)
  
  if (nrow(pop_dat) == 0) {
    stop(sprintf("No population data found for %s %s %s", 
                 species, metadata$site_orig, metadata$sex))
  }
  
  # Calculate energy at PBT using the same functions as in energy_functions.R
  # Extract parameters
  Tmin <- pop_dat$Tmin[1]
  Topt <- pop_dat$Topt[1]
  Above <- pop_dat$Above[1]
  Tmax <- Topt + Above
  Ropt <- pop_dat$Ropt[1]
  mass <- pop_dat$mass[1]
  
  # Get sex effect coefficient
  sex_effect <- ifelse(species == "MS", -0.415, 0.123)
  
  # Calculate sex multiplier
  sex_multiplier <- if(metadata$sex == "M") {
    exp(sex_effect * 0.5)
  } else {
    exp(sex_effect * -0.5)
  }
  
  # Calculate feeding rate at PBT
  base_rate <- lrf_1991_correct(pbt, Ropt, Topt, Tmin, Tmax)
  adjusted_rate <- base_rate * sex_multiplier * mass
  
  # Convert to energy gain (same conversion as in get_energy_gains)
  dry_fec_mg_hr <- adjusted_rate
  dry_wga_mg_hr <- dry_fec_mg_hr * 12.8/14.8 * 0.40/(1-0.40)
  kcal_hr <- dry_wga_mg_hr * 14.8/1000
  gain_at_pbt <- kcal_hr * 4.184  # kJ/hr
  
  # Calculate metabolic rate at PBT
  mass <- pop_dat$mass[1]
  elev <- pop_dat$elev[1]
  rmr_b0 <- pop_dat$rmr_b0[1]
  rmr_b1 <- pop_dat$rmr_b1[1]
  rmr_b2 <- pop_dat$rmr_b2[1]
  rmr_b3 <- pop_dat$rmr_b3[1]
  
  mrs_at_pbt <- get_mrs(pbt, mass, elev, rmr_b0, rmr_b1, rmr_b2, rmr_b3)
  loss_at_pbt <- get_mr_losses(mrs_at_pbt)
  net_gain_at_pbt <- gain_at_pbt - loss_at_pbt
  
  cat(sprintf("At PBT: gain=%.4f, loss=%.4f, net=%.4f kJ/hr\n", 
              gain_at_pbt, loss_at_pbt, net_gain_at_pbt))
  
  # Process results to add PBT thermoregulation
  cat("Adding PBT thermoregulation condition...\n")
  
  # Group by timestamp to find available microclimates
  results_pbt <- results %>%
    group_by(dtuse) %>%
    mutate(
      # Find if PBT is within range of available Tbs
      tb_min = min(Tb, na.rm = TRUE),
      tb_max = max(Tb, na.rm = TRUE),
      pbt_in_range = pbt >= tb_min & pbt <= tb_max,
      
      # Find closest Tb to PBT
      tb_diff = abs(Tb - pbt),
      closest_to_pbt = tb_diff == min(tb_diff, na.rm = TRUE)
    ) %>%
    group_by(dtuse) %>%
    mutate(
      # Select the microclimate for PBT condition
      # If PBT is in range, use PBT values; otherwise use closest Tb values
      gains_pbt = ifelse(pbt_in_range[1], gain_at_pbt, gains[closest_to_pbt][1]),
      losses_pbt = ifelse(pbt_in_range[1], loss_at_pbt, losses[closest_to_pbt][1]),
      net_gains_pbt = ifelse(pbt_in_range[1], net_gain_at_pbt, net_gains[closest_to_pbt][1]),
      tb_pbt = ifelse(pbt_in_range[1], pbt, Tb[closest_to_pbt][1]),
      
      # Record which microclimate was selected
      selected_height = height[closest_to_pbt][1],
      selected_shade = shade[closest_to_pbt][1]
    ) %>%
    ungroup()
  
  # Add a summary row for each timestamp with PBT thermoregulation
  pbt_summary <- results_pbt %>%
    group_by(dtuse) %>%
    slice(1) %>%  # Take one row per timestamp as template
    mutate(
      # Override with PBT values
      height = NA,  # Indicates this is the PBT row
      shade = NA,
      condition = "PBT",
      Tb = tb_pbt,
      gains = gains_pbt,
      losses = losses_pbt,
      net_gains = net_gains_pbt,
      
      # Keep track of what was selected
      pbt_value = pbt,
      pbt_was_in_range = pbt_in_range,
      pbt_selected_height = selected_height,
      pbt_selected_shade = selected_shade
    ) %>%
    select(-c(tb_min, tb_max, pbt_in_range, tb_diff, closest_to_pbt, 
              gains_pbt, losses_pbt, net_gains_pbt, tb_pbt,
              selected_height, selected_shade))
  
  # Combine original results with PBT rows
  results_with_pbt <- bind_rows(
    results %>% mutate(condition = "microclimate"),
    pbt_summary
  )
  
  # Sort by timestamp and condition
  results_with_pbt <- results_with_pbt %>%
    arrange(dtuse, condition, height, shade)
  
  # Save results
  if (is.null(output_file)) {
    output_file <- results_file
  }
  
  saveRDS(results_with_pbt, output_file)
  cat(sprintf("Results with PBT thermoregulation saved to %s\n", output_file))
  
  # Print summary statistics
  pbt_stats <- pbt_summary %>%
    summarize(
      total_hours = n(),
      hours_pbt_achievable = sum(pbt_was_in_range),
      pct_achievable = hours_pbt_achievable / total_hours * 100,
      total_net_gain = sum(net_gains),
      mean_tb = mean(Tb),
      mean_daily_gain = total_net_gain / n_distinct(as.Date(dtuse))
    )
  
  cat("\nPBT Thermoregulation Summary:\n")
  cat(sprintf("  Total hours: %d\n", pbt_stats$total_hours))
  cat(sprintf("  Hours PBT achievable: %d (%.1f%%)\n", 
              pbt_stats$hours_pbt_achievable, pbt_stats$pct_achievable))
  cat(sprintf("  Total net energy gain: %.1f kJ\n", pbt_stats$total_net_gain))
  cat(sprintf("  Mean body temperature: %.1f°C\n", pbt_stats$mean_tb))
  cat(sprintf("  Mean daily energy gain: %.1f kJ/day\n", pbt_stats$mean_daily_gain))
  
  return(results_with_pbt)
}

# Function to process all results in a directory
process_directory <- function(results_dir, pops, pattern = "eb_results_.*\\.rds$") {
  # Find all result files
  result_files <- list.files(results_dir, pattern = pattern, 
                             recursive = TRUE, full.names = TRUE)
  
  cat(sprintf("Found %d result files to process\n", length(result_files)))
  
  # Process each file
  for (file in result_files) {
    # Extract species from filename
    filename <- basename(file)
    species <- str_extract(filename, "MB|MS")
    
    if (is.na(species)) {
      warning(sprintf("Could not extract species from filename: %s", filename))
      next
    }
    
    cat(sprintf("\nProcessing %s (species: %s)\n", filename, species))
    
    tryCatch({
      add_pbt_thermoregulation(file, species, pops)
    }, error = function(e) {
      warning(sprintf("Error processing %s: %s", file, e$message))
    })
  }
}

# Main execution when run as script
if (!interactive()) {
  # Load required functions
  source("R/energy_functions.R")
  
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript add_pbt_thermoregulation.R <results_file_or_directory> [pops_file]\n")
    cat("  results_file_or_directory: Path to a single .rds file or directory containing results\n")
    cat("  pops_file: Path to pops.rds file (default: data/pops.rds)\n")
    quit(status = 1)
  }
  
  input_path <- args[1]
  pops_file <- ifelse(length(args) > 1, args[2], "data/pops.rds")
  
  # Load population data
  cat(sprintf("Loading population data from %s\n", pops_file))
  pops <- readRDS(pops_file)
  
  # Process input
  if (dir.exists(input_path)) {
    # Process all files in directory
    process_directory(input_path, pops)
  } else if (file.exists(input_path)) {
    # Process single file
    # Extract species from filename
    species <- str_extract(basename(input_path), "MB|MS")
    if (is.na(species)) {
      stop(sprintf("Could not extract species from filename: %s", input_path))
    }
    add_pbt_thermoregulation(input_path, species, pops)
  } else {
    stop(sprintf("Input path does not exist: %s", input_path))
  }
}