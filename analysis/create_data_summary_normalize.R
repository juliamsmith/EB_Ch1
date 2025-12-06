library(tidyverse)

# Lobry-Rosso-Flandrois TPC function
lrf_1991_exact <- function(temp, rmax, topt, tmin, tmax) {
  if (temp <= tmin || temp >= tmax) {
    return(0)
  }
  
  numerator <- (temp - tmax) * (temp - tmin)^2
  denominator <- (topt - tmin) * ((topt - tmin) * (temp - topt) - (topt - tmax) * (topt + tmin - 2 * temp))
  
  rate <- rmax * numerator / denominator
  
  if (is.na(rate) || rate < 0) {
    return(0)
  }
  
  return(rate)
}

# Function to calculate area under TPC curve using LRF model
calculate_tpc_area <- function(Tmin, Topt, Above, Ropt) {
  Tmax <- Topt + Above
  T_range <- seq(Tmin, Tmax, by = 0.1)
  
  tpc_values <- sapply(T_range, function(T) {
    lrf_1991_exact(temp = T, rmax = Ropt, topt = Topt, tmin = Tmin, tmax = Tmax)
  })
  
  area <- sum(diff(T_range) * (head(tpc_values, -1) + tail(tpc_values, -1)) / 2)
  
  if (is.na(area) || area == 0 || is.infinite(area)) {
    warning(paste("Invalid TPC area:", area, 
                  "- Parameters: Tmin=", Tmin, "Topt=", Topt, 
                  "Above=", Above, "Ropt=", Ropt))
    return(NA)
  }
  
  return(area)
}

# Function to calculate species-specific average TPC areas
calculate_species_avg_tpc <- function(pops_df) {
  species_avgs <- pops_df %>%
    rowwise() %>%
    mutate(tpc_area = calculate_tpc_area(Tmin, Topt, Above, Ropt)) %>%
    ungroup() %>%
    group_by(spp) %>%
    summarize(avg_tpc_area = mean(tpc_area, na.rm = TRUE), .groups = "drop")
  
  return(species_avgs)
}

# Function to extract info from file path
extract_filename_info <- function(filepath) {
  filename <- basename(filepath)
  parts <- strsplit(tools::file_path_sans_ext(filename), "_")[[1]]
  
  species <- parts[3]
  year <- parts[4]
  site_orig <- parts[5]
  site_clim <- parts[6]
  
  dir_parts <- strsplit(dirname(filepath), "/")[[1]]
  sex <- dir_parts[length(dir_parts)]
  
  year_period <- case_when(
    as.numeric(year) < 1965 ~ "historical",
    as.numeric(year) > 2000 ~ "contemporary",
    TRUE ~ "other"
  )
  
  return(list(
    species = species,
    year = year,
    site_orig = site_orig,
    site_clim = site_clim,
    sex = sex,
    year_period = year_period
  ))
}

# Modified function to process a single RDS file with TPC area normalization
process_eb_file <- function(filepath, pops_df, species_avg_tpc) {
  file_info <- extract_filename_info(filepath)
  
  tryCatch({
    data <- readRDS(filepath)
    
    # Get TPC parameters for this species/site/sex combination
    tpc_params <- pops_df %>%
      filter(spp == file_info$species,
             site == file_info$site_orig,
             sex == file_info$sex)
    
    if (nrow(tpc_params) == 0) {
      warning(paste("No TPC parameters found for:", 
                    file_info$species, file_info$site_orig, file_info$sex))
      return(NULL)
    }
    
    # Calculate area under TPC for this population
    tpc_area <- calculate_tpc_area(
      Tmin = tpc_params$Tmin[1],
      Topt = tpc_params$Topt[1],
      Above = tpc_params$Above[1],
      Ropt = tpc_params$Ropt[1]
    )
    
    if (is.na(tpc_area)) {
      warning(paste("Skipping file due to invalid TPC area:", filepath))
      return(NULL)
    }
    
    # Get species average TPC area
    species_avg <- species_avg_tpc %>%
      filter(spp == file_info$species) %>%
      pull(avg_tpc_area)
    
    if (length(species_avg) == 0) {
      warning(paste("No species average TPC found for:", file_info$species))
      return(NULL)
    }
    
    # Normalize: divide by population TPC area, multiply by species average
    # This removes population variation but preserves species-level scaling
    normalization_factor <- species_avg / tpc_area
    
    data <- data %>%
      mutate(
        gains_normalized = gains * normalization_factor,
        net_gains_normalized = gains_normalized - losses
      )
    
    # Summarize with normalized values
    summary <- data %>%
      group_by(shade, height) %>%
      summarize(
        total_net_gains = sum(net_gains_normalized, na.rm = TRUE),
        avg_gains = mean(gains_normalized, na.rm = TRUE),
        avg_losses = mean(losses, na.rm = TRUE),
        hours = n(),
        .groups = "drop"
      ) %>%
      mutate(
        species = file_info$species,
        year = file_info$year,
        site_orig = file_info$site_orig,
        site_clim = file_info$site_clim,
        sex = file_info$sex,
        year_period = file_info$year_period,
        filepath = filepath
      )
    
    return(summary)
  }, error = function(e) {
    warning(paste("Error processing file:", filepath, "\nError:", e$message))
    return(NULL)
  })
}

# Function to create a wide format of the results
create_wide_format <- function(summary_df) {
  summary_df <- summary_df %>%
    mutate(condition = paste0(
      case_when(
        shade == 0 ~ "sun",
        shade == 0.3 ~ "partial_shade",
        shade == 0.9 ~ "full_shade",
        TRUE ~ paste0("shade_", shade)
      ),
      "_",
      case_when(
        height == 0.005 ~ "ground",
        height == 0.03 ~ "low_veg",
        height == 0.15 ~ "high_veg",
        TRUE ~ paste0("h_", height)
      )
    ))
  
  wide_df <- summary_df %>%
    select(species, year, site_orig, site_clim, sex, year_period, 
           condition, total_net_gains) %>%
    pivot_wider(
      names_from = condition,
      values_from = total_net_gains
    )
  
  return(wide_df)
}

# Main function to run the analysis with TPC normalization
summarize_eb_results <- function(pops_df, 
                                 base_dir = "C:/Users/jmsmi/OneDrive/Documents/GitHub/EB_Ch1/output/results") {
  # Calculate species-specific average TPC areas
  cat("Calculating species-specific average TPC areas...\n")
  species_avg_tpc <- calculate_species_avg_tpc(pops_df)
  print(species_avg_tpc)
  
  rds_files <- list.files(base_dir, pattern = "^eb_results_.*\\.rds$", 
                          recursive = TRUE, full.names = TRUE)
  
  cat("Found", length(rds_files), "RDS files to process.\n")
  
  all_results <- list()
  for (i in seq_along(rds_files)) {
    file <- rds_files[i]
    cat("Processing file", i, "of", length(rds_files), ":", basename(file), "\n")
    result <- process_eb_file(file, pops_df, species_avg_tpc)
    if (!is.null(result)) {
      all_results[[i]] <- result
    }
  }
  
  combined_results <- bind_rows(all_results)
  wide_results <- create_wide_format(combined_results)
  
  return(list(
    long = combined_results,
    wide = wide_results
  ))
}

# Run the analysis
results <- summarize_eb_results(pops_df = pops)
write_csv(results$wide, "eb_results_summary_wide_tpc_normalized.csv")










# library(tidyverse)
# 
# # Lobry-Rosso-Flandrois TPC function
# lrf_1991_exact <- function(temp, rmax, topt, tmin, tmax) {
#   # Return 0 for temps outside valid range
#   if (temp <= tmin || temp >= tmax) {
#     return(0)
#   }
#   
#   # Implement the exact equation
#   numerator <- (temp - tmax) * (temp - tmin)^2
#   denominator <- (topt - tmin) * ((topt - tmin) * (temp - topt) - (topt - tmax) * (topt + tmin - 2 * temp))
#   
#   rate <- rmax * numerator / denominator
#   
#   # Handle NaN or negative values
#   if (is.na(rate) || rate < 0) {
#     return(0)
#   }
#   
#   return(rate)
# }
# 
# # Function to calculate area under TPC curve using LRF model
# calculate_tpc_area <- function(Tmin, Topt, Above, Ropt) {
#   # Tmax is Topt + Above
#   Tmax <- Topt + Above
#   
#   # Create temperature range from Tmin to Tmax
#   T_range <- seq(Tmin, Tmax, by = 0.1)
#   
#   # Calculate TPC values using the LRF model
#   tpc_values <- sapply(T_range, function(T) {
#     lrf_1991_exact(temp = T, rmax = Ropt, topt = Topt, tmin = Tmin, tmax = Tmax)
#   })
#   
#   # Calculate area using trapezoidal integration
#   area <- sum(diff(T_range) * (head(tpc_values, -1) + tail(tpc_values, -1)) / 2)
#   
#   # Check for invalid area
#   if (is.na(area) || area == 0 || is.infinite(area)) {
#     warning(paste("Invalid TPC area:", area, 
#                   "- Parameters: Tmin=", Tmin, "Topt=", Topt, 
#                   "Above=", Above, "Ropt=", Ropt))
#     return(NA)
#   }
#   
#   return(area)
# }
# 
# # Function to extract info from file path
# extract_filename_info <- function(filepath) {
#   filename <- basename(filepath)
#   parts <- strsplit(tools::file_path_sans_ext(filename), "_")[[1]]
#   
#   species <- parts[3]
#   year <- parts[4]
#   site_orig <- parts[5]
#   site_clim <- parts[6]
#   
#   dir_parts <- strsplit(dirname(filepath), "/")[[1]]
#   sex <- dir_parts[length(dir_parts)]
#   
#   year_period <- case_when(
#     as.numeric(year) < 1965 ~ "historical",
#     as.numeric(year) > 2000 ~ "contemporary",
#     TRUE ~ "other"
#   )
#   
#   return(list(
#     species = species,
#     year = year,
#     site_orig = site_orig,
#     site_clim = site_clim,
#     sex = sex,
#     year_period = year_period
#   ))
# }
# 
# # Modified function to process a single RDS file with TPC area normalization
# process_eb_file <- function(filepath, pops_df) {
#   file_info <- extract_filename_info(filepath)
#   
#   tryCatch({
#     data <- readRDS(filepath)
#     
#     # Get TPC parameters for this species/site/sex combination
#     tpc_params <- pops_df %>%
#       filter(spp == file_info$species,
#              site == file_info$site_orig,
#              sex == file_info$sex)
#     
#     if (nrow(tpc_params) == 0) {
#       warning(paste("No TPC parameters found for:", 
#                     file_info$species, file_info$site_orig, file_info$sex))
#       return(NULL)
#     }
#     
#     # Calculate area under TPC using LRF model
#     tpc_area <- calculate_tpc_area(
#       Tmin = tpc_params$Tmin[1],
#       Topt = tpc_params$Topt[1],
#       Above = tpc_params$Above[1],
#       Ropt = tpc_params$Ropt[1]
#     )
#     
#     # Check if TPC area calculation failed
#     if (is.na(tpc_area)) {
#       warning(paste("Skipping file due to invalid TPC area:", filepath))
#       return(NULL)
#     }
#     
#     # Normalize gains by TPC area and recalculate net gains
#     data <- data %>%
#       mutate(
#         gains_normalized = gains / tpc_area,
#         net_gains_normalized = gains_normalized - losses
#       )
#     
#     # Summarize with normalized values
#     summary <- data %>%
#       group_by(shade, height) %>%
#       summarize(
#         total_net_gains = sum(net_gains_normalized, na.rm = TRUE),
#         avg_gains = mean(gains_normalized, na.rm = TRUE),
#         avg_losses = mean(losses, na.rm = TRUE),
#         hours = n(),
#         .groups = "drop"
#       ) %>%
#       mutate(
#         species = file_info$species,
#         year = file_info$year,
#         site_orig = file_info$site_orig,
#         site_clim = file_info$site_clim,
#         sex = file_info$sex,
#         year_period = file_info$year_period,
#         filepath = filepath
#       )
#     
#     return(summary)
#   }, error = function(e) {
#     warning(paste("Error processing file:", filepath, "\nError:", e$message))
#     return(NULL)
#   })
# }
# 
# # Function to create a wide format of the results
# create_wide_format <- function(summary_df) {
#   summary_df <- summary_df %>%
#     mutate(condition = paste0(
#       case_when(
#         shade == 0 ~ "sun",
#         shade == 0.3 ~ "partial_shade",
#         shade == 0.9 ~ "full_shade",
#         TRUE ~ paste0("shade_", shade)
#       ),
#       "_",
#       case_when(
#         height == 0.005 ~ "ground",
#         height == 0.03 ~ "low_veg",
#         height == 0.15 ~ "high_veg",
#         TRUE ~ paste0("h_", height)
#       )
#     ))
#   
#   wide_df <- summary_df %>%
#     select(species, year, site_orig, site_clim, sex, year_period, 
#            condition, total_net_gains) %>%
#     pivot_wider(
#       names_from = condition,
#       values_from = total_net_gains
#     )
#   
#   return(wide_df)
# }
# 
# # Main function to run the analysis with TPC normalization
# summarize_eb_results <- function(pops_df, 
#                                  base_dir = "C:/Users/jmsmi/OneDrive/Documents/GitHub/EB_Ch1/output/results") {
#   rds_files <- list.files(base_dir, pattern = "^eb_results_.*\\.rds$", 
#                           recursive = TRUE, full.names = TRUE)
#   
#   cat("Found", length(rds_files), "RDS files to process.\n")
#   
#   all_results <- list()
#   for (i in seq_along(rds_files)) {
#     file <- rds_files[i]
#     cat("Processing file", i, "of", length(rds_files), ":", basename(file), "\n")
#     result <- process_eb_file(file, pops_df)
#     if (!is.null(result)) {
#       all_results[[i]] <- result
#     }
#   }
#   
#   combined_results <- bind_rows(all_results)
#   wide_results <- create_wide_format(combined_results)
#   
#   return(list(
#     long = combined_results,
#     wide = wide_results
#   ))
# }
# 
# # Run the analysis
# results <- summarize_eb_results(pops_df = pops)
# write_csv(results$wide, "eb_results_summary_wide_tpc_normalized.csv")
# 
# 
# 
# # Test with one of your populations (MB, A1, F)
# test_area <- calculate_tpc_area(Tmin = 10.73, Topt = 39.52, Above = 16.73, Ropt = 2.99)
# print(paste("TPC area:", test_area))
# 
# # Or visualize the curve
# T_test <- seq(10.73, 39.52 + 16.73, by = 0.1)
# rates <- sapply(T_test, function(t) lrf_1991_exact(t, rmax = 2.99, topt = 39.52, tmin = 10.73, tmax = 56.25))
# plot(T_test, rates, type = "l", xlab = "Temperature (°C)", ylab = "Rate", 
#      main = "LRF TPC Curve")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(tidyverse)
# 
# # Function to calculate area under TPC curve
# # Using the Sharpe-Schoolfield model based on your parameters
# calculate_tpc_area <- function(rmr_b0, rmr_b1, rmr_b2, rmr_b3, Tmin, Topt) {
#   # Define temperature range for integration (reasonable biological range)
#   # Using Tmin to a bit above Topt to capture the full curve
#   T_range <- seq(Tmin, Topt + 10, by = 0.1)
#   
#   # Sharpe-Schoolfield TPC function
#   # R(T) = rmr_b0 + rmr_b1*T + rmr_b2*T^2 + rmr_b3*T^3
#   tpc_values <- rmr_b0 + rmr_b1*T_range + rmr_b2*T_range^2 + rmr_b3*T_range^3
#   
#   # Set negative values to 0 (metabolic rate can't be negative)
#   tpc_values[tpc_values < 0] <- 0
#   
#   # Calculate area using trapezoidal integration
#   area <- sum(diff(T_range) * (head(tpc_values, -1) + tail(tpc_values, -1)) / 2)
#   
#   return(area)
# }
# 
# # Function to extract info from file path
# extract_filename_info <- function(filepath) {
#   filename <- basename(filepath)
#   parts <- strsplit(tools::file_path_sans_ext(filename), "_")[[1]]
#   
#   species <- parts[3]
#   year <- parts[4]
#   site_orig <- parts[5]
#   site_clim <- parts[6]
#   
#   dir_parts <- strsplit(dirname(filepath), "/")[[1]]
#   sex <- dir_parts[length(dir_parts)]
#   
#   year_period <- case_when(
#     as.numeric(year) < 1965 ~ "historical",
#     as.numeric(year) > 2000 ~ "contemporary",
#     TRUE ~ "other"
#   )
#   
#   return(list(
#     species = species,
#     year = year,
#     site_orig = site_orig,
#     site_clim = site_clim,
#     sex = sex,
#     year_period = year_period
#   ))
# }
# 
# # Modified function to process a single RDS file with TPC area normalization
# process_eb_file <- function(filepath, pops_df) {
#   file_info <- extract_filename_info(filepath)
#   
#   tryCatch({
#     data <- readRDS(filepath)
#     
#     # Get TPC parameters for this species/site/sex combination
#     tpc_params <- pops_df %>%
#       filter(spp == file_info$species,
#              site == file_info$site_orig,
#              sex == file_info$sex)
#     
#     if (nrow(tpc_params) == 0) {
#       warning(paste("No TPC parameters found for:", 
#                     file_info$species, file_info$site_orig, file_info$sex))
#       return(NULL)
#     }
#     
#     # Calculate area under TPC
#     tpc_area <- calculate_tpc_area(
#       rmr_b0 = tpc_params$rmr_b0[1],
#       rmr_b1 = tpc_params$rmr_b1[1],
#       rmr_b2 = tpc_params$rmr_b2[1],
#       rmr_b3 = tpc_params$rmr_b3[1],
#       Tmin = tpc_params$Tmin[1],
#       Topt = tpc_params$Topt[1]
#     )
#     
#     # Normalize gains by TPC area and recalculate net gains
#     data <- data %>%
#       mutate(
#         gains_normalized = gains / tpc_area,
#         net_gains_normalized = gains_normalized - losses
#       )
#     
#     # Summarize with normalized values
#     summary <- data %>%
#       group_by(shade, height) %>%
#       summarize(
#         total_net_gains = sum(net_gains_normalized, na.rm = TRUE),
#         avg_gains = mean(gains_normalized, na.rm = TRUE),
#         avg_losses = mean(losses, na.rm = TRUE),
#         hours = n(),
#         .groups = "drop"
#       ) %>%
#       mutate(
#         species = file_info$species,
#         year = file_info$year,
#         site_orig = file_info$site_orig,
#         site_clim = file_info$site_clim,
#         sex = file_info$sex,
#         year_period = file_info$year_period,
#         filepath = filepath,
#         # Add TPC parameters and area for reference
#         tpc_area = tpc_area,
#         elev = tpc_params$elev[1],
#         mass = tpc_params$mass[1],
#         rmr_b0 = tpc_params$rmr_b0[1],
#         rmr_b1 = tpc_params$rmr_b1[1],
#         rmr_b2 = tpc_params$rmr_b2[1],
#         rmr_b3 = tpc_params$rmr_b3[1],
#         Tmin = tpc_params$Tmin[1],
#         Topt = tpc_params$Topt[1],
#         Above = tpc_params$Above[1],
#         Ropt = tpc_params$Ropt[1]
#       )
#     
#     return(summary)
#   }, error = function(e) {
#     warning(paste("Error processing file:", filepath, "\nError:", e$message))
#     return(NULL)
#   })
# }
# 
# # Function to create a wide format of the results
# create_wide_format <- function(summary_df) {
#   summary_df <- summary_df %>%
#     mutate(condition = paste0(
#       case_when(
#         shade == 0 ~ "sun",
#         shade == 0.3 ~ "partial_shade",
#         shade == 0.9 ~ "full_shade",
#         TRUE ~ paste0("shade_", shade)
#       ),
#       "_",
#       case_when(
#         height == 0.005 ~ "ground",
#         height == 0.03 ~ "low_veg",
#         height == 0.15 ~ "high_veg",
#         TRUE ~ paste0("h_", height)
#       )
#     ))
#   
#   wide_df <- summary_df %>%
#     select(species, year, site_orig, site_clim, sex, year_period, 
#            condition, total_net_gains, 
#            # Include TPC info (will be same for all conditions, but useful to have)
#            tpc_area, elev, mass, rmr_b0, rmr_b1, rmr_b2, rmr_b3, 
#            Tmin, Topt, Above, Ropt) %>%
#     pivot_wider(
#       names_from = condition,
#       values_from = total_net_gains
#     )
#   
#   return(wide_df)
# }
# 
# # Main function to run the analysis with TPC normalization
# summarize_eb_results <- function(pops_df, 
#                                  base_dir = "C:/Users/jmsmi/OneDrive/Documents/GitHub/EB_Ch1/output/results") {
#   rds_files <- list.files(base_dir, pattern = "^eb_results_.*\\.rds$", 
#                           recursive = TRUE, full.names = TRUE)
#   
#   cat("Found", length(rds_files), "RDS files to process.\n")
#   
#   all_results <- list()
#   for (i in seq_along(rds_files)) {
#     file <- rds_files[i]
#     cat("Processing file", i, "of", length(rds_files), ":", basename(file), "\n")
#     result <- process_eb_file(file, pops_df)
#     if (!is.null(result)) {
#       all_results[[i]] <- result
#     }
#   }
#   
#   combined_results <- bind_rows(all_results)
#   wide_results <- create_wide_format(combined_results)
#   
#   return(list(
#     long = combined_results,
#     wide = wide_results
#   ))
# }
# 
# # Run the analysis (assuming your pops dataframe is loaded)
# results <- summarize_eb_results(pops_df = pops)
# write_csv(results$wide, "eb_results_summary_wide_tpc_normalized.csv")