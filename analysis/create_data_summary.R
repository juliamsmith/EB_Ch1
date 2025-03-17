library(tidyverse)

#might need to set wd to source file location

# Function to extract info from file path
extract_filename_info <- function(filepath) {
  # Extract components from filename (assumes filename pattern: eb_results_SPECIES_YEAR_SITE-ORIG_SITE-CLIM.rds)
  filename <- basename(filepath)
  parts <- strsplit(tools::file_path_sans_ext(filename), "_")[[1]]
  
  # Extract species, year, site_orig, site_clim
  species <- parts[3]
  year <- parts[4]
  site_orig <- parts[5]
  site_clim <- parts[6]
  
  # Extract sex from directory structure
  dir_parts <- strsplit(dirname(filepath), "/")[[1]]
  sex <- dir_parts[length(dir_parts)]
  
  # Determine period based on year
  # Adjust these thresholds based on your actual year groupings
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

# Function to process a single RDS file
process_eb_file <- function(filepath) {
  # Extract info from filename
  file_info <- extract_filename_info(filepath)
  
  # Read RDS file
  tryCatch({
    data <- readRDS(filepath)
    
    # Summarize net gains by shade and height
    summary <- data %>%
      group_by(shade, height) %>%
      summarize(total_net_gains = sum(net_gains, na.rm = TRUE),
                avg_gains = mean(gains, na.rm = TRUE),
                avg_losses = mean(losses, na.rm = TRUE),
                hours = n(),
                .groups = "drop") %>%
      # Add file info columns
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
  # Create a descriptive condition name for each shade/height combination
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
  
  # Pivot to wide format
  wide_df <- summary_df %>%
    select(species, year, site_orig, site_clim, sex, year_period, condition, total_net_gains) %>%
    pivot_wider(
      names_from = condition,
      values_from = total_net_gains
    )
  
  return(wide_df)
}

#sorry for hardcoding of my own folder locations

# Main function to run the analysis
summarize_eb_results <- function(base_dir = "C:/Users/jmsmi/OneDrive/Documents/GitHub/EB_Ch1/output/results") {
  # Find all RDS files
  rds_files <- list.files(base_dir, pattern = "^eb_results_.*\\.rds$", 
                          recursive = TRUE, full.names = TRUE)
  
  cat("Found", length(rds_files), "RDS files to process.\n")
  
  # Process each file
  all_results <- list()
  for (i in seq_along(rds_files)) {
    file <- rds_files[i]
    cat("Processing file", i, "of", length(rds_files), ":", basename(file), "\n")
    result <- process_eb_file(file)
    if (!is.null(result)) {
      all_results[[i]] <- result
    }
  }
  
  # Combine all results
  combined_results <- bind_rows(all_results)
  
  # Create wide format version
  wide_results <- create_wide_format(combined_results)
  
  # Return both formats
  return(list(
    long = combined_results,
    wide = wide_results
  ))
}

# Run the analysis
results <- summarize_eb_results()

# # Save the outputs -- FIRST SET WD
# write_csv(results$long, "eb_results_summary_long.csv")
# write_csv(results$wide, "eb_results_summary_wide.csv") #could uncomment

# Preview the results
cat("\nSummary of long format results:\n")
print(dim(results$long))
print(head(results$long))

cat("\nSummary of wide format results:\n")
print(dim(results$wide))
print(head(results$wide))

# Example analysis: Compare average net gains across different conditions
cat("\nAverage net gains by species and condition:\n")
results$long %>%
  group_by(species, condition = paste(shade, height)) %>%
  summarize(
    avg_net_gains = mean(total_net_gains, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(avg_net_gains)) %>%
  print(n = 10)

# Example analysis: Compare historical vs contemporary periods
cat("\nComparison of historical vs contemporary periods:\n")
results$long %>%
  filter(year_period %in% c("historical", "contemporary")) %>%
  group_by(year_period, species, condition = paste(shade, height)) %>%
  summarize(
    avg_net_gains = mean(total_net_gains, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = year_period,
    values_from = avg_net_gains
  ) %>%
  mutate(change = contemporary - historical) %>%
  arrange(desc(change)) %>%
  print(n = 10)