# Required libraries
library(tidyverse)
library(ggridges)
library(lubridate)
library(viridis)

#-------------------------------------------------------------
# Function to load and process data from simplified RDS files
#-------------------------------------------------------------

load_simplified_data <- function(base_path, sites, years, 
                                 height_filter = 0.03, 
                                 shade_filter = 0.3,
                                 daytime_only = TRUE) {
  
  # Initialize an empty dataframe to store all data
  all_data <- data.frame()
  
  # Track which files were processed
  processed_files <- character()
  
  # Loop through sites and years
  for (site in sites) {
    for (year in years) {
      # Construct file path
      file_path <- file.path(base_path, site, paste0("simplified_", site, "_", year, ".rds"))
      
      # Check if file exists
      if (file.exists(file_path)) {
        tryCatch({
          # Read the data
          data <- readRDS(file_path)
          
          # Apply filters
          filtered_data <- data %>%
            filter(
              # Filter for specific height
              near(height, height_filter, tol = 0.001),
              # Filter for specific shade
              near(shade, shade_filter, tol = 0.01)
            )
          
          # Additional filter for daytime hours if requested
          if (daytime_only) {
            filtered_data <- filtered_data %>%
              filter(hour >= 8 & hour <= 18)
          }
          
          # Add to the combined dataset
          all_data <- bind_rows(all_data, filtered_data)
          
          # Add to list of processed files
          processed_files <- c(processed_files, file_path)
          
        }, error = function(e) {
          warning(paste("Error processing file:", file_path, "-", e$message))
        })
      } else {
        message(paste("File not found:", file_path))
      }
    }
  }
  
  # Return both the combined data and the list of processed files
  return(list(
    data = all_data,
    processed_files = processed_files
  ))
}

#-------------------------------------------------------------
# Function to create the ridge plot comparing sites
#-------------------------------------------------------------

create_site_temperature_plot <- function(data, 
                                         site_order = c("Eldo", "A1", "B1", "C1", "D1"),
                                         include_stats = TRUE) {
  
  # Prepare data for plotting
  plot_data <- data %>%
    mutate(
      # Ensure sites are ordered correctly
      site_factor = factor(site, levels = site_order)
    )
  
  # Calculate summary statistics if requested
  if (include_stats) {
    stats_data <- plot_data %>%
      group_by(site) %>%
      summarise(
        mean_temp = mean(Tair, na.rm = TRUE),
        median_temp = median(Tair, na.rm = TRUE),
        min_temp = min(Tair, na.rm = TRUE),
        max_temp = max(Tair, na.rm = TRUE),
        q25 = quantile(Tair, 0.25, na.rm = TRUE),
        q75 = quantile(Tair, 0.75, na.rm = TRUE),
        sd_temp = sd(Tair, na.rm = TRUE),
        n = n()
      )
    
    # Print the stats
    print(stats_data)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Tair, y = site_factor, fill = site_factor)) +
    geom_density_ridges(
      alpha = 0.7, 
      scale = 1.5, 
      rel_min_height = 0.01,
      quantile_lines = TRUE,
      quantiles = 2  # Show median line
    ) +
    scale_fill_viridis_d(
      name = "Site",
      option = "plasma",
      direction = -1
    ) +
    labs(
      title = "Air Temperature Distribution Across Sites",
      subtitle = paste0(
        "Years (", min(plot_data$year), "-", max(plot_data$year), 
        ") combined, ",
        ifelse(all(plot_data$hour >= 8 & plot_data$hour <= 18), "daytime hours only", "all hours"),
        ", height = ", round(mean(plot_data$height), 3),
        ", shade = ", round(mean(plot_data$shade), 2)
      ),
      x = "Air Temperature (°C)",
      y = "Site (ordered by increasing elevation)",
      caption = "Median temperatures shown with darker lines"
    ) +
    theme_ridges(center_axis_labels = TRUE) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    scale_x_continuous(
      expand = c(0.01, 0),
      breaks = seq(0, 50, by = 5)
    )
  
  return(p)
}

#-------------------------------------------------------------
# Function to create a monthly comparison
#-------------------------------------------------------------

create_monthly_site_plot <- function(data, site_order = c("Eldo", "A1", "B1", "C1", "D1")) {
  
  # Prepare data
  plot_data <- data %>%
    mutate(
      site_factor = factor(site, levels = site_order),
      month_name = factor(month.name[month], levels = month.name[6:8])
    )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Tair, y = site_factor, fill = site_factor)) +
    geom_density_ridges(
      alpha = 0.7, 
      scale = 1.5, 
      rel_min_height = 0.01,
      quantile_lines = TRUE,
      quantiles = 2
    ) +
    scale_fill_viridis_d(
      name = "Site",
      option = "plasma",
      direction = -1
    ) +
    facet_wrap(~ month_name, ncol = 3) +
    labs(
      title = "Monthly Air Temperature Distribution Across Sites",
      subtitle = paste0(
        "Years (", min(plot_data$year), "-", max(plot_data$year), 
        ") combined, ",
        ifelse(all(plot_data$hour >= 8 & plot_data$hour <= 18), "daytime hours only", "all hours")
      ),
      x = "Air Temperature (°C)",
      y = "Site (ordered by increasing elevation)"
    ) +
    theme_ridges(center_axis_labels = TRUE) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    scale_x_continuous(
      expand = c(0.01, 0),
      breaks = seq(0, 50, by = 5)
    )
  
  return(p)
}

#-------------------------------------------------------------
# Function to create year comparison plot for a specific site
#-------------------------------------------------------------

create_year_comparison_plot <- function(data, target_site, year_range = NULL) {
  
  # Filter data for the specific site
  site_data <- data %>%
    filter(site == target_site)
  
  # If year range is specified, filter for those years
  if (!is.null(year_range)) {
    site_data <- site_data %>%
      filter(year >= year_range[1] & year <= year_range[2])
  }
  
  # Ensure we have data
  if (nrow(site_data) == 0) {
    warning(paste("No data for site:", target_site))
    return(NULL)
  }
  
  # Create year labels
  site_data <- site_data %>%
    mutate(year_factor = factor(year))
  
  # Create the plot
  p <- ggplot(site_data, aes(x = Tair, y = year_factor, fill = year_factor)) +
    geom_density_ridges(
      alpha = 0.7, 
      scale = 0.9, 
      rel_min_height = 0.01,
      quantile_lines = TRUE,
      quantiles = 2
    ) +
    scale_fill_viridis_d(
      name = "Year",
      option = "viridis"
    ) +
    labs(
      title = paste("Air Temperature Distribution at", target_site, "Over Years"),
      subtitle = paste0(
        "Daytime hours only, ",
        "height = ", round(mean(site_data$height), 3),
        ", shade = ", round(mean(site_data$shade), 2)
      ),
      x = "Air Temperature (°C)",
      y = "Year"
    ) +
    theme_ridges(center_axis_labels = TRUE) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    scale_x_continuous(
      expand = c(0.01, 0),
      breaks = seq(0, 50, by = 5)
    )
  
  return(p)
}

#-------------------------------------------------------------
# Main code to execute the analysis
#-------------------------------------------------------------

# Define the path to your data
# base_path <- "G:/Shared drives/RoL_FitnessConstraints/projects/Julia_Energy_Budget/processed_era"

base_path <- "C:/Users/jmsmi/Downloads/ERA5_land_data/processed_era5land"

# Define sites and years to process
sites <- c("Eldo", "A1", "B1", "C1", "D1")
years <- 2001:2024  # Years after 2000

# Load the data
result <- load_simplified_data(
  base_path = base_path,
  sites = sites,
  years = years,
  height_filter = 0.03,
  shade_filter = 0.3,
  daytime_only = TRUE
)

# Display how many files were processed
cat("Processed", length(result$processed_files), "files\n")

# Create and display the site comparison plot
site_plot <- create_site_temperature_plot(result$data)
print(site_plot)

# Save the plot if desired
#ggsave("temperature_distribution_by_site.png", site_plot, width = 10, height = 8, dpi = 300)

# Create monthly comparison
monthly_plot <- create_monthly_site_plot(result$data)
print(monthly_plot)
#ggsave("monthly_temperature_by_site.png", monthly_plot, width = 12, height = 8, dpi = 300)

# Create year comparison for a specific site (e.g., B1)
b1_years_plot <- create_year_comparison_plot(result$data, "B1")
print(b1_years_plot)
#ggsave("B1_temperature_by_year.png", b1_years_plot, width = 10, height = 10, dpi = 300)

#-------------------------------------------------------------
# Diagnostic function to help troubleshoot site data
#-------------------------------------------------------------

diagnose_site_data <- function(data, site_to_check) {
  # Filter for the specific site
  site_data <- data %>% filter(site == site_to_check)
  
  if (nrow(site_data) == 0) {
    cat("No data found for site:", site_to_check, "\n")
    return(NULL)
  }
  
  # Print summary statistics
  cat("\n==== Summary Statistics for", site_to_check, "====\n")
  cat("Number of observations:", nrow(site_data), "\n")
  cat("Year range:", min(site_data$year), "to", max(site_data$year), "\n")
  cat("Temperature range:", min(site_data$Tair, na.rm = TRUE), "to", 
      max(site_data$Tair, na.rm = TRUE), "°C\n")
  cat("Mean temperature:", mean(site_data$Tair, na.rm = TRUE), "°C\n")
  cat("Median temperature:", median(site_data$Tair, na.rm = TRUE), "°C\n")
  
  # Check for bimodality with a histogram
  hist_plot <- ggplot(site_data, aes(x = Tair)) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
    labs(
      title = paste("Temperature Histogram for", site_to_check),
      x = "Air Temperature (°C)",
      y = "Count"
    ) +
    theme_minimal()
  
  # Check for unusual patterns by year
  yearly_boxplot <- ggplot(site_data, aes(x = factor(year), y = Tair)) +
    geom_boxplot(fill = "lightblue") +
    labs(
      title = paste("Yearly Temperature Distribution for", site_to_check),
      x = "Year",
      y = "Air Temperature (°C)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Return diagnostic plots
  return(list(
    histogram = hist_plot,
    yearly_boxplot = yearly_boxplot
  ))
}

# Example usage:
# diagnostic_plots <- diagnose_site_data(result$data, "B1")
# print(diagnostic_plots$histogram)
# print(diagnostic_plots$yearly_boxplot)