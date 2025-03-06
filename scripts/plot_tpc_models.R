#!/usr/bin/env Rscript

# Plot Thermal Performance Curve Models
# This script creates visualizations of the TPC models used for feeding rates

# Load required packages
library(tidyverse)
library(cowplot)

# Set working directory to the script location
script_dir <- dirname(sys.frame(1)$ofile)
setwd(file.path(script_dir, ".."))

# Source the energy functions
source("R/energy_functions.R")

# Function to generate data for plotting
generate_tpc_data <- function(temp_range=10:50) {
  # Create a dataframe for all combinations
  sites <- c("Eldo", "A1", "B1", "C1", "D1")
  species <- c("MB", "MS")
  sexes <- c("F", "M")
  
  # Create all combinations
  combinations <- expand.grid(
    species = species,
    site = sites,
    sex = sexes,
    temp = temp_range
  )
  
  # Filter out combinations that don't exist
  combinations <- combinations %>%
    filter(!(species == "MB" & site == "Eldo")) %>%
    filter(!(species == "MS" & (site == "C1" | site == "D1")))
  
  # Calculate feeding rates
  results <- combinations %>%
    rowwise() %>%
    mutate(
      feeding_rate = if(species == "MS") {
        get_ms_feeding_rate(temp, site, sex)
      } else {
        get_mb_feeding_rate(temp, site, sex)
      }
    )
  
  return(results)
}

# Generate data
tpc_data <- generate_tpc_data()

# Plot by species and site
plot_by_species_site <- function(data) {
  plots <- list()
  
  for (sp in unique(data$species)) {
    sp_data <- data %>% filter(species == sp)
    
    # Create plot
    p <- ggplot(sp_data, aes(x = temp, y = feeding_rate, color = site, linetype = sex)) +
      geom_line(size = 1) +
      labs(
        title = paste("Thermal Performance Curve for", sp),
        x = "Temperature (°C)",
        y = "Feeding Rate",
        color = "Site",
        linetype = "Sex"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
      )
    
    plots[[sp]] <- p
  }
  
  return(plots)
}

# Create plots
tpc_plots <- plot_by_species_site(tpc_data)

# Save plots
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

for (sp in names(tpc_plots)) {
  ggsave(
    filename = file.path("output/figures", paste0("tpc_", sp, ".pdf")),
    plot = tpc_plots[[sp]],
    width = 8,
    height = 6
  )
  
  ggsave(
    filename = file.path("output/figures", paste0("tpc_", sp, ".png")),
    plot = tpc_plots[[sp]],
    width = 8,
    height = 6,
    dpi = 300
  )
}

# Create a combined plot
combined_plot <- plot_grid(plotlist = tpc_plots, ncol = 1, labels = "AUTO")

ggsave(
  filename = file.path("output/figures", "tpc_combined.pdf"),
  plot = combined_plot,
  width = 8,
  height = 10
)

ggsave(
  filename = file.path("output/figures", "tpc_combined.png"),
  plot = combined_plot,
  width = 8,
  height = 10,
  dpi = 300
)

cat("TPC model plots have been saved to output/figures/\n")