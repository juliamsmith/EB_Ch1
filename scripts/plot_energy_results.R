#!/usr/bin/env Rscript

# Plot Energy Budget Results
# This script creates visualizations of the energy budget calculations

# Load required packages
library(tidyverse)
library(cowplot)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
results_file <- ifelse(length(args) > 0, args[1], "output/combined_results.rds")
output_dir <- ifelse(length(args) > 1, args[2], "output/figures")

# Check if results file exists
if (!file.exists(results_file)) {
  stop(paste("Results file not found:", results_file))
}

# Load results
results <- readRDS(results_file)

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Plot 1: Total energy gain by species, site origin, and climate site for each period
plot_total_energy <- function(data) {
  # Summarize data
  energy_summary <- data %>%
    group_by(species, sex, site_orig, site_clim, year_period) %>%
    summarize(
      total_net_gain = sum(net_gains, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create plot
  p <- ggplot(energy_summary, aes(x = site_clim, y = total_net_gain, fill = year_period)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(sex ~ species + site_orig, scales = "free_y") +
    labs(
      title = "Total Net Energy Gain by Site and Period",
      x = "Climate Site",
      y = "Total Net Energy Gain (kJ)",
      fill = "Period"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# Plot 2: Body temperature distributions by species, site, and period
plot_temperature_distributions <- function(data) {
  # Create plot
  p <- ggplot(data, aes(x = Tb, fill = year_period)) +
    geom_density(alpha = 0.5) +
    facet_grid(site_clim ~ species + site_orig) +
    labs(
      title = "Body Temperature Distributions by Site and Period",
      x = "Body Temperature (°C)",
      y = "Density",
      fill = "Period"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p)
}

# Plot 3: Change in energy gain from historical to recent periods
plot_energy_change <- function(data) {
  # Summarize data
  energy_change <- data %>%
    group_by(species, sex, site_orig, site_clim, year_period) %>%
    summarize(
      total_net_gain = sum(net_gains, na.rm = TRUE),
      n_days = n_distinct(as.Date(dtuse)),
      daily_net_gain = total_net_gain / n_days,
      .groups = "drop"
    ) %>%
    pivot_wider(
      id_cols = c(species, sex, site_orig, site_clim),
      names_from = year_period,
      values_from = daily_net_gain
    ) %>%
    mutate(
      change = recent - historical,
      percent_change = (change / historical) * 100
    )
  
  # Create plot
  p <- ggplot(energy_change, aes(x = site_clim, y = percent_change, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(sex ~ site_orig) +
    labs(
      title = "Percent Change in Daily Energy Gain (Historical to Recent)",
      x = "Climate Site",
      y = "Percent Change",
      fill = "Species"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  
  return(p)
}

# Create plots
total_energy_plot <- plot_total_energy(results)
temp_dist_plot <- plot_temperature_distributions(results)
energy_change_plot <- plot_energy_change(results)

# Save plots
ggsave(
  filename = file.path(output_dir, "total_energy.pdf"),
  plot = total_energy_plot,
  width = 12,
  height = 8
)

ggsave(
  filename = file.path(output_dir, "temp_distributions.pdf"),
  plot = temp_dist_plot,
  width = 12,
  height = 8
)

ggsave(
  filename = file.path(output_dir, "energy_change.pdf"),
  plot = energy_change_plot,
  width = 12,
  height = 8
)

# Also save PNG versions
ggsave(
  filename = file.path(output_dir, "total_energy.png"),
  plot = total_energy_plot,
  width = 12,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "temp_distributions.png"),
  plot = temp_dist_plot,
  width = 12,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "energy_change.png"),
  plot = energy_change_plot,
  width = 12,
  height = 8,
  dpi = 300
)

cat("Energy budget result plots have been saved to", output_dir, "\n")