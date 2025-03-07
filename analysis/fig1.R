# Load required libraries
library(tidyverse)
library(ggplot2)
library(patchwork)  # For combining plots

# Implement LRF function exactly as shown in the pasted code
lrf_1991_exact <- function(temp, rmax, topt, tmin, tmax) {
  # Return 0 for temps outside valid range
  if (temp <= tmin || temp >= tmax) {
    return(0)
  }
  
  # Implement the exact equation as shown in the image
  numerator <- (temp - tmax) * (temp - tmin)^2
  denominator <- (topt - tmin) * ((topt - tmin) * (temp - topt) - (topt - tmax) * (topt + tmin - 2 * temp))
  
  rate <- rmax * numerator / denominator
  
  # Handle NaN or negative values
  if (is.na(rate) || rate < 0) {
    return(0)
  }
  
  return(rate)
}

# Function to calculate feeding rate for MS using exact LRF model (averaged across sexes)
get_ms_feeding_rate <- function(temp, site) {
  # Parameters for MS from the model output
  site_params <- list(
    Eldo = list(Tmin = 14.74, Topt = 39.88, Above = 9.56, Ropt = 3.50),
    A1 = list(Tmin = 11.67, Topt = 39.89, Above = 6.68, Ropt = 3.63),
    B1 = list(Tmin = 10.96, Topt = 39.78, Above = 8.20, Ropt = 4.30)
  )
  
  # For D1, use A1 parameters as a placeholder
  if (site == "D1") {
    site_params$D1 <- site_params$A1
  }
  
  # Get parameters for this site
  params <- site_params[[site]]
  
  # Calculate Tmax (Topt + Above)
  Tmax <- params$Topt + params$Above
  
  # Calculate base feeding rate using exact LRF
  base_rate <- tryCatch({
    lrf_1991_exact(temp, params$Ropt, params$Topt, params$Tmin, Tmax)
  }, error = function(e) {
    # Return zero for any errors
    return(0)
  })
  
  # We're not applying sex effects anymore as we're averaging across sexes
  return(base_rate)
}

# Function to calculate feeding rate for MB using exact LRF model (averaged across sexes)
get_mb_feeding_rate <- function(temp, site) {
  # Parameters for MB from the model output
  site_params <- list(
    A1 = list(Tmin = 13.57, Topt = 37.54, Above = 11.70, Ropt = 2.64),
    B1 = list(Tmin = 12.80, Topt = 41.76, Above = 10.74, Ropt = 3.29),
    C1 = list(Tmin = 13.79, Topt = 41.03, Above = 11.77, Ropt = 3.93)
  )
  
  # For Eldo and D1, use A1 and C1 parameters respectively as placeholders
  site_params$Eldo <- site_params$A1
  site_params$D1 <- list(Tmin = 14.0, Topt = 42.0, Above = 12.0, Ropt = 4.0)  # Made-up values for D1
  
  # Get parameters for this site
  params <- site_params[[site]]
  
  # Calculate Tmax (Topt + Above)
  Tmax <- params$Topt + params$Above
  
  # Calculate base feeding rate using exact LRF
  base_rate <- tryCatch({
    lrf_1991_exact(temp, params$Ropt, params$Topt, params$Tmin, Tmax)
  }, error = function(e) {
    # Return zero for any errors
    return(0)
  })
  
  # We're not applying sex effects anymore as we're averaging across sexes
  return(base_rate)
}

# Read the data
# Assuming the ad.csv file is in the working directory
#ad_data <- read.csv("ad.csv")

# Generate curve data
generate_curve_data <- function() {
  # Temperature range for curves - use fine-grained steps for smoother curves
  temp_range <- seq(15, 45, by = 0.1)
  
  # Define which sites to use for each species
  mb_sites <- c("A1", "B1", "C1", "D1")
  ms_sites <- c("Eldo", "A1", "B1")
  
  # Create empty dataframe to store results
  curve_data <- data.frame()
  
  # Generate data for MB
  for (site in mb_sites) {
    site_curve <- data.frame(
      spp = "MB",
      site = site,
      temp = temp_range,
      rate = sapply(temp_range, function(t) get_mb_feeding_rate(t, site)),
      type = "curve"
    )
    curve_data <- rbind(curve_data, site_curve)
  }
  
  # Generate data for MS
  for (site in ms_sites) {
    site_curve <- data.frame(
      spp = "MS",
      site = site,
      temp = temp_range,
      rate = sapply(temp_range, function(t) get_ms_feeding_rate(t, site)),
      type = "curve"
    )
    curve_data <- rbind(curve_data, site_curve)
  }
  
  return(curve_data)
}

# Generate the curve data
curve_data <- generate_curve_data()

# Add type column to the original data
ad_data$type <- "point"

# Combine original data with curve data (dropping sex from curve data since we're ignoring it)
combined_data <- rbind(
  ad_data[, c("spp", "site", "sex", "temp", "rate", "type")],
  transform(curve_data, sex = NA)  # Add NA for sex in curve data
)

# Set the site order based on elevation (low to high)
site_order <- c("Eldo", "A1", "B1", "C1", "D1")
combined_data$site <- factor(combined_data$site, levels = site_order)

# Define colors for sites (cooler as elevation increases)
site_colors <- c(
  "Eldo" = "#FF4500",  # Red-orange
  "A1" = "#FF8C00",    # Dark orange
  "B1" = "#FFD700",    # Gold/yellow
  "C1" = "#4682B4",    # Steel blue
  "D1" = "#0000CD"     # Medium blue
)

# Create the plot
create_tpc_plot <- function(data) {
  # Filter for just the point data for the plot
  point_data <- data %>% 
    filter(type == "point")
  
  # Filter for just the curve data
  curve_data <- data %>% filter(type == "curve")
  
  # Create a plot for each species and site
  p <- ggplot() +
    # Add the curves first (so points appear on top)
    geom_line(
      data = curve_data,
      aes(x = temp, y = rate, color = site),
      size = 1.2
    ) +
    # Add the scatter points with jitter
    geom_jitter(
      data = point_data,
      aes(x = temp, y = rate, color = site),
      width = 0.5, height = 0, alpha = 0.6, size = 1.8, shape = 16
    ) +
    # Set up faceting by species and site
    facet_grid(site ~ spp) +
    # Set axis labels
    labs(
      x = "Temperature (°C)",
      y = "Mass-adjusted feces production rate (mg/g hopper/hr)",
      title = "Thermal Performance Curves for Grasshopper Species"
    ) +
    # Customize the color scheme
    scale_color_manual(
      values = site_colors,
      name = "Site"
    ) +
    # Set the x-axis range
    scale_x_continuous(limits = c(15, 45), breaks = seq(15, 45, by = 5)) +
    # Set the y-axis range
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
    # Theme customization
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      strip.background = element_rect(fill = "lightgrey"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Create the plot
tpc_plot <- create_tpc_plot(combined_data)

# Display the plot
print(tpc_plot)

# Save the plot
#ggsave("grasshopper_tpc_plot_simplified.png", tpc_plot, width = 12, height = 8, dpi = 300)

# Function to create a more compact multi-panel plot (alternative layout)
create_compact_tpc_plot <- function(data) {
  # Filter points and curves
  point_data <- data %>% filter(type == "point")
  curve_data <- data %>% filter(type == "curve")
  
  # Create a separate plot for each species
  plots <- list()
  
  for (species in c("MB", "MS")) {
    # Determine which sites to include based on species
    if (species == "MB") {
      relevant_sites <- c("A1", "B1", "C1", "D1")
    } else {
      relevant_sites <- c("Eldo", "A1", "B1")
    }
    
    # Filter data for this species
    species_points <- point_data %>% 
      filter(spp == species, site %in% relevant_sites)
    
    species_curves <- curve_data %>% 
      filter(spp == species, site %in% relevant_sites)
    
    # Create the plot
    p <- ggplot() +
      # Add curves
      geom_line(
        data = species_curves,
        aes(x = temp, y = rate, color = site),
        size = 1.2
      ) +
      # Add points with jitter
      geom_jitter(
        data = species_points,
        aes(x = temp, y = rate, color = site),
        width = 0.5, height = 0, alpha = 0.6, size = 1.8, shape = 16
      ) +
      # Facet by site only (since we're making separate plots for each species)
      facet_wrap(~ site, ncol = 1) +
      # Labels
      labs(
        x = "Temperature (°C)",
        y = "Mass-adjusted feces production rate",
        title = species
      ) +
      # Colors
      scale_color_manual(values = site_colors) +
      # Axis ranges
      scale_x_continuous(limits = c(15, 45), breaks = seq(15, 45, by = 10)) +
      scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 5)) +
      # Theme
      theme_bw() +
      theme(
        strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14)
      )
    
    plots[[species]] <- p
  }
  
  # Create a legend-only plot
  legend_data <- data.frame(
    site = site_order,
    x = 1,
    y = 1
  )
  
  legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = site)) +
    geom_point() +
    scale_color_manual(values = site_colors, name = "Site") +
    theme_void() +
    theme(legend.position = "bottom")
  
  # Extract the legend
  legend <- cowplot::get_legend(legend_plot)
  
  # Combine plots with patchwork
  combined_plot <- plots$MB + plots$MS + 
    plot_layout(ncol = 2, widths = c(1, 1)) +
    plot_annotation(
      title = "Thermal Performance Curves for Grasshopper Species",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
    )
  
  # Add legend at the bottom
  final_plot <- combined_plot / legend
  
  return(final_plot)
}

# Create the alternative compact plot
compact_plot <- create_compact_tpc_plot(combined_data)

# Save the compact plot
#ggsave("grasshopper_tpc_compact_plot_simplified.png", compact_plot, width = 10, height = 12, dpi = 300)

# Test function to visualize a single curve
test_curve <- function() {
  # Parameters for one example
  tmin <- 11.67   # A1 MS
  topt <- 39.89
  tmax <- topt + 6.68
  rmax <- 3.63
  
  # Generate curve
  temps <- seq(5, 50, by=0.1)
  rates <- sapply(temps, function(t) lrf_1991_exact(t, rmax, topt, tmin, tmax))
  
  # Plot
  df <- data.frame(temp = temps, rate = rates)
  p <- ggplot(df, aes(x = temp, y = rate)) +
    geom_line(color = "#FF8C00", size = 1.2) +
    labs(title = "Example LRF Curve (MS at A1)", 
         x = "Temperature (°C)", 
         y = "Rate") +
    theme_bw()
  
  print(p)
  #ggsave("test_curve.png", p, width = 8, height = 6, dpi = 300)
}

# Run the test curve visualization
test_curve()














# Load required libraries
library(tidyverse)
library(ggplot2)
library(patchwork)  # For combining plots

# Implement LRF function
lrf_1991_exact <- function(temp, rmax, topt, tmin, tmax) {
  # Return 0 for temps outside valid range
  if (temp <= tmin || temp >= tmax) {
    return(0)
  }
  
  # Implement the exact equation
  numerator <- (temp - tmax) * (temp - tmin)^2
  denominator <- (topt - tmin) * ((topt - tmin) * (temp - topt) - (topt - tmax) * (topt + tmin - 2 * temp))
  
  rate <- rmax * numerator / denominator
  
  # Handle NaN or negative values
  if (is.na(rate) || rate < 0) {
    return(0)
  }
  
  return(rate)
}

# Define the parameters and credible intervals from the models
# For MS (from feed_mod3)
ms_params <- list(
  Eldo = list(
    Tmin = 14.74, Tmin_lower = 14.17, Tmin_upper = 15.20,
    Topt = 39.88, Topt_lower = 37.60, Topt_upper = 43.36,
    Above = 9.56, Above_lower = 4.15, Above_upper = 17.40,
    Ropt = 3.50, Ropt_lower = 2.73, Ropt_upper = 4.42
  ),
  A1 = list(
    Tmin = 11.67, Tmin_lower = 9.69, Tmin_upper = 13.22,
    Topt = 39.89, Topt_lower = 37.63, Topt_upper = 42.79,
    Above = 6.68, Above_lower = 2.67, Above_upper = 13.54,
    Ropt = 3.63, Ropt_lower = 2.84, Ropt_upper = 4.57
  ),
  B1 = list(
    Tmin = 10.96, Tmin_lower = 8.32, Tmin_upper = 13.09,
    Topt = 39.78, Topt_lower = 37.20, Topt_upper = 43.21,
    Above = 8.20, Above_lower = 3.06, Above_upper = 16.26,
    Ropt = 4.30, Ropt_lower = 3.41, Ropt_upper = 5.29
  )
)

# For MB (from feed_modMB)
mb_params <- list(
  A1 = list(
    Tmin = 13.57, Tmin_lower = 11.28, Tmin_upper = 15.10,
    Topt = 37.54, Topt_lower = 33.39, Topt_upper = 45.12,
    Above = 11.70, Above_lower = 4.09, Above_upper = 23.70,
    Ropt = 2.64, Ropt_lower = 1.61, Ropt_upper = 4.12
  ),
  B1 = list(
    Tmin = 12.80, Tmin_lower = 10.29, Tmin_upper = 14.84,
    Topt = 41.76, Topt_lower = 35.73, Topt_upper = 51.37,
    Above = 10.74, Above_lower = 1.56, Above_upper = 23.92,
    Ropt = 3.29, Ropt_lower = 2.02, Ropt_upper = 5.93
  ),
  C1 = list(
    Tmin = 13.79, Tmin_lower = 12.03, Tmin_upper = 15.08,
    Topt = 41.03, Topt_lower = 35.81, Topt_upper = 50.23,
    Above = 11.77, Above_lower = 2.44, Above_upper = 25.60,
    Ropt = 3.93, Ropt_lower = 2.52, Ropt_upper = 6.12
  )
)

# Add D1 parameters (made-up with reasonable credible intervals)
# Using values similar to our previous script but with credible intervals
mb_params$D1 <- list(
  Tmin = 14.00, Tmin_lower = 12.50, Tmin_upper = 15.50,
  Topt = 42.00, Topt_lower = 37.00, Topt_upper = 47.00,
  Above = 12.00, Above_lower = 5.00, Above_upper = 20.00,
  Ropt = 4.00, Ropt_lower = 2.80, Ropt_upper = 5.20
)

# Set site order and colors
site_order <- c("Eldo", "A1", "B1", "C1", "D1")
site_colors <- c(
  "Eldo" = "#FF4500",  # Red-orange
  "A1" = "#FF8C00",    # Dark orange
  "B1" = "#FFD700",    # Gold/yellow
  "C1" = "#4682B4",    # Steel blue
  "D1" = "#0000CD"     # Medium blue
)

# Function to calculate the curve
calculate_curve <- function(params, temp_range) {
  # For each temperature, calculate the rate
  result <- data.frame(
    temp = temp_range,
    rate = sapply(temp_range, function(t) {
      lrf_1991_exact(t, params$Ropt, params$Topt, params$Tmin, params$Topt + params$Above)
    })
  )
  
  return(result)
}

# Generate curve data for MS
generate_ms_curves <- function() {
  # Temperature range
  temp_range <- seq(5, 60, by = 0.1)
  
  # Create empty list to store results
  curves <- list()
  
  # Generate curves for each site
  for (site in names(ms_params)) {
    curves[[site]] <- calculate_curve(ms_params[[site]], temp_range)
    curves[[site]]$site <- site
  }
  
  # Combine all curves
  result <- bind_rows(curves)
  result$site <- factor(result$site, levels = site_order)
  
  return(result)
}

# Generate curve data for MB
generate_mb_curves <- function() {
  # Temperature range
  temp_range <- seq(5, 60, by = 0.1)
  
  # Create empty list to store results
  curves <- list()
  
  # Generate curves for each site
  for (site in names(mb_params)) {
    curves[[site]] <- calculate_curve(mb_params[[site]], temp_range)
    curves[[site]]$site <- site
  }
  
  # Combine all curves
  result <- bind_rows(curves)
  result$site <- factor(result$site, levels = site_order)
  
  return(result)
}

# Generate the curves
ms_curves <- generate_ms_curves()
mb_curves <- generate_mb_curves()

# Calculate the maximum rate value in both datasets for consistent scaling
max_rate_ms <- max(ms_curves$rate)
max_rate_mb <- max(mb_curves$rate)
max_rate <- max(max_rate_ms, max_rate_mb)

# Create parameter summary points for plot
create_param_summary <- function(params, species) {
  # Create a dataframe with all necessary columns
  summary_points <- data.frame(
    parameter = character(),
    site = character(),
    x = numeric(),
    x_lower = numeric(),
    x_upper = numeric(),
    y = numeric(),
    y_lower = numeric(),
    y_upper = numeric(),
    species = character(),
    stringsAsFactors = FALSE
  )
  
  for (site_name in names(params)) {
    site <- params[[site_name]]
    
    # Add Tmin point with error bars
    tmin_point <- data.frame(
      parameter = "Tmin",
      site = site_name,
      x = site$Tmin,
      x_lower = site$Tmin_lower,
      x_upper = site$Tmin_upper,
      y = 0,  # Base position (will be offset later)
      y_lower = NA,  # Add these for consistency
      y_upper = NA,  # Add these for consistency
      species = species,
      stringsAsFactors = FALSE
    )
    
    # Add Topt point with error bars
    topt_point <- data.frame(
      parameter = "Topt",
      site = site_name,
      x = site$Topt,
      x_lower = site$Topt_lower,
      x_upper = site$Topt_upper,
      y = site$Ropt,  # At y=Ropt for Topt
      y_lower = NA,  # Add these for consistency
      y_upper = NA,  # Add these for consistency
      species = species,
      stringsAsFactors = FALSE
    )
    
    # Add Tmax point with error bars
    tmax_point <- data.frame(
      parameter = "Tmax",
      site = site_name,
      x = site$Topt + site$Above,
      x_lower = site$Topt_lower + site$Above_lower,
      x_upper = site$Topt_upper + site$Above_upper,
      y = 0,  # Base position (will be offset later)
      y_lower = NA,  # Add these for consistency
      y_upper = NA,  # Add these for consistency
      species = species,
      stringsAsFactors = FALSE
    )
    
    # Add Ropt point with vertical error bars
    ropt_point <- data.frame(
      parameter = "Ropt",
      site = site_name,
      x = 65,  # Fixed position on the right side
      x_lower = NA,  # Add these for consistency
      x_upper = NA,  # Add these for consistency
      y = site$Ropt,
      y_lower = site$Ropt_lower,
      y_upper = site$Ropt_upper,
      species = species,
      stringsAsFactors = FALSE
    )
    
    # Bind rows to the main dataframe
    summary_points <- rbind(summary_points, tmin_point, topt_point, tmax_point, ropt_point)
  }
  
  summary_points$site <- factor(summary_points$site, levels = site_order)
  return(summary_points)
}

# Generate parameter summary points
ms_param_summary <- create_param_summary(ms_params, "MS")
mb_param_summary <- create_param_summary(mb_params, "MB")

# Function to create a TPC plot with improved offset credible intervals
create_tpc_plot <- function(curves, param_summary, title, max_y) {
  # Get numerical site indices (to help with ordered offsets)
  site_indices <- setNames(1:length(site_order), site_order)
  
  # Find the maximum rate in the curves to place Topt points above
  curve_max <- max(curves$rate)
  
  # Create a modified parameter summary with better offset y-positions
  modified_param_summary <- param_summary %>%
    mutate(
      # Set position index by site (ordered by elevation)
      site_idx = as.integer(site_indices[as.character(site)]),
      
      # Create different offsets for different parameter types
      y_display = case_when(
        # Tmin parameters go below the x-axis, ordered by site elevation with more spacing
        parameter == "Tmin" ~ -2.0 - (site_idx * 0.8),
        
        # Topt parameters now all above the curves with proper spacing
        parameter == "Topt" ~ curve_max + 1 + (site_idx * 0.8),
        
        # Tmax parameters go below the x-axis, but at different positions than Tmin
        parameter == "Tmax" ~ -2.0 - (site_idx * 0.8),
        
        # Ropt parameters go on the right side with horizontal spacing
        parameter == "Ropt" ~ y,
        
        TRUE ~ y
      ),
      
      # For Ropt, use a horizontal offset
      x_display = case_when(
        parameter == "Ropt" ~ 60 + (site_idx * 1.0),
        TRUE ~ x
      )
    )
  
  # Split out Ropt points for special handling
  ropt_points <- modified_param_summary %>% 
    filter(parameter == "Ropt")
  
  # Other parameter points
  other_points <- modified_param_summary %>% 
    filter(parameter != "Ropt")
  
  p <- ggplot() +
    # Add the curve
    geom_line(
      data = curves,
      aes(x = temp, y = rate, color = site),
      size = 1.5
    ) +
    # Add points for parameters (except Ropt)
    geom_point(
      data = other_points,
      aes(x = x, y = y_display, color = site),
      size = 3
    ) +
    # Add horizontal error bars for parameters (except Ropt)
    geom_errorbarh(
      data = other_points,
      aes(x = x, y = y_display, xmin = x_lower, xmax = x_upper, color = site),
      height = 0.2,
      size = 1
    ) +
    # Add Ropt points at the side with vertical error bars
    geom_point(
      data = ropt_points,
      aes(x = x_display, y = y, color = site),
      size = 3
    ) +
    # Add vertical error bars for Ropt
    geom_errorbar(
      data = ropt_points,
      aes(x = x_display, ymin = y_lower, ymax = y_upper, color = site),
      width = 0.5,
      size = 1
    ) +
    # Add dotted line at y=0
    geom_hline(yintercept = 0, linetype = "dotted", color = "darkgray") +
    # Set axis labels
    labs(
      x = "Temperature (°C)",
      y = "Feces mass / hopper mass / hour (mg/g/hr)",
      title = title
    ) +
    # Customize the color scheme
    scale_color_manual(values = site_colors) +
    # Set the x-axis range
    scale_x_continuous(limits = c(5, 65), breaks = seq(10, 60, by = 10)) +
    # Set the y-axis range, allowing space for the offset parameter points
    scale_y_continuous(limits = c(-6, max_y), breaks = seq(0, max_y, by = 2)) +
    # Theme customization
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Create the plots
ms_plot <- create_tpc_plot(
  ms_curves, 
  ms_param_summary, 
  "MS Thermal Performance Curves",
  max_y = 13  # Maximum y-axis value
)

mb_plot <- create_tpc_plot(
  mb_curves, 
  mb_param_summary, 
  "MB Thermal Performance Curves",
  max_y = 13  # Maximum y-axis value
)

# Combine plots
combined_plot <- ms_plot + mb_plot + 
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Thermal Performance Curves for Grasshopper Species",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18))
  )

# Create a legend-only plot for documentation
legend_data <- data.frame(
  site = site_order,
  x = rep(1, length(site_order)),
  y = 1:length(site_order)
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = site)) +
  geom_point(size = 3) +
  scale_color_manual(values = site_colors, name = "Site") +
  theme_void() +
  theme(legend.position = "bottom")

# Save just the legend
legend_only <- cowplot::get_legend(legend_plot)
ggsave("site_legend.png", plot = cowplot::ggdraw(legend_only), width = 6, height = 2, dpi = 300)

# Save the plots
#ggsave("ms_thermal_performance_final.png", ms_plot, width = 8, height = 6, dpi = 300)
#ggsave("mb_thermal_performance_final.png", mb_plot, width = 8, height = 6, dpi = 300)
#ggsave("combined_thermal_performance_final.png", combined_plot, width = 16, height = 8, dpi = 300)

# Display the plots
print(ms_plot)
print(mb_plot)
print(combined_plot)


