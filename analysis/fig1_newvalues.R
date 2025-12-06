# Load required libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
library(grid)
library(gtable)

ad_data <- read_csv("C:/Users/jmsmi/OneDrive/Documents/GitHub/thermal_perf/data/ad.csv")

#look at sample size in various ways
ad_data %>% group_by(spp, site, sex, temp) %>% summarize(n())

#and get the number of individuals by site and sex
ad_data %>% group_by(spp, site, sex) %>% summarize(length(unique(full_ID)))

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

# Site coordinates with elevations
site_coords <- data.frame(
  name = c("Eldo", "A1", "B1", "C1", "D1"),
  lat = c(39.944, 40.015, 40.023, 40.036, 40.059),
  lon = c(-105.262, -105.377, -105.430, -105.547, -105.617),
  elev = c(1740, 2195, 2591, 3048, 3515)
)

# Create elevation labels
site_coords$elev_label <- paste0(site_coords$elev, "m")

# Function to calculate feeding rate for MS using exact LRF model (averaged across sexes)
get_ms_feeding_rate <- function(temp, site) {
  # Parameters for MS from the model output
  site_params <- list(
    Eldo = list(Tmin = 12.63, Topt = 40.69, Above = 10.07, Ropt = 3.40),
    A1 = list(Tmin = 12.37, Topt = 42.17, Above = 9.02, Ropt = 4.13),
    B1 = list(Tmin = 11.57, Topt = 39.68, Above = 7.56, Ropt = 4.89)
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
    return(0)
  })
  
  return(base_rate)
}

# Function to calculate feeding rate for MB using exact LRF model (averaged across sexes)
get_mb_feeding_rate <- function(temp, site) {
  # Parameters for MB from the model output
  site_params <- list(
    A1 = list(Tmin = 10.73, Topt = 39.52, Above = 16.73, Ropt = 2.99),
    B1 = list(Tmin = 11.31, Topt = 40.86, Above = 12.3, Ropt = 3.93),
    C1 = list(Tmin = 11.07, Topt = 40.76, Above = 12.91, Ropt = 4.09),
    D1 = list(Tmin = 10.86, Topt = 40.08, Above = 13.84, Ropt = 4.20)
  )
  
  # For Eldo, use A1 parameters as placeholder
  site_params$Eldo <- site_params$A1
  
  # Get parameters for this site
  params <- site_params[[site]]
  
  # Calculate Tmax (Topt + Above)
  Tmax <- params$Topt + params$Above
  
  # Calculate base feeding rate using exact LRF
  base_rate <- tryCatch({
    lrf_1991_exact(temp, params$Ropt, params$Topt, params$Tmin, Tmax)
  }, error = function(e) {
    return(0)
  })
  
  return(base_rate)
}

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

# Combine original data with curve data
combined_data <- rbind(
  ad_data[, c("spp", "site", "sex", "temp", "rate", "type")],
  transform(curve_data, sex = NA)
)

# Map site names to elevation labels
site_name_to_elev <- setNames(site_coords$elev_label, site_coords$name)
combined_data$site_label <- site_name_to_elev[as.character(combined_data$site)]

# Set the site order based on elevation (low to high)
site_order <- c("Eldo", "A1", "B1", "C1", "D1")
elev_order <- site_coords$elev_label[match(site_order, site_coords$name)]
combined_data$site <- factor(combined_data$site, levels = site_order)
combined_data$site_label <- factor(combined_data$site_label, levels = elev_order)

# Define colors for sites (cooler as elevation increases)
site_colors <- c(
  "Eldo" = "#FF4500",  # Red-orange
  "A1" = "#FF8C00",    # Dark orange
  "B1" = "#FFD700",    # Gold/yellow
  "C1" = "#4682B4",    # Steel blue
  "D1" = "#0000CD"     # Medium blue
)

# Create color mapping for elevation labels
elev_colors <- setNames(site_colors, site_name_to_elev[names(site_colors)])

# Filter points and curves
point_data <- combined_data %>% filter(type == "point")
curve_data <- combined_data %>% filter(type == "curve")

# Create the plot with faceting
p <- ggplot() +
  # Add curves
  geom_line(
    data = curve_data,
    aes(x = temp, y = rate, color = site_label),
    size = 1.2
  ) +
  # Add points with jitter
  geom_jitter(
    data = point_data,
    aes(x = temp, y = rate, color = site_label),
    width = 0.5, height = 0, alpha = 0.6, size = 1.8, shape = 16
  ) +
  # Set up faceting by site (rows) and species (columns)
  facet_grid(site_label ~ spp,
             labeller = labeller(
               spp = as_labeller(c("MB" = "bolditalic('M. boulderensis')", 
                                       "MS" = "bolditalic('M. sanguinipes')"), 
                                     label_parsed))) +
  # Set axis labels
  labs(
    x = "Temperature (°C)",
    y = "Mass-adjusted feces production rate (mg/g hopper/hr)"
  ) +
  # Customize the color scheme
  scale_color_manual(
    values = elev_colors,
    name = "Site"
  ) +
  # Set the x-axis range
  scale_x_continuous(limits = c(15, 45), breaks = seq(15, 45, by = 5)) +
  # Set the y-axis range
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  # Theme customization
  theme_minimal() +
  theme(
    strip.background.x = element_rect(fill = "#D2B48C", color = NA),  # Brown for species
    strip.background.y = element_rect(fill = "#D2B48C", color = NA),  # Brown for sites
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none",  # Remove legend
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

# Convert to grob for adding labels
g <- ggplotGrob(p)

# Add "Species" label at the top
species_label <- textGrob("Species", gp = gpar(fontsize = 12, fontface = "bold", col = "#8B4513"))
g <- gtable_add_rows(g, heights = unit(0.6, "cm"), pos = 0)
panel_cols <- which(grepl("panel", g$layout$name))
g <- gtable_add_grob(g, species_label, t = 1, 
                     l = min(g$layout$l[panel_cols]),
                     r = max(g$layout$r[panel_cols]))

# Add "Site" label on the right
site_label_grob <- textGrob("Population", rot = 270, gp = gpar(fontsize = 12, fontface = "bold", col = "#8B4513"))
g <- gtable_add_cols(g, widths = unit(0.6, "cm"), pos = -1)
panel_rows <- which(grepl("panel", g$layout$name))
g <- gtable_add_grob(g, site_label_grob, 
                     t = min(g$layout$t[panel_rows]),
                     b = max(g$layout$b[panel_rows]),
                     l = ncol(g), r = ncol(g))

# Display the plot
grid.draw(g)

# Save the plot
# ggsave("grasshopper_tpc_faceted.png", plot = g, width = 10, height = 12, dpi = 300)





#better B

library(tidyverse)
library(ggplot2)
library(patchwork)
library(grid)
library(gtable)

# ===== STEP 1: Hard-coded 90% ETI parameters =====

# For MB model
mb_params <- list(
  A1 = list(
    Tmin = 10.727636, Tmin_lower = 7.683487, Tmin_upper = 13.465337,
    Topt = 39.524259, Topt_lower = 36.128778, Topt_upper = 44.445992,
    Above = 16.727712,
    Above_lower = 11.070841,
    Above_upper = 33.106123,
    Ropt = 2.990098, Ropt_lower = 2.484207, Ropt_upper = 3.572438,
    Tmax = 56.252348, Tmax_lower = 47.199619, Tmax_upper = 69.234901,
    Tb = 45.524712
  ),
  B1 = list(
    Tmin = 11.314475, Tmin_lower = 9.406580, Tmin_upper = 13.026377,
    Topt = 40.856685, Topt_lower = 37.824315, Topt_upper = 45.531373,
    Above = 12.298740,
    Above_lower = 8.875827,
    Above_upper = 27.494798,
    Ropt = 3.927178, Ropt_lower = 3.220651, Ropt_upper = 4.734666,
    Tmax = 53.155425, Tmax_lower = 45.700142, Tmax_upper = 65.319113,
    Tb = 41.840950
  ),
  C1 = list(
    Tmin = 11.066096, Tmin_lower = 9.210138, Tmin_upper = 12.796826,
    Topt = 40.758789, Topt_lower = 37.831284, Topt_upper = 45.413410,
    Above = 12.911405,
    Above_lower = 8.050695,
    Above_upper = 28.357692,
    Ropt = 4.093169, Ropt_lower = 3.435337, Ropt_upper = 4.858143,
    Tmax = 53.670194, Tmax_lower = 45.881979, Tmax_upper = 66.188976,
    Tb = 42.604099
  ),
  D1 = list(
    Tmin = 10.855164, Tmin_lower = 7.486004, Tmin_upper = 13.611458,
    Topt = 40.083176, Topt_lower = 35.980420, Topt_upper = 45.607234,
    Above = 13.836828,
    Above_lower = 9.631193,
    Above_upper = 31.087635,
    Ropt = 4.203449, Ropt_lower = 3.230990, Ropt_upper = 5.401352,
    Tmax = 53.920004, Tmax_lower = 45.611613, Tmax_upper = 67.068055,
    Tb = 43.064840
  )
)

# For MS model
ms_params <- list(
  Eldo = list(
    Tmin = 12.632975, Tmin_lower = 11.245301, Tmin_upper = 13.785139,
    Topt = 40.688372, Topt_lower = 38.271171, Topt_upper = 44.655864,
    Above = 10.067734,
    Above_lower = 6.293769,
    Above_upper = 23.398960,
    Ropt = 3.400167, Ropt_lower = 2.768452, Ropt_upper = 4.119330,
    Tmax = 50.756106, Tmax_lower = 44.564940, Tmax_upper = 61.670131,
    Tb = 38.123131
  ),
  A1 = list(
    Tmin = 12.373525, Tmin_lower = 11.215186, Tmin_upper = 13.383920,
    Topt = 42.166751, Topt_lower = 39.205276, Topt_upper = 47.115379,
    Above = 9.023474,
    Above_lower = 4.532040,
    Above_upper = 22.362056,
    Ropt = 4.126447, Ropt_lower = 3.364652, Ropt_upper = 5.044721,
    Tmax = 51.190225, Tmax_lower = 44.737316, Tmax_upper = 61.567332,
    Tb = 38.816700
  ),
  B1 = list(
    Tmin = 11.567755, Tmin_lower = 10.120242, Tmin_upper = 12.763509,
    Topt = 39.678300, Topt_lower = 37.900113, Topt_upper = 42.099225,
    Above = 7.556771,
    Above_lower = 5.875222,
    Above_upper = 17.607547,
    Ropt = 4.887275, Ropt_lower = 3.939816, Ropt_upper = 6.005675,
    Tmax = 47.235071, Tmax_lower = 43.775347, Tmax_upper = 55.507660,
    Tb = 35.667317
  )
)

# ===== STEP 2: Rest of plotting code =====

# Implement LRF function
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

# Set site order and colors
site_order <- c("Eldo", "A1", "B1", "C1", "D1")
site_colors <- c(
  "Eldo" = "#FF4500",
  "A1" = "#FF8C00",
  "B1" = "#FFD700",
  "C1" = "#4682B4",
  "D1" = "#0000CD"
)

# Site elevations for legend
site_elevations <- c(
  "Eldo" = "1740m",
  "A1" = "2195m",
  "B1" = "2591m",
  "C1" = "3048m",
  "D1" = "3515m"
)

# Function to calculate the curve
calculate_curve <- function(params, temp_range) {
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
  temp_range <- seq(5, 70, by = 0.1)
  curves <- list()
  
  for (site in names(ms_params)) {
    curves[[site]] <- calculate_curve(ms_params[[site]], temp_range)
    curves[[site]]$site <- site
    curves[[site]]$species <- "MS"
  }
  
  result <- bind_rows(curves)
  result$site <- factor(result$site, levels = site_order)
  return(result)
}

# Generate curve data for MB
generate_mb_curves <- function() {
  temp_range <- seq(5, 70, by = 0.1)
  curves <- list()
  
  for (site in names(mb_params)) {
    curves[[site]] <- calculate_curve(mb_params[[site]], temp_range)
    curves[[site]]$site <- site
    curves[[site]]$species <- "MB"
  }
  
  result <- bind_rows(curves)
  result$site <- factor(result$site, levels = site_order)
  return(result)
}

# Generate the curves
ms_curves <- generate_ms_curves()
mb_curves <- generate_mb_curves()

# Combine both datasets
all_curves <- bind_rows(mb_curves, ms_curves)
all_curves$species <- factor(all_curves$species, levels = c("MB", "MS"))

# Create parameter summary points for plot
create_param_summary <- function(params, species) {
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
    
    tmin_point <- data.frame(
      parameter = "Tmin",
      site = site_name,
      x = site$Tmin,
      x_lower = site$Tmin_lower,
      x_upper = site$Tmin_upper,
      y = 0,
      y_lower = NA,
      y_upper = NA,
      species = species,
      stringsAsFactors = FALSE
    )
    
    topt_point <- data.frame(
      parameter = "Topt",
      site = site_name,
      x = site$Topt,
      x_lower = site$Topt_lower,
      x_upper = site$Topt_upper,
      y = site$Ropt,
      y_lower = NA,
      y_upper = NA,
      species = species,
      stringsAsFactors = FALSE
    )
    
    tmax_point <- data.frame(
      parameter = "Tmax",
      site = site_name,
      x = site$Tmax,
      x_lower = site$Tmax_lower,
      x_upper = site$Tmax_upper,
      y = 0,
      y_lower = NA,
      y_upper = NA,
      species = species,
      stringsAsFactors = FALSE
    )
    
    ropt_point <- data.frame(
      parameter = "Ropt",
      site = site_name,
      x = 65,
      x_lower = NA,
      x_upper = NA,
      y = site$Ropt,
      y_lower = site$Ropt_lower,
      y_upper = site$Ropt_upper,
      species = species,
      stringsAsFactors = FALSE
    )
    
    summary_points <- rbind(summary_points, tmin_point, topt_point, tmax_point, ropt_point)
  }
  
  summary_points$site <- factor(summary_points$site, levels = site_order)
  return(summary_points)
}

# Generate parameter summary points
ms_param_summary <- create_param_summary(ms_params, "MS")
mb_param_summary <- create_param_summary(mb_params, "MB")

# Combine parameter summaries
all_param_summary <- bind_rows(mb_param_summary, ms_param_summary)
all_param_summary$species <- factor(all_param_summary$species, levels = c("MB", "MS"))

# Calculate max rate for each species separately
max_rate_mb <- max(mb_curves$rate)
max_rate_ms <- max(ms_curves$rate)

# Create modified parameter summary with offsets PER SPECIES
modified_param_summary <- all_param_summary %>%
  group_by(species) %>%
  mutate(
    # Get sites present in this species
    sites_present = list(unique(site[!is.na(x)])),
    # Assign indices based on order within this species only
    site_idx = match(site, unique(site[!is.na(x)])),
    # Get max rate for this species
    species_max_rate = ifelse(species == "MB", max_rate_mb, max_rate_ms)
  ) %>%
  ungroup() %>%
  mutate(
    y_display = case_when(
      parameter == "Tmin" ~ -0.5 - (site_idx * 0.3),
      parameter == "Topt" ~ species_max_rate + 0.5 + (site_idx * 0.3),
      parameter == "Tmax" ~ -0.5 - (site_idx * 0.3),
      parameter == "Ropt" ~ y,
      TRUE ~ y
    ),
    
    # Increase spacing for Ropt points
    x_display = case_when(
      parameter == "Ropt" ~ 55 + (site_idx * 2.5),
      TRUE ~ x
    )
  )

# Split out Ropt points
ropt_points <- modified_param_summary %>% 
  filter(parameter == "Ropt")

# Other parameter points
other_points <- modified_param_summary %>% 
  filter(parameter != "Ropt")

# Create the faceted plot
p <- ggplot() +
  geom_line(
    data = all_curves,
    aes(x = temp, y = rate, color = site),
    size = 1.5
  ) +
  geom_point(
    data = other_points,
    aes(x = x, y = y_display, color = site),
    size = 3
  ) +
  geom_errorbarh(
    data = other_points,
    aes(x = x, y = y_display, xmin = x_lower, xmax = x_upper, color = site),
    height = 0.2,
    size = 1
  ) +
  geom_point(
    data = ropt_points,
    aes(x = x_display, y = y, color = site),
    size = 3
  ) +
  geom_errorbar(
    data = ropt_points,
    aes(x = x_display, ymin = y_lower, ymax = y_upper, color = site),
    width = 2,
    size = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "darkgray") +
  facet_wrap(~ species, ncol = 2,
             labeller = labeller(
               species = as_labeller(c("MB" = "bolditalic('M. boulderensis')", 
                                       "MS" = "bolditalic('M. sanguinipes')"), 
                                     label_parsed))) +
  labs(
    x = "Temperature (°C)",
    y = "Feces mass / hopper mass / hour (mg/g/hr)",
    color = "Population"
  ) +
  scale_color_manual(
    values = site_colors,
    labels = site_elevations,
    breaks = names(site_elevations)
  ) +
  scale_x_continuous(limits = c(5, 70), breaks = seq(10, 70, by = 10)) +
  scale_y_continuous(limits = c(-2, 8), breaks = seq(0, 10, by = 2)) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "#D2B48C", color = NA),
    strip.text = element_text(face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Add custom legend in top left of MS panel only
# Create legend data for MS panel only
legend_data <- data.frame(
  site = factor(names(site_elevations), levels = names(site_elevations)),
  elev = site_elevations,
  y = seq(6.75, 4.75, length.out = 5),  # Lower position with more spacing
  x = 9.75,  # Left side of plot
  species = "MS"  # Only in MS panel
)

# Legend box coordinates
legend_box <- data.frame(
  xmin = 5.5,
  xmax = 24,
  ymin = 4.25,
  ymax = 7.5,
  species = "MS"
)

# Legend title data (only for MS panel)
legend_title <- data.frame(
  x = 7,
  y = 7.25,
  label = "Population",
  species = "MS"
)

# Add legend to the plot (only appears in MS facet due to species filter)
p <- p +
  # Legend box
  geom_rect(data = legend_box,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", color = "black", linewidth = 0.5,
            inherit.aes = FALSE) +
  # Legend title (only in MS panel)
  geom_text(data = legend_title,
            aes(x = x, y = y, label = label),
            fontface = "bold", size = 4, hjust = 0,
            inherit.aes = FALSE) +
  # Legend items - only in MS panel
  geom_point(data = legend_data, 
             aes(x = x, y = y, color = site),
             size = 3, inherit.aes = FALSE) +
  geom_text(data = legend_data,
            aes(x = x + 1, y = y, label = elev),
            hjust = 0, size = 3.5, inherit.aes = FALSE)

# Convert to grob for adding "Species" label
g <- ggplotGrob(p)

# Add "Species" label at the top
species_label <- textGrob("Species", gp = gpar(fontsize = 12, fontface = "bold", col = "#8B4513"))
g <- gtable_add_rows(g, heights = unit(0.6, "cm"), pos = 0)
panel_cols <- which(grepl("panel", g$layout$name))
g <- gtable_add_grob(g, species_label, t = 1, 
                     l = min(g$layout$l[panel_cols]),
                     r = max(g$layout$r[panel_cols]))

# Display the plot
grid.draw(g)

# Or save it
# ggsave("thermal_performance_faceted.png", plot = g, width = 12, height = 6, dpi = 300)

















#new
# plot_tpc_figures.R
# Script to create TPC figures using MAP+HDI parameters from pops.rds

library(tidyverse)
library(ggplot2)
library(patchwork)
library(grid)
library(gtable)

# Load the pops data
pops <- readRDS("data/pops.rds")

# Load the raw data for plotting points
ad_data <- read_csv("C:/Users/jmsmi/OneDrive/Documents/GitHub/thermal_perf/data/ad.csv")

# ===== HELPER FUNCTIONS =====

# Implement LRF function
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

# Site coordinates with elevations
site_coords <- data.frame(
  name = c("Eldo", "A1", "B1", "C1", "D1"),
  lat = c(39.944, 40.015, 40.023, 40.036, 40.059),
  lon = c(-105.262, -105.377, -105.430, -105.547, -105.617),
  elev = c(1740, 2195, 2591, 3048, 3515)
)
site_coords$elev_label <- paste0(site_coords$elev, "m")

# Site order and colors
site_order <- c("Eldo", "A1", "B1", "C1", "D1")
site_colors <- c(
  "Eldo" = "#FF4500",
  "A1" = "#FF8C00",
  "B1" = "#FFD700",
  "C1" = "#4682B4",
  "D1" = "#0000CD"
)
site_elevations <- c(
  "Eldo" = "1740m",
  "A1" = "2195m",
  "B1" = "2591m",
  "C1" = "3048m",
  "D1" = "3515m"
)

# ===== FIGURE A: DATA WITH FITTED CURVES =====

# Extract unique TPC parameters from pops (one row per site/species combination)
tpc_params <- pops %>%
  select(spp, site, Tmin, Topt, Above, Ropt, Tmax) %>%
  distinct()

# Function to calculate feeding rate using parameters from pops
get_feeding_rate <- function(temp, site, spp, params_df) {
  params <- params_df %>% filter(site == !!site, spp == !!spp)
  
  if (nrow(params) == 0) return(0)
  
  base_rate <- tryCatch({
    lrf_1991_exact(temp, params$Ropt, params$Topt, params$Tmin, params$Topt + params$Above)
  }, error = function(e) {
    return(0)
  })
  
  return(base_rate)
}

# Generate curve data
generate_curve_data <- function(params_df) {
  temp_range <- seq(15, 45, by = 0.1)
  
  mb_sites <- c("A1", "B1", "C1", "D1")
  ms_sites <- c("Eldo", "A1", "B1")
  
  curve_data <- data.frame()
  
  # Generate data for MB
  for (site in mb_sites) {
    site_curve <- data.frame(
      spp = "MB",
      site = site,
      temp = temp_range,
      rate = sapply(temp_range, function(t) get_feeding_rate(t, site, "MB", params_df)),
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
      rate = sapply(temp_range, function(t) get_feeding_rate(t, site, "MS", params_df)),
      type = "curve"
    )
    curve_data <- rbind(curve_data, site_curve)
  }
  
  return(curve_data)
}

# Generate the curve data
curve_data <- generate_curve_data(tpc_params)

# Add type column to the original data
ad_data$type <- "point"

# Combine original data with curve data
combined_data <- rbind(
  ad_data[, c("spp", "site", "sex", "temp", "rate", "type")],
  transform(curve_data, sex = NA)
)

# Map site names to elevation labels
site_name_to_elev <- setNames(site_coords$elev_label, site_coords$name)
combined_data$site_label <- site_name_to_elev[as.character(combined_data$site)]

combined_data$site <- factor(combined_data$site, levels = site_order)
elev_order <- site_coords$elev_label[match(site_order, site_coords$name)]
combined_data$site_label <- factor(combined_data$site_label, levels = elev_order)

# Create color mapping for elevation labels
elev_colors <- setNames(site_colors, site_name_to_elev[names(site_colors)])

# Filter points and curves
point_data <- combined_data %>% filter(type == "point")
curve_data_plot <- combined_data %>% filter(type == "curve")

# Create Figure A
fig_a <- ggplot() +
  geom_line(
    data = curve_data_plot,
    aes(x = temp, y = rate, color = site_label),
    size = 1.2
  ) +
  geom_jitter(
    data = point_data,
    aes(x = temp, y = rate, color = site_label),
    width = 0.5, height = 0, alpha = 0.6, size = 1.8, shape = 16
  ) +
  facet_grid(site_label ~ spp,
             labeller = labeller(
               spp = as_labeller(c("MB" = "bolditalic('M. boulderensis')", 
                                   "MS" = "bolditalic('M. sanguinipes')"), 
                                 label_parsed))) +
  labs(
    x = "Temperature (°C)",
    y = "Mass-adjusted feces production rate (mg/g hopper/hr)"
  ) +
  scale_color_manual(
    values = elev_colors,
    name = "Population"
  ) +
  scale_x_continuous(limits = c(15, 45), breaks = seq(15, 45, by = 5)) +
  scale_y_continuous(limits = c(0, 13), breaks = seq(0, 12, by = 2), 
                     labels = c("0", "", "4", "", "8", "", "12")) +
  theme_minimal() +
  theme(
    strip.background.x = element_rect(fill = "#D2B48C", color = NA),
    strip.background.y = element_rect(fill = "#D2B48C", color = NA),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

# Convert to grob for adding labels
g_a <- ggplotGrob(fig_a)

# Add "Species" label at the top
species_label <- textGrob("Species", gp = gpar(fontface = "bold", col = "#8B4513")) #fontsize = 12, 
g_a <- gtable_add_rows(g_a, heights = unit(0.6, "cm"), pos = 0)
panel_cols <- which(grepl("panel", g_a$layout$name))
g_a <- gtable_add_grob(g_a, species_label, t = 1, 
                       l = min(g_a$layout$l[panel_cols]),
                       r = max(g_a$layout$r[panel_cols]))

# Add "Population" label on the right
site_label_grob <- textGrob("Population", rot = 270, gp = gpar(fontface = "bold", col = "#8B4513")) #fontsize = 12, 
g_a <- gtable_add_cols(g_a, widths = unit(0.6, "cm"), pos = -1)
panel_rows <- which(grepl("panel", g_a$layout$name))
g_a <- gtable_add_grob(g_a, site_label_grob, 
                       t = min(g_a$layout$t[panel_rows]),
                       b = max(g_a$layout$b[panel_rows]),
                       l = ncol(g_a), r = ncol(g_a))

# ===== FIGURE B: PARAMETER ESTIMATES WITH CREDIBLE INTERVALS =====

# Extract parameters with credible intervals from pops
# Create a simpler structure - list of lists
mb_params_df <- pops %>%
  filter(spp == "MB") %>%
  select(site, Tmin, Tmin_lower, Tmin_upper, Topt, Topt_lower, Topt_upper,
         Above, Above_lower, Above_upper, Ropt, Ropt_lower, Ropt_upper,
         Tmax, Tmax_lower, Tmax_upper) %>%
  distinct()

ms_params_df <- pops %>%
  filter(spp == "MS") %>%
  select(site, Tmin, Tmin_lower, Tmin_upper, Topt, Topt_lower, Topt_upper,
         Above, Above_lower, Above_upper, Ropt, Ropt_lower, Ropt_upper,
         Tmax, Tmax_lower, Tmax_upper) %>%
  distinct()

# Convert to list format expected by plotting functions
mb_params <- list()
for (i in 1:nrow(mb_params_df)) {
  site_name <- mb_params_df$site[i]
  mb_params[[site_name]] <- as.list(mb_params_df[i, ])
}

ms_params <- list()
for (i in 1:nrow(ms_params_df)) {
  site_name <- ms_params_df$site[i]
  ms_params[[site_name]] <- as.list(ms_params_df[i, ])
}

# Function to calculate the curve
calculate_curve <- function(params, temp_range) {
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
  temp_range <- seq(5, 70, by = 0.1)
  curves <- list()
  
  for (site in names(ms_params)) {
    curves[[site]] <- calculate_curve(ms_params[[site]], temp_range)
    curves[[site]]$site <- site
    curves[[site]]$species <- "MS"
  }
  
  result <- bind_rows(curves)
  result$site <- factor(result$site, levels = site_order)
  return(result)
}

# Generate curve data for MB
generate_mb_curves <- function() {
  temp_range <- seq(5, 70, by = 0.1)
  curves <- list()
  
  for (site in names(mb_params)) {
    curves[[site]] <- calculate_curve(mb_params[[site]], temp_range)
    curves[[site]]$site <- site
    curves[[site]]$species <- "MB"
  }
  
  result <- bind_rows(curves)
  result$site <- factor(result$site, levels = site_order)
  return(result)
}

# Generate the curves
ms_curves <- generate_ms_curves()
mb_curves <- generate_mb_curves()

# Combine both datasets
all_curves <- bind_rows(mb_curves, ms_curves)
all_curves$species <- factor(all_curves$species, levels = c("MB", "MS"))

# Create parameter summary points for plot
create_param_summary <- function(params, species) {
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
    
    tmin_point <- data.frame(
      parameter = "Tmin",
      site = site_name,
      x = as.numeric(site$Tmin),
      x_lower = as.numeric(site$Tmin_lower),
      x_upper = as.numeric(site$Tmin_upper),
      y = 0,
      y_lower = NA,
      y_upper = NA,
      species = species,
      stringsAsFactors = FALSE
    )
    
    topt_point <- data.frame(
      parameter = "Topt",
      site = site_name,
      x = as.numeric(site$Topt),
      x_lower = as.numeric(site$Topt_lower),
      x_upper = as.numeric(site$Topt_upper),
      y = as.numeric(site$Ropt),
      y_lower = NA,
      y_upper = NA,
      species = species,
      stringsAsFactors = FALSE
    )
    
    tmax_point <- data.frame(
      parameter = "Tmax",
      site = site_name,
      x = as.numeric(site$Tmax),
      x_lower = as.numeric(site$Tmax_lower),
      x_upper = as.numeric(site$Tmax_upper),
      y = 0,
      y_lower = NA,
      y_upper = NA,
      species = species,
      stringsAsFactors = FALSE
    )
    
    ropt_point <- data.frame(
      parameter = "Ropt",
      site = site_name,
      x = 65,
      x_lower = NA,
      x_upper = NA,
      y = as.numeric(site$Ropt),
      y_lower = as.numeric(site$Ropt_lower),
      y_upper = as.numeric(site$Ropt_upper),
      species = species,
      stringsAsFactors = FALSE
    )
    
    summary_points <- rbind(summary_points, tmin_point, topt_point, tmax_point, ropt_point)
  }
  
  summary_points$site <- factor(summary_points$site, levels = site_order)
  return(summary_points)
}

# Generate parameter summary points
ms_param_summary <- create_param_summary(ms_params, "MS")
mb_param_summary <- create_param_summary(mb_params, "MB")

# Combine parameter summaries
all_param_summary <- bind_rows(mb_param_summary, ms_param_summary)
all_param_summary$species <- factor(all_param_summary$species, levels = c("MB", "MS"))

# Calculate max rate for each species separately
max_rate_mb <- max(mb_curves$rate)
max_rate_ms <- max(ms_curves$rate)

# Create modified parameter summary with offsets PER SPECIES
modified_param_summary <- all_param_summary %>%
  group_by(species) %>%
  mutate(
    sites_present = list(unique(site[!is.na(x)])),
    l_sites_pres = length(unique(site[!is.na(x)])),
    site_idx = match(site, unique(site[!is.na(x)])),
    species_max_rate = ifelse(species == "MB", max_rate_mb, max_rate_ms)
  ) %>%
  ungroup() %>%
  mutate(
    y_display = case_when(
      parameter == "Tmin" ~ -0.2 - (site_idx * 0.4),
      parameter == "Topt" ~ species_max_rate + 0.3 + 0.4*(l_sites_pres+1) - (site_idx * 0.4),
      parameter == "Tmax" ~ -0.2 - (site_idx * 0.4),
      parameter == "Ropt" ~ y,
      TRUE ~ y
    ),
    x_display = case_when(
      parameter == "Ropt" ~ 55 + (site_idx * 2.5),
      TRUE ~ x
    )
  )

# Split out Ropt points
ropt_points <- modified_param_summary %>% 
  filter(parameter == "Ropt")

# Other parameter points
other_points <- modified_param_summary %>% 
  filter(parameter != "Ropt")

# Create Figure B
fig_b <- ggplot() +
  geom_line(
    data = all_curves,
    aes(x = temp, y = rate, color = site),
    size = 1.5
  ) +
  geom_point(
    data = other_points,
    aes(x = x, y = y_display, color = site),
    size = 3
  ) +
  geom_errorbarh(
    data = other_points,
    aes(x = x, y = y_display, xmin = x_lower, xmax = x_upper, color = site),
    height = 0.2,
    size = 1
  ) +
  geom_point(
    data = ropt_points,
    aes(x = x_display, y = y, color = site),
    size = 3
  ) +
  geom_errorbar(
    data = ropt_points,
    aes(x = x_display, ymin = y_lower, ymax = y_upper, color = site),
    width = 2,
    size = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "darkgray") +
  facet_wrap(~ species, ncol = 2,
             labeller = labeller(
               species = as_labeller(c("MB" = "bolditalic('M. boulderensis')", 
                                       "MS" = "bolditalic('M. sanguinipes')"), 
                                     label_parsed))) +
  labs(
    x = "Temperature (°C)",
    y = "Feces mass / hopper mass / hour (mg/g/hr)",
    color = "Population"
  ) +
  scale_color_manual(
    values = site_colors,
    labels = site_elevations,
    breaks = names(site_elevations)
  ) +
  scale_x_continuous(limits = c(5, 70), breaks = seq(10, 70, by = 10)) +
  scale_y_continuous(limits = c(-2, 6.5), breaks = seq(0, 10, by = 2)) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "#D2B48C", color = NA),
    strip.text = element_text(size=12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Add custom legend in top left of MS panel only
legend_data <- data.frame(
  site = factor(names(site_elevations), levels = names(site_elevations)),
  elev = site_elevations,
  y = seq(5.25, 3, length.out = 5),
  x = 9.75,
  species = "MS"
)

legend_box <- data.frame(
  xmin = 6,
  xmax = 24,
  ymin = 2.5,
  ymax = 6.25,
  species = "MS"
)

legend_title <- data.frame(
  x = 7,
  y = 6,
  label = "Population",
  species = "MS"
)

fig_b <- fig_b +
  geom_rect(data = legend_box,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", color = "black", linewidth = 0.5,
            inherit.aes = FALSE) +
  geom_text(data = legend_title,
            aes(x = x, y = y, label = label),
            fontface = "bold", size = 4, hjust = 0,
            inherit.aes = FALSE) +
  geom_point(data = legend_data, 
             aes(x = x, y = y, color = site),
             size = 3, inherit.aes = FALSE) +
  geom_text(data = legend_data,
            aes(x = x + 1, y = y, label = elev),
            hjust = 0, size = 3.5, inherit.aes = FALSE)

# Convert to grob for adding "Species" label
g_b <- ggplotGrob(fig_b)

# Add "Species" label at the top
species_label_b <- textGrob("Species", gp = gpar(fontface = "bold", col = "#8B4513")) #fontsize = 12, 
g_b <- gtable_add_rows(g_b, heights = unit(0.6, "cm"), pos = 0)
panel_cols_b <- which(grepl("panel", g_b$layout$name))
g_b <- gtable_add_grob(g_b, species_label_b, t = 1, 
                       l = min(g_b$layout$l[panel_cols_b]),
                       r = max(g_b$layout$r[panel_cols_b]))

# ===== COMBINE FIGURES A AND B =====

# Convert grobs to ggplot objects for patchwork
fig_a_gg <- wrap_elements(g_a)
fig_b_gg <- wrap_elements(g_b)

# Combine with patchwork - make A twice as tall as B
combined_figure <- fig_a_gg / fig_b_gg +
  plot_layout(heights = c(4, 3)) +  # A is 2 units tall, B is 1 unit tall
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Display the combined figure
print(combined_figure)

# Save the combined figure
#ggsave("figures/combined_tpc_figure.png", combined_figure, width = 12, height = 18, dpi = 300)

#cat("Combined figure created and saved to figures/combined_tpc_figure.png\n")