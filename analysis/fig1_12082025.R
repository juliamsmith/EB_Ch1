# plot_tpc_figures_modified.R
# Script to create TPC figures using MAP+HDI parameters from pops.rds
# Modified version with:
# 1) Reversed site order in panel A (highest elevation at top)
# 2) Legend in empty facet of panel A
# 3) Flipped credible interval order in B (except Pmax)
# 4) Parameter labels (Tmin, Tmax, Topt, Pmax) in panel B
# 5) Consistent y-axis labels
# 6) Appropriate ggsave dimensions

library(tidyverse)
library(ggplot2)
library(patchwork)
library(grid)
library(gtable)

setwd("C:/Users/jmsmi/OneDrive/Documents/GitHub/EB_Ch1/analysis")

# Load the pops data
pops <- readRDS("C:/Users/jmsmi/OneDrive/Documents/GitHub/EB_Ch1/data/pops.rds")

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

# Site order: low to high for data/legend, but we'll reverse for facet display in A
site_order <- c("Eldo", "A1", "B1", "C1", "D1")
site_order_reversed <- rev(site_order)  # High to low for panel A facets

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

# CHANGE 1: Reverse order for panel A - highest elevation at top
combined_data$site <- factor(combined_data$site, levels = site_order_reversed)
elev_order_reversed <- site_coords$elev_label[match(site_order_reversed, site_coords$name)]
combined_data$site_label <- factor(combined_data$site_label, levels = elev_order_reversed)

# Create color mapping for elevation labels
elev_colors <- setNames(site_colors, site_name_to_elev[names(site_colors)])

# Filter points and curves
point_data <- combined_data %>% filter(type == "point")
curve_data_plot <- combined_data %>% filter(type == "curve")

# CHANGE 2: Create legend data for panel A (horizontal layout in 1740m/MB facet - truly empty)
# Each item: point then text next to it
legend_data_a <- data.frame(
  site = factor(names(site_elevations), levels = site_order_reversed),
  elev = site_elevations,
  y = 5.5,  # Moved down to create space from title
  x = c(17, 22, 27, 32, 37),  # Spaced out to allow room for text
  spp = "MB",
  site_label = factor("1740m", levels = elev_order_reversed)  # 1740m is empty for MB
)

legend_box_a <- data.frame(
  xmin = 15.5,
  xmax = 44.5,  # Wide box
  ymin = 3.5,
  ymax = 10,
  spp = "MB",
  site_label = factor("1740m", levels = elev_order_reversed)
)

legend_title_a <- data.frame(
  x = 30,  # Centered
  y = 8.5,
  label = "Population",
  spp = "MB",
  site_label = factor("1740m", levels = elev_order_reversed)
)

# CHANGE 5: Use consistent y-axis label
y_axis_label <- "Feces mass / hopper mass / hour (mg/g/hr)"

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
  # Add legend box (single facet, horizontal layout)
  geom_rect(data = legend_box_a,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", color = "black", linewidth = 0.5,
            inherit.aes = FALSE) +
  # Legend title (centered above)
  geom_text(data = legend_title_a,
            aes(x = x, y = y, label = label),
            fontface = "bold", size = 3, hjust = 0.5,
            inherit.aes = FALSE) +
  # Legend points (horizontal)
  geom_point(data = legend_data_a, 
             aes(x = x, y = y, color = site),
             size = 2, inherit.aes = FALSE) +
  # Legend text (to the right of each point)
  geom_text(data = legend_data_a,
            aes(x = x + 0.8, y = y, label = elev),
            hjust = 0, size = 2.2, inherit.aes = FALSE) +
  facet_grid(site_label ~ spp,
             labeller = labeller(
               spp = as_labeller(c("MB" = "bolditalic('M. boulderensis')", 
                                   "MS" = "bolditalic('M. sanguinipes')"), 
                                 label_parsed))) +
  labs(
    x = "Temperature (°C)",
    y = y_axis_label
  ) +
  scale_color_manual(
    values = c(elev_colors, site_colors),  # Include both mappings
    name = "Population"
  ) +
  scale_x_continuous(limits = c(15, 45), breaks = seq(15, 45, by = 5)) +
  scale_y_continuous(limits = c(0, 13), breaks = seq(0, 12, by = 2), 
                     labels = c("0", "", "4", "", "8", "", "12")) +
  theme_minimal() +
  theme(
    strip.background.x = element_rect(fill = "#D2B48C", color = NA),
    strip.background.y = element_rect(fill = "#D2B48C", color = NA),
    strip.text.x = element_text(size = 12, face = "bold", color = "black"),
    strip.text.y = element_text(size = 9, face = "bold", color = "black"),  # CHANGE 3: smaller elevation labels
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(5, 5, 5, 25, "pt")  # CHANGE 4: add left margin for balance
  )

# Convert to grob for adding labels
g_a <- ggplotGrob(fig_a)

# Add "Species" label at the top
species_label <- textGrob("Species", gp = gpar(fontface = "bold", col = "#8B4513"))
g_a <- gtable_add_rows(g_a, heights = unit(0.6, "cm"), pos = 0)
panel_cols <- which(grepl("panel", g_a$layout$name))
g_a <- gtable_add_grob(g_a, species_label, t = 1, 
                       l = min(g_a$layout$l[panel_cols]),
                       r = max(g_a$layout$r[panel_cols]))

# Add "Population" label on the right
site_label_grob <- textGrob("Population", rot = 270, gp = gpar(fontface = "bold", col = "#8B4513"))
g_a <- gtable_add_cols(g_a, widths = unit(0.6, "cm"), pos = -1)
panel_rows <- which(grepl("panel", g_a$layout$name))
g_a <- gtable_add_grob(g_a, site_label_grob, 
                       t = min(g_a$layout$t[panel_rows]),
                       b = max(g_a$layout$b[panel_rows]),
                       l = ncol(g_a), r = ncol(g_a))

# ===== FIGURE B: PARAMETER ESTIMATES WITH CREDIBLE INTERVALS =====

# Extract parameters with credible intervals from pops
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

# CHANGE 3: Flip the order for Tmin, Tmax, Topt (but NOT Pmax/Ropt)
# Now high elevation (D1) should be at top (lowest site_idx)
# CHANGE 2: Increased spacing between CrIs
modified_param_summary <- all_param_summary %>%
  group_by(species) %>%
  mutate(
    sites_present = list(unique(site[!is.na(x)])),
    l_sites_pres = length(unique(site[!is.na(x)])),
    # Original index (low elevation first)
    site_idx_orig = match(site, unique(site[!is.na(x)])),
    # Reversed index for Tmin, Tmax, Topt (high elevation first)
    site_idx_rev = l_sites_pres - site_idx_orig + 1,
    species_max_rate = ifelse(species == "MB", max_rate_mb, max_rate_ms)
  ) %>%
  ungroup() %>%
  mutate(
    # Use reversed index for Tmin, Tmax, Topt; original for Ropt
    # Increased spacing (0.75 to better match Pmax spacing)
    y_display = case_when(
      parameter == "Tmin" ~ -0.5 - (site_idx_rev * 0.75),
      parameter == "Topt" ~ species_max_rate + 0.7 + 0.75*(l_sites_pres+1) - (site_idx_rev * 0.75),
      parameter == "Tmax" ~ -0.5 - (site_idx_rev * 0.75),
      parameter == "Ropt" ~ y,
      TRUE ~ y
    ),
    x_display = case_when(
      parameter == "Ropt" ~ 55 + (site_idx_orig * 2.5),  # Keep original order for Ropt
      TRUE ~ x
    )
  )

# Split out Ropt points
ropt_points <- modified_param_summary %>% 
  filter(parameter == "Ropt")

# Other parameter points
other_points <- modified_param_summary %>% 
  filter(parameter != "Ropt")

# CHANGE 4: Create labels for parameter groups
# We need one label per parameter per species panel
# Calculate positions for labels based on the CrI positions

# For Tmin labels (to the RIGHT of the CrIs) - CHANGE 5
# Add extra buffer from the maximum x_upper value
tmin_labels <- modified_param_summary %>%
  filter(parameter == "Tmin") %>%
  group_by(species) %>%
  summarize(
    x = max(x_upper, na.rm = TRUE) + 4,  # More buffer space
    y = mean(y_display, na.rm = TRUE),
    .groups = "drop"
  )

# For Tmax labels (to the LEFT of the CrIs) - CHANGE 5
# Add extra buffer from the minimum x_lower value
tmax_labels <- modified_param_summary %>%
  filter(parameter == "Tmax") %>%
  group_by(species) %>%
  summarize(
    x = min(x_lower, na.rm = TRUE) - 4,  # More buffer space
    y = mean(y_display, na.rm = TRUE),
    .groups = "drop"
  )

# For Topt labels (to the LEFT of the CrIs)
# Add extra buffer from the minimum x_lower value
topt_labels <- modified_param_summary %>%
  filter(parameter == "Topt") %>%
  group_by(species) %>%
  summarize(
    x = min(x_lower, na.rm = TRUE) - 4,  # More buffer space
    y = mean(y_display, na.rm = TRUE),
    .groups = "drop"
  )

# For Pmax labels (ABOVE the CrIs) - CHANGE 5
# Add extra buffer from the maximum y_upper value
pmax_labels <- modified_param_summary %>%
  filter(parameter == "Ropt") %>%
  group_by(species) %>%
  summarize(
    x = mean(x_display, na.rm = TRUE),
    y = max(y_upper, na.rm = TRUE) + 0.7,  # More buffer space
    .groups = "drop"
  )

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
    height = 0.5,  # Larger feet/brackets
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
  # CHANGE 4: Add parameter labels without subscripts
  geom_text(data = tmin_labels,
            aes(x = x, y = y),
            label = "Tmin",
            size = 3.5, color = "black") +
  geom_text(data = tmax_labels,
            aes(x = x, y = y),
            label = "Tmax",
            size = 3.5, color = "black") +
  geom_text(data = topt_labels,
            aes(x = x, y = y),
            label = "Topt",
            size = 3.5, color = "black") +
  geom_text(data = pmax_labels,
            aes(x = x, y = y),
            label = "Pmax",
            size = 3.5, color = "black") +
  facet_wrap(~ species, ncol = 2,
             labeller = labeller(
               species = as_labeller(c("MB" = "bolditalic('M. boulderensis')", 
                                       "MS" = "bolditalic('M. sanguinipes')"), 
                                     label_parsed))) +
  labs(
    x = "Temperature (°C)",
    y = y_axis_label,  # CHANGE 5: Consistent label
    color = "Population"
  ) +
  scale_color_manual(
    values = site_colors,
    labels = site_elevations,
    breaks = names(site_elevations)
  ) +
  scale_x_continuous(limits = c(5, 70), breaks = seq(10, 70, by = 10)) +
  scale_y_continuous(limits = c(-4, 8.5), breaks = seq(0, 8, by = 2)) +  # Expanded y limits
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "#D2B48C", color = NA),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.position = "none",
    axis.title = element_text(size = 10),  # Smaller axis title
    axis.title.y = element_text(margin = margin(r = 2)),  # Reduce y-axis title margin
    axis.text = element_text(size = 8),  # Smaller axis text
    plot.margin = margin(5, 5, 5, 5, "pt")
  )

# Add custom legend in top left of MS panel only - with better spacing
legend_data <- data.frame(
  site = factor(names(site_elevations), levels = names(site_elevations)),
  elev = site_elevations,
  y = seq(7, 3.5, length.out = 5),  # Moved down to create space from title
  x = 9.75,
  species = "MS"
)

legend_box <- data.frame(
  xmin = 6,
  xmax = 22,
  ymin = 2.8,
  ymax = 8.3,
  species = "MS"
)

legend_title <- data.frame(
  x = 7,
  y = 8,
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
            fontface = "bold", size = 3.5, hjust = 0,
            inherit.aes = FALSE) +
  geom_point(data = legend_data, 
             aes(x = x, y = y, color = site),
             size = 2.5, inherit.aes = FALSE) +
  geom_text(data = legend_data,
            aes(x = x + 1, y = y, label = elev),
            hjust = 0, size = 3, inherit.aes = FALSE)

# Convert to grob for adding "Species" label
g_b <- ggplotGrob(fig_b)

# Add "Species" label at the top
species_label_b <- textGrob("Species", gp = gpar(fontface = "bold", col = "#8B4513"))
g_b <- gtable_add_rows(g_b, heights = unit(0.6, "cm"), pos = 0)
panel_cols_b <- which(grepl("panel", g_b$layout$name))
g_b <- gtable_add_grob(g_b, species_label_b, t = 1, 
                       l = min(g_b$layout$l[panel_cols_b]),
                       r = max(g_b$layout$r[panel_cols_b]))

# ===== COMBINE FIGURES A AND B =====

# Convert grobs to ggplot objects for patchwork
fig_a_gg <- wrap_elements(g_a)
fig_b_gg <- wrap_elements(g_b)

# Combine with patchwork
combined_figure <- fig_a_gg / fig_b_gg +
  plot_layout(heights = c(4, 3)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Display the combined figure
print(combined_figure)

# CHANGE 6: Save with appropriate dimensions
# Based on image dimensions ~975x1010 pixels at screen resolution
# For publication quality at 300 dpi, using similar aspect ratio
ggsave("figure1_tpc.png", 
       plot = combined_figure, 
       width = 8,      # inches
       height = 8.3,   # inches (slightly taller than wide, matching original)
       dpi = 300,
       bg = "white")

# Also save as PDF for publication
ggsave("figure1_tpc.pdf", 
       plot = combined_figure, 
       width = 8, 
       height = 8.3,
       bg = "white")