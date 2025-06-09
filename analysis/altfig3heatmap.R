library(tidyverse)
library(ggplot2)
library(viridis)  # For better color scales

# Function to create sample data (you'll replace this with your actual data)
create_sample_real_data <- function() {
  # Same function as before
  # Define site constraints for each species
  mb_sites <- c("A1", "B1", "C1", "D1")
  ms_sites <- c("Eldo", "A1", "B1")
  all_sites <- c("Eldo", "A1", "B1", "C1", "D1")
  
  # Create years vectors
  recent_years <- c(2005:2024)  # 20 years of "recent" data
  
  # Set a consistent seed for reproducibility
  set.seed(42)
  
  # Base values that decrease with elevation
  base_values <- data.frame(
    site = c("Eldo", "A1", "B1", "C1", "D1"),
    base_value = c(80, 60, 50, 40, 30)
  )
  
  # Create MB combinations - transplant only
  mb_transplant <- expand.grid(
    species = "MB",
    year = recent_years,
    site_orig = mb_sites,
    site_clim = all_sites,
    sex = c("F"),  # Females only
    stringsAsFactors = FALSE
  ) %>%
    mutate(year_period = "recent")
  
  # Create MS combinations - transplant only
  ms_transplant <- expand.grid(
    species = "MS",
    year = recent_years,
    site_orig = ms_sites,
    site_clim = all_sites,
    sex = c("F"),  # Females only
    stringsAsFactors = FALSE
  ) %>%
    mutate(year_period = "recent")
  
  # Combine transplant data
  all_data <- bind_rows(mb_transplant, ms_transplant)
  
  # Join with base values
  all_data <- all_data %>%
    left_join(base_values, by = c("site_orig" = "site")) %>%
    rename(orig_base = base_value) %>%
    left_join(base_values, by = c("site_clim" = "site")) %>%
    rename(clim_base = base_value)
  
  # Generate consistent energy values
  all_data <- all_data %>%
    mutate(
      # Species factor
      species_factor = case_when(
        species == "MB" ~ 0.9,
        species == "MS" ~ 1.1
      ),
      
      # Climate effect (organisms perform better in their home climate)
      climate_match = abs(as.numeric(factor(site_orig, levels = all_sites)) - 
                            as.numeric(factor(site_clim, levels = all_sites))),
      climate_penalty = climate_match * 5,
      
      # Add home climate bonus
      home_bonus = ifelse(site_orig == site_clim, 15, 0),
      
      # Calculate energy values consistently
      partial_shade_low_veg = (orig_base * species_factor - climate_penalty + home_bonus) + 
        rnorm(n(), mean = 0, sd = 3),
      
      # Ensure no negative values
      partial_shade_low_veg = pmax(partial_shade_low_veg, 10)
    ) %>%
    select(-orig_base, -clim_base, -species_factor, -climate_match, 
           -climate_penalty, -home_bonus)
  
  return(all_data)
}

# Create the sample data with a structure matching your real data
# You would replace this with how you load your actual data
#thing <- create_sample_real_data()

# Function to create a heatmap representation of the reciprocal transplant data
create_reciprocal_transplant_heatmap <- function(data) {
  # Define site order
  site_order <- c("Eldo", "A1", "B1", "C1", "D1")
  
  # Calculate means for each site_orig x site_clim combination
  # Filter for females only
  mb_data <- data %>%
    filter(species == "MB", sex == "F") %>%
    group_by(site_orig, site_clim) %>%
    summarize(
      mean_energy = mean(partial_shade_low_veg, na.rm = TRUE),
      sd_energy = sd(partial_shade_low_veg, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      # Factor both variables for proper ordering
      site_orig = factor(site_orig, levels = site_order),
      site_clim = factor(site_clim, levels = site_order),
      # Add a flag for home climate
      home_climate = ifelse(site_orig == site_clim, "Home", "Transplant")
    )
  
  ms_data <- data %>%
    filter(species == "MS", sex == "F") %>%
    group_by(site_orig, site_clim) %>%
    summarize(
      mean_energy = mean(partial_shade_low_veg, na.rm = TRUE),
      sd_energy = sd(partial_shade_low_veg, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      # Factor both variables for proper ordering
      site_orig = factor(site_orig, levels = site_order),
      site_clim = factor(site_clim, levels = site_order),
      # Add a flag for home climate
      home_climate = ifelse(site_orig == site_clim, "Home", "Transplant")
    )
  
  # Find combined range for consistent color scale
  min_energy <- min(c(mb_data$mean_energy, ms_data$mean_energy), na.rm = TRUE)
  max_energy <- max(c(mb_data$mean_energy, ms_data$mean_energy), na.rm = TRUE)
  
  # Create the heatmap for MB
  mb_heatmap <- ggplot(mb_data, aes(x = site_clim, y = site_orig, fill = mean_energy)) +
    geom_tile(color = "white", size = 0.5) +
    # Add values as text
    geom_text(aes(label = sprintf("%.1f", mean_energy), 
                  color = ifelse(mean_energy > (min_energy + max_energy)/2, "white", "black")), 
              size = 3.5) +
    # Use viridis color scale for better perception
    scale_fill_viridis_c(
      name = "Energy (kJ)",
      limits = c(min_energy, max_energy),
      option = "D"  # Use the default "D" option
    ) +
    # Simple black/white text that's visible on all backgrounds
    scale_color_identity() +
    # Add special border for home climate cells
    geom_tile(data = subset(mb_data, home_climate == "Home"),
              aes(x = site_clim, y = site_orig),
              fill = NA, color = "black", size = 1.5) +
    # Better labels
    labs(
      title = "MB - Mean Energy in Virtual Reciprocal Transplants (Females)",
      subtitle = "MB occurs only at A1, B1, C1, and D1 sites",
      x = "Climate Site",
      y = "Population Origin"
    ) +
    # Adjust theme elements
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold")
    )
  
  # Create the heatmap for MS
  ms_heatmap <- ggplot(ms_data, aes(x = site_clim, y = site_orig, fill = mean_energy)) +
    geom_tile(color = "white", size = 0.5) +
    # Add values as text
    geom_text(aes(label = sprintf("%.1f", mean_energy), 
                  color = ifelse(mean_energy > (min_energy + max_energy)/2, "white", "black")), 
              size = 3.5) +
    # Use viridis color scale for better perception
    scale_fill_viridis_c(
      name = "Energy (kJ)",
      limits = c(min_energy, max_energy),
      option = "D"  # Use the default "D" option
    ) +
    # Simple black/white text that's visible on all backgrounds
    scale_color_identity() +
    # Add special border for home climate cells
    geom_tile(data = subset(ms_data, home_climate == "Home"),
              aes(x = site_clim, y = site_orig),
              fill = NA, color = "black", size = 1.5) +
    # Better labels
    labs(
      title = "MS - Mean Energy in Virtual Reciprocal Transplants (Females)",
      subtitle = "MS occurs only at Eldo, A1, and B1 sites",
      x = "Climate Site",
      y = "Population Origin"
    ) +
    # Adjust theme elements
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold")
    )
  
  # Return both heatmaps
  return(list(mb_heatmap = mb_heatmap, ms_heatmap = ms_heatmap))
}

# Create the heatmaps
heatmaps <- create_reciprocal_transplant_heatmap(thing)

# Display the heatmaps
print(heatmaps$mb_heatmap)
print(heatmaps$ms_heatmap)

# If you want to save the heatmaps, uncomment these lines:
# ggsave("MB_reciprocal_transplant_heatmap.png", heatmaps$mb_heatmap, width = 8, height = 6)
# ggsave("MS_reciprocal_transplant_heatmap.png", heatmaps$ms_heatmap, width = 8, height = 6)

# BONUS: Combined heatmap with facets for both species
create_combined_heatmap <- function(data) {
  # Define site order
  site_order <- c("Eldo", "A1", "B1", "C1", "D1")
  
  # Calculate means for each species, site_orig x site_clim combination
  # Filter for females only
  combined_data <- data %>%
    filter(sex == "F") %>%
    group_by(species, site_orig, site_clim) %>%
    summarize(
      mean_energy = mean(partial_shade_low_veg, na.rm = TRUE),
      sd_energy = sd(partial_shade_low_veg, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      # Factor both variables for proper ordering
      site_orig = factor(site_orig, levels = site_order),
      site_clim = factor(site_clim, levels = site_order),
      # Add a flag for home climate
      home_climate = ifelse(site_orig == site_clim, "Home", "Transplant")
    )
  
  # Create combined heatmap
  combined_heatmap <- ggplot(combined_data, aes(x = site_clim, y = site_orig, fill = mean_energy)) +
    facet_wrap(~ species) +
    geom_tile(color = "white", size = 0.5) +
    # Add values as text
    geom_text(aes(label = sprintf("%.1f", mean_energy), 
                  color = ifelse(mean_energy > mean(range(combined_data$mean_energy, na.rm = TRUE)), 
                                 "white", "black")), 
              size = 3.5) +
    # Use viridis color scale for better perception
    scale_fill_viridis_c(name = "Energy (kJ)", option = "D") +
    # Simple black/white text that's visible on all backgrounds
    scale_color_identity() +
    # Add special border for home climate cells
    geom_tile(data = subset(combined_data, home_climate == "Home"),
              aes(x = site_clim, y = site_orig),
              fill = NA, color = "black", size = 1.5) +
    # Better labels
    labs(
      title = "Mean Energy in Virtual Reciprocal Transplants (Females)",
      subtitle = "MB: A1-D1 sites, MS: Eldo-B1 sites",
      x = "Climate Site",
      y = "Population Origin"
    ) +
    # Adjust theme elements
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  return(combined_heatmap)
}

# Create and display the combined heatmap
combined_heatmap <- create_combined_heatmap(thing)
print(combined_heatmap)
# ggsave("combined_reciprocal_transplant_heatmap.png", combined_heatmap, width = 10, height = 6)