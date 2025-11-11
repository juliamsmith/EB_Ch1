library(tidyverse)
library(patchwork)
library(TrenchR)
library(grid)
library(gtable)

#directory is output/results/

# LRF function for thermal performance curves
lrf_1991_exact <- function(temp, rmax, topt, tmin, tmax) {
  # Return 0 for temps outside valid range
  if (temp <= tmin || temp >= tmax) {
    return(0)
  }
  
  # Implement the exact equation
  tryCatch({
    numerator <- (temp - tmax) * (temp - tmin)^2
    denominator <- (topt - tmin) * ((topt - tmin) * (temp - topt) - (topt - tmax) * (topt + tmin - 2 * temp))
    
    rate <- rmax * numerator / denominator
    
    # Handle NaN or negative values
    if (is.na(rate) || rate < 0) {
      return(0)
    }
    
    return(rate)
  }, error = function(e) {
    warning(paste("Error in LRF calculation:", e$message))
    return(0)
  })
}

# Function to calculate TPC gains with energy conversion
get_tpc_gains <- function(temp, rmax, topt, tmin, tmax, assim.rate = 0.4) {
  # Calculate LRF-based dry feces production
  dry.fec.mg.hr <- lrf_1991_exact(temp, rmax, topt, tmin, tmax)
  
  # Ensure non-negative values
  dry.fec.mg.hr[dry.fec.mg.hr < 0] <- 0
  
  # Convert to wet grass-adjusted (WGA) production
  # Using Harrison & Fewell conversion
  dry.wga.mg.hr <- dry.fec.mg.hr * 12.8/14.8 * assim.rate / (1 - assim.rate)
  
  # Convert to kilocalories per hour
  kcal.hr <- dry.wga.mg.hr * 14.8 / 1000
  
  # Convert to kilojoules per hour
  kJ.hr <- kcal.hr * 4.184
  
  return(kJ.hr)
}

# Site coordinates
site_coords <- data.frame(
  name = c("Eldo", "A1", "B1", "C1", "D1"),
  lat = c(39.944, 40.015, 40.023, 40.036, 40.059),
  lon = c(-105.262, -105.377, -105.430, -105.547, -105.617),
  elev = c(1740, 2195, 2591, 3048, 3515)
)

# Valid species-site combinations
valid_combinations <- list(
  MS = c("Eldo", "A1", "B1"),
  MB = c("A1", "B1", "C1", "D1")
)

# Site order for plotting and elevation mapping
site_order <- c("Eldo", "A1", "B1", "C1", "D1")
site_elev <- c("Eldo" = "1740m", "A1" = "2195m", "B1" = "2591m", 
               "C1" = "3048m", "D1" = "3515m")

# Function to load RDS files and filter to specific microclimate condition
load_and_process_direct <- function(base_dir = ".", pattern = "^eb_results_", 
                                    use_condition = "partial_shade_low_veg") {
  # Find all RDS files matching the pattern
  all_files <- list.files(path = base_dir, 
                          pattern = paste0(pattern, ".*\\.rds$"), 
                          recursive = TRUE, 
                          full.names = TRUE)
  
  if (length(all_files) == 0) {
    stop("No result files found! Check the base directory and pattern.")
  }
  
  message("Found ", length(all_files), " files to process")
  
  # Read and combine all files
  temp_data <- lapply(all_files, function(file) {
    tryCatch({
      data <- readRDS(file)
      
      # Filter to the specific microclimate condition
      if (use_condition == "partial_shade_low_veg") {
        # Filter to shade = 0.3, height = 0.03 (partial_shade_low_veg)
        filtered_data <- data %>%
          filter(shade == 0.3 & height == 0.03)
      } else if (use_condition == "PBT") {
        # Use thermoregulated body temperature (PBT condition)
        filtered_data <- data %>%
          filter(condition == "PBT")
      } else if (use_condition == "sun_low_veg") {
        # Filter to shade = 0.0, height = 0.03 (full sun, low vegetation)
        filtered_data <- data %>%
          filter(shade == 0.0 & height == 0.03)
      } else if (use_condition == "shade_low_veg") {
        # Filter to shade = 0.8, height = 0.03 (heavy shade, low vegetation)
        filtered_data <- data %>%
          filter(shade == 0.8 & height == 0.03)
      } else {
        warning("Unknown condition: ", use_condition, ". Using all data.")
        filtered_data <- data
      }
      
      if (nrow(filtered_data) == 0) {
        warning("No data found for condition ", use_condition, " in file ", basename(file))
        return(NULL)
      }
      
      message("Loaded ", nrow(filtered_data), " rows from ", basename(file), 
              " for condition: ", use_condition)
      
      return(filtered_data)
    }, error = function(e) {
      warning(paste("Error reading file", file, ":", e$message))
      return(NULL)
    })
  })
  
  # Combine all data frames, filtering out NULL results
  combined_data <- bind_rows(compact(temp_data))
  
  if (nrow(combined_data) == 0) {
    stop("No data loaded from result files for condition: ", use_condition)
  }
  
  message("Total rows loaded: ", nrow(combined_data))
  return(combined_data)
}

# Function to calculate metabolic rate
get_mrs <- function(tb, mass, elev, b0, b1, b2, b3, k=8.62*10^-5) {
  rmr <- exp(b0 + b1*log(mass) + b2*(1/(k*(tb+273.15))) + b3*elev)
  return(rmr)
}

# Function to convert metabolic rates to energy losses
get_mr_losses <- function(mrs) {
  o2_ml <- 1/.7 * mrs
  lipid_mg <- o2_ml/2
  loss_kJ <- lipid_mg*39/1000
  return(loss_kJ)
}

# Function to generate energy curves
generate_energy_curves <- function(species, site, pop_data, year_period) {
  # Get parameters for this site-species combination
  p <- pop_data %>% 
    filter(site == !!site, spp == species) %>% 
    slice(1)  # Take first row (assumes consistent params for site-species)
  
  # Skip if no matching data found
  if (nrow(p) == 0) {
    return(NULL)
  }
  
  # Define temperature range
  temp_range <- seq(p$Tmin, p$Topt + p$Above, by = 0.1)
  
  # Calculate TPC (gains) using LRF - already mass standardized
  tpc_values <- sapply(temp_range, function(temp) {
    get_tpc_gains(temp, p$Ropt, p$Topt, p$Tmin, p$Topt + p$Above)
  })
  
  # Calculate metabolic rate (losses)
  mr_values <- sapply(temp_range, function(temp) {
    get_mrs(temp, p$mass, p$elev, p$rmr_b0, p$rmr_b1, p$rmr_b2, p$rmr_b3)
  })
  
  # Only standardize metabolic losses by mass since TPCs are already standardized
  mr_loss_values <- sapply(mr_values, get_mr_losses)
  mr_loss_values_per_g <- mr_loss_values / p$mass
  
  # Calculate net energy gain (per gram)
  net_values <- tpc_values - mr_loss_values_per_g
  
  # Create data frame
  data.frame(
    temperature = temp_range,
    Gains = tpc_values,
    Losses = mr_loss_values_per_g,
    Net = net_values,
    site = site,
    spp = species,
    year_period = year_period
  )
}

# Function to filter daylight hours using TrenchR with 1-hour buffer
filter_daylight_hours <- function(data, lat, lon, buffer_hours = 1) {
  # First, let's check the hour format
  message("  Hour range in data: ", min(data$hour, na.rm = TRUE), " to ", 
          max(data$hour, na.rm = TRUE))
  
  # Calculate zenith angle for each observation using vectorized operations
  # Note: hour should be 0-24 format. If your hour is in a different format, adjust here
  data <- data %>%
    mutate(
      zenith = zenith_angle(doy = doy, lat = lat, lon = lon, hour = hour)
    )
  
  # For each unique day, find sunrise and sunset times
  # Create a lookup table of sunrise/sunset by day of year
  unique_doys <- unique(data$doy)
  
  sunrise_sunset <- map_dfr(unique_doys, function(d) {
    # Test hours from 0-24 to find when zenith crosses 90 degrees
    test_hours <- seq(0, 24, by = 0.1)
    zeniths <- zenith_angle(doy = d, lat = lat, lon = lon, hour = test_hours)
    
    # Find sunrise (first time zenith goes below 90)
    sunrise_idx <- which(zeniths < 90)[1]
    sunrise_hour <- if(!is.na(sunrise_idx)) test_hours[sunrise_idx] else NA
    
    # Find sunset (last time zenith is below 90)
    sunset_idx <- tail(which(zeniths < 90), 1)
    sunset_hour <- if(length(sunset_idx) > 0) test_hours[sunset_idx] else NA
    
    tibble(
      doy = d,
      sunrise = sunrise_hour,
      sunset = sunset_hour
    )
  })
  
  # Join sunrise/sunset times to data
  data <- data %>%
    left_join(sunrise_sunset, by = "doy")
  
  # Filter to daylight hours with buffer
  daylight_data <- data %>%
    filter(
      !is.na(sunrise) & !is.na(sunset),
      hour >= (sunrise + buffer_hours),
      hour <= (sunset - buffer_hours)
    )
  
  message("  Sunrise/sunset range: ", 
          round(min(sunrise_sunset$sunrise, na.rm = TRUE), 1), "-",
          round(max(sunrise_sunset$sunrise, na.rm = TRUE), 1), " to ",
          round(min(sunrise_sunset$sunset, na.rm = TRUE), 1), "-",
          round(max(sunrise_sunset$sunset, na.rm = TRUE), 1))
  message("  Filtered from ", nrow(data), " to ", nrow(daylight_data), 
          " daylight observations (with ", buffer_hours, " hour buffer)")
  
  return(daylight_data)
}

# Main function to create the plot
process_energy_plot <- function(base_dir = ".", pops_data, use_condition = "partial_shade_low_veg") {
  # Load temperature data and filter to specific condition
  message("Loading temperature data from RDS files for condition: ", use_condition)
  all_temps <- load_and_process_direct(base_dir, use_condition = use_condition)
  
  message("Loaded temperature data. Number of rows:", nrow(all_temps))
  message("Column names: ", paste(colnames(all_temps), collapse = ", "))
  
  # Check what condition was actually loaded
  if ("condition" %in% colnames(all_temps)) {
    message("Conditions in data: ", paste(unique(all_temps$condition), collapse = ", "))
  }
  if ("shade" %in% colnames(all_temps) && "height" %in% colnames(all_temps)) {
    message("Shade/height combinations: ")
    print(table(all_temps$shade, all_temps$height, useNA = "ifany"))
  }
  
  # Generate energy curves for valid site-species combinations
  message("Generating energy curves...")
  energy_curves <- map_dfr(c("MB", "MS"), function(species) {
    # Filter valid sites for this species
    valid_sites <- valid_combinations[[species]]
    
    # Process each valid site
    map_dfr(valid_sites, function(site) {
      # Skip if site-species combination doesn't exist in pops_data
      if (nrow(filter(pops_data, spp == species, site == !!site)) == 0) {
        message(paste("Skipping", species, site, "- not found in population data"))
        return(NULL)
      }
      
      message(paste("Processing energy curves for", species, "at site", site))
      bind_rows(
        generate_energy_curves(
          species = species, 
          site = site, 
          pop_data = pops_data, 
          year_period = "historical"
        ),
        generate_energy_curves(
          species = species, 
          site = site, 
          pop_data = pops_data, 
          year_period = "recent"
        )
      )
    })
  })
  
  # Print summary to verify we have unique curves
  message("Energy curves summary:")
  energy_summary <- energy_curves %>%
    group_by(spp, site) %>%
    summarize(
      min_temp = min(temperature),
      max_temp = max(temperature),
      max_gain = max(Gains),
      max_loss = max(Losses),
      .groups = "drop"
    )
  print(energy_summary)
  
  message("Filtering temperature data...")
  # Filter temperature data to match valid combinations
  filtered_temps <- all_temps %>%
    filter(
      # Ensure we have all required columns
      !is.na(Tb),
      !is.na(species),
      !is.na(site_clim),
      !is.na(year_period),
      
      # Filter to temperatures <= 65°C
      Tb <= 65,
      
      # Filter to valid species-site combinations
      (species == "MS" & site_clim %in% valid_combinations[["MS"]]) |
        (species == "MB" & site_clim %in% valid_combinations[["MB"]])
    )
  
  message("Before daylight filtering: ", nrow(filtered_temps), " rows")
  
  # Filter to daylight hours using site-specific lat/lon
  message("Applying daylight filtering by site...")
  filtered_temps <- filtered_temps %>%
    group_split(site_clim) %>%
    map_dfr(~ {
      # Get the site name for this group
      current_site <- as.character(first(.x$site_clim))
      
      # Get lat/lon for this site from site_coords
      site_match <- site_coords %>% filter(name == current_site)
      
      if (nrow(site_match) > 0) {
        lat_val <- site_match$lat[1]
        lon_val <- site_match$lon[1]
        
        message("Filtering daylight for site ", current_site, 
                " (lat: ", lat_val, ", lon: ", lon_val, ")")
        filter_daylight_hours(.x, lat = lat_val, lon = lon_val)
      } else {
        message("No coordinates found for site ", current_site, ", keeping all data")
        .x
      }
    })
  
  message("After daylight filtering: ", nrow(filtered_temps), " rows")
  
  # Continue with remaining processing
  filtered_temps <- filtered_temps %>%
    mutate(
      site_clim = factor(site_clim, levels = site_order),
      species = factor(species, levels = c("MB", "MS")),
      # Add elevation labels
      site_label = site_elev[as.character(site_clim)],
      site_label = factor(site_label, levels = site_elev)
    )
  
  message("Filtered temperature data. Number of rows:", nrow(filtered_temps))
  message("Temperature range: ", round(min(filtered_temps$Tb), 1), " to ", round(max(filtered_temps$Tb), 1), "°C")
  
  # Print summary of what we're plotting
  temp_summary <- filtered_temps %>%
    group_by(species, site_clim, year_period) %>%
    summarize(
      n_obs = n(),
      mean_temp = mean(Tb),
      min_temp = min(Tb),
      max_temp = max(Tb),
      .groups = "drop"
    )
  message("Temperature data summary by group:")
  print(temp_summary)
  
  # Filter energy curves by temperature
  energy_curves_long <- energy_curves %>%
    filter(temperature <= 65) %>%
    pivot_longer(cols = c(Gains, Losses, Net), 
                 names_to = "type", 
                 values_to = "value") %>%
    mutate(
      type = factor(type, levels = c("Gains", "Losses", "Net")),
      site = factor(site, levels = site_order),
      spp = factor(spp, levels = c("MB", "MS"))
    )
  
  message("Calculating scaling factor...")
  # Calculate scaling factor
  max_density <- max(sapply(split(filtered_temps$Tb, 
                                  paste(filtered_temps$site_clim, filtered_temps$species)), 
                            function(x) {
                              if(length(x) > 1) return(max(stats::density(x)$y))
                              else return(0)
                            }))
  max_rate <- max(abs(energy_curves_long$value))
  scaling_factor <- max_density / max_rate
  
  message("Max density: ", round(max_density, 4))
  message("Max rate: ", round(max_rate, 4))
  message("Scaling factor: ", round(scaling_factor, 4))
  
  # Modify energy_curves_long to match faceting structure with elevation labels
  energy_curves_long <- energy_curves_long %>%
    rename(site_clim = site) %>%
    mutate(
      site_clim = factor(site_clim, levels = site_order),
      species = factor(spp, levels = c("MB", "MS")),
      site_label = site_elev[as.character(site_clim)],
      site_label = factor(site_label, levels = site_elev)
    )
  
  # Verify we have the correct combinations
  message("Energy curves by site and species:")
  print(table(energy_curves_long$site_clim, energy_curves_long$species))
  
  # Create plot with rotated layout (sites as rows)
  p <- ggplot() +
    # Temperature distributions - updated colors
    geom_density(data = filtered_temps, 
                 aes(x = Tb, fill = year_period),
                 alpha = 0.5) +
    # Energy curves - now green to match axis
    geom_line(data = energy_curves_long,
              aes(x = temperature, 
                  y = value * scaling_factor, 
                  linetype = type,
                  group = interaction(type, site_clim, species)),
              color = "darkgreen",
              linewidth = 1) +
    scale_fill_manual(
      values = c("historical" = "#8DA0CB", "recent" = "#FC8D62"),
      labels = c("historical" = "Historical (1950-1959)", "recent" = "Contemporary (2015-2024)"),
      name = "Period",
      drop = FALSE
    ) +
    scale_linetype_manual(values = c("Gains" = "solid", 
                                     "Losses" = "dashed", 
                                     "Net" = "dotted")) +
    scale_y_continuous(
      name = "Temperature Density",
      sec.axis = sec_axis(~./scaling_factor, name = "Energy Rate (kJ/hr/g)")
    ) +
    scale_x_continuous(limits = c(-5, 65)) +
    facet_grid(site_label ~ species, scales = "free_y") +
    labs(x = "Body Temperature (°C)") +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "#D2B48C", color = NA),
      strip.text = element_text(face = "bold", color = "black"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines"),
      axis.title.y.right = element_text(color = "darkgreen"),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.7)))
  
  # Convert to grob for adding labels
  g <- ggplotGrob(p)
  
  # Add "Species" label at the top
  species_label <- textGrob("Species", gp = gpar(fontsize = 12, fontface = "bold", col = "#8B4513"))
  g <- gtable_add_rows(g, heights = unit(0.6, "cm"), pos = 0)
  panel_cols <- which(grepl("panel", g$layout$name))  # Match all panels
  g <- gtable_add_grob(g, species_label, t = 1, 
                       l = min(g$layout$l[panel_cols]),  # Leftmost panel column
                       r = max(g$layout$r[panel_cols]))  # Rightmost panel column
  
  # Add "Site" label on the right - centered vertically across all panels
  site_label_grob <- textGrob("Site", rot = 270, gp = gpar(fontsize = 12, fontface = "bold", col = "#8B4513"))
  g <- gtable_add_cols(g, widths = unit(0.6, "cm"), pos = -1)
  panel_rows <- which(grepl("panel", g$layout$name))
  g <- gtable_add_grob(g, site_label_grob, 
                       t = min(g$layout$t[panel_rows]),  # Top of first panel
                       b = max(g$layout$b[panel_rows]),  # Bottom of last panel
                       l = ncol(g), r = ncol(g))
  
  return(g)
}

# Example usage:
plot_partial_shade <- process_energy_plot(".", pops, use_condition = "partial_shade_low_veg")
grid.newpage()
grid.draw(plot_partial_shade)




#now doing some diagnostics:


library(tidyverse)

# Function to examine extreme temperatures
diagnose_extreme_temps <- function(data, temp_col = "Tb", 
                                   high_threshold = 65, 
                                   low_threshold = 0) {
  
  # Create datetime column if it doesn't exist
  if (!("datetime" %in% colnames(data))) {
    data <- data %>%
      mutate(datetime = as.POSIXct(paste(year, month, day, hour, sep = "-"),
                                   format = "%Y-%m-%d-%H"))
  }
  
  # Identify high temperature observations
  high_temps <- data %>%
    filter(!!sym(temp_col) > high_threshold) %>%
    arrange(desc(!!sym(temp_col))) %>%
    select(datetime, year, month, day, hour, doy, 
           species, site_orig, site_clim, year_period,
           Tb, Tair, Tsoil, rad, wind, shade, height,
           everything())
  
  # Identify low temperature observations
  low_temps <- data %>%
    filter(!!sym(temp_col) < low_threshold) %>%
    arrange(!!sym(temp_col)) %>%
    select(datetime, year, month, day, hour, doy,
           species, site_orig, site_clim, year_period,
           Tb, Tair, Tsoil, rad, wind, shade, height,
           everything())
  
  # Summary statistics
  message("\n=== EXTREME TEMPERATURE DIAGNOSTICS ===\n")
  
  message("HIGH TEMPERATURES (>", high_threshold, "°C):")
  message("  Total observations: ", nrow(high_temps))
  if (nrow(high_temps) > 0) {
    message("  Temperature range: ", round(min(high_temps$Tb), 1), " to ", 
            round(max(high_temps$Tb), 1), "°C")
    message("\n  By species:")
    print(table(high_temps$species))
    message("\n  By site:")
    print(table(high_temps$site_clim))
    message("\n  By year:")
    print(table(high_temps$year))
    message("\n  By month:")
    print(table(high_temps$month))
    message("\n  By hour:")
    print(table(high_temps$hour))
  }
  
  message("\n\nLOW TEMPERATURES (<", low_threshold, "°C):")
  message("  Total observations: ", nrow(low_temps))
  if (nrow(low_temps) > 0) {
    message("  Temperature range: ", round(min(low_temps$Tb), 1), " to ", 
            round(max(low_temps$Tb), 1), "°C")
    message("\n  By species:")
    print(table(low_temps$species))
    message("\n  By site:")
    print(table(low_temps$site_clim))
    message("\n  By year:")
    print(table(low_temps$year))
  }
  
  # Create summary plots
  if (nrow(high_temps) > 0) {
    p_high <- ggplot(high_temps, aes(x = hour, y = Tb)) +
      geom_point(aes(color = site_clim), alpha = 0.5) +
      facet_wrap(~species) +
      labs(title = paste("High Temperatures (>", high_threshold, "°C) by Hour"),
           x = "Hour of Day", y = "Body Temperature (°C)",
           color = "Site") +
      theme_minimal()
    print(p_high)
  }
  
  # Return list of data frames
  return(list(
    high_temps = high_temps,
    low_temps = low_temps,
    summary = list(
      n_high = nrow(high_temps),
      n_low = nrow(low_temps),
      high_range = if(nrow(high_temps) > 0) range(high_temps$Tb) else NA,
      low_range = if(nrow(low_temps) > 0) range(low_temps$Tb) else NA
    )
  ))
}

# Example usage:

all_temps <- load_and_process_direct(base_dir=".", use_condition = "partial_shade_low_veg")
# Load your data first, then:
extreme_temps <- diagnose_extreme_temps(all_temps, high_threshold = 65, low_threshold = 0)
# 
# View the high temperature observations:
# View(extreme_temps$high_temps)
# 
# View the low temperature observations:
# View(extreme_temps$low_temps)
#
# Export to CSV if needed:
# write_csv(extreme_temps$high_temps, "high_temps_diagnosis.csv")
# write_csv(extreme_temps$low_temps, "low_temps_diagnosis.csv")





#slightly older

library(tidyverse)
library(patchwork)

#directory is output/results/

# LRF function for thermal performance curves
lrf_1991_exact <- function(temp, rmax, topt, tmin, tmax) {
  # Return 0 for temps outside valid range
  if (temp <= tmin || temp >= tmax) {
    return(0)
  }
  
  # Implement the exact equation
  tryCatch({
    numerator <- (temp - tmax) * (temp - tmin)^2
    denominator <- (topt - tmin) * ((topt - tmin) * (temp - topt) - (topt - tmax) * (topt + tmin - 2 * temp))
    
    rate <- rmax * numerator / denominator
    
    # Handle NaN or negative values
    if (is.na(rate) || rate < 0) {
      return(0)
    }
    
    return(rate)
  }, error = function(e) {
    warning(paste("Error in LRF calculation:", e$message))
    return(0)
  })
}

# Function to calculate TPC gains with energy conversion
get_tpc_gains <- function(temp, rmax, topt, tmin, tmax, assim.rate = 0.4) {
  # Calculate LRF-based dry feces production
  dry.fec.mg.hr <- lrf_1991_exact(temp, rmax, topt, tmin, tmax)
  
  # Ensure non-negative values
  dry.fec.mg.hr[dry.fec.mg.hr < 0] <- 0
  
  # Convert to wet grass-adjusted (WGA) production
  # Using Harrison & Fewell conversion
  dry.wga.mg.hr <- dry.fec.mg.hr * 12.8/14.8 * assim.rate / (1 - assim.rate)
  
  # Convert to kilocalories per hour
  kcal.hr <- dry.wga.mg.hr * 14.8 / 1000
  
  # Convert to kilojoules per hour
  kJ.hr <- kcal.hr * 4.184
  
  return(kJ.hr)
}

# Valid species-site combinations
valid_combinations <- list(
  MS = c("Eldo", "A1", "B1"),
  MB = c("A1", "B1", "C1", "D1")
)

# Site order for plotting
site_order <- c("Eldo", "A1", "B1", "C1", "D1")

# Function to load RDS files and filter to specific microclimate condition
load_and_process_direct <- function(base_dir = ".", pattern = "^eb_results_", 
                                    use_condition = "partial_shade_low_veg") {
  # Find all RDS files matching the pattern
  all_files <- list.files(path = base_dir, 
                          pattern = paste0(pattern, ".*\\.rds$"), 
                          recursive = TRUE, 
                          full.names = TRUE)
  
  if (length(all_files) == 0) {
    stop("No result files found! Check the base directory and pattern.")
  }
  
  message("Found ", length(all_files), " files to process")
  
  # Read and combine all files
  temp_data <- lapply(all_files, function(file) {
    tryCatch({
      data <- readRDS(file)
      
      # Filter to the specific microclimate condition
      if (use_condition == "partial_shade_low_veg") {
        # Filter to shade = 0.3, height = 0.03 (partial_shade_low_veg)
        filtered_data <- data %>%
          filter(shade == 0.3 & height == 0.03)
      } else if (use_condition == "PBT") {
        # Use thermoregulated body temperature (PBT condition)
        filtered_data <- data %>%
          filter(condition == "PBT")
      } else if (use_condition == "sun_low_veg") {
        # Filter to shade = 0.0, height = 0.03 (full sun, low vegetation)
        filtered_data <- data %>%
          filter(shade == 0.0 & height == 0.03)
      } else if (use_condition == "shade_low_veg") {
        # Filter to shade = 0.8, height = 0.03 (heavy shade, low vegetation)
        filtered_data <- data %>%
          filter(shade == 0.8 & height == 0.03)
      } else {
        warning("Unknown condition: ", use_condition, ". Using all data.")
        filtered_data <- data
      }
      
      if (nrow(filtered_data) == 0) {
        warning("No data found for condition ", use_condition, " in file ", basename(file))
        return(NULL)
      }
      
      message("Loaded ", nrow(filtered_data), " rows from ", basename(file), 
              " for condition: ", use_condition)
      
      return(filtered_data)
    }, error = function(e) {
      warning(paste("Error reading file", file, ":", e$message))
      return(NULL)
    })
  })
  
  # Combine all data frames, filtering out NULL results
  combined_data <- bind_rows(compact(temp_data))
  
  if (nrow(combined_data) == 0) {
    stop("No data loaded from result files for condition: ", use_condition)
  }
  
  message("Total rows loaded: ", nrow(combined_data))
  return(combined_data)
}

# Function to calculate metabolic rate
get_mrs <- function(tb, mass, elev, b0, b1, b2, b3, k=8.62*10^-5) {
  rmr <- exp(b0 + b1*log(mass) + b2*(1/(k*(tb+273.15))) + b3*elev)
  return(rmr)
}

# Function to convert metabolic rates to energy losses
get_mr_losses <- function(mrs) {
  o2_ml <- 1/.7 * mrs
  lipid_mg <- o2_ml/2
  loss_kJ <- lipid_mg*39/1000
  return(loss_kJ)
}

# Function to generate energy curves
generate_energy_curves <- function(species, site, pop_data, year_period) {
  # Get parameters for this site-species combination
  p <- pop_data %>% 
    filter(site == !!site, spp == species) %>% 
    slice(1)  # Take first row (assumes consistent params for site-species)
  
  # Skip if no matching data found
  if (nrow(p) == 0) {
    return(NULL)
  }
  
  # Define temperature range
  temp_range <- seq(p$Tmin, p$Topt + p$Above, by = 0.1)
  
  # Calculate TPC (gains) using LRF - already mass standardized
  tpc_values <- sapply(temp_range, function(temp) {
    get_tpc_gains(temp, p$Ropt, p$Topt, p$Tmin, p$Topt + p$Above)
  })
  
  # Calculate metabolic rate (losses)
  mr_values <- sapply(temp_range, function(temp) {
    get_mrs(temp, p$mass, p$elev, p$rmr_b0, p$rmr_b1, p$rmr_b2, p$rmr_b3)
  })
  
  # Only standardize metabolic losses by mass since TPCs are already standardized
  mr_loss_values <- sapply(mr_values, get_mr_losses)
  mr_loss_values_per_g <- mr_loss_values / p$mass
  
  # Calculate net energy gain (per gram)
  net_values <- tpc_values - mr_loss_values_per_g
  
  # Create data frame
  data.frame(
    temperature = temp_range,
    Gains = tpc_values,
    Losses = mr_loss_values_per_g,
    Net = net_values,
    site = site,
    spp = species,
    year_period = year_period
  )
}

# Main function to create the plot
process_energy_plot <- function(base_dir = ".", pops_data, use_condition = "partial_shade_low_veg") {
  # Load temperature data and filter to specific condition
  message("Loading temperature data from RDS files for condition: ", use_condition)
  all_temps <- load_and_process_direct(base_dir, use_condition = use_condition)
  
  message("Loaded temperature data. Number of rows:", nrow(all_temps))
  message("Column names: ", paste(colnames(all_temps), collapse = ", "))
  
  # Check what condition was actually loaded
  if ("condition" %in% colnames(all_temps)) {
    message("Conditions in data: ", paste(unique(all_temps$condition), collapse = ", "))
  }
  if ("shade" %in% colnames(all_temps) && "height" %in% colnames(all_temps)) {
    message("Shade/height combinations: ")
    print(table(all_temps$shade, all_temps$height, useNA = "ifany"))
  }
  
  # Generate energy curves for valid site-species combinations
  message("Generating energy curves...")
  energy_curves <- map_dfr(c("MB", "MS"), function(species) {
    # Filter valid sites for this species
    valid_sites <- valid_combinations[[species]]
    
    # Process each valid site
    map_dfr(valid_sites, function(site) {
      # Skip if site-species combination doesn't exist in pops_data
      if (nrow(filter(pops_data, spp == species, site == !!site)) == 0) {
        message(paste("Skipping", species, site, "- not found in population data"))
        return(NULL)
      }
      
      message(paste("Processing energy curves for", species, "at site", site))
      bind_rows(
        generate_energy_curves(
          species = species, 
          site = site, 
          pop_data = pops_data, 
          year_period = "historical"
        ),
        generate_energy_curves(
          species = species, 
          site = site, 
          pop_data = pops_data, 
          year_period = "current"
        )
      )
    })
  })
  
  # Print summary to verify we have unique curves
  message("Energy curves summary:")
  energy_summary <- energy_curves %>%
    group_by(spp, site) %>%
    summarize(
      min_temp = min(temperature),
      max_temp = max(temperature),
      max_gain = max(Gains),
      max_loss = max(Losses),
      .groups = "drop"
    )
  print(energy_summary)
  
  message("Filtering temperature data...")
  # Filter temperature data to match valid combinations
  filtered_temps <- all_temps %>%
    filter(
      # Ensure we have all required columns
      !is.na(Tb),
      !is.na(species),
      !is.na(site_clim),
      !is.na(year_period),
      
      # Filter to temperatures <= 65°C
      Tb <= 65,
      
      # Filter to valid species-site combinations
      (species == "MS" & site_clim %in% valid_combinations[["MS"]]) |
        (species == "MB" & site_clim %in% valid_combinations[["MB"]])
    ) %>%
    mutate(
      site_clim = factor(site_clim, levels = site_order),
      species = factor(species, levels = c("MB", "MS"))
    )
  
  message("Filtered temperature data. Number of rows:", nrow(filtered_temps))
  message("Temperature range: ", round(min(filtered_temps$Tb), 1), " to ", round(max(filtered_temps$Tb), 1), "°C")
  
  # Print summary of what we're plotting
  temp_summary <- filtered_temps %>%
    group_by(species, site_clim, year_period) %>%
    summarize(
      n_obs = n(),
      mean_temp = mean(Tb),
      min_temp = min(Tb),
      max_temp = max(Tb),
      .groups = "drop"
    )
  message("Temperature data summary by group:")
  print(temp_summary)
  
  # Filter energy curves by temperature
  energy_curves_long <- energy_curves %>%
    filter(temperature <= 65) %>%
    pivot_longer(cols = c(Gains, Losses, Net), 
                 names_to = "type", 
                 values_to = "value") %>%
    mutate(
      type = factor(type, levels = c("Gains", "Losses", "Net")),
      site = factor(site, levels = site_order),
      spp = factor(spp, levels = c("MB", "MS"))
    )
  
  message("Calculating scaling factor...")
  # Calculate scaling factor
  max_density <- max(sapply(split(filtered_temps$Tb, 
                                  paste(filtered_temps$site_clim, filtered_temps$species)), 
                            function(x) {
                              if(length(x) > 1) return(max(stats::density(x)$y))
                              else return(0)
                            }))
  max_rate <- max(abs(energy_curves_long$value))
  scaling_factor <- max_density / max_rate
  
  message("Max density: ", round(max_density, 4))
  message("Max rate: ", round(max_rate, 4))
  message("Scaling factor: ", round(scaling_factor, 4))
  
  # Modify energy_curves_long to match faceting structure
  energy_curves_long <- energy_curves_long %>%
    rename(site_clim = site) %>%  # Rename to match filtered_temps
    mutate(
      site_clim = factor(site_clim, levels = site_order),
      species = factor(spp, levels = c("MB", "MS"))
    )
  
  # Verify we have the correct combinations
  message("Energy curves by site and species:")
  print(table(energy_curves_long$site_clim, energy_curves_long$species))
  
  # Create plot with rotated layout (sites as rows)
  p <- ggplot() +
    # Temperature distributions
    geom_density(data = filtered_temps, 
                 aes(x = Tb, fill = year_period),
                 alpha = 0.5) +
    # Energy curves - now using site_clim and species for faceting
    geom_line(data = energy_curves_long,
              aes(x = temperature, 
                  y = value * scaling_factor, 
                  linetype = type,
                  group = interaction(type, site_clim, species)), # Add group to ensure proper separation
              color = "black",
              linewidth = 1) +
    scale_fill_manual(values = c("historical" = "#FFB3BA", "current" = "#BAFFC9")) +
    scale_linetype_manual(values = c("Gains" = "solid", 
                                     "Losses" = "dashed", 
                                     "Net" = "dotted")) +
    scale_y_continuous(
      name = "Temperature Density",
      sec.axis = sec_axis(~./scaling_factor, name = "Energy Rate (kJ/hr/g)")
    ) +
    scale_x_continuous(limits = c(0, 65)) +  # Limit x-axis to 65°C
    facet_grid(site_clim ~ species, scales = "free_y") +  # Sites as rows
    labs(x = "Temperature (°C)",
         title = paste("Temperature Distributions and Energy Rates -", use_condition),
         fill = "Period") +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines"),
      axis.title.y.right = element_text(color = "darkgreen"),
      legend.position = "bottom"
    )
  
  return(p)
}

# Example usage with different conditions:

# Use partial shade, low vegetation (shade = 0.3, height = 0.03) - YOUR PREFERRED CONDITION
plot_partial_shade <- process_energy_plot(".", pops, use_condition = "partial_shade_low_veg")
print(plot_partial_shade)

# Alternative: Use thermoregulated body temperature (PBT)
# plot_pbt <- process_energy_plot(".", pops, use_condition = "PBT")
# print(plot_pbt)

# Alternative: Use full sun, low vegetation (shade = 0.0, height = 0.03)
# plot_sun <- process_energy_plot(".", pops, use_condition = "sun_low_veg")
# print(plot_sun)

# Save the plot
# ggsave("temperature_energy_rates_partial_shade.png", plot_partial_shade, 
#        width = 12, height = 10, dpi = 300)



















#Old version






library(tidyverse)
library(patchwork)

#directory is output/results/

# LRF function for thermal performance curves
lrf_1991_exact <- function(temp, rmax, topt, tmin, tmax) {
  # Return 0 for temps outside valid range
  if (temp <= tmin || temp >= tmax) {
    return(0)
  }
  
  # Implement the exact equation
  tryCatch({
    numerator <- (temp - tmax) * (temp - tmin)^2
    denominator <- (topt - tmin) * ((topt - tmin) * (temp - topt) - (topt - tmax) * (topt + tmin - 2 * temp))
    
    rate <- rmax * numerator / denominator
    
    # Handle NaN or negative values
    if (is.na(rate) || rate < 0) {
      return(0)
    }
    
    return(rate)
  }, error = function(e) {
    warning(paste("Error in LRF calculation:", e$message))
    return(0)
  })
}

# Function to calculate TPC gains with energy conversion
get_tpc_gains <- function(temp, rmax, topt, tmin, tmax, assim.rate = 0.4) {
  # Calculate LRF-based dry feces production
  dry.fec.mg.hr <- lrf_1991_exact(temp, rmax, topt, tmin, tmax)
  
  # Ensure non-negative values
  dry.fec.mg.hr[dry.fec.mg.hr < 0] <- 0
  
  # Convert to wet grass-adjusted (WGA) production
  # Using Harrison & Fewell conversion
  dry.wga.mg.hr <- dry.fec.mg.hr * 12.8/14.8 * assim.rate / (1 - assim.rate)
  
  # Convert to kilocalories per hour
  kcal.hr <- dry.wga.mg.hr * 14.8 / 1000
  
  # Convert to kilojoules per hour
  kJ.hr <- kcal.hr * 4.184
  
  return(kJ.hr)
}

# Valid species-site combinations
valid_combinations <- list(
  MS = c("Eldo", "A1", "B1"),
  MB = c("A1", "B1", "C1", "D1")
)

# Site order for plotting
site_order <- c("Eldo", "A1", "B1", "C1", "D1")

# Function to load RDS files directly based on your file structure
load_and_process_direct <- function(base_dir = ".", pattern = "^eb_results_") {
  # Find all RDS files matching the pattern
  all_files <- list.files(path = base_dir, 
                          pattern = paste0(pattern, ".*\\.rds$"), 
                          recursive = TRUE, 
                          full.names = TRUE)
  
  if (length(all_files) == 0) {
    stop("No result files found! Check the base directory and pattern.")
  }
  
  # Read and combine all files
  temp_data <- lapply(all_files, function(file) {
    tryCatch({
      data <- readRDS(file)
      # Extract file components from filename
      file_parts <- basename(file) %>% 
        gsub("\\.rds$", "", .) %>%
        strsplit("_") %>% 
        unlist()
      
      # For debugging
      # print(paste("Processing file:", basename(file)))
      # print(paste("Parts:", paste(file_parts, collapse = ", ")))
      
      # Return the data
      data
    }, error = function(e) {
      warning(paste("Error reading file", file, ":", e$message))
      return(NULL)
    })
  })
  
  # Combine all data frames, filtering out NULL results
  combined_data <- bind_rows(compact(temp_data))
  
  if (nrow(combined_data) == 0) {
    stop("No data loaded from result files!")
  }
  
  return(combined_data)
}

# Function to calculate metabolic rate
get_mrs <- function(tb, mass, elev, b0, b1, b2, b3, k=8.62*10^-5) {
  rmr <- exp(b0 + b1*log(mass) + b2*(1/(k*(tb+273.15))) + b3*elev)
  return(rmr)
}

# Function to convert metabolic rates to energy losses
get_mr_losses <- function(mrs) {
  o2_ml <- 1/.7 * mrs
  lipid_mg <- o2_ml/2
  loss_kJ <- lipid_mg*39/1000
  return(loss_kJ)
}

# Function to generate energy curves
generate_energy_curves <- function(species, site, pop_data, year_period) {
  # Get parameters for this site-species combination
  p <- pop_data %>% 
    filter(site == !!site, spp == species) %>% 
    slice(1)  # Take first row (assumes consistent params for site-species)
  
  # Skip if no matching data found
  if (nrow(p) == 0) {
    return(NULL)
  }
  
  # Define temperature range
  temp_range <- seq(p$Tmin, p$Topt + p$Above, by = 0.1)
  
  # Calculate TPC (gains) using LRF - already mass standardized
  tpc_values <- sapply(temp_range, function(temp) {
    get_tpc_gains(temp, p$Ropt, p$Topt, p$Tmin, p$Topt + p$Above)
  })
  
  # Calculate metabolic rate (losses)
  mr_values <- sapply(temp_range, function(temp) {
    get_mrs(temp, p$mass, p$elev, p$rmr_b0, p$rmr_b1, p$rmr_b2, p$rmr_b3)
  })
  
  # Only standardize metabolic losses by mass since TPCs are already standardized
  mr_loss_values <- sapply(mr_values, get_mr_losses)
  mr_loss_values_per_g <- mr_loss_values / p$mass
  
  # Calculate net energy gain (per gram)
  net_values <- tpc_values - mr_loss_values_per_g
  
  # Create data frame
  data.frame(
    temperature = temp_range,
    Gains = tpc_values,
    Losses = mr_loss_values_per_g,
    Net = net_values,
    site = site,
    spp = species,
    year_period = year_period
  )
}

# Main function to create the plot
process_energy_plot <- function(base_dir = ".", pops_data) {
  # Load temperature data directly from RDS files
  message("Loading temperature data from RDS files...")
  all_temps <- load_and_process_direct(base_dir)
  
  message("Loaded temperature data. Number of rows:", nrow(all_temps))
  message("Column names: ", paste(colnames(all_temps), collapse = ", "))
  
  # Generate energy curves for valid site-species combinations
  message("Generating energy curves...")
  energy_curves <- map_dfr(c("MB", "MS"), function(species) {
    # Filter valid sites for this species
    valid_sites <- valid_combinations[[species]]
    
    # Process each valid site
    map_dfr(valid_sites, function(site) {
      # Skip if site-species combination doesn't exist in pops_data
      if (nrow(filter(pops_data, spp == species, site == !!site)) == 0) {
        message(paste("Skipping", species, site, "- not found in population data"))
        return(NULL)
      }
      
      message(paste("Processing energy curves for", species, "at site", site))
      bind_rows(
        generate_energy_curves(
          species = species, 
          site = site, 
          pop_data = pops_data, 
          year_period = "historical"
        ),
        generate_energy_curves(
          species = species, 
          site = site, 
          pop_data = pops_data, 
          year_period = "current"
        )
      )
    })
  })
  
  # Print summary to verify we have unique curves
  message("Energy curves summary:")
  energy_summary <- energy_curves %>%
    group_by(spp, site) %>%
    summarize(
      min_temp = min(temperature),
      max_temp = max(temperature),
      max_gain = max(Gains),
      max_loss = max(Losses),
      .groups = "drop"
    )
  print(energy_summary)
  
  message("Filtering temperature data...")
  # Filter temperature data to match valid combinations
  filtered_temps <- all_temps %>%
    filter(
      # Ensure we have all required columns
      !is.na(Tb),
      !is.na(species),
      !is.na(site_clim),
      !is.na(year_period),
      
      # Filter to temperatures <= 65°C
      Tb <= 65,
      
      # Filter to valid species-site combinations
      (species == "MS" & site_clim %in% valid_combinations[["MS"]]) |
        (species == "MB" & site_clim %in% valid_combinations[["MB"]])
    ) %>%
    mutate(
      site_clim = factor(site_clim, levels = site_order),
      species = factor(species, levels = c("MB", "MS"))
    )
  
  message("Filtered temperature data. Number of rows:", nrow(filtered_temps))
  
  # Filter energy curves by temperature
  energy_curves_long <- energy_curves %>%
    filter(temperature <= 65) %>%
    pivot_longer(cols = c(Gains, Losses, Net), 
                 names_to = "type", 
                 values_to = "value") %>%
    mutate(
      type = factor(type, levels = c("Gains", "Losses", "Net")),
      site = factor(site, levels = site_order),
      spp = factor(spp, levels = c("MB", "MS"))
    )
  
  message("Calculating scaling factor...")
  # Calculate scaling factor
  max_density <- max(sapply(split(filtered_temps$Tb, 
                                  paste(filtered_temps$site_clim, filtered_temps$species)), 
                            function(x) {
                              if(length(x) > 1) return(max(stats::density(x)$y))
                              else return(0)
                            }))
  max_rate <- max(abs(energy_curves_long$value))
  scaling_factor <- max_density / max_rate
  
  message("Creating plot...")
  # Join energy curves with temperature data to ensure correct faceting
  message("Preparing data for plotting...")
  
  # Create a mapping dataframe to connect site and species to the right facet
  site_species_map <- expand.grid(
    site = site_order,
    spp = c("MB", "MS")
  ) %>%
    filter(
      (spp == "MB" & site %in% valid_combinations[["MB"]]) |
        (spp == "MS" & site %in% valid_combinations[["MS"]])
    )
  
  # Modify energy_curves_long to match faceting structure
  energy_curves_long <- energy_curves_long %>%
    rename(site_clim = site) %>%  # Rename to match filtered_temps
    mutate(
      site_clim = factor(site_clim, levels = site_order),
      species = factor(spp, levels = c("MB", "MS"))
    )
  
  # Verify we have the correct combinations
  message("Energy curves by site and species:")
  print(table(energy_curves_long$site_clim, energy_curves_long$species))
  
  # Create plot with rotated layout (sites as rows)
  ggplot() +
    # Temperature distributions
    geom_density(data = filtered_temps, 
                 aes(x = Tb, fill = year_period),
                 alpha = 0.5) +
    # Energy curves - now using site_clim and species for faceting
    geom_line(data = energy_curves_long,
              aes(x = temperature, 
                  y = value * scaling_factor, 
                  linetype = type,
                  group = interaction(type, site_clim, species)), # Add group to ensure proper separation
              color = "black",
              linewidth = 1) +
    scale_fill_manual(values = c("historical" = "#FFB3BA", "current" = "#BAFFC9")) +
    scale_linetype_manual(values = c("Gains" = "solid", 
                                     "Losses" = "dashed", 
                                     "Net" = "dotted")) +
    scale_y_continuous(
      name = "Temperature Density",
      sec.axis = sec_axis(~./scaling_factor, name = "Energy Rate (kJ/hr/g)")
    ) +
    scale_x_continuous(limits = c(0, 65)) +  # Limit x-axis to 65°C
    facet_grid(site_clim ~ species, scales = "free_y") +  # Sites as rows
    labs(x = "Temperature (°C)",
         title = "Temperature Distributions and Energy Rates",
         fill = "Period") +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines"),
      axis.title.y.right = element_text(color = "darkgreen"),
      legend.position = "bottom"
    )
}

# Example usage
#pops <- read.csv("pops.csv") # Load your population data
plot <- process_energy_plot(".", pops)
print(plot)
# ggsave("temperature_energy_rates.png", width = 12, height = 10, dpi = 300)