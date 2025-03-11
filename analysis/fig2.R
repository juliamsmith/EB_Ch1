library(tidyverse)
library(patchwork)

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
# pops <- read.csv("pops.csv") # Load your population data
plot <- process_energy_plot(".", pops)
print(plot)
# ggsave("temperature_energy_rates.png", width = 12, height = 10, dpi = 300)