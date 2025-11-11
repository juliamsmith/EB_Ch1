library(brms)
library(tidyverse)
library(ggplot2)

ms_fit <- readRDS("~/GitHub/thermal_perf/ms_improved7_betterprior.rds")
mb_fit <- readRDS("~/GitHub/thermal_perf/mb_both_years_with_yeareffect3betterpriors.rds")

mb_both_years_fit <- mb_fit

# Prepare data for MB both years (2023 and 2024)
prepare_mb_both_years <- function() {
  cat("Preparing MB data for both 2023 and 2024...\n")
  
  # Read data
  data <- read.csv("data/ad.csv")
  cat("Total rows in data:", nrow(data), "\n")
  
  # Filter for MB both years
  data_subset <- data %>% 
    filter(spp == "MB", yr %in% c(2023, 2024))
  
  cat("Rows after filtering for MB 2023-2024:", nrow(data_subset), "\n")
  cat("Data by year:\n")
  print(table(data_subset$yr))
  
  # Calculate min_nonzero for adjustment (across both years)
  min_nonzero <- min(data_subset$rate[data_subset$rate > 0], na.rm = TRUE)
  cat("Min nonzero rate (both years):", min_nonzero, "\n")
  
  # Process data
  model_data <- data_subset %>%
    mutate(
      # Site as unordered factor
      pop = factor(site, levels = c("A1", "B1", "C1", "D1")),
      # Sex as -0.5/0.5 coding (M = 0.5, F = -0.5)
      sex = if_else(sex == "M", 0.5, -0.5),
      # Year effect: 0 = 2024 (reference), 1 = 2023
      yeareffect = if_else(yr == 2024, 0, 1),
      # Individual ID (make unique across years)
      full_ID = factor(paste(full_ID, yr, sep = "_")),
      # Adjusted response variable
      fec_double_adj = rate + min_nonzero/100
    ) %>%
    select(fec_double_adj, temp, sex, pop, full_ID, yeareffect, yr)
  
  cat("Final data dimensions:", nrow(model_data), "x", ncol(model_data), "\n")
  cat("Number of unique individuals:", length(unique(model_data$full_ID)), "\n")
  cat("Temperature range:", min(model_data$temp), "to", max(model_data$temp), "\n")
  cat("Year effect coding:\n")
  print(table(model_data$yr, model_data$yeareffect))
  
  return(model_data)
}

# Prepare data for ms (simplified - no standardization)
prepare_ms <- function() {
  cat("Preparing ms data...\n")
  
  # Read data
  data <- read.csv("data/ad.csv")
  cat("Total rows in data:", nrow(data), "\n")
  
  # Filter for ms only
  data_subset <- data %>% 
    filter(spp == "MS")
  
  cat("Rows after filtering for ms:", nrow(data_subset), "\n")
  
  # Calculate min_nonzero for adjustment
  min_nonzero <- min(data_subset$rate[data_subset$rate > 0], na.rm = TRUE)
  cat("Min nonzero rate:", min_nonzero, "\n")
  
  # Process data - simplified, no temperature standardization
  model_data <- data_subset %>%
    mutate(
      # Site as unordered factor
      pop = factor(site, levels = c("Eldo", "A1", "B1")),
      # Sex as -0.5/0.5 coding (M = 0.5, F = -0.5)
      sex = if_else(sex == "M", 0.5, -0.5),
      # Individual ID
      full_ID = factor(full_ID),
      # Adjusted response variable
      fec_double_adj = rate + min_nonzero/100
    ) %>%
    select(fec_double_adj, temp, sex, pop, full_ID)
  
  cat("Final data dimensions:", nrow(model_data), "x", ncol(model_data), "\n")
  cat("Number of unique individuals:", length(unique(model_data$full_ID)), "\n")
  cat("Temperature range:", min(model_data$temp), "to", max(model_data$temp), "\n")
  
  return(model_data)
}

mb_data <- prepare_mb_both_years()
ms_data <- prepare_ms()

library(brms)
library(tidyverse)
library(ggplot2)

# Define LRF function
LRF <- function(temp, Tmin, Tmax, Topt, Ropt){
  ifelse(temp > Tmax | temp < Tmin, 0, 
         Ropt * ((temp + 273.15) - (Tmax + 273.15)) * 
           ((temp + 273.15) - (Tmin + 273.15)) ^ 2 / 
           (((Topt + 273.15) - (Tmin + 273.15)) * 
              (((Topt + 273.15) - (Tmin + 273.15)) * 
                 ((temp + 273.15) - (Topt + 273.15)) - 
                 ((Topt + 273.15) - (Tmax + 273.15)) * 
                 ((Topt + 273.15) + (Tmin + 273.15) - 2 * (temp + 273.15)))))
}

cat("=== SIMPLE TPC PLOTS USING PARAMETER ESTIMATES ===\n")
cat("Using mean parameter estimates for fast, clean plotting\n\n")

# =====================================
# EXTRACT PARAMETER ESTIMATES
# =====================================

# Function to extract TPC parameters from model
extract_tpc_params <- function(model, populations) {
  
  # Get fixed effects estimates (means)
  fixed_effects <- fixef(model)[, "Estimate"]
  
  # Extract parameters for each population
  param_list <- list()
  
  for(pop in populations) {
    param_list[[pop]] <- list(
      Tmin = fixed_effects[paste0("Tmin_pop", pop)],
      Topt = fixed_effects[paste0("Topt_pop", pop)],
      Above = fixed_effects[paste0("Above_pop", pop)],
      Ropt = fixed_effects[paste0("Ropt_pop", pop)],
      Tmax = fixed_effects[paste0("Topt_pop", pop)] + fixed_effects[paste0("Above_pop", pop)]
    )
  }
  
  # Extract other effects if they exist
  other_effects <- list()
  if("yeareffectheight_Intercept" %in% names(fixed_effects)) {
    other_effects$year_effect <- fixed_effects[["yeareffectheight_Intercept"]]
  }
  if("sexeffect_Intercept" %in% names(fixed_effects)) {
    other_effects$sex_effect <- fixed_effects[["sexeffect_Intercept"]]
  }
  
  return(list(populations = param_list, effects = other_effects))
}

# =====================================
# CREATE TPC PREDICTION CURVES WITH SEX EFFECTS
# =====================================

create_simple_tpc_curves <- function(params, temp_range = c(5, 60), temp_step = 0.1) {
  
  # Finer resolution to avoid wobbly lines
  temp_seq <- seq(temp_range[1], temp_range[2], by = temp_step)
  curves <- data.frame()
  
  for(pop_name in names(params$populations)) {
    pop_params <- params$populations[[pop_name]]
    
    # Calculate base TPC curve
    tpc_values <- LRF(temp_seq, pop_params$Tmin, pop_params$Tmax, 
                      pop_params$Topt, pop_params$Ropt)
    
    # Sex effects: Male (sex = 0.5), Female (sex = -0.5)
    sex_effect <- params$effects$sex_effect
    male_multiplier <- exp(sex_effect * 0.5)
    female_multiplier <- exp(sex_effect * -0.5)
    
    # Create curves for different years and sexes
    if("year_effect" %in% names(params$effects)) {
      # MB species: has year effect
      year_multiplier_2024 <- 1.0
      year_multiplier_2023 <- exp(params$effects$year_effect)
      
      # 2024 curves
      curves <- rbind(curves, data.frame(
        temperature = temp_seq,
        performance = tpc_values * year_multiplier_2024 * male_multiplier,
        population = pop_name,
        year = "2024",
        sex = "Male"
      ))
      
      curves <- rbind(curves, data.frame(
        temperature = temp_seq,
        performance = tpc_values * year_multiplier_2024 * female_multiplier,
        population = pop_name,
        year = "2024", 
        sex = "Female"
      ))
      
      # 2023 curves
      curves <- rbind(curves, data.frame(
        temperature = temp_seq,
        performance = tpc_values * year_multiplier_2023 * male_multiplier,
        population = pop_name,
        year = "2023",
        sex = "Male"
      ))
      
      curves <- rbind(curves, data.frame(
        temperature = temp_seq,
        performance = tpc_values * year_multiplier_2023 * female_multiplier,
        population = pop_name,
        year = "2023",
        sex = "Female"
      ))
      
    } else {
      # MS species: single year
      curves <- rbind(curves, data.frame(
        temperature = temp_seq,
        performance = tpc_values * male_multiplier,
        population = pop_name,
        year = "single",
        sex = "Male"
      ))
      
      curves <- rbind(curves, data.frame(
        temperature = temp_seq,
        performance = tpc_values * female_multiplier,
        population = pop_name,
        year = "single",
        sex = "Female"
      ))
    }
  }
  
  return(curves)
}

# =====================================
# PLOTTING FUNCTIONS
# =====================================

# MS Species Plot
plot_ms_simple <- function(ms_model, ms_data) {
  
  cat("Creating MS plot with parameter estimates and sex-specific curves...\n")
  
  # Extract parameters
  ms_populations <- c("Eldo", "A1", "B1")
  ms_params <- extract_tpc_params(ms_model, ms_populations)
  
  # Create TPC curves
  ms_curves <- create_simple_tpc_curves(ms_params)
  
  # Prepare data - FIXED: add population column to match curves
  ms_plot_data <- ms_data %>%
    mutate(
      sex_label = if_else(sex == 0.5, "Male", "Female"),
      population = pop  # Add population column to match curves
    )
  
  # Create plot
  ggplot() +
    # TPC curves by sex - Now properly matched to facets
    geom_line(data = ms_curves,
              aes(x = temperature, y = performance, color = sex),
              size = 1.2, alpha = 0.8) +
    
    # Data points
    geom_point(data = ms_plot_data,
               aes(x = temp, y = fec_double_adj, color = sex_label),
               size = 2, alpha = 0.7) +
    
    # Facet by population - now both data and curves have 'population' column
    facet_wrap(~population, ncol = 3,
               labeller = labeller(population = function(x) paste("Population", x))) +
    
    # Styling
    scale_color_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +
    scale_y_continuous(limits = c(0, 15)) +  # Consistent y-axis
    labs(
      title = "MS Species: Thermal Performance Curves by Sex",
      subtitle = paste("Sex effect:", round(ms_params$effects$sex_effect, 3), 
                       "| Ropt range:", round(min(sapply(ms_params$populations, function(x) x$Ropt)), 2),
                       "-", round(max(sapply(ms_params$populations, function(x) x$Ropt)), 2)),
      x = "Temperature (°C)",
      y = "Performance Rate",
      color = "Sex"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40")
    )
}

# MB Species Plot
plot_mb_simple <- function(mb_model, mb_data) {
  
  cat("Creating MB plot with parameter estimates and sex-specific curves...\n")
  
  # Extract parameters
  mb_populations <- c("A1", "B1", "C1", "D1")
  mb_params <- extract_tpc_params(mb_model, mb_populations)
  
  # Create TPC curves
  mb_curves <- create_simple_tpc_curves(mb_params)
  
  # Prepare data - FIXED: use consistent population column name
  mb_plot_data <- mb_data %>%
    mutate(
      sex_label = if_else(sex == 0.5, "Male", "Female"),
      year_label = as.character(yr),
      population = pop  # Add population column to match curves
    )
  
  # Create plot
  ggplot() +
    # TPC curves by sex and year - FIXED: will now match to correct facets
    geom_line(data = mb_curves,
              aes(x = temperature, y = performance, 
                  color = sex, linetype = year),
              size = 1.2, alpha = 0.8) +
    
    # Data points
    geom_point(data = mb_plot_data,
               aes(x = temp, y = fec_double_adj, 
                   color = sex_label, shape = factor(yr)),
               size = 2.5, alpha = 0.7) +
    
    # Facet by population - this will now properly match
    facet_wrap(~population, ncol = 2,
               labeller = labeller(population = function(x) paste("Population", x))) +
    
    # Styling
    scale_color_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +
    scale_linetype_manual(values = c("2024" = "solid", "2023" = "dashed")) +
    scale_shape_manual(values = c("2023" = 17, "2024" = 16),
                       labels = c("2023", "2024")) +
    scale_y_continuous(limits = c(0, 15)) +  # Consistent y-axis
    
    labs(
      title = "MB Species: Thermal Performance Curves by Sex and Year",
      subtitle = paste("Year effect:", round(mb_params$effects$year_effect, 3), 
                       "| Sex effect:", round(mb_params$effects$sex_effect, 3),
                       "| Ropt range:", round(min(sapply(mb_params$populations, function(x) x$Ropt)), 2),
                       "-", round(max(sapply(mb_params$populations, function(x) x$Ropt)), 2)),
      x = "Temperature (°C)",
      y = "Performance Rate",
      color = "Sex",
      shape = "Year",
      linetype = "Year"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40")
    )
}

# =====================================
# GENERATE PLOTS
# =====================================

cat("Generating simplified TPC plots...\n\n")

# Create plots (assuming your models are loaded)
ms_plot_simple <- plot_ms_simple(ms_fit, ms_data)
mb_plot_simple <- plot_mb_simple(mb_both_years_fit, mb_data)

# Display plots
print(ms_plot_simple)
print(mb_plot_simple)

# Save plots
#ggsave("MS_TPC_simple.png", ms_plot_simple, width = 12, height = 6, dpi = 300)
#ggsave("MB_TPC_simple.png", mb_plot_simple, width = 12, height = 8, dpi = 300)

cat("=== RESULTS SUMMARY ===\n")
#cat("=" %R% 25, "\n")

# Extract and display key results
ms_params <- extract_tpc_params(ms_fit, c("Eldo", "A1", "B1"))
mb_params <- extract_tpc_params(mb_both_years_fit, c("A1", "B1", "C1", "D1"))

cat("MS SPECIES THERMAL PARAMETERS:\n")
for(pop in names(ms_params$populations)) {
  p <- ms_params$populations[[pop]]
  cat(sprintf("Population %s: Tmin=%.1f, Topt=%.1f, Tmax=%.1f, Ropt=%.2f\n", 
              pop, p$Tmin, p$Topt, p$Tmax, p$Ropt))
}

cat("\nMB SPECIES THERMAL PARAMETERS:\n")
for(pop in names(mb_params$populations)) {
  p <- mb_params$populations[[pop]]
  cat(sprintf("Population %s: Tmin=%.1f, Topt=%.1f, Tmax=%.1f, Ropt=%.2f\n", 
              pop, p$Tmin, p$Topt, p$Tmax, p$Ropt))
}

cat("\nMB YEAR EFFECT:", round(mb_params$effects$year_effect, 3))
cat(" (2023 performance = ", round(exp(mb_params$effects$year_effect), 3), " × 2024 performance)\n")

cat("\nSEX EFFECTS:\n")
cat("MS sex effect:", round(ms_params$effects$sex_effect, 3), "\n")
cat("MB sex effect:", round(mb_params$effects$sex_effect, 3), "\n")

