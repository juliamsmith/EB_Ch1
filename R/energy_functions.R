# R/energy_functions.R

# Calculate metabolic rate based on temperature and other factors
get_mrs <- function(tb, mass, elev, b0, b1, b2, b3, k=8.62*10^-5) {
  rmr <- exp(b0+b1*log(mass) + b2*(1/(k*(tb+273.15))) + b3*elev)
  return(rmr)
}

# Convert metabolic rates to energy loss
get_mr_losses <- function(mrs){
  o2.ml <- 1/.7*mrs
  lipid.mg <- o2.ml/2
  loss.kJ <- lipid.mg*39/1000 # assume 39kJ/g
  return(loss.kJ)
}

# # LRF function for thermal performance curves
# lrf_function <- function(temp, Tmin, Tmax, Topt, Ropt) {
#   if (temp <= Tmin || temp >= Tmax || Topt <= Tmin || Topt >= Tmax) {
#     return(1e-10)  # Return a small number for invalid temperatures
#   }
#   
#   # Calculate phi
#   phi <- (Tmax - Topt) / (Topt - Tmin)
#   
#   # Calculate the main function
#   numerator <- Ropt * ((Tmax - temp) / (Tmax - Topt)) * exp(phi * (temp - Topt) / (Tmax - Topt))
#   denominator <- 1 + exp(phi * (temp - Topt) / (Tmax - Topt))
#   
#   return(numerator / denominator)
# }
# 
# # Function to calculate feeding rate for MS using LRF model
# get_ms_feeding_rate <- function(temp, site, sex) {
#   # Parameters for MS from the model output
#   site_params <- list(
#     Eldo = list(Tmin = 14.74, Topt = 39.88, Above = 9.56, Ropt = 3.50),
#     A1 = list(Tmin = 11.67, Topt = 39.89, Above = 6.68, Ropt = 3.63),
#     B1 = list(Tmin = 10.96, Topt = 39.78, Above = 8.20, Ropt = 4.30)
#   )
#   
#   # For D1, use A1 parameters as a placeholder
#   if (site == "D1") {
#     site_params$D1 <- site_params$A1
#   }
#   
#   # Sex effect
#   sex_effect <- -0.37
#   
#   # Get parameters for this site
#   params <- site_params[[site]]
#   
#   # Calculate Tmax (Topt + Above)
#   Tmax <- params$Topt + params$Above
#   
#   # Calculate base feeding rate using LRF
#   base_rate <- lrf_function(temp, params$Tmin, Tmax, params$Topt, params$Ropt)
#   
#   # Apply sex effect if male
#   if (sex == "M") {
#     base_rate <- exp(log(base_rate) + sex_effect)
#   }
#   
#   return(base_rate)
# }
# 
# # Function to calculate feeding rate for MB using LRF model
# get_mb_feeding_rate <- function(temp, site, sex) {
#   # Parameters for MB from the model output
#   site_params <- list(
#     A1 = list(Tmin = 13.57, Topt = 37.54, Above = 11.70, Ropt = 2.64),
#     B1 = list(Tmin = 12.80, Topt = 41.76, Above = 10.74, Ropt = 3.29),
#     C1 = list(Tmin = 13.79, Topt = 41.03, Above = 11.77, Ropt = 3.93)
#   )
#   
#   # For Eldo and D1, use A1 and C1 parameters respectively as placeholders
#   site_params$Eldo <- site_params$A1
#   site_params$D1 <- list(Tmin = 14.0, Topt = 42.0, Above = 12.0, Ropt = 4.0)  # Made-up values for D1
#   
#   # Sex effect
#   sex_effect <- -0.06
#   
#   # Get parameters for this site
#   params <- site_params[[site]]
#   
#   # Calculate Tmax (Topt + Above)
#   Tmax <- params$Topt + params$Above
#   
#   # Calculate base feeding rate using LRF
#   base_rate <- lrf_function(temp, params$Tmin, Tmax, params$Topt, params$Ropt)
#   
#   # Apply sex effect if male
#   if (sex == "M") {
#     base_rate <- exp(log(base_rate) + sex_effect)
#   }
#   
#   return(base_rate)
# }
# 
# # Calculate feeding gains for a temperature series
# calculate_feeding_gains <- function(tbs, species, site, sex, assim.rate=0.40) {
#   # Choose the appropriate model based on species
#   if (species == "MS") {
#     dry_fec_mg_hr <- sapply(tbs, function(tb) get_ms_feeding_rate(tb, site, sex))
#   } else if (species == "MB") {
#     dry_fec_mg_hr <- sapply(tbs, function(tb) get_mb_feeding_rate(tb, site, sex))
#   } else {
#     stop(paste("Unknown species:", species))
#   }
#   
#   # Convert to energy gains (adapted from get_tpc_gains)
#   dry_fec_mg_hr[dry_fec_mg_hr < 0] = 0
#   dry_wga_mg_hr <- dry_fec_mg_hr * 12.8/14.8 * assim.rate/(1-assim.rate) # see Harrison & Fewell
#   kcal_hr <- dry_wga_mg_hr * 14.8/1000 # from Fewell & Harrison
#   kJ_hr <- kcal_hr * 4.184
#   
#   return(kJ_hr)
# }

# Update to get_energy_gains function to use the correct LRF formula
# Update to get_energy_gains function with correct sex effect application
get_energy_gains <- function(species, site_orig, sex, tbs, pops, dts){
  # Calculate time intervals in hours
  dt_ints <- as.numeric(difftime(dts$dtuse, lag(dts$dtuse), units = "hours"))
  dt_ints[1] <- dt_ints[2]  # Handle the first interval
  
  # Find population data
  pop_dat <- pops %>% 
    filter(spp == species & site == site_orig & sex == sex) %>%
    slice(1)
  
  if (nrow(pop_dat) == 0) {
    stop(paste("No population data found for species =", species, 
               ", site =", site_orig, ", sex =", sex))
  }
  
  # Extract parameters
  Tmin <- pop_dat$Tmin[1]
  Topt <- pop_dat$Topt[1]
  Above <- pop_dat$Above[1]
  Tmax <- Topt + Above
  Ropt <- pop_dat$Ropt[1]
  
  # Get sex effect coefficient
  sex_effect <- ifelse(species == "MS", -0.37, -0.06)
  
  # Calculate sex multiplier based on coding (M=0.5, F=-0.5)
  sex_multiplier <- if (sex == "M") {
    exp(sex_effect * 0.5)  # For males (sex=0.5)
  } else {
    exp(sex_effect * -0.5)  # For females (sex=-0.5)
  }
  
  # Define correct LRF function
  lrf_1991_correct <- function(temp, rmax, topt, tmin, tmax) {
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
  
  # Get hourly rates using the correct LRF function
  gains_per_hour <- sapply(tbs, function(tb) {
    # Calculate base rate with correct LRF function
    base_rate <- lrf_1991_correct(tb, Ropt, Topt, Tmin, Tmax)
    
    # Apply sex effect
    adjusted_rate <- base_rate * sex_multiplier
    
    # Convert to energy gain
    # Assimilation rate of 0.40
    dry_fec_mg_hr <- adjusted_rate
    dry_wga_mg_hr <- dry_fec_mg_hr * 12.8/14.8 * 0.40/(1-0.40)
    kcal_hr <- dry_wga_mg_hr * 14.8/1000
    kJ_hr <- kcal_hr * 4.184
    
    return(kJ_hr)
  })
  
  # Calculate metabolic rates and losses
  mass <- pop_dat$mass[1]
  elev <- pop_dat$elev[1]
  rmr_b0 <- pop_dat$rmr_b0[1]
  rmr_b1 <- pop_dat$rmr_b1[1]
  rmr_b2 <- pop_dat$rmr_b2[1]
  rmr_b3 <- pop_dat$rmr_b3[1]
  
  mrs <- get_mrs(tbs, mass, elev, rmr_b0, rmr_b1, rmr_b2, rmr_b3)
  losses_per_hour <- get_mr_losses(mrs)
  
  # Convert to interval amounts
  gains <- gains_per_hour * dt_ints
  losses <- losses_per_hour * dt_ints
  net_gains <- gains - losses
  
  # Return as data frame with clear column names
  df <- data.frame(
    gains = gains,
    losses = losses,
    net_gains = net_gains
  )
  return(df) 
}

# Calculate body temperatures and energy budget for a given scenario
calculate_energy_budget <- function(climate_data, pops, 
                                    species, sex, site_orig, site_clim, year) {
  
  # Add zenith angle if not already present
  if (!("psi" %in% colnames(climate_data))) {
    climate_data <- add_psi_to_climate(climate_data)
  }
  
  # Calculate body temperature for each row
  climate_data <- climate_data %>%
    mutate(Tb = Tb_grasshopper2.5(
      T_a = Tair,              # Use the height-specific air temperature directly
      T_g = Tsoil,             # Ground temperature
      u = ifelse(wind == 0, 0.001, wind),  # Use height-specific wind directly
      S = (1 - shade) * ifelse(psi != 90 & psi != -90, rad, 0),  # Convert shade to sun level
      K_t = 0.7,               # Clearness index
      psi = psi,               # Solar zenith angle
      l = 0.03,                # Grasshopper length (3cm)
      z = height,              # Height above ground
      Acondfact = 0.25         # Area conduction factor
    ))
  
  # Calculate energy budget
  energy_data <- get_energy_gains(species, site_orig, sex, 
                                  climate_data$Tb, pops, 
                                  climate_data)
  
  # Combine results
  result <- cbind(climate_data, energy_data) %>%
    mutate(
      species = species,
      sex = sex,
      site_orig = site_orig,
      site_clim = site_clim,
      year_period = ifelse(year >= 2005, "recent", "historical")
    )
  
  return(result)
}
