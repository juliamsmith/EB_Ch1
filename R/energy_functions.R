# Improved energy_functions.R

# Calculate metabolic rate based on temperature and other factors
get_mrs <- function(tb, mass, elev, b0, b1, b2, b3, k=8.62*10^-5) {
  # Convert temperatures to Kelvin for Arrhenius equation
  tb_K <- tb + 273.15
  
  # Ensure parameters are numeric
  mass <- as.numeric(mass)
  elev <- as.numeric(elev)
  b0 <- as.numeric(b0)
  b1 <- as.numeric(b1)
  b2 <- as.numeric(b2)
  b3 <- as.numeric(b3)
  
  # Compute metabolic rate with input validation
  rmr <- exp(b0 + b1*log(mass) + b2*(1/(k*tb_K)) + b3*elev)
  
  # Handle extreme values
  rmr[is.na(rmr) | is.infinite(rmr) | rmr < 0] <- NA
  
  return(rmr)
}

# Convert metabolic rates to energy loss
get_mr_losses <- function(mrs) {
  # Handle missing values
  if(all(is.na(mrs))) {
    warning("All metabolic rates are NA")
    return(rep(NA, length(mrs)))
  }
  
  # Convert oxygen consumption to energy loss
  o2.ml <- 1/.7 * mrs
  lipid.mg <- o2.ml/2
  loss.kJ <- lipid.mg*39/1000 # assume 39kJ/g
  
  # Ensure losses are non-negative
  loss.kJ[loss.kJ < 0] <- 0
  
  return(loss.kJ)
}

# Improved LRF function for thermal performance curves
lrf_1991_correct <- function(temp, rmax, topt, tmin, tmax) {
  # Input validation
  if(is.na(temp) || is.na(rmax) || is.na(topt) || is.na(tmin) || is.na(tmax)) {
    return(0)
  }
  
  # Return 0 for temps outside valid range
  if(temp <= tmin || temp >= tmax || topt <= tmin || topt >= tmax) {
    return(0)
  }
  
  # Implement the exact equation with safeguards
  numerator <- (temp - tmax) * (temp - tmin)^2
  denominator <- (topt - tmin) * ((topt - tmin) * (temp - topt) - (topt - tmax) * (topt + tmin - 2 * temp))
  
  # Avoid division by zero
  if(denominator == 0) {
    return(0)
  }
  
  rate <- rmax * numerator / denominator
  
  # Handle NaN or negative values
  if(is.na(rate) || rate < 0) {
    return(0)
  }
  
  return(rate)
}

# Improved energy gains function with safer time interval handling
get_energy_gains <- function(species, site_orig, sex, tbs, pops, dts) {
  # Input validation
  if(length(tbs) == 0 || length(dts$dtuse) == 0) {
    warning("Empty input data in get_energy_gains")
    return(data.frame(gains=numeric(0), losses=numeric(0), net_gains=numeric(0)))
  }
  
  # Sort data by time first to ensure chronological order
  time_order <- order(dts$dtuse)
  dts <- dts[time_order,]
  tbs <- tbs[time_order]
  
  # Calculate time intervals in hours - with safer handling
  dt_ints <- as.numeric(difftime(dts$dtuse, lag(dts$dtuse), units = "hours"))
  
  # For the first interval, use a reasonable default (1 hour is typical)
  dt_ints[1] <- 1
  
  # Handle any negative or unreasonable intervals
  dt_ints <- pmin(pmax(dt_ints, 0), 24)  # Limit to 0-24 hour range
  
  # Find population data
  pop_dat <- pops %>% 
    filter(spp == species & site == site_orig & sex == sex) %>%
    slice(1)
  
  if(nrow(pop_dat) == 0) {
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
  sex_multiplier <- if(sex == "M") {
    exp(sex_effect * 0.5)  # For males (sex=0.5)
  } else {
    exp(sex_effect * -0.5)  # For females (sex=-0.5)
  }
  
  # Calculate hourly rates using the corrected LRF function
  gains_per_hour <- sapply(tbs, function(tb) {
    # Skip calculation for NA temperatures
    if(is.na(tb)) {
      return(0)
    }
    
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
  
  # Ensure no negative values
  gains[is.na(gains) | gains < 0] <- 0
  losses[is.na(losses) | losses < 0] <- 0
  net_gains <- gains - losses
  
  # Return as data frame with clear column names
  df <- data.frame(
    gains = gains,
    losses = losses,
    net_gains = net_gains
  )
  return(df) 
}

# Improved function to calculate body temperatures and energy budget
calculate_energy_budget <- function(climate_data, pops, 
                                    species, sex, site_orig, site_clim, year) {
  
  # Check if climate data exists
  if(nrow(climate_data) == 0) {
    stop("Empty climate data provided")
  }
  
  # Add zenith angle if not already present
  if(!("psi" %in% colnames(climate_data))) {
    climate_data <- add_psi_to_climate(climate_data)
  }
  
  # Calculate body temperature for each row with tryCatch for error handling
  climate_data <- climate_data %>%
    mutate(Tb = NA)  # Initialize Tb column
  
  # Process rows in batches for better error handling
  batch_size <- 1000
  num_batches <- ceiling(nrow(climate_data) / batch_size)
  
  for(batch in 1:num_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(climate_data))
    batch_rows <- start_idx:end_idx
    
    tryCatch({
      climate_data$Tb[batch_rows] <- Tb_grasshopper2.5(
        T_a = climate_data$Tair[batch_rows],              # Use air temperature
        T_g = climate_data$Tsoil[batch_rows],             # Ground temperature
        u = ifelse(climate_data$wind[batch_rows] == 0, 0.001, climate_data$wind[batch_rows]),  # Ensure positive wind
        S = ifelse(psi != 90 & psi != -90, rad, 0),
        K_t = 0.7,               # Clearness index
        psi = climate_data$psi[batch_rows],               # Solar zenith angle
        l = 0.03,                # Grasshopper length (3cm)
        z = climate_data$height[batch_rows],              # Height above ground
        Acondfact = 0.25         # Area conduction factor
      )
    }, error = function(e) {
      warning(paste("Error calculating body temperatures for batch", batch, ":", e$message))
      climate_data$Tb[batch_rows] <- NA
    })
  }
  
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
  
  # Final data cleaning for extreme values
  result <- result %>%
    mutate(
      # Cap unreasonable body temperatures
      Tb = ifelse(Tb > 60, 60, ifelse(Tb < -20, -20, Tb)),
      # Ensure energy values are reasonable
      gains = ifelse(gains < 0, 0, gains),
      losses = ifelse(losses < 0, 0, losses),
      net_gains = gains - losses
    )
  
  return(result)
}
