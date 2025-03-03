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

# Calculate thermal performance curve for feeding
get_tpc_gains <- function(tbs, q10, a, b, c, assim.rate=.40){ # assimilation assumption
  dry.fec.mg.hr <- rezende_2019(tbs, q10, a, b, c)
  dry.fec.mg.hr[dry.fec.mg.hr<0] = 0
  dry.wga.mg.hr <- dry.fec.mg.hr*12.8/14.8*assim.rate/(1-assim.rate) # see harrison & fewell
  kcal.hr <- dry.wga.mg.hr*14.8/1000 # from Fewell & Harrison
  kJ.hr <- kcal.hr*4.184
  return(kJ.hr)
}

# Calculate energy gains based on TPC and temperatures
get_energy_gains <- function(sppi, sitei, sexi, tbs, pops, dts){
  # Calculate time intervals in hours
  dt_ints <- as.numeric(difftime(dts$dtuse, lag(dts$dtuse), units = "hours"))
  dt_ints[1] <- dt_ints[2]  # Handle the first interval
  
  pop_dat <- pops %>% filter(spp==sppi & site==sitei & sex==sexi)
  
  # Get hourly rates
  gains_per_hour <- get_tpc_gains(tbs, 
                                  pop_dat$tpc_q10, 
                                  pop_dat$tpc_a, 
                                  pop_dat$tpc_b, 
                                  pop_dat$tpc_c)
  
  mrs <- get_mrs(tbs, pop_dat$mass[1], pop_dat$elev[1], pop_dat$rmr_b0[1], 
                 pop_dat$rmr_b1[1], pop_dat$rmr_b2[1], pop_dat$rmr_b3[1])
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
calculate_energy_budget <- function(climate_data, height, shade_level, surface_roughness, pops, 
                                    species, sex, site_orig, site_clim, year) {
  
  # Filter climate data for the specific site
  site_climate <- climate_data %>%
    filter(site == site_clim, 
           height == height, 
           shade == shade_level)
  
  # Add zenith angle if not already present
  if(!("psi" %in% colnames(site_climate))) {
    site_climate <- add_psi_to_climate(site_climate)
  }
  
  # Convert shade level to sun level (1 - shade)
  sun_level <- 1 - shade_level
  
  # Calculate body temperature
  site_climate <- site_climate %>%
    mutate(Tb = Tb_grasshopper2.5(
      T_a = T_specified_height,
      T_g = T_soilest,
      u = wind_speed_profile_neutral(ifelse(wsuse==0, .001, wsuse), 1, surface_roughness, height),
      S = sun_level * ifelse(psi!=90 & psi!= -90, sruse, 0),
      K_t = 0.7,
      psi = psi,
      l = 0.03,
      z = height,
      Acondfact = 0.25
    ))
  
  # Get population parameters for the specific species and site
  pop_idx <- which(pops$spp == species & pops$site == site_orig & pops$sex == sex)[1]
  
  if(length(pop_idx) == 0) {
    stop(paste("No population data found for species =", species, 
               ", site =", site_orig, ", sex =", sex))
  }
  
  # Calculate energy budget
  energy_data <- get_energy_gains(species, site_orig, sex, 
                                  site_climate$Tb, pops, 
                                  site_climate)
  
  # Combine results
  result <- cbind(
    site_climate,
    energy_data
  ) %>%
    mutate(
      species = species,
      sex = sex,
      site_orig = site_orig,
      site_clim = site_clim,
      year = year,
      height = height,
      shade = shade_level
    )
  
  return(result)
}