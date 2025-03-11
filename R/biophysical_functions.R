# Improved biophysical_functions.R

# Fixed Grasshopper body temperature model with safeguards
Tb_grasshopper2.5 <- function (T_a, T_g, u, S, K_t, psi, l, Acondfact = 0.25, z = 0.001, 
                               abs = 0.7, r_g = 0.3) 
{
  # Input validation with informative warnings
  if (any(u < 0)) {
    warning("Wind speed cannot be negative. Setting negative values to 0.")
    u[u < 0] <- 0
  }
  if (any(S < 0)) {
    warning("Solar radiation cannot be negative. Setting negative values to 0.")
    S[S < 0] <- 0
  }
  if (any(K_t < 0) || any(K_t > 1)) {
    warning("Clearness index must be between 0 and 1. Clamping values.")
    K_t[K_t < 0] <- 0
    K_t[K_t > 1] <- 1
  }
  if (any(psi < -90) || any(psi > 90)) {
    warning("Zenith angle must be between -90 and 90 degrees. Clamping values.")
    psi[psi < -90] <- -90
    psi[psi > 90] <- 90
  }
  if (any(l <= 0)) {
    warning("Length must be positive. Setting to default 0.03m.")
    l[l <= 0] <- 0.03
  }
  if (any(Acondfact < 0) || any(Acondfact > 1)) {
    warning("Area conduction factor must be between 0 and 1. Clamping values.")
    Acondfact[Acondfact < 0] <- 0
    Acondfact[Acondfact > 1] <- 1
  }
  if (any(z < 0)) {
    warning("Height above ground cannot be negative. Setting to 0.")
    z[z < 0] <- 0
  }
  if (any(abs < 0) || any(abs > 1)) {
    warning("Absorptivity must be between 0 and 1. Clamping values.")
    abs[abs < 0] <- 0
    abs[abs > 1] <- 1
  }
  if (any(r_g < 0) || any(r_g > 1)) {
    warning("Ground reflectivity must be between 0 and 1. Clamping values.")
    r_g[r_g < 0] <- 0
    r_g[r_g > 1] <- 1
  }
  
  # Ensure positive wind speed (avoid division by zero)
  u[u == 0] <- 0.001
  
  # Convert temperatures to Kelvin
  T_a <- celsius_to_kelvin(T_a)
  T_g <- celsius_to_kelvin(T_g)
  
  # Constants
  sigma <- stefan_boltzmann_constant()
  epsilon <- 1
  Kf <- 0.025
  v <- 15.68 * 10^-6
  
  # Morphological calculations
  c <- l/2
  a <- (0.365 + 0.241 * l * 1000)/1000
  e <- sqrt(1 - a^2/c^2)
  A <- 2 * pi * a^2 + 2 * pi * a * c/e * asin(e)
  
  # Calculate diffuse fraction (kd) based on clearness index (K_t)
  kd <- 1 - 0.09 * K_t
  kd[K_t > 0.22 & K_t <= 0.8] <- 0.9511 - 0.1604 * K_t[K_t > 0.22 & K_t <= 0.8] + 
    4.388 * K_t[K_t > 0.22 & K_t <= 0.8]^2 - 
    16.638 * K_t[K_t > 0.22 & K_t <= 0.8]^3 + 
    12.336 * K_t[K_t > 0.22 & K_t <= 0.8]^4
  kd[K_t > 0.8] <- 0.165
  
  # Split total solar radiation into direct and diffuse
  Sttl <- S
  Sdir <- Sttl * (1 - kd)
  Sdif <- Sttl * kd
  
  # Convert zenith angle to radians
  psi_r <- degrees_to_radians(psi)
  
  # Heat transfer calculations
  Re <- u * l/v
  Nu <- 0.41 * Re^0.5
  h_c <- Nu * Kf/l
  hc_s <- h_c * (-0.007 * z/l + 1.71)
  Thick <- .025
  hcut <- 0.15
  Acond <- A * Acondfact
  
  # Calculate solar projection area
  sa <- 0.19 - 0.00173 * psi
  Adir <- A * sa
  Aref <- Adir
  
  # Calculate radiative heat fluxes - WITH SAFEGUARDS FOR HIGH ZENITH ANGLES
  # For direct radiation, add safeguard when sun is near horizon (cos(psi_r) near zero)
  cos_psi_r <- cos(psi_r)
  # Set a minimum threshold to avoid division by very small numbers
  cos_psi_r_safe <- pmax(cos_psi_r, 0.01)
  
  # Set Qdir to zero when cos(psi_r) is very small or negative (sun below/at horizon)
  Qdir <- ifelse(cos_psi_r > 0.01, 
                 abs * Adir * Sdir / cos_psi_r_safe, 
                 0)
  
  Qdif <- abs * Aref * Sdif
  Qref <- r_g * Aref * Sttl
  Qabs <- Qdir + Qdif + Qref
  
  # Sky temperature calculation
  T_sky <- 0.0552 * (T_a)^1.5
  
  # Heat balance equation coefficients
  a <- A * epsilon * sigma
  b <- hc_s * A + hcut * Acond/Thick
  d <- hc_s * A * T_a + 0.5 * A * epsilon * sigma * (T_sky^4 + T_g^4) + 
    hcut * Acond * T_g/Thick + Qabs
  
  # Body temperature calculation (complex equation from physics model)
  # With safeguards to catch numerical issues
  T_b <- 1/2 * sqrt((2 * b)/(a * sqrt((sqrt(3) * sqrt(256 * 
                                                        a^3 * d^3 + 27 * a^2 * b^4) + 9 * a * b^2)^(1/3)/(2^(1/3) * 
                                                                                                            3^(2/3) * a) - (4 * (2/3)^(1/3) * d)/(sqrt(3) * sqrt(256 * 
                                                                                                                                                                   a^3 * d^3 + 27 * a^2 * b^4) + 9 * a * b^2)^(1/3))) - 
                      (sqrt(3) * sqrt(256 * a^3 * d^3 + 27 * a^2 * b^4) + 
                         9 * a * b^2)^(1/3)/(2^(1/3) * 3^(2/3) * a) + (4 * 
                                                                         (2/3)^(1/3) * d)/(sqrt(3) * sqrt(256 * a^3 * d^3 + 27 * 
                                                                                                            a^2 * b^4) + 9 * a * b^2)^(1/3)) - 1/2 * sqrt((sqrt(3) * 
                                                                                                                                                             sqrt(256 * a^3 * d^3 + 27 * a^2 * b^4) + 9 * a * b^2)^(1/3)/(2^(1/3) * 
                                                                                                                                                                                                                            3^(2/3) * a) - (4 * (2/3)^(1/3) * d)/(sqrt(3) * sqrt(256 * 
                                                                                                                                                                                                                                                                                   a^3 * d^3 + 27 * a^2 * b^4) + 9 * a * b^2)^(1/3))
  
  # Handle NA or extreme values
  T_b[which(is.na(T_b))] <- NA
  
  # Convert back to Celsius
  T_b_celsius <- kelvin_to_celsius(T_b)
  
  # Limit to physically plausible range (-20 to 60°C)
  T_b_celsius <- pmin(pmax(T_b_celsius, -20), 60)
  
  return(T_b_celsius)
}

# Add psi (zenith angle) to climate data - with error handling
add_psi_to_climate <- function(climate_data) {
  # Check if required columns exist
  required_cols <- c("dtuse", "site")
  missing_cols <- required_cols[!required_cols %in% names(climate_data)]
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse=", ")))
  }
  
  # Handle time zones consistently
  if (!inherits(climate_data$dtuse, "POSIXct")) {
    warning("dtuse is not a POSIXct object, attempting to convert")
    climate_data$dtuse <- as.POSIXct(climate_data$dtuse)
  }
  
  # Use tryCatch to handle errors in get_psi calls
  climate_data <- climate_data %>% 
    rowwise() %>% 
    mutate(psi = tryCatch(
      get_psi(dtuse, site),
      error = function(e) {
        warning(paste("Error calculating psi for site", site, "at time", dtuse, ":", e$message))
        return(90) # Default to horizon (no direct radiation)
      }
    ))
  
  return(climate_data)
}

# Improved function for zenith angle calculation
get_psi <- function(dt, site) {
  # Site coordinates lookup with error handling
  sites <- list(
    Eldo = list(lat = 39.9436, lon = -105.262),
    A1 = list(lat = 40.015, lon = -105.376),
    B1 = list(lat = 40.019, lon = -105.455),
    C1 = list(lat = 40.0301, lon = -105.541),
    D1 = list(lat = 40.0401, lon = -105.600)
  )
  
  if (!site %in% names(sites)) {
    stop(paste("Unknown site:", site))
  }
  
  lat <- sites[[site]]$lat
  lon <- sites[[site]]$lon
  
  # Robust date/time handling
  if (inherits(dt, "character")) {
    date <- as.POSIXct(dt, format = "%Y-%m-%d")
    if (is.na(date)) {
      stop(paste("Could not parse date string:", dt))
    }
  } else if (inherits(dt, "POSIXct")) {
    date <- dt
  } else {
    stop("dt must be a POSIXct object or a character string")
  }
  
  # Extract hour
  hour <- as.numeric(format(date, format="%H")) + 
    as.numeric(format(date, format="%M"))/60
  
  # Get day of year
  doy <- day_of_year(day=date, format="%Y-%m-%d")
  
  # Calculate zenith angle
  zenith <- zenith_angle(doy, lat, lon, hour)
  
  # Ensure zenith is within valid range
  zenith <- max(min(zenith, 90), -90)
  
  return(zenith)
}