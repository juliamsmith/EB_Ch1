# R/biophysical_functions.R

# Grasshopper body temperature model
Tb_grasshopper2.5 <- function (T_a, T_g, u, S, K_t, psi, l, Acondfact = 0.25, z = 0.001, 
                               abs = 0.7, r_g = 0.3) 
{
  stopifnot(u >= 0, S >= 0, K_t >= 0, K_t <= 1, psi >= -90, 
            psi <= 90, l >= 0, Acondfact >= 0, Acondfact <= 1, z >= 
              0, abs >= 0, abs <= 1, r_g >= 0, r_g <= 1)
  T_a <- celsius_to_kelvin(T_a)
  T_g <- celsius_to_kelvin(T_g)
  sigma <- stefan_boltzmann_constant()
  epsilon <- 1
  Kf <- 0.025
  v <- 15.68 * 10^-6
  c <- l/2
  a <- (0.365 + 0.241 * l * 1000)/1000
  e <- sqrt(1 - a^2/c^2)
  A <- 2 * pi * a^2 + 2 * pi * a * c/e * asin(e)
  kd <- 1 - 0.09 * K_t
  kd[K_t > 0.22 & K_t <= 0.8] <- 0.9511 - 0.1604 * K_t + 4.388 * 
    K_t^2 - 16.638 * K_t^3 + 12.336 * K_t^4
  kd[K_t > 0.8] <- 0.165
  Sttl <- S
  Sdir <- Sttl * (1 - kd)
  Sdif <- Sttl * kd
  psi_r <- degrees_to_radians(psi)
  Re <- u * l/v
  Nu <- 0.41 * Re^0.5
  h_c <- Nu * Kf/l
  hc_s <- h_c * (-0.007 * z/l + 1.71)
  Thick <- .025
  hcut <- 0.15
  Acond <- A * Acondfact
  sa <- 0.19 - 0.00173 * psi
  Adir <- A * sa
  Aref <- Adir
  Qdir <- abs * Adir * Sdir/cos(psi_r)
  Qdif <- abs * Aref * Sdif
  Qref <- r_g * Aref * Sttl
  Qabs <- Qdir + Qdif + Qref
  T_sky <- 0.0552 * (T_a)^1.5
  a <- A * epsilon * sigma
  b <- hc_s * A + hcut * Acond/Thick
  d <- hc_s * A * T_a + 0.5 * A * epsilon * sigma * (T_sky^4 + 
                                                       T_g^4) + hcut * Acond * T_g/Thick + Qabs
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
  T_b[which(is.na(T_b))] <- NA
  kelvin_to_celsius(T_b)
}

# Calculate zenith angle of sun
get_psi <- function(dt, site){
  if(site=="Eldo"){
    lat <- 39.9436
    lon <- -105.262
  } else if(site=="A1") {
    lat <- 40.015
    lon <- -105.376
  } else if(site=="B1") {
    lat <- 40.019
    lon <- -105.455
  } else {
    lat <- 40.0301
    lon <- -105.541
  }
  date <- as.POSIXct(dt, format= "%Y-%m-%d")
  hour <- as.numeric(format(as.POSIXct(dt), format="%H"))
  doy <- day_of_year(day=as.POSIXct(dt, format= "%Y-%m-%d"), format ="%Y-%m-%d")
  zenith <- zenith_angle(doy, 
                         lat, 
                         lon, 
                         hour)
  return(zenith)
}

# Add psi (zenith angle) to climate data
add_psi_to_climate <- function(climate_data) {
  climate_data <- climate_data %>% 
    rowwise() %>% 
    mutate(psi = get_psi(dtuse, site))
  
  return(climate_data)
}