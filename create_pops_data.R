# create_pops_data.R
# Script to create the population parameters file (pops.rds) with LRF parameters

library(tidyverse)

# Create population data frame with LRF parameters
pops <- data.frame(
  spp = c(rep("MB", 8), 
          rep("MS", 6)),
  site = c(rep("A1", 2), 
           rep("B1", 2),
           rep("C1", 2),
           rep("D1", 2),
           rep("Eldo", 2),
           rep("A1", 2), 
           rep("B1", 2)),
  elev = c(rep(2195, 2), rep(2591, 2), rep(3014, 2), rep(3300, 2), # Used estimated elevation for D1
           rep(1740, 2), rep(2195, 2), rep(2591, 2)),
  sex = rep(c("M", "F"), 7),
  
  # Mass values (in g) - adjust these based on your data
  mass = c(0.33, 0.75,   # MB A1 M, F
           0.28, 0.62,   # MB B1 M, F
           0.24, 0.52,   # MB C1 M, F
           0.22, 0.45,   # MB D1 M, F (estimated)
           0.31, 0.45,   # MS Eldo M, F
           0.33, 0.46,   # MS A1 M, F
           0.29, 0.38),  # MS B1 M, F
  
  # Metabolic rate parameters
  rmr_b0 = c(rep(17.1, 8), # MB 
             rep(16.5, 6)), # MS
  rmr_b1 = c(rep(0.98, 8), 
             rep(0.85, 6)),
  rmr_b2 = c(rep(-0.48, 8), 
             rep(-0.47, 6)),
  rmr_b3 = c(rep(1.24*10^(-4), 8), 
             rep(1.08*10^(-4), 6))
)

# Add LRF parameters based on your model results
# MB parameters
mb_params <- data.frame(
  site = c("A1", "B1", "C1", "D1"),
  Tmin = c(13.57, 12.80, 13.79, 14.0),   # D1 is estimated
  Topt = c(37.54, 41.76, 41.03, 42.0),   # D1 is estimated
  Above = c(11.70, 10.74, 11.77, 12.0),  # D1 is estimated
  Ropt = c(2.64, 3.29, 3.93, 4.0)        # D1 is estimated
)

# MS parameters
ms_params <- data.frame(
  site = c("Eldo", "A1", "B1"),
  Tmin = c(14.74, 11.67, 10.96),
  Topt = c(39.88, 39.89, 39.78),
  Above = c(9.56, 6.68, 8.20),
  Ropt = c(3.50, 3.63, 4.30)
)

# Merge LRF parameters into pops dataframe
pops <- pops %>%
  left_join(
    bind_rows(
      mb_params %>% mutate(spp = "MB"),
      ms_params %>% mutate(spp = "MS")
    ),
    by = c("spp", "site")
  )

# Create directory if it doesn't exist
dir.create("data", showWarnings = FALSE)

# Save to RDS file
saveRDS(pops, "data/pops.rds")

# Also save as CSV for easier inspection
write_csv(pops, "data/pops.csv")

cat("Population data created and saved to data/pops.rds and data/pops.csv\n")

# Print a summary of the data
cat("\nPopulation data summary:\n")
cat(sprintf("Number of populations: %d\n", nrow(pops)))
cat(sprintf("Species: %s\n", paste(unique(pops$spp), collapse=", ")))
cat(sprintf("Sites: %s\n", paste(unique(pops$site), collapse=", ")))