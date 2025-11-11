# create_pops_data.R
# Script to create the population parameters file (pops.rds) with LRF parameters
# Using MAP (Maximum A Posteriori) + HDI (Highest Density Interval) at 90%
library(tidyverse)
library(bayestestR)
library(coda)

# Load the model fits
ms_fit <- readRDS("~/GitHub/thermal_perf/ms_improved7_betterprior.rds")
mb_fit <- readRDS("~/GitHub/thermal_perf/mb_both_years_with_yeareffect3betterpriors.rds")

# ===== SET METHOD =====
method <- "MAP"
CI <- 0.9

# ===== EXTRACT PARAMETERS FROM MB MODEL =====

n_chains <- mb_fit$fit@sim$chains
mb_posteriors <- c()
for(i in 1:n_chains){
  b <- as.data.frame(as.mcmc(mb_fit)[[i]])
  b$chain <- i
  mb_posteriors <- rbind(b, mb_posteriors)
}

pops_mb <- c("A1", "B1", "C1", "D1")
mb_params <- data.frame()

for(x in pops_mb){
  Tmin_col <- parse(text = paste0('mb_posteriors$`b_Tmin_pop', x, '`'), keep.source = FALSE)[[1]]
  Topt_col <- parse(text = paste0('mb_posteriors$`b_Topt_pop', x, '`'), keep.source = FALSE)[[1]]
  Above_col <- parse(text = paste0('mb_posteriors$`b_Above_pop', x, '`'), keep.source = FALSE)[[1]]
  Ropt_col <- parse(text = paste0('mb_posteriors$`b_Ropt_pop', x, '`'), keep.source = FALSE)[[1]]
  Tmax_col <- parse(text = paste0('mb_posteriors$`b_Above_pop', x, '` + mb_posteriors$`b_Topt_pop', x, '`'), keep.source = FALSE)[[1]]
  
  mb_params <- rbind(mb_params, data.frame(
    site = x,
    Tmin = as.numeric(point_estimate(eval(Tmin_col), centrality = method)),
    Tmin_lower = as.numeric(ci(eval(Tmin_col), ci = CI, method = "HDI")[2]),
    Tmin_upper = as.numeric(ci(eval(Tmin_col), ci = CI, method = "HDI")[3]),
    Topt = as.numeric(point_estimate(eval(Topt_col), centrality = method)),
    Topt_lower = as.numeric(ci(eval(Topt_col), ci = CI, method = "HDI")[2]),
    Topt_upper = as.numeric(ci(eval(Topt_col), ci = CI, method = "HDI")[3]),
    Above = as.numeric(point_estimate(eval(Above_col), centrality = method)),
    Above_lower = as.numeric(ci(eval(Above_col), ci = CI, method = "HDI")[2]),
    Above_upper = as.numeric(ci(eval(Above_col), ci = CI, method = "HDI")[3]),
    Ropt = as.numeric(point_estimate(eval(Ropt_col), centrality = method)),
    Ropt_lower = as.numeric(ci(eval(Ropt_col), ci = CI, method = "HDI")[2]),
    Ropt_upper = as.numeric(ci(eval(Ropt_col), ci = CI, method = "HDI")[3]),
    Tmax = as.numeric(point_estimate(eval(Tmax_col), centrality = method)),
    Tmax_lower = as.numeric(ci(eval(Tmax_col), ci = CI, method = "HDI")[2]),
    Tmax_upper = as.numeric(ci(eval(Tmax_col), ci = CI, method = "HDI")[3])
  ))
}

# ===== EXTRACT PARAMETERS FROM MS MODEL =====

n_chains <- ms_fit$fit@sim$chains
ms_posteriors <- c()
for(i in 1:n_chains){
  b <- as.data.frame(as.mcmc(ms_fit)[[i]])
  b$chain <- i
  ms_posteriors <- rbind(b, ms_posteriors)
}

pops_ms <- c("Eldo", "A1", "B1")
ms_params <- data.frame()

for(x in pops_ms){
  Tmin_col <- parse(text = paste0('ms_posteriors$`b_Tmin_pop', x, '`'), keep.source = FALSE)[[1]]
  Topt_col <- parse(text = paste0('ms_posteriors$`b_Topt_pop', x, '`'), keep.source = FALSE)[[1]]
  Above_col <- parse(text = paste0('ms_posteriors$`b_Above_pop', x, '`'), keep.source = FALSE)[[1]]
  Ropt_col <- parse(text = paste0('ms_posteriors$`b_Ropt_pop', x, '`'), keep.source = FALSE)[[1]]
  Tmax_col <- parse(text = paste0('ms_posteriors$`b_Above_pop', x, '` + ms_posteriors$`b_Topt_pop', x, '`'), keep.source = FALSE)[[1]]
  
  ms_params <- rbind(ms_params, data.frame(
    site = x,
    Tmin = as.numeric(point_estimate(eval(Tmin_col), centrality = method)),
    Tmin_lower = as.numeric(ci(eval(Tmin_col), ci = CI, method = "HDI")[2]),
    Tmin_upper = as.numeric(ci(eval(Tmin_col), ci = CI, method = "HDI")[3]),
    Topt = as.numeric(point_estimate(eval(Topt_col), centrality = method)),
    Topt_lower = as.numeric(ci(eval(Topt_col), ci = CI, method = "HDI")[2]),
    Topt_upper = as.numeric(ci(eval(Topt_col), ci = CI, method = "HDI")[3]),
    Above = as.numeric(point_estimate(eval(Above_col), centrality = method)),
    Above_lower = as.numeric(ci(eval(Above_col), ci = CI, method = "HDI")[2]),
    Above_upper = as.numeric(ci(eval(Above_col), ci = CI, method = "HDI")[3]),
    Ropt = as.numeric(point_estimate(eval(Ropt_col), centrality = method)),
    Ropt_lower = as.numeric(ci(eval(Ropt_col), ci = CI, method = "HDI")[2]),
    Ropt_upper = as.numeric(ci(eval(Ropt_col), ci = CI, method = "HDI")[3]),
    Tmax = as.numeric(point_estimate(eval(Tmax_col), centrality = method)),
    Tmax_lower = as.numeric(ci(eval(Tmax_col), ci = CI, method = "HDI")[2]),
    Tmax_upper = as.numeric(ci(eval(Tmax_col), ci = CI, method = "HDI")[3])
  ))
}

# ===== CREATE POPS DATAFRAME =====

# Create population data frame with base information
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

# Merge LRF parameters into pops dataframe
pops <- pops %>%
  left_join(
    bind_rows(
      mb_params %>% mutate(spp = "MB"),
      ms_params %>% mutate(spp = "MS")
    ),
    by = c("spp", "site")
  )

# ===== SAVE FILES =====

# Create directory if it doesn't exist
dir.create("data", showWarnings = FALSE)

# Save to RDS file
saveRDS(pops, "data/pops.rds")

# Also save as CSV for easier inspection
write_csv(pops, "data/pops.csv")

cat("Population data created and saved to data/pops.rds and data/pops.csv\n")
cat(sprintf("Method used: MAP + HDI with %d%% credible intervals\n", CI * 100))

# Print a summary of the data
cat("\nPopulation data summary:\n")
cat(sprintf("Number of populations: %d\n", nrow(pops)))
cat(sprintf("Species: %s\n", paste(unique(pops$spp), collapse=", ")))
cat(sprintf("Sites: %s\n", paste(unique(pops$site), collapse=", ")))

# Show sample of the TPC parameters
cat("\nSample of TPC parameters (first 3 rows):\n")
print(pops[1:3, c("spp", "site", "Tmin", "Topt", "Above", "Ropt", "Tmax", 
                  "Tmin_lower", "Tmin_upper", "Tmax_lower", "Tmax_upper")])

# Also print the full parameter table for each species
cat("\n=== MB TPC Parameters (MAP + 90% HDI) ===\n")
print(mb_params)

cat("\n=== MS TPC Parameters (MAP + 90% HDI) ===\n")
print(ms_params)