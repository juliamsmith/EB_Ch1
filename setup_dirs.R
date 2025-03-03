# setup_dirs.R
# Create required directories
dirs <- c("R", "scripts", "slurm", "data", "output/logs", "output/results")
sapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# Move existing surface roughness data if available
if (file.exists("C1_2022_surfroughness_mod.csv")) {
  # Calculate surface roughness from existing data
  library(TrenchR)
  library(tidyverse)
  
  C1wsprofile <- read_csv("C1_2022_surfroughness_mod.csv")
  wslow <- mean(C1wsprofile$windspeed1)
  wsmed <- mean(C1wsprofile$windspeed2)
  wshi <- mean(C1wsprofile$windspeed3)
  surf <- surface_roughness(u_r=c(wslow, wsmed, wshi), zr=c(.57, .82, 1.05))
  
  # Save surface roughness for later use
  saveRDS(surf, "data/surface_roughness.rds")
}

cat("Directory structure created successfully.\n")