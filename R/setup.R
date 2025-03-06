# R/setup.R
setup_packages <- function() {
  # Try to create a user library if it doesn't exist
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (user_lib == "") {
    user_lib <- "~/R/library"
  }
  
  if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Add the library path
  .libPaths(c(user_lib, .libPaths()))
  
  # Also try to use the specified gscratch library if available
  gscratch_lib <- "/gscratch/biology/jmsmith/R"
  if (dir.exists(gscratch_lib)) {
    .libPaths(c(gscratch_lib, .libPaths()))
  }
  
  # Print the library paths for debugging
  cat("Library paths:\n")
  print(.libPaths())
  
  # Check if packages are installed, and load them
  packages <- c("tidyverse", "TrenchR", "rTPC", "lubridate")
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("Package %s not found, attempting to install...\n", pkg))
      install.packages(pkg, lib = user_lib)
      library(pkg, character.only = TRUE)
    } else {
      cat(sprintf("Package %s loaded successfully\n", pkg))
    }
  }
}
