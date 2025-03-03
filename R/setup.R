# R/setup.R
setup_packages <- function() {
  packages <- c("tidyverse", "TrenchR", "rTPC", "lubridate")
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}