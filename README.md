# EB_Ch1 - Grasshopper Energy Budget Analysis
This repository contains code for calculating the energy budget of grasshoppers under different environmental conditions. The code is designed to run on the Hyak computing cluster, using climate data from multiple sites and years to determine how various factors affect grasshopper energy availability.
## Overview
This project calculates discretionary energy for two grasshopper species (MB and MS) across different sites at varying elevations. The energy budget simulations account for:

Different populations (site/elevation origins)
Different climate conditions (site variables)
Different years
Different thermoregulation scenarios (varying height and shade levels)

The code uses biophysical models to calculate body temperatures based on climate data, and then uses these temperatures to estimate energy gains from feeding and energy losses from metabolism.
## Directory Structure

├── R/                     # Core functions
│   ├── biophysical_functions.R  # Body temperature calculations
│   ├── energy_functions.R       # Energy budget calculations
│   ├── data_functions.R         # Data loading and saving 
│   └── setup.R                  # Package loading & setup
├── scripts/
│   ├── calculate_energy_budget.R  # Main script for calculations
│   ├── submit_jobs.R              # Script to generate batch jobs
│   └── combine_results.R          # Script to combine results
├── slurm/
│   └── energy_budget.sh           # Slurm job template
├── data/
│   ├── pops.rds                   # Population parameters
│   └── surface_roughness.rds      # Surface roughness data
└── output/
    ├── logs/                      # Slurm job logs
    └── results/                   # Calculation results

## Climate Data
The code expects climate data in a specific format, located in a directory outside this repository:
Copy../multi_microclimate/
├── A1/
│   ├── climateuse_combined_A1_2011_6_7_8.rds
│   ├── climateuse_combined_A1_2022_6_7_8.rds
│   └── ...
├── B1/
├── C1/
└── Eldo/
Each climate data file contains observations with different shade levels and heights for each time point.
## Usage
### Setup

1. Clone this repository to your Hyak directory:

git clone https://github.com/yourusername/EB_Ch1.git
cd EB_Ch1

2. Create the necessary directory structure:

Rscript setup_dirs.R

3. Prepare your population data file (data/pops.rds).

### Running Jobs
#### Option 1: Generate and submit jobs in one step
Rscript scripts/submit_jobs.R
This will create and submit jobs for all combinations of species, years, sites, heights, and shade levels.
#### Option 2: Generate job scripts without submitting
Rscript scripts/submit_jobs.R --dry-run
This will create all job scripts without submitting them, plus a submit_all_jobs.sh script that you can use to submit them later.
#### Option 3: Run a single combination
Rscript scripts/calculate_energy_budget.R MB 2022 A1 B1 0.03 0.25 output/results
This will calculate the energy budget for MB grasshoppers from A1, using climate data from B1 in 2022, at height 0.03m with 25% shade.
### Combining Results
After all jobs have completed, combine the results:
Rscript scripts/combine_results.R
This will create:

 - output/combined_results.rds - The complete dataset
 - output/results_summary.csv - A summary of key metrics by combination

## Parameters

Species: MB or MS
Sites: A1, B1, C1, Eldo
Heights: Different heights above ground (e.g., 0.03, 0.1, 0.3)
Shade levels: Amount of shade from 0 (full sun) to 0.9 (90% shade)

## Modifying Parameters
To change the set of parameters used in simulations, edit scripts/submit_jobs.R and modify the parameter lists:

species_list <- c("MB", "MS")
years <- c(2011, 2022, 2023)
site_list <- c("Eldo", "A1", "B1", "C1")
heights <- c(0.03, 0.1, 0.3) #NOTE: only with values in the .rds
shade_levels <- c(0, 0.25, 0.5, 0.75, 0.9) #NOTE: see above

## Dependencies

 - R packages: tidyverse, TrenchR, rTPC (not for long), lubridate
 - External: Apptainer/Singularity with a tidyverse container

## Notes

The calculations are computationally intensive. Plan cluster resource usage accordingly.
The total number of jobs can be large based on the combinations of parameters.
Check job logs in output/logs/ for any errors.
