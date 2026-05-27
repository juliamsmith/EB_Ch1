# EB_Ch1 - The code associated with the "How do climate and population differences affect grasshopper energy gain across an elevational gradient?" manuscript

## Project overview ##

We measure the digestion thermal performance curve and estimate the energy budgets for populations of two grasshopper species -- *Melanoplus sanguinipes* and *Melanoplus boulderensis* -- along an elevation gradient in Boulder County, Colorado. We want to know:

1. Do the digestion thermal performance curves of different populations differ?
2. How do physiology and climate affect the discretionary energy of the populations? 
3. How have energy budgets shifted from the 1950s to today?


## Guide to folders ##

**climate_and_microclimate_scripts/** contains the script used to download and merge ERA5 and ERA5-Land data and the script used to get microclimate data using mcera5's micro_era5() function (due to package updates and such the latter may not be quite reproducible).

**data/** contains data about the different populations (a .rds created by create_pops_data.R) and our data on the thermal dependence of digestion (ad.csv).

**R/, scripts/, and slurm/** are used to run the energy budget for different populations and conditions on a computing cluster (UW's hyak).

**output_old/** contains the outputs of the energy budget runs, including climate variables and grasshopper position, as well as energy gains, losses, and net gains. 

**analysis/** contains the scripts that produce the figures, tables, and analyses for the manuscript. create_data_summary.R summarizes that data in output_old/, creating the csv that many of the other scripts work with. tpc_brms_fits/ contains the brms objects (they were fit using the ad.csv data). Most other scripts are as labelled. figs3n4.R produces figures 3 and 4 but can be modified to produce the figures in the appendix.

## Guide to thermal performance curve (TPC) fitting ##

Following the methods of Siemers et al. 2024 (https://doi.org/10.1002/ecs2.4842), we fit our data to the using the Lobry-Rosso-Flandrois equation with a lognormal error distribution. We used the brms package. We let the parameter fits (e.g., Topt) vary by population, we included a fixed effect (essentially on the height/magnitude of the curve but not the shape) for sex and a random effect for individual, and (in the case of *M. boulderensis*) a fixed effect for year of data collection. Our cleaned data is data/ad.csv and the brms fits for each species are analysis/tpc_brms_fits/mb_both_years_with_yeareffect3betterpriors.rds and analysis/tpc_brms_fits/ms_improved7_betterprior.rds. You can examine the fits of these brms objects (df_name will give you a summary of the brms call and fit, df_name$prior will tell you about the priors, etc.).

## Guide to running the energy budget ##

This code calculates discretionary energy for two grasshopper species (*M. boulderensis* = MB and *M. sanguinipes* = MS) across different sites at varying elevations. The energy budget simulations account for:
1) Different populations (site/elevation origins)
2) Different climate conditions (site variables)
3) Different years
4) Different thermoregulation scenarios (varying height and shade levels)
5) Different sexes (though they are ultimately combined in the analysis)

The code uses biophysical models to calculate body temperatures based on climate data, and then uses these temperatures to estimate energy gains from feeding and energy losses from metabolism.

Here is the approach that worked for me on UW hyak. It involved using apptainer to have access to tidyverse and other R packages on the computing cluster.

**Step 1:** Generate .sh scripts for different combinations

In the terminal:
```
apptainer exec --bind $PWD:$PWD --bind /gscratch/biology/jmsmith/R:/gscratch/biology/jmsmith/R --bind /mmfs1/gscratch/biology/jmsmith/targeted_microclimate:/mmfs1/gscratch/biology/jmsmith/targeted_microclimate tidyverse_latest.sif Rscript scripts/submit_bundled_jobs.R both --dry-run
```

I'll print one for you to show you what they look like:

here is bundle_job_MB_recent_D1_D1_M.sh:
```
#!/bin/bash
#SBATCH --job-name=eb_MB_recent_D1_D1_M
#SBATCH --output=output/logs/%j.out
#SBATCH --error=output/logs/%j.err
#SBATCH --time=00:10:00
#SBATCH -p ckpt
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#Set working directory
WORKDIR=/mmfs1/gscratch/biology/jmsmith/EB_Ch1
cd $WORKDIR
module load apptainer/1.1.5
apptainer exec
--bind $WORKDIR:$WORKDIR
--bind /gscratch/biology/jmsmith/R:/gscratch/biology/jmsmith/R
--bind /mmfs1/gscratch/biology/jmsmith/targeted_microclimate:/mmfs1/gscratch/biology/jmsmith/targeted_microclimate
tidyverse_latest.sif
Rscript $WORKDIR/scripts/calculate_energy_bundle.R MB "2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024" D1 D1 M output/results/D1
```
**Step 2**: Submit jobs to the cluster (not really important if you are not using the cluster):
This terminal call will run all of your .sh scripts that you just generated:
```
./submit_all_MS_MB_bundles.sh
```

**Step 3** (optional):
Run the R script R/add_pbt_thermoregulation.R to calculate the energy gains and losses for a grasshopper thermoregulating to get as close to its preferred body temperature as possible (ultimately allowing us to generate Figures A3 and A4). This script adds rows to the data in output_old/results. 
```
Rscript.exe R/add_pbt_thermoregulation.R output_old/results/ 
```

One more note is that the microclimate data that is this code uses to get Tbs exists in a folder called targeted_microclimate just outside of this repo (because these are large files). You can see how it was generated by looking at the scripts in the climate_and_microclimate_scripts/ folder.

You can find the results of these runs in output_old/results. I also summarize them (just total discretionary energy for the season rather than each timestep) analysis/eb_results_summary_wide.csv in (generated by analysis/create_data_summary.R).

## Guide to figures and tables ##

The code needed to run to generate figures, tables, and analyses is in the analysis/ folder. 
First run create_data_summary.R to produce eb_results_summary_wide_old.R (and to produce data for Figures A1 and A2 which show energy budget outputs with a normalized area under the TPC, run create_data_summary_normalize.R to produce eb_results_summary_wide_normalized_old.R), which is needed for the code that visualizes/analyses energy budget results. Run fig1.R, fig2.R, and fig3n4.R to produce the figures in the main text. Modifications to fig3n4.R can produce supplementary figures A1-A4. Run make_tpc_tables.R to produce tables A3 and A4. statsandtables.R produces tables and statistics related to the energy budgets, namely tables A5-A9 and the hierarchical partitioning results referenced in the manuscript.

