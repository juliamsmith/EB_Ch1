# EB_Ch1 - The code associated with the "How do climate and population differences affect grasshopper energy gain across an elevational gradient?" manuscript

**climate_and_microclimate_scripts/** contains the script used to download and merge ERA5 and ERA5-Land data and the script used to get microclimate data using mcera5's micro_era5() function (due to package updates and such the latter may not be quite reproducible).

**data/** contains data about the different populations (a .rds created by create_pops_data.R) and our data on the thermal dependence of digestion (ad.csv).

**R/, scripts/, and slurm/** are used to run the energy budget for different populations and conditions on a computing cluster (UW's hyak).

**output_old/** contains the outputs of the energy budget runs, including climate variables and grasshopper position, as well as energy gains, losses, and net gains. 

**analysis/** contains the scripts that produce the figures, tables, and analyses for the manuscript. create_data_summary.R summarizes that data in output_old/, creating the csv that many of the other scripts work with. tpc_brms_fits/ contains the brms objects (they were fit using the ad.csv data). Most other scripts are as labelled. figs3n4.R produces figures 3 and 4 but can be modified to produce the figures in the appendix.
