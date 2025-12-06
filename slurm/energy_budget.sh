#!/bin/bash
#SBATCH --job-name=eb_{JOB_NAME}
#SBATCH --output=output/logs/%j.out
#SBATCH --error=output/logs/%j.err
#SBATCH --time=1:00:00
#SBATCH -p compute
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Set working directory
WORKDIR=/mmfs1/gscratch/biology/jmsmith/EB_Ch1
cd $WORKDIR

module load apptainer/1.1.5
apptainer exec \
  --bind $WORKDIR:$WORKDIR \
  --bind /mmfs1/gscratch/biology/jmsmith/targeted_microclimate:/mmfs1/gscratch/biology/jmsmith/targeted_microclimate3 \
  tidyverse_latest.sif \
  Rscript $WORKDIR/scripts/calculate_energy_budget.R {SPECIES} {YEAR} {SITE_ORIG} {SITE_CLIM} {SEX} {OUTPUT_DIR}
