#!/bin/bash
#SBATCH --job-name=eb_{JOB_NAME}
#SBATCH --output=output/logs/%j.out
#SBATCH --error=output/logs/%j.err
#SBATCH --time=6:00:00
#SBATCH -p ckpt
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Set working directory
WORKDIR=/mmfs1/gscratch/biology/jmsmith/EB_Ch1
cd $WORKDIR

module load apptainer/1.1.5
apptainer exec \
  --bind $WORKDIR:$WORKDIR \
  --bind /gscratch/biology/jmsmith/R:/gscratch/biology/jmsmith/R \
  --bind /mmfs1/gscratch/biology/jmsmith/targeted_microclimate:/mmfs1/gscratch/biology/jmsmith/targeted_microclimate \
  tidyverse_latest.sif \
  Rscript $WORKDIR/scripts/calculate_energy_bundle.R {SPECIES} "{YEARS}" {SITE_ORIG} {SITE_CLIM} {SEX} {OUTPUT_DIR}
