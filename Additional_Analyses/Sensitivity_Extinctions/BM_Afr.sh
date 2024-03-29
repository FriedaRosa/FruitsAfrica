#!/bin/bash
 
#SBATCH --job-name=BM_Afr
#SBATCH --chdir=/work/woelke
#SBATCH --output=/work/%u/%x-%A-%a.log
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=FAIL,END
 
module load foss/2019b R/4.0.0-2

output_dir=/work/$USER/$SLURM_JOB_NAME
mkdir -p "$output_dir"


Rscript --vanilla $HOME/BM_Afr_Sensitivity.r \
  $HOME/tree$SLURM_ARRAY_TASK_ID.nex \
  "$output_dir"