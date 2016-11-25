#!/bin/bash -l

##SBATCH --job-name = "GMM_Estim"

#SBATCH --ntasks=1
#SBATCH --time=01:00:00

#SBATCH --output=slurm_estimation.out
#SBATCH --error=slurm_estimation.err

module add matlab/R2015b

sleep $[($RANDOM % 240) + 120]

matlab -nodisplay -r Main_Boot_p0_array





