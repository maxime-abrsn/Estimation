#!/bin/bash -l

#SBATCH -J Main_Boot_p0_array

#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --error=Main_Boot_p0_array.err
#SBATCH --output=Main_Boot_p0_array.out

module add matlab/R2016b

matlab -nodisplay -nosplash -nodesktop < Main_Boot_p0_array.m



