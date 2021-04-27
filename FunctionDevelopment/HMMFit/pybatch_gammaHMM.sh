#!/bin/bash

#SBATCH --time=12:00:00       # walltime, >24 for full fit
#SBATCH --partition=cpu_short  #_medium for full fit
#SBATCH -J "gammaHMM"       # job name
#SBATCH --mem-per-cpu=5GB
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --array=1-250      # job array of input size


file=$1

module load python/cpu/3.6.5
# Slurm array must be max number of pairs per session
# If there are fewer pairs in session, these jobs will
# finish gracefully within MATLAB
srun python fitTransMat.py $file $SLURM_ARRAY_TASK_ID
