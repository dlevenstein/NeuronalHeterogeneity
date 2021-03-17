#!/bin/bash

#SBATCH --time=4:00:00       # walltime
#SBATCH --partition=cpu_short
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
srun python fitGammaHMM.py $file $SLURM_ARRAY_TASK_ID
