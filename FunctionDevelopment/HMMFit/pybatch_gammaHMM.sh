#!/bin/bash

#SBATCH --time=2:00:00       # walltime
#SBATCH --partition=cpu_dev
#SBATCH -J "gammaHMM"       # job name
#SBATCH --mem-per-cpu=5GB
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --array=1-10      # job array of input size


file=$1

module load python/gcc/3.6.5
# Slurm array must be max number of pairs per session
# If there are fewer pairs in session, these jobs will
# finish gracefully within MATLAB
srun python fitTransMat.py '$file' $SLURM_ARRAY_TASK_ID