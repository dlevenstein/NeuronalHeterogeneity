#!/bin/bash
#SBATCH -p cpu_medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=6
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=24G

savepath=/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/Modeling/Simulation_Data/Recurrence2

module load matlab/R2018a

export SCRATCH=/gpfs/scratch/dl2820
mkdir -p $SCRATCH/$SLURM_JOB_ID

matlab -nodisplay -nodesktop -singleCompThread -r "cd('/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/ClusterCodes');addPaths;tic;iSTDPRecurrence('$savepath');toc;exit;"

rm -rf $SCRATCH/$SLURM_JOB_ID
