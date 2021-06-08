#!/bin/bash
#SBATCH -p fn_medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=10
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=12G

savepath=/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/Modeling/Simulation_Data/RecurrenceJun03

module load matlab/R2018a

export SCRATCH=/gpfs/scratch/dl2820
mkdir -p $SCRATCH/$SLURM_JOB_ID

matlab -nodisplay -nodesktop -singleCompThread -r "cd('/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/ClusterCodes');addPaths;tic;iSTDPRecurrence('$savepath');toc;exit;"

rm -rf $SCRATCH/$SLURM_JOB_ID
