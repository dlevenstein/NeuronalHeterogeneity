#!/bin/bash
#SBATCH -p fn_medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=13
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=24G
echo $1
module load matlab/R2018a

export SCRATCH=/gpfs/scratch/dl2820
mkdir -p $SCRATCH/$SLURM_JOB_ID


matlab -nodisplay -nodesktop -singleCompThread -r "cd('/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/ClusterCodes');addPaths;bz_RunAnalysis('ISIModeLFPAnalysis_vCTX','$1','basePath',true);exit;"


rm -rf $SCRATCH/$SLURM_JOB_ID
