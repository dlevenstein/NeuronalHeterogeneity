#!/bin/bash
#SBATCH -p cpu_medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=48
echo $1
module load matlab/R2018a


matlab -nodisplay -nodesktop -singleCompThread -r "cd('/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/ClusterCodes');addPaths;bz_RunAnalysis('PopActivityModulationAnalysis','$1','basePath',true,'savein','basePath');exit;"

