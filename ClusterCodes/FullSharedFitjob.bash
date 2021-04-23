#!/bin/bash
#SBATCH -p cpu_medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=48G
echo $1
module load matlab/R2018a


matlab -nodisplay -nodesktop -singleCompThread -r "cd('/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/ClusterCodes');addPaths;MultiRecordingGammaModes('/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/SharedGammaModeFitAnalysis','saveName','CA1');exit;"
