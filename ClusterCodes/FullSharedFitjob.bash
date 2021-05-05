#!/bin/bash
#SBATCH -p cpu_long
#SBATCH --nodes=1
#SBATCH --tasks-per-node=36
#SBATCH --time=15-00:00:00
#SBATCH --mem-per-cpu=6G
echo $1
module load matlab/R2018a

export SCRATCH=/gpfs/scratch/dl2820
mkdir -p $SCRATCH/$SLURM_JOB_ID


matlab -nodisplay -nodesktop -singleCompThread -r "cd('/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/ClusterCodes');addPaths;MultiRecordingGammaModes('/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/SharedGammaModeFitAnalysis','saveName','CA1','clusterpar',true);exit;"

rm -rf $SCRATCH/$SLURM_JOB_ID
