#!/bin/bash
#SBATCH -p cpu_medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=41
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=6G
echo $1
echo $2
echo $3
module load matlab/R2018a

export SCRATCH=/gpfs/scratch/dl2820
mkdir -p $SCRATCH/$SLURM_JOB_ID


matlab -nodisplay -nodesktop -singleCompThread -r "cd('/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/ClusterCodes');addPaths;MultiRecordingGammaModes($2,'saveName',$1,'region',$3,'saveFolder','/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/SharedGammaModeFitAnalysis','clusterpar',true);exit;"

rm -rf $SCRATCH/$SLURM_JOB_ID
