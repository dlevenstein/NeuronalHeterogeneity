function [] = RescoreYS(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%% Set up for parallel in cluster
% pc = parcluster('local');
% % store temporary files in the 'scratch' drive on the cluster, labeled by job ID
% pc.JobStorageLocation = strcat(getenv('SCRATCH'), '/', getenv('SLURM_JOB_ID'));
% % enable MATLAB to utilize the multiple cores allocated in the job script
% % SLURM_NTASKS_PER_NODE is a variable set in the job script by the flag --tasks-per-node
% % we use SLURM_NTASKS_PER_NODE - 1, because one of these tasks is the original MATLAB script itself
% parpool(pc, str2num(getenv('SLURM_NTASKS_PER_NODE'))-1);

%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC/Achilles_10252013';
%basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX/YMV08_170922';
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%% Rename StateScoreFigures Folder
ssffolder = fullfile(basePath,'StateScoreFigures');
ssffolder_new = fullfile(basePath,'StateScoreFigures_old');
if exist(ssffolder)
    copyfile(ssffolder,ssffolder_new)
end
%%
lfpfile = fullfile(basePath,[baseName,'.SleepScoreLFP.LFP.mat']);
if exist(lfpfile,'file')
    load(lfpfile)

    %Keep yuta's time window
    scoretime = [SleepScoreLFP.t(1) SleepScoreLFP.t(end)];
else
    scoretime = [0 Inf];
end
%%
SleepScoreMaster(basePath,'overwrite',false,'scoretime',scoretime)

