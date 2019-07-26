function [] = APGetCellClass(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%


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

%% Make the spikes file
SpikeDataFile = fullfile(basePath,'Analysis','SpikeData.mat');
load(SpikeDataFile);
loadcells = [shank cellIx];
%%
spikes = bz_GetSpikes('basepath',basePath,'saveMat',true,'forceReload',true,...
    'onlyLoad',loadcells);

if length(spikes.UID) ~= length(spikes.region)
    keyboard
end

%%
thalcells = spikes.UID(strcmp(spikes.region,'THAL'));

ignorecells = spikes.UID(~strcmp(spikes.region,'THAL')); %spikes.UID(strcmp(spikes.region,'SUB') | strcmp(spikes.region,'HPC'));

%%
bz_CellClassification(basePath,'keepKnown',true,...
    'knownE',thalcells,...
    'ignorecells',ignorecells,'forceReload',true)
