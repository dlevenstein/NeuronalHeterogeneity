function [ISIStats,CellClass] = SpikeStatsAnalysis_Load(basePath,figfolder)
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
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
%basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
%baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
%sessionInfo = bz_getSessionInfo(basePath);
%spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
%SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
%states = fieldnames(SleepState.ints);
% states{4} = 'ALL';
% SleepState.ints.ALL = [0 Inf];
% statecolors = {'k','b','r',[0.6 0.6 0.6]};
% 
% try
%     celltypes = CellClass.celltypes;
% catch
%     celltypes = unique(CellClass.label);
% end
% cellcolor = {'k','r'};
if isfield(CellClass,'Waveforms')
    CellClass = rmfield(CellClass,'Waveforms');
end
ISIStats = rmfield(ISIStats,'allspikes');
