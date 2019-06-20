function [] = YSGetCellClass(basePath,figfolder)
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

%% Load the old format
load(fullfile(basePath,[baseName,'_UnitFeature.mat']))
load(fullfile(basePath,['depthsort_parameter_1.mat']))

spikes = bz_GetSpikes('basepath',basePath);

%Make sure there's the same number of cells
checkcells = length(spikes.UID) == length(UnitFeature.Depth);
if ~checkcells
   display('Wrong Number of Cells in UnitFeature.mat and spikes.cellinfo.mat!') 
   keyboard
   return
end
%%
[~,depthsort] = sort(depth);
newUID = spikes.UID(depthsort);

%%
figure
subplot(2,2,1)
plot(UnitFeature.SpkWaveFormNormMax)
subplot(2,2,2)
plot(UnitFeature.SpkWaveFormNormMax)
%%
%sum(~UnitFeature.IDwaveE & ~UnitFeature.IDwaveI & ~UnitFeature.IDposwav)
%%
bz_CellClassification(basePath,'keepKnown',true,...
    'knownE',newUID(UnitFeature.IDwaveE),'knownI',newUID(UnitFeature.IDwaveI),...
    'ignorecells',newUID(UnitFeature.IDposwav),'forceReload',true)
