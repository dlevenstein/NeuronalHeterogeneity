function [ ] = AnalysisXXXXXXXX(basePath,figfolder)
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
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
%spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
%CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
%ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
statenames = fieldnames(SleepState.ints);
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

% try
%     celltypes = CellClass.celltypes;
% catch
%     celltypes = unique(CellClass.label);
% end
cellcolor = {'k','r'};


%%
for ss = 1:3
    SleepState.dur.(statenames{ss}) = diff(SleepState.ints.(statenames{ss}),[],2);
end

MAthresh = 60;
SleepState.ints.MA = SleepState.ints.WAKEstate(SleepState.dur.WAKEstate<=MAthresh,:);
SleepState.ints.WAKE = SleepState.ints.WAKEstate(SleepState.dur.(statenames{ss})>MAthresh,:);

%%


[ transprob,normprob ] = IntTransitionProbabilities({SleepState.ints.WAKE,...
    SleepState.ints.(statenames{2}),...
    SleepState.ints.MA,SleepState.ints.(statenames{3})});

%%
figure
subplot(2,2,1)
imagesc(transprob)
set(gca,'ytick',[1 2 3 4]);set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabels',{'W','N','M','R'});
set(gca,'yticklabels',{'W','N','M','R'});
crameri bilbao
colorbar

subplot(2,2,2)
imagesc(normprob)
set(gca,'ytick',[1 2 3 4]);set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabels',{'W','N','M','R'});
set(gca,'yticklabels',{'W','N','M','R'});
crameri('vik','pivot',1)
colorbar
%%
figure
for ss = 1:3
    subplot(3,3,ss)
        hist(log10(SleepState.dur.(statenames{ss})),20)
end
