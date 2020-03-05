function [ transprob,normprob,durhist] = BehaviorTransitionsAnalysis(basePath,figfolder)
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
%basePath = pwd;
%basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
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
totnumstates = sum([size(SleepState.ints.(statenames{1}),1),size(SleepState.ints.(statenames{2}),1),...
    size(SleepState.ints.(statenames{3}),1)]);
durhist.bins = linspace(0.5,4,25);
for ss = 1:3
    SleepState.dur.(statenames{ss}) = diff(SleepState.ints.(statenames{ss}),[],2);
    durhist.(statenames{ss}) = hist(log10(SleepState.dur.(statenames{ss})),durhist.bins);
    durhist.(statenames{ss}) = durhist.(statenames{ss})./totnumstates;
end

MAthresh = 180;
SleepState.ints.MA = SleepState.ints.WAKEstate(SleepState.dur.WAKEstate<=MAthresh,:);
SleepState.ints.WAKE = SleepState.ints.WAKEstate(SleepState.dur.WAKEstate>MAthresh,:);

%%


[ transprob,normprob ] = IntTransitionProbabilities({SleepState.ints.WAKE,...
    SleepState.ints.(statenames{2}),...
    SleepState.ints.(statenames{3}),SleepState.ints.MA});

%%
figure
subplot(2,2,1)
imagesc(transprob)
set(gca,'ytick',[1 2 3 4]);set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabels',{'W','N','R','M'});
set(gca,'yticklabels',{'W','N','R','M'});
crameri bilbao
colorbar

subplot(2,2,2)
imagesc(normprob)
set(gca,'ytick',[1 2 3 4]);set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabels',{'W','N','R','M'});
set(gca,'yticklabels',{'W','N','R','M'});
crameri('vik','pivot',1)
colorbar

subplot(2,2,3)
hold on
for ss = 1:3
    plot(durhist.bins,durhist.(statenames{ss}),statecolors{ss})

end
        hold on
        axis tight
        plot(log10(MAthresh).*[1 1],ylim(gca),'r--')
        title(statenames{ss})
        LogScale('x',10)
NiceSave('BehaviorTransition',figfolder,baseName);
