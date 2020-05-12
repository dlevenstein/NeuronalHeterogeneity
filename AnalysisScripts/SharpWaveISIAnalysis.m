function [PeriSWISIDist_next,PeriSWISIDist ] = SharpWaveISIAnalysis(basePath,figfolder)
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
%reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
% basePath = '/Users/dl2820/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
% %basePath = pwd;
% %basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
%ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
%states = fieldnames(SleepState.ints);
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
%statecolors = {'k','b','r',[0.6 0.6 0.6]};

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};

%%
SharpWaves = bz_LoadEvents(basePath,'SWR');
%FindbestRippleChannel
%bz_GetBestRippleChan

%%

%ripples = bz_FindRipples(basePath,rpchan)

%%
eventimes = SharpWaves.peaktimes;


%%
[PeriSWISIDist_next] = bz_PeriEventISIDist(spikes.times,eventimes,...
    'numXbins',160,'minX',40,'whichISIs','next',...
    'cellclass','load','basePath',basePath);

[PeriSWISIDist] = bz_PeriEventISIDist(spikes.times,eventimes,...
    'numXbins',160,'minX',40,'whichISIs','both',...
    'cellclass','load','basePath',basePath);
%%
figure
for tt = 1:length(celltypes)
subplot(2,2,tt)
    imagesc(PeriSWISIDist.pop.(celltypes{tt}).Xbins,...
        PeriSWISIDist.pop.(celltypes{tt}).Ybins,PeriSWISIDist.pop.(celltypes{tt}).pYX')
    hold on
    plot(PeriSWISIDist.pop.(celltypes{tt}).Xbins,...
        log10(1./PeriSWISIDist.pop.(celltypes{tt}).rate),cellcolor{tt},'linewidth',1)
    axis tight
    plot([0 0],ylim(gca),'w--')
    LogScale('y',10,'nohalf',true)
    xlabel('t (s) - relative to SW');ylabel('ISI (s)')
    title(celltypes{tt})
    
    
subplot(2,2,tt+2)
    imagesc(PeriSWISIDist_next.pop.(celltypes{tt}).Xbins,...
        PeriSWISIDist_next.pop.(celltypes{tt}).Ybins,PeriSWISIDist_next.pop.(celltypes{tt}).pYX')
    hold on
    plot(PeriSWISIDist_next.pop.(celltypes{tt}).Xbins,...
        log10(1./PeriSWISIDist_next.pop.(celltypes{tt}).rate),cellcolor{tt},'linewidth',1)
    axis tight
    plot([0 0],ylim(gca),'w--')
    LogScale('y',10,'nohalf',true)
    xlabel('t (s) - relative to SW');ylabel('ISI (s)')
    title(celltypes{tt})
end
NiceSave('PeriSWISI_next',figfolder,baseName)


