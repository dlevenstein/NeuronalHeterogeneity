function [PeriSWISIDist_next,PeriSWISIDist,SW_ISIstats ] = SharpWaveISIAnalysis(basePath,figfolder)
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

%% ISI dist/return map: in and out of SWR
SWRints.SWR = SharpWaves.times;
SWRints.iSWR = [SharpWaves.times(1:end-1,2) SharpWaves.times(2:end,1)];
SWRints.iSWR = RestrictInts(SWRints.iSWR,SleepState.ints.NREMstate);
SWRints.SWR = RestrictInts(SWRints.SWR,SleepState.ints.NREMstate);
SW_ISIstats = bz_ISIStats(spikes,'ints',SWRints,'showfig',true,'cellclass',CellClass.label);
SW_ISIstats = rmfield(SW_ISIstats,'allspikes');


%%
swrlabels = {'SWR','iSWR'};
%%
figure
subplot(2,2,1)
    hold on
    for tt = 1:length(celltypes)
        plot(log10(SW_ISIstats.summstats.SWR.meanrate(CellClass.(celltypes{tt}))),...
            log10(SW_ISIstats.summstats.iSWR.meanrate(CellClass.(celltypes{tt}))),...
            '.','color',cellcolor{tt})
    end
    hold on
    UnityLine
    xlabel('SWR Rate');ylabel('iSWR Rate')
    

for tt = 1:length(celltypes) 
    
    subplot(4,4,tt+6)
        plot(SW_ISIstats.ISIhist.logbins,SW_ISIstats.meandists.SWR.(celltypes{tt}).ISIdist,'k')
        hold on
        plot(SW_ISIstats.ISIhist.logbins,SW_ISIstats.meandists.iSWR.(celltypes{tt}).ISIdist,'r')
        LogScale('x',10,'exp',true)
        title(celltypes{tt})
        box off 
        
    for ss = 1:length(swrlabels)
        subplot(4,4,10+tt+(ss-1)*4)
            imagesc(SW_ISIstats.ISIhist.logbins,SW_ISIstats.ISIhist.logbins,...
            SW_ISIstats.meandists.(swrlabels{ss}).(celltypes{tt}).Return)
            LogScale('xy',10,'exp',true)
            axis xy
    end
end

NiceSave('SW_ISIStats',figfolder,baseName)


%%
[PeriSWISIDist_next] = bz_PeriEventISIDist(spikes.times,eventimes,...
    'numXbins',160,'minX',40,'whichISIs','next','winsize',[-0.5 0.5],...
    'cellclass','load','basePath',basePath);

[PeriSWISIDist] = bz_PeriEventISIDist(spikes.times,eventimes,...
    'numXbins',160,'minX',40,'whichISIs','both','winsize',[-0.5 0.5],...
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
NiceSave('PeriSWISI',figfolder,baseName)


