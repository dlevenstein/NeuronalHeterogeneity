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
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
basePath = '/Users/dl2820/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
%basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};
%%
ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end);nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);

%% Parms
winsize = [-1 1]; %window around event time

%%
SlowWaves = bz_LoadEvents(basePath,'SlowWaves');
eventimes = SlowWaves.timestamps;


%%
[PeriSWISIDist_next] = bz_PeriEventISIDist(spikes.times,eventimes,...
    'numXbins',80,'minX',40,'whichISIs','next',...
    'cellclass','load','basePath',basePath);

[PeriSWISIDist] = bz_PeriEventISIDist(spikes.times,eventimes,...
    'numXbins',80,'minX',40,'whichISIs','both',...
    'cellclass','load','basePath',basePath);
%%
figure
for tt = 1:length(celltypes)
subplot(2,2,tt)
imagesc(PeriSWISIDist.pop.(celltypes{tt}).Xbins,...
    PeriSWISIDist.pop.(celltypes{tt}).Ybins,PeriSWISIDist.pop.(celltypes{tt}).pYX')
hold on
plot(PeriSWISIDist.pop.(celltypes{tt}).Xbins,...
    log10(1./PeriSWISIDist.pop.(celltypes{tt}).rate),cellcolor{tt})
plot([0 0],ylim(gca),'w--')
LogScale('y',10)
end
%%
excell = 5;
spiketimes = ISIStats.allspikes.times{excell};
ISIs = ISIStats.allspikes.ISIs{excell};
ISInp1 = ISIStats.allspikes.ISInp1{excell};

%%
ISIs_rel = [];
ISInp1_rel = [];
spiketimes_rel = [];
for ee = 1:length(eventimes)
    reltime = spiketimes-eventimes(ee);
    inwin = reltime >= winsize(1) & reltime <= winsize(2);
    
    ISIs_rel = [ISIs_rel;ISIs(inwin)];
    ISInp1_rel = [ISInp1_rel;ISInp1(inwin)];
    spiketimes_rel = [spiketimes_rel;reltime(inwin)];
end

%%
[ PeriSWISIDist ] = ConditionalHist( [spiketimes_rel;spiketimes_rel],log10([ISIs_rel;ISInp1_rel]),...
    'Xbounds',winsize,'numXbins',50,'Ybounds',[-3 2],'numYbins',125,'minX',50);

PeriSWISIDist.rate = PeriSWISIDist.Xhist./(diff(PeriSWISIDist.Xbins([1 2])).*length(eventimes));
%%
figure
imagesc(PeriSWISIDist.Xbins,(PeriSWISIDist.Ybins),PeriSWISIDist.pYX')
hold on
plot(PeriSWISIDist.Xbins,log10(1./PeriSWISIDist.rate),'r')
LogScale('y',10)
%% For each spike, get it's event-relative time 
%(note, this is bad for spikes that are near multiple events)
ISIStats.allspikes.eventtime = cellfun(@(X) interp1(headdir.timestamps,headdir.data,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);

%%
[ PeriSWISIDist ] = cellfun(@(X,Y,Z) ConditionalHist( [Z;Z],log10([X;Y]),...
    'Xbounds',winsize,'numXbins',20,'Ybounds',[-3 2],'numYbins',125,'minX',25),...%,'conditionby',headdir.data(headdir.instate)),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.eventtime,...
    'UniformOutput',false);

PeriSWISIDist = cat(1,PeriSWISIDist{:});
PeriSWISIDist = CollapseStruct( PeriSWISIDist,3);
%ISIbyPos.pYX = ISIbyPos.pYX./mode(diff(T));
PeriSWISIDist.rate = sum(PeriSWISIDist.pYX,2);


