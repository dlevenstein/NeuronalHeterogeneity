function [ ] = Analysis20190209(basePath,figfolder)
% Date 02/09/2019
%
%Question: What is the time occupancy of different ISIs? How does this
%relate to other conditional variables (i.e. power, synchrony) of interest?
%Define R(T) = probability that time t is during an ISI of interval of
%duration T.
%
%Plots
%-R(T) for different states
%-R(T|Synchrony)/E/I?
%
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r'};%,[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};


%% Load the LFP if needed

lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID;
downsamplefactor = 2;
lfp = bz_GetLFP(lfpchan,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%Noralize the LFP
lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

%% Interpolate containing ISI at each ms

lfp.ISI = cellfun(@(X,Y) interp1(X,Y,lfp.timestamps,'next'),...
    ISIStats.allspikes.times,ISIStats.allspikes.ISIs,'UniformOutput',false);

%% Restrict to state

state = states{2};
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);
lfp.instate = InIntervals(lfp.timestamps,SleepState.ints.(state));

%% Calculate occupancy histogram
ISIoccupancy.bins = linspace(0,20,100);
ISIoccupancy.logbins = linspace(-2.5,2.5,75);

ISIoccupancy.(state).hist = cellfun(@(X) hist(X(lfp.instate),ISIoccupancy.bins),...
    lfp.ISI,'UniformOutput',false);
ISIoccupancy.(state).hist = cellfun(@(X) X./sum(X),ISIoccupancy.(state).hist,...
    'UniformOutput',false);
ISIoccupancy.(state).hist = cat(1,ISIoccupancy.(state).hist{:});

ISIoccupancy.(state).loghist = cellfun(@(X) hist(log10(X(lfp.instate)),ISIoccupancy.logbins),...
    lfp.ISI,'UniformOutput',false);
ISIoccupancy.(state).loghist = cellfun(@(X) X./sum(X),ISIoccupancy.(state).loghist,...
    'UniformOutput',false);
ISIoccupancy.(state).loghist = cat(1,ISIoccupancy.(state).loghist{:});


%%
%cmap = [1 1 1;colormap(parula)];

figure
%colormap(cmap)
subplot(2,1,1)
    imagesc(ISIoccupancy.logbins,[1 spikes.numcells],...
        ISIoccupancy.(state).loghist(ISIStats.sorts.(state).ratebyclass,:))
    hold on
    plot(log10(1./ISIStats.summstats.(state).meanrate(ISIStats.sorts.(state).ratebyclass)),...
        [1:spikes.numcells],'.')
    LogScale('x',10)
    ColorbarWithAxis([0 0.1],'P_t(log(ISI))')
    xlabel('ISI')
subplot(2,1,2)
    imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
        ISIoccupancy.(state).hist(ISIStats.sorts.(state).ratebyclass,:))
    
NiceSave(['ISIoccupancy_',state],figfolder,baseName,'includeDate',true)
%%
excell = randsample(spikes.numcells,1);
viewwin = bz_RandomWindowInIntervals(SleepState.ints.(state),30);

figure
    subplot(4,1,1)
        plot(lfp.timestamps,log10(lfp.ISI{excell}),'k')
        axis tight
        hold on
        StateScorePlot(SleepState.ints,statecolors)
        box off 
        LogScale('y',10)
        %hold on
        %plot(ISIStats.allspikes.times{excell},zeros(size(ISIStats.allspikes.times{excell})),'+')
        %xlim(viewwin)
    subplot(4,1,2)
        plot(lfp.timestamps,log10(1./lfp.ISI{excell}),'k')
        %hold on
        axis tight
        box off
        LogScale('y',10)
        %plot(ISIStats.allspikes.times{excell},zeros(size(ISIStats.allspikes.times{excell})),'+')
        %xlim(viewwin)
    subplot(2,2,3)
        plot(lfp.timestamps,lfp.ISI{excell},'k')
        hold on
        plot(ISIStats.allspikes.times{excell},zeros(size(ISIStats.allspikes.times{excell})),'+')
        xlim(viewwin)
    subplot(4,2,6)
        bar(ISIoccupancy.bins,ISIoccupancy.(state).hist(excell,:))
        axis tight
        box off
    subplot(4,2,8)
        bar(ISIoccupancy.logbins,ISIoccupancy.(state).loghist(excell,:))
        axis tight
        box off
        LogScale('x',10)




