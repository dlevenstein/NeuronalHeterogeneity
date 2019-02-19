function [ ] = Analysis20190215(basePath,figfolder)
% Date 02/15/2019
%
%
%Question: What does the LFP power look like in a given frequency band
%when a cell spikes with a given ISI - can we show that increased power in
%some frequency bands is associated with specific (activated) ISIs? (or
%certain phase coupling?...)
%
%
%Plots
%-power distrubtion given ISI (spike pre/postceded, time during)
%-mean power (relative to median) given ISI
%-start with easy: theta, WAKE. then gamma/ripple
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = [reporoot,'Datasets/onDesktop/AG_HPC/Achilles_10252013'];
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
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};


%% Load the LFP if needed

 lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID;
 downsamplefactor = 5;
 lfp = bz_GetLFP(lfpchan,...
     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%Noralize the LFP
%lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);


%% Filter in theta band
thetalfp = bz_Filter(lfp,'passband',[6 10]);
thetalfp.amp = NormToInt((thetalfp.amp),'median',SleepState.ints.WAKEstate);
%% Get Theta Power at each spike
ISIStats.allspikes.thetapower = cellfun(@(X) ...
    interp1(thetalfp.timestamps,thetalfp.amp,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);

ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end); nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);

%% Restrict spikes to state
state = states{1}; %was using REM?! oops....
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);

%% Conditional THeta power at spike given prev/next ISI
[ CONDXY ] = cellfun(@(X,Y,Z,W) ConditionalHist( log10([X(W);Y(W)]),[Z(W);Z(W)],...
    'Xbounds',[-3 1.5],'Ybounds',[0 2.5]),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.thetapower,ISIStats.allspikes.instate,...
    'UniformOutput',false);
CONDXY = cat(1,CONDXY{:});
CONDXY = CollapseStruct( CONDXY,3);


for tt = 1:length(celltypes)
    thetadistbyPOP.(celltypes{tt}) = nanmean(CONDXY.pYX(:,:,CellClass.(celltypes{tt})),3);
    meanthetabyPOP.(celltypes{tt}) = nanmean(CONDXY.meanYX(:,:,CellClass.(celltypes{tt})),3);
end
%%
excell = randsample(spikes.numcells,1);
figure
subplot(2,2,1)
    plot(log10(ISIStats.allspikes.ISIs{excell}),ISIStats.allspikes.thetapower{excell},'k.')
    hold on
    plot(log10(ISIStats.allspikes.ISIs{excell}),ISIStats.allspikes.thetapower{excell},'k.')

    LogScale('x',10)
    xlabel('ISI (s)');ylabel('Theta Power')
subplot(2,2,2)
    imagesc(CONDXY.Xbins(1,:,1),CONDXY.Ybins(1,:,1), CONDXY.pYX(:,:,excell)')
    hold on
    plot(CONDXY.Xbins(1,:,1),CONDXY.meanYX(:,:,excell),'w')
    axis xy
    LogScale('x',10)
    xlabel('ISI (s)');ylabel('Theta Power')
%%
figure
for tt = 1:length(celltypes)
subplot(4,3,tt)
    imagesc(CONDXY.Xbins(1,:,1),CONDXY.Ybins(1,:,1), thetadistbyPOP.(celltypes{tt})')
    hold on
    plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('x',10)
    xlabel('ISI (s)');ylabel('Theta Power')
    title((celltypes{tt}))
end   
    
subplot(4,3,3)
    imagesc(CONDXY.Xbins(1,:,1),[1 spikes.numcells],squeeze(CONDXY.meanYX(:,:,ISIStats.sorts.WAKEstate.ratebyclass))')
    ColorbarWithAxis([0.5 1.5],'Theta Power (med^-^1)')   
    
NiceSave('ThetaActivation',figfolder,baseName,'includeDate',true)


%% Conditional THeta power at time point given ISI its in
%% Interpolate containing ISI at each ms

thetalfp.ISI = cellfun(@(X,Y) interp1(X,Y,thetalfp.timestamps,'next'),...
    ISIStats.allspikes.times,ISIStats.allspikes.ISIs,'UniformOutput',false);
%thetalfp.ISI = cat(2,thetalfp.ISI{:});

%%
thetalfp.instate = InIntervals(thetalfp.timestamps,double(SleepState.ints.(state)));
%%

[ CONDXY ] = cellfun(@(X) ConditionalHist( log10(X(thetalfp.instate)),thetalfp.amp(thetalfp.instate),...
    'Xbounds',[-2.5 1.5],'Ybounds',[0 2.5]),...
    thetalfp.ISI,...
    'UniformOutput',false);
CONDXY = cat(1,CONDXY{:});
CONDXY = CollapseStruct( CONDXY,3);


for tt = 1:length(celltypes)
    thetadistbyPOP.(celltypes{tt}) = nanmean(CONDXY.pYX(:,:,CellClass.(celltypes{tt})),3);
    meanthetabyPOP.(celltypes{tt}) = nanmean(CONDXY.meanYX(:,:,CellClass.(celltypes{tt})),3);
end
%%
excell = randsample(spikes.numcells,1);
figure
subplot(2,2,1)
    plot(log10(ISIStats.allspikes.ISIs{excell}),ISIStats.allspikes.thetapower{excell},'k.')
    hold on
    plot(log10(ISIStats.allspikes.ISIs{excell}),ISIStats.allspikes.thetapower{excell},'k.')

    LogScale('x',10)
    xlabel('ISI (s)');ylabel('Theta Power')
subplot(2,2,2)
    imagesc(CONDXY.Xbins(1,:,1),CONDXY.Ybins(1,:,1), CONDXY.pYX(:,:,excell)')
    hold on
    plot(CONDXY.Xbins(1,:,1),CONDXY.meanYX(:,:,excell),'w')
    axis xy
    LogScale('x',10)
    xlabel('ISI (s)');ylabel('Theta Power')
%%
figure
for tt = 1:length(celltypes)
subplot(4,3,tt)
    imagesc(CONDXY.Xbins(1,:,1),CONDXY.Ybins(1,:,1), thetadistbyPOP.(celltypes{tt})')
    hold on
    plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('x',10)
    xlabel('ISI (s)');ylabel('Theta Power')
    title((celltypes{tt}))
end   
    
subplot(4,3,3)
    s = imagesc(CONDXY.Xbins(1,:,1),[1 spikes.numcells],squeeze(CONDXY.meanYX(:,:,ISIStats.sorts.WAKEstate.ratebyclass))');
    alpha(s,single(~isnan(squeeze(CONDXY.meanYX(:,:,ISIStats.sorts.WAKEstate.ratebyclass))')))
    xlabel('ISI (s)')
    LogScale('x',10)
    ylabel('Cell, Sorted by Rate/Type')

    ColorbarWithAxis([0.5 1.5],'Theta Power (med^-^1)')   
    
NiceSave('ThetaActivationOcc',figfolder,baseName,'includeDate',true)


