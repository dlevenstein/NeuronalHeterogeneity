function [ ] = Analysis20190219(basePath,figfolder)
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
%-ISIdistirbution given theta
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


%% Filter in theta/gamma bands
thetalfp = bz_Filter(lfp,'passband',[6 10]);
thetalfp.amp = NormToInt((thetalfp.amp),'mean',SleepState.ints.WAKEstate);

gammalfp = bz_Filter(lfp,'passband',[40 100]);
gammalfp.amp = NormToInt((gammalfp.amp),'mean',SleepState.ints.WAKEstate);
%% Get Theta Power at each spike
ISIStats.allspikes.thetapower =cellfun(@(X) ...
    interp1(thetalfp.timestamps,thetalfp.amp,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);

ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end); nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);

%% Restrict spikes to state
state = states{1}; %was using REM?! oops....
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);

%% Conditional (prev/next) ISI given THeta power
[ thetaXY ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(W);Z(W)],log10([X(W);Y(W)]),...
    'Ybounds',[-3 2],'Xbounds',[0 2],'numYbins',100,'numXbins',30),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.thetapower,ISIStats.allspikes.instate,...
    'UniformOutput',false);
thetaXY = cat(1,thetaXY{:});
thetaXY = CollapseStruct( thetaXY,3);


for tt = 1:length(celltypes)
    thetadistbyPOP.(celltypes{tt}) = nanmean(thetaXY.pYX(:,:,CellClass.(celltypes{tt})),3);
    meanthetabyPOP.(celltypes{tt}) = nanmean(thetaXY.meanYX(:,:,CellClass.(celltypes{tt})),3);
end




%% Get Gamma Power at each spike
ISIStats.allspikes.gammapower =cellfun(@(X) ...
    interp1(gammalfp.timestamps,gammalfp.amp,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);

%% Conditional (prev/next) ISI given Gamma power
[ gammaXY ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(W);Z(W)],log10([X(W);Y(W)]),...
    'Ybounds',[-3 2],'Xbounds',[0 2],'numYbins',100,'numXbins',30),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.gammapower,ISIStats.allspikes.instate,...
    'UniformOutput',false);
gammaXY = cat(1,gammaXY{:});
gammaXY = CollapseStruct( gammaXY,3);


for tt = 1:length(celltypes)
    gammadistbyPOP.(celltypes{tt}) = nanmean(gammaXY.pYX(:,:,CellClass.(celltypes{tt})),3);
   % meanthetabyPOP.(celltypes{tt}) = nanmean(CONDXY.meanYX(:,:,CellClass.(celltypes{tt})),3);
end

%%
figure
for tt = 1:length(celltypes)
subplot(3,3,tt*3-1)
    imagesc(gammaXY.Xbins(1,:,1),gammaXY.Ybins(1,:,1), gammadistbyPOP.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('y',10)
    ylabel('ISI (s)');xlabel('Gamma Power')
    title((celltypes{tt}))
    colorbar
    if tt ==1 
        caxis([0 0.025])
    elseif tt==2
         caxis([0 0.035])
    end

end   
      
   
for tt = 1:length(celltypes)
subplot(3,3,tt*3-2)
    imagesc(thetaXY.Xbins(1,:,1),thetaXY.Ybins(1,:,1), thetadistbyPOP.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('y',10)
    ylabel('ISI (s)');xlabel('Theta Power')
    title((celltypes{tt}))
    colorbar
    if tt ==1 
        caxis([0 0.025])
    elseif tt==2
         caxis([0 0.035])
    end

end   
      
    
NiceSave('ISIGivenThetaGamma',figfolder,baseName,'includeDate',true)

