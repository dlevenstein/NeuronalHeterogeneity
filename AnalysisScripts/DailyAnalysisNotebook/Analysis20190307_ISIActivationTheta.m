function [ ] = Analysis20190307(basePath,figfolder)
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
    'Ybounds',[-3 2],'Xbounds',[0 2],'numYbins',150,'numXbins',30),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.thetapower,ISIStats.allspikes.instate,...
    'UniformOutput',false);
thetaXY = cat(1,thetaXY{:});
thetaXY = CollapseStruct( thetaXY,3);


for tt = 1:length(celltypes)
    thetadistbyPOP.(celltypes{tt}) = nanmean(thetaXY.pYX(:,:,CellClass.(celltypes{tt})),3);
    meanthetabyPOP.(celltypes{tt}) = nanmean(thetaXY.meanYX(:,:,CellClass.(celltypes{tt})),3);
end


%% ISI dist Low/High Theta



returnmap.lowthetathresh = 0.75;
returnmap.highthetathresh = 1.25;
returnmap.logISIbins = linspace(-3,2,60);
returnmap.allcells.hightheta = cellfun(@(n,np1,instate,theta) ...
    hist3(log10([n(instate&theta>returnmap.highthetathresh),np1(instate&theta>returnmap.highthetathresh)]),...
    {returnmap.logISIbins,returnmap.logISIbins}),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.instate,ISIStats.allspikes.thetapower,...
    'UniformOutput',false);
returnmap.allcells.hightheta = cellfun(@(X) X./sum(X(:)),...
    returnmap.allcells.hightheta,'UniformOutput',false);

returnmap.allcells.lowtheta = cellfun(@(n,np1,instate,theta) ...
    hist3(log10([n(instate&theta<returnmap.lowthetathresh),np1(instate&theta<returnmap.lowthetathresh)]),...
    {returnmap.logISIbins,returnmap.logISIbins}),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.instate,ISIStats.allspikes.thetapower,...
    'UniformOutput',false);
returnmap.allcells.lowtheta = cellfun(@(X) X./sum(X(:)),...
    returnmap.allcells.lowtheta,'UniformOutput',false);
    
for tt = 1:length(celltypes)
    returnmap.(celltypes{tt}).hightheta = nanmean(cat(3,returnmap.allcells.hightheta{CellClass.(celltypes{tt})}),3);
    returnmap.(celltypes{tt}).lowtheta = nanmean(cat(3,returnmap.allcells.lowtheta{CellClass.(celltypes{tt})}),3);
end

isidist.lowthetathresh = returnmap.lowthetathresh;
isidist.highthetathresh = returnmap.highthetathresh;
isidist.logISIbins = linspace(-3,2,100);
isidist.allcells.hightheta = cellfun(@(n,instate,theta) ...
    hist(log10(n(instate&theta>isidist.highthetathresh)),isidist.logISIbins),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.instate,...
    ISIStats.allspikes.thetapower,'UniformOutput',false);
isidist.allcells.hightheta = cellfun(@(X) X./sum(X(:)),...
    isidist.allcells.hightheta,'UniformOutput',false);

isidist.allcells.lowtheta = cellfun(@(n,instate,theta) ...
    hist(log10(n(instate&theta<isidist.lowthetathresh)),isidist.logISIbins),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.instate,...
    ISIStats.allspikes.thetapower,'UniformOutput',false);
isidist.allcells.lowtheta = cellfun(@(X) X./sum(X(:)),...
    isidist.allcells.lowtheta,'UniformOutput',false);
    
for tt = 1:length(celltypes)
    isidist.(celltypes{tt}).hightheta = nanmean(cat(3,isidist.allcells.hightheta{CellClass.(celltypes{tt})}),3);
    isidist.(celltypes{tt}).lowtheta = nanmean(cat(3,isidist.allcells.lowtheta{CellClass.(celltypes{tt})}),3);
end


%%
figure
   
for tt = 1:length(celltypes)
subplot(3,2,tt*2-1)
    imagesc(thetaXY.Xbins(1,:,1),thetaXY.Ybins(1,:,1), thetadistbyPOP.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('y',10)
    ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Power')
    title((celltypes{tt}))
    colorbar
    if tt ==1 
        caxis([0 0.0175])
    elseif tt==2
         caxis([0 0.025])
    end
    
end   


for tt=1:length(celltypes)
    subplot(3,4,4.*tt)
    imagesc(returnmap.logISIbins,returnmap.logISIbins,returnmap.(celltypes{tt}).hightheta)
    axis xy
    LogScale('xy',10)
    if tt==1
        title('High Theta')
    end
    xlabel('ISI n')
end

for tt=1:length(celltypes)
    subplot(3,4,4.*tt-1)
    imagesc(returnmap.logISIbins,returnmap.logISIbins,returnmap.(celltypes{tt}).lowtheta)
    axis xy
    LogScale('xy',10)
    if tt==1
        title('Low Theta')
    end
    xlabel('ISI n');ylabel('ISI n+1')
end

for tt=1:length(celltypes)
    subplot(3,2,4+tt)
    plot(isidist.logISIbins,isidist.(celltypes{tt}).lowtheta,...
        'color',cellcolor{tt},'linewidth',2)
    hold on
    plot(isidist.logISIbins,isidist.(celltypes{tt}).hightheta,...
        ':','color',cellcolor{tt},'linewidth',2)
    LogScale('x',10)
    legend('Low Theta','High Theta')
    xlabel('ISI (s)')
    set(gca,'ytick',[])
    box off
    ylabel('P(ISI)')
end
    
NiceSave('ISIGivenTheta',figfolder,baseName,'includeDate',true)

