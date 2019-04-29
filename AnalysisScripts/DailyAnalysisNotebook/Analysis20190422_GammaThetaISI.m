function [ ] = Analysis20190422(basePath,figfolder)
% Date 03/27/2019
%
%Goal: invesitgate the relationship between theta and high gamma/ripple
%oscillation, and specifically how it effects the ISI distribution
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC/Achilles_10252013';
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
statecolors = {'k','b','r',[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};

%% The necessary inputs
states = fieldnames(SleepState.ints);
numstates = length(states);


%% Get theta power at each time/spike

 thchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID;
 downsamplefactor = 5;
 th_lfp = bz_GetLFP(thchan,...
     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%%
thetalfp = bz_Filter(th_lfp,'passband',[6 9]);
deltalfp = bz_Filter(th_lfp,'passband',[2 20]);

thetalfp.thetadelta = NormToInt(thetalfp.amp./deltalfp.amp,'mean',SleepState.ints.WAKEstate);
thetalfp.amp = NormToInt((thetalfp.amp),'mean',SleepState.ints.WAKEstate);


ISIStats.allspikes.thetapower =cellfun(@(X) ...
    interp1(thetalfp.timestamps,thetalfp.amp,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);
ISIStats.allspikes.thetarat =cellfun(@(X) ...
    interp1(thetalfp.timestamps,thetalfp.thetadelta,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);
ISIStats.allspikes.thetarat_log = cellfun(@(X) log2(X),ISIStats.allspikes.thetarat,'UniformOutput',false);


ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end); nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);

%% Restrict spikes to state
state = states{1}; %was using REM?! oops....
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);

%% Conditional (prev/next) ISI given THeta ratio
[ thetaXY ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(W);Z(W)],log10([X(W);Y(W)]),...
    'Xbounds',[-2 2],'numXbins',15,'Ybounds',[-3 2],'numYbins',150),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.thetarat_log,ISIStats.allspikes.instate,...
    'UniformOutput',false);
thetaXY = cat(1,thetaXY{:});
thetaXY = CollapseStruct( thetaXY,3);


for tt = 1:length(celltypes)
    thetadistbyPOP.(celltypes{tt}) = nanmean(thetaXY.pYX(:,:,CellClass.(celltypes{tt})),3);
    meanthetabyPOP.(celltypes{tt}) = nanmean(thetaXY.meanYX(:,:,CellClass.(celltypes{tt})),3);
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

    NiceSave('ISIbyThetaRat',figfolder,baseName,'includeDate',true)

end
