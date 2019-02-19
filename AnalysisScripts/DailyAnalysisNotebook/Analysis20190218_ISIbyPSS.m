function [ ] = Analysis20190218(basePath,figfolder)
% Date 02/18/2019
%
%Question: ISI distribution conditioned on PSS - is it a sharp transition?
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
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

%% Load the PSS
specslope = bz_PowerSpectrumSlope([],[],[],'showfig',true,...
    'saveMat',basePath);

%%
numbins = 40;
PSShist.bins = linspace(-2,0,numbins);
for ss = 1:length(states)
    specslope.timeidx.(states{ss}) = InIntervals(specslope.timestamps,SleepState.ints.(states{ss}));
    
    PSShist.(states{ss}) = hist(specslope.data(specslope.timeidx.(states{ss})),PSShist.bins);
end

%% Get PSS at each spike

ISIStats.allspikes.PSS = cellfun(@(X) interp1(specslope.timestamps,specslope.data,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);



%% Conditional prev/next ISI at each spike given PSS

% Next ISI each spike
ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end); nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);


[ CONDXY ] = cellfun(@(X,Y,Z) ConditionalHist( [Z;Z],log10([X;Y]),...
    'Xbounds',[-1.5 -0.4],'Ybounds',[-3 1.5],'minX',25),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.PSS,...
    'UniformOutput',false);
CONDXY = cat(1,CONDXY{:});
CONDXY = CollapseStruct( CONDXY,3);


for tt = 1:length(celltypes)
    ISIdistbyPSS.(celltypes{tt}) = nanmean(CONDXY.pYX(:,:,CellClass.(celltypes{tt})),3);
    %meanthetabyPOP.(celltypes{tt}) = nanmean(CONDXY.meanYX(:,:,CellClass.(celltypes{tt})),3);
end

%%
figure
for tt = 1:length(celltypes)
subplot(5,4,tt*4)
    imagesc(CONDXY.Xbins(1,:,1),CONDXY.Ybins(1,:,1), ISIdistbyPSS.(celltypes{tt})')
    axis xy
    LogScale('y',10,'exp',false)
    xlabel('PSS');ylabel('ISI (s)')
    title((celltypes{tt}))
    %colorbar
    %caxis([0 0.06])
    xlim([-1.5 -0.4])

end   
  
subplot(10,4,20)
    for ss = 1:3
    plot(PSShist.bins,PSShist.(states{ss}),'color',statecolors{ss},'linewidth',2)
    hold on
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
    set(gca,'ytick',[])
    %legend(states{:},'location','eastoutside')
    axis tight
    box off
    xlim([-1.5 -0.4])
    
    
NiceSave('PSSISI',figfolder,baseName,'includeDate',true)

