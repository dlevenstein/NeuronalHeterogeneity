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
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
basePath = pwd;
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
% 
lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID;
downsamplefactor = 5;
lfp = bz_GetLFP(lfpchan,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%Noralize the LFP
lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

%% PlaceFields
for tt = 1:2
firingMaps.mean{tt} = squeeze(mean(firingMaps.rateMaps{tt},2));
for cc = 1:size(firingMaps.mean{tt},1)
    [maxpeak, firingMaps.peak{tt}(cc)] = max(firingMaps.mean{tt}(cc,5:end-5));
    firingMaps.norm{tt}(cc,:) = firingMaps.mean{tt}(cc,:)./maxpeak;
    
end
[~,sorts.placepeak{tt}] = sort(firingMaps.peak{tt});
end


%%
figure
imagesc(firingMaps.norm{2}(sorts.placepeak{tt},:))
clim([0 1])
%% Restrict to state

% state = states{3};
% ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
%     ISIStats.allspikes.times,'UniformOutput',false);

%%
ISIrate = bz_LoadCellinfo(basePath,'ISIrate');


%%
for ss = 1:3
%ss = 1;
[~,sorts.MTO] = sort(ISIrate.OccupancyStats.(states{ss}).median);
exwin = bz_RandomWindowInIntervals(SleepState.ints.(states{ss}),25);
inwin = InIntervals(ISIrate.timestamps,exwin);

sorts.MTO_E = intersect(sorts.MTO,find(CellClass.pE),'stable');
sorts.MTO_I = intersect(sorts.MTO,find(CellClass.pI),'stable');


figure
subplot(3,1,1)
bz_MultiLFPPlot(lfp,'timewin',exwin,'spikes',spikes,...
    'sortmetric',ISIrate.OccupancyStats.(states{ss}).median,...
    'cellgroups',{CellClass.pE,CellClass.pI})
subplot(3,1,2)
imagesc(ISIrate.timestamps(inwin),[1 spikes.numcells],log10(1./ISIrate.ISI(inwin,sorts.MTO_E))')
colorbar('location','east')
%xlim(exwin)

subplot(3,1,3)
imagesc(ISIrate.timestamps(inwin),[1 spikes.numcells],log10(1./ISIrate.ISI(inwin,sorts.MTO_I))')
colorbar('location','east')
%xlim(exwin)

NiceSave(['ISIraster_',(states{ss})],figfolder,baseName,'includeDate',true)
end
