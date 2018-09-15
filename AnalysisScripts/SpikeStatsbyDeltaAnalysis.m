function [  ] = SpikeStatsbyDeltaAnalysis(basePath,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyDeltaAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
SleepState.ints.ALL = [0 Inf];
statenames = fieldnames(SleepState.ints);
numstates = length(statenames);
statecolors = {'k','b','r'};
classnames = unique(CellClass.label);
numclasses = length(classnames);
classcolors = {'k','r'};

downsamplefactor = 2;
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
     'basepath',basePath,'downsample',downsamplefactor);
sessionInfo = bz_getSessionInfo(basePath);

%% Filter the LFP for delta activity

deLFP = bz_Filter(lfp,'passband',[0.5 8],'order',1,'filter','fir1');

%% Calculate the power-phase ratemap, cv2map

[PowerPhaseRatemap,spikebinIDs] = bz_PowerPhaseRatemap(spikes,deLFP,...
    'ints',SleepState.ints.NREMstate);

%% GLM for coupling
excell = 9;
deLFP.data = deLFP.hilb;
deLFP.freqs = 4;
[ GLMFP ] = GLMLFP(spikes.times(excell),deLFP,...
    'intervals',SleepState.ints.NREMstate );


%% Simulate spikes from GLM
for tt = 1:length(GLMFP.timestamps)
    GLMFP.spkmat(tt) = rand(1)<=GLMFP.predRate(tt);
end
simspikes.times = {GLMFP.timestamps(GLMFP.spkmat)};


[PowerPhaseRatemap_sim] = bz_PowerPhaseRatemap(simspikes,deLFP,...
    'ints',SleepState.ints.NREMstate);
%% Figure
figure
subplot(3,3,1)
    imagesc(PowerPhaseRatemap.phasebins,PowerPhaseRatemap.powerbins,...
        PowerPhaseRatemap.ratemap{excell})
    hold on
    imagesc(PowerPhaseRatemap.phasebins+2*pi,PowerPhaseRatemap.powerbins,...
        PowerPhaseRatemap.ratemap{excell})
    xlim([-pi 3*pi])
    caxis([0 15])
    axis xy
    colorbar
    
subplot(3,3,2)
    imagesc(PowerPhaseRatemap_sim.phasebins,PowerPhaseRatemap_sim.powerbins,...
        PowerPhaseRatemap_sim.ratemap{1})
    hold on
    imagesc(PowerPhaseRatemap_sim.phasebins+2*pi,PowerPhaseRatemap_sim.powerbins,...
        PowerPhaseRatemap_sim.ratemap{1})
    xlim([-pi 3*pi])
    caxis([0 15])
    axis xy
    colorbar

subplot(3,3,4)
    plot(GLMFP.powerbins,GLMFP.Rpower,'k','linewidth',2)
    axis tight
    box off
    xlabel('Power');ylabel('Phase Coupling')


end

