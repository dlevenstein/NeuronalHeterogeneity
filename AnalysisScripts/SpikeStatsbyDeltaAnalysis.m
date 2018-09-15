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

deLFP = bz_Filter(lfp,'passband',[0.5 6],'order',1,'filter','fir1');

%% Calculate the power-phase ratemap, cv2map

[PowerPhaseRatemap,spikebinIDs] = bz_PowerPhaseRatemap(spikes,deLFP,...
    'ints',SleepState.ints.NREMstate);

%% GLM for coupling
excell = 9;

predLFP = deLFP;
predLFP.data = predLFP.hilb;
predLFP.freqs = 4;
[ GLMFP ] = GLMLFP(spikes.times(excell),deLFP,...
    'intervals',SleepState.ints.NREMstate );


%% Simulate spikes from GLM
for tt = 1:length(GLMFP.timestamps)
    GLMFP.spkmat(tt) = rand(1)<=GLMFP.predRate(tt);
end
simspikes.times = {GLMFP.timestamps(GLMFP.spkmat)};
simspikes.ISIs = diff(simspikes.times{1});
simspikes.ISIhist = hist(log10(simspikes.ISIs),ISIStats.ISIhist.logbins);
simspikes.ISIhist = simspikes.ISIhist./sum(simspikes.ISIhist);
[PowerPhaseRatemap_sim] = bz_PowerPhaseRatemap(simspikes,deLFP,...
    'ints',SleepState.ints.NREMstate);
%% Figure
viewwin = bz_RandomWindowInIntervals(SleepState.ints.NREMstate,5);

figure
subplot(4,1,1)
    plot(lfp.timestamps,lfp.data,'k')
    hold on
    plot(deLFP.timestamps,deLFP.data,'b')
    plot(spikes.times{excell},max(double(lfp.data)).*ones(size(spikes.times{excell})),'k.')
    xlim(viewwin)
    ylabel({'LFP','Real Spikes'});
    box off
    set(gca,'xticklabels',[]);set(gca,'yticklabels',[])
    
subplot(4,1,2)
    plot(GLMFP.timestamps,GLMFP.predRate,'k')
    hold on
    plot(simspikes.times{1},max(GLMFP.predRate).*ones(size(simspikes.times{1})),'k.')
    xlim(viewwin)
    ylabel({'Predicted Rate','Simulated Spikes'});
    box off
    set(gca,'xticklabels',[])

subplot(4,3,7)
    imagesc(PowerPhaseRatemap.phasebins,PowerPhaseRatemap.powerbins,...
        PowerPhaseRatemap.ratemap{excell})
    hold on
    imagesc(PowerPhaseRatemap.phasebins+2*pi,PowerPhaseRatemap.powerbins,...
        PowerPhaseRatemap.ratemap{excell})
    xlim([-pi 3*pi])
    caxis([0 15])
    axis xy
    colorbar
    xlabel('Phase');ylabel('Power')
    title('Real Delta Ratemap')
    
subplot(4,3,8)
    imagesc(PowerPhaseRatemap_sim.phasebins,PowerPhaseRatemap_sim.powerbins,...
        PowerPhaseRatemap_sim.ratemap{1})
    hold on
    imagesc(PowerPhaseRatemap_sim.phasebins+2*pi,PowerPhaseRatemap_sim.powerbins,...
        PowerPhaseRatemap_sim.ratemap{1})
    xlim([-pi 3*pi])
    caxis([0 15])
    axis xy
    colorbar
    xlabel('Phase');ylabel('Power')
    title('Simulated Delta Ratemap')

subplot(4,3,11)
    plot(GLMFP.powerbins,GLMFP.Rpower,'k','linewidth',2)
    hold on
    text(-1.5,1,['R0 = ',num2str(round(GLMFP.R0*lfp.samplingRate,1)),' Hz'])
    axis tight
    box off
    xlabel('Power');ylabel('GLM: Phase Coupling')
    %title('GLM Kernel')

subplot(4,3,9)
    plot(ISIStats.ISIhist.logbins,ISIStats.ISIhist.NREMstate.log(excell,:),'k','linewidth',2)
    hold on
    plot(ISIStats.ISIhist.logbins,simspikes.ISIhist,'r','linewidth',2)
    axis tight
    box off
    xlabel('ISI (s)')
    LogScale('x',10)
    legend('location','northeast','Real','Sim.')
    
NiceSave('ExCellDelta',figfolder,baseName)
end

