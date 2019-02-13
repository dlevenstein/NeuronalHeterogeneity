function [ ] = Analysis20190201(basePath,figfolder)
% Friday 02/01/2019
%
%Question: How is spiking statistics coupled to high frequency oscillations
%(Gamma/Spindle/Ripple)?
%-Gamma power in return map
%-Phase/Power Ratemap with gamma
%-Phase/Power coupling by cell
%-GLM with gamma phase/power couplings
%% Load Header
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/Dino_mPFC/Dino_061814';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/Analysis20190201'];
baseName = bz_BasenameFromBasepath(basePath);

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


%%
%Use slow wave channel for the LFP. Eventually will want to see if channel
%selection matters - calculate coupling for each channel for each cell and
%to determine the spatial extent of oscillations

downsamplefactor = 1;
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

%% Restrict to state
%Start with NREM.

state = states{3};
instatespiketimes = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);
%% Spike coupling with wavelets
%
bz_GenSpikeLFPCoupling(spikes.times,lfp.data,'sf_LFP',lfp.samplingRate,...
    'int',SleepState.ints.(state),'DOWNSAMPLE',1,'frange',[16 156],'ncyc',15,...
    'subpop',CellClass.pE+2.*CellClass.pI,'saveFig',figfolder);


%% Filter in gamma band and calculate phase-power rate coupling
%Which gamma frequencies? Can pick by calculating phase/power coupling to 
%wavelet spectrogram.
gammaband = [50 100];
gamma = bz_Filter(lfp,'passband',gammaband);

bz_GenSpikeLFPCoupling(spikes.times,lfp.data,'sf_LFP',lfp.samplingRate,...
    'int',SleepState.ints.(state),'DOWNSAMPLE',1,'frange',gammaband,'ncyc',10,...
    'subpop',CellClass.pE+2.*CellClass.pI,'nfreqs',1)

bz_PowerPhaseRatemap(spikes,gamma,'ints',SleepState.ints.(state));
%% Calculate gamma power/phase at every spike -
%show phase and power for all spikes in the return map and calculate weighted mrl 
%every point in the return map, using gaussian weighting.

gamma.amp = NormToInt(gamma.amp,'modZ',SleepState.ints.(state),gamma.samplingRate);

ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end);nan],ISIStats.allspikes.ISIs,...
    'UniformOutput',false);

ISIStats.allspikes.gammapower = ...
    cellfun(@(Y) interp1(gamma.timestamps,gamma.amp,Y,'nearest'),...
    ISIStats.allspikes.times,...
    'UniformOutput',false);


ISIStats.allspikes.gammahilb = ...
    cellfun(@(Y) interp1(gamma.timestamps,gamma.hilb,Y,'nearest'),...
    ISIStats.allspikes.times,...
    'UniformOutput',false);

ISIStats.allspikes.gammaphase = ...
    cellfun(@(Y) interp1(gamma.timestamps,gamma.phase,Y,'nearest'),...
    ISIStats.allspikes.times,...
    'UniformOutput',false);
%%
%Gamma power in the return map
ISIbounds = [-3 1.5];
[meanZ,countmap ] = cellfun(@(X,Y,Z,Q) ConditionalHist3(log10(X(Q)),...
    log10(Y(Q)),(Z(Q)),...
    'numXbins',80,'numYbins',80,'Xbounds',ISIbounds,'Ybounds',ISIbounds,...
    'minXY',30,'sigma',0.25),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,ISIStats.allspikes.gammahilb,...
    instatespiketimes,'UniformOutput',false);


for tt = 1:length(celltypes)
    meangammapowerreturn.(celltypes{tt}) = nanmean(cat(3,meanZ{CellClass.(celltypes{tt})}),3);
    meanreturn.(celltypes{tt}) = nanmean(cat(3,countmap{CellClass.(celltypes{tt})}),3);

end
%% Figure: an example cell
excell = randsample(spikes.numcells,1);

figure
%subplot(2,2,1)
    %plot(ISIStats.allspikes.cellrate{excell},ISIStats.allspikes.gammapower{excell},'.')
subplot(2,2,3)
    scatter(log10(ISIStats.allspikes.ISIs{excell}(instatespiketimes{excell})),...
        log10(ISIStats.allspikes.ISInp1{excell}(instatespiketimes{excell})),1,...
        (ISIStats.allspikes.gammapower{excell}(instatespiketimes{excell})))
    colorbar
    
subplot(2,2,1)
    imagesc(ISIbounds,ISIbounds,meanZ{excell}')
    axis xy
    LogScale('xy',10)
    
    
%% Figure: population return map
figure
for tt = 1:length(celltypes)
    subplot(2,2,tt)
        s = imagesc(ISIbounds,ISIbounds,meangammapowerreturn.(celltypes{tt})');
        alpha(s,meanreturn.(celltypes{tt})'*25)
        axis xy
        LogScale('xy',10)
        colorbar
end
