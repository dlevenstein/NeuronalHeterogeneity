function [ ] = Analysis20190208(basePath,figfolder)
% Date 02/08/2019
%
%Questions: How does gamma coupling look during NREM,WAKE? Do some cells have 
%better coupling/CV2 relation than others and does it relate to sFR?
%Does population CV2 at finer timescales (100ms?) relate to gamma bursts?
%Need a good way to compare popCV2 binning timescale
%Or can I identify transition periods in the population or in individual 
%cells with CV2? How frequently are bursts during gamma/ripples? 
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


%% Load the LFP if needed

lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID;
ctxchan = sessionInfo.channels(ismember(sessionInfo.region,'CTX'));
downsamplefactor = 1;
lfp = bz_GetLFP(ctxchan,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%Noralize the LFP
lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

%% Channel Sorting by Spike group
%Find the channels that are used and order them by their spike group
spikegroupordering = [sessionInfo.spikeGroups.groups{:}];
[spikegroupordering,~,spikegroupsort] = intersect(spikegroupordering,ctxchan,'stable');

%Find the cells that are used and order them by their channel (spike group
%ordering)
%cellchannelIDX(i) is the channelIDX (ordered) corresponding to the i'th
%channel in spikes
%cellchannel sort will sort cells by their channel
%celltypechannelsort will sort cells by their channel and type
[~,cellchannelIDX] = ismember(spikes.maxWaveformCh,spikegroupordering);
[~,cellchanneltypesort] = sort(cellchannelIDX+max(cellchannelIDX).*CellClass.pI);
[cellchannelIDX,cellchannelsort] = sort(cellchannelIDX);

%% Restrict to state

state = states{1}; %was using REM?! oops....
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);


%% Spike coupling with wavelets

%Take only subset of time (random intervals) so wavelets doesn't break
%computer (total 625s)
if sum(diff(SleepState.ints.(state),1,2))>625
    nwin = 25;
    winsize = 25; %s
    windows = bz_RandomWindowInIntervals( SleepState.ints.(state),winsize,nwin );
else
    windows = SleepState.ints.(state);
end
%Calculate pop-phase coupling for all channels
nfreqs = 100;
%synchphasecoupling = zeros(length(ctxchan),nfreqs);
%synchpowercorr = zeros(length(ctxchan),nfreqs);
for cc = 1:length(ctxchan)
    cc
[freqs,synchcoupling] = ...
    bz_GenSpikeLFPCoupling(spikes.times,lfp.data(:,cc),'sf_LFP',lfp.samplingRate,...
    'int',windows,'DOWNSAMPLE',1,'frange',[20 200],'ncyc',15,...
    'subpop',CellClass.pE+2.*CellClass.pI,'synchwin',0.002,'synchdt',0.002,...
    'nfreqs',nfreqs);
    for tt = 1:length(celltypes)
        synchphasecoupling.(celltypes{tt})(cc,:) = synchcoupling(tt).phasemag;
        synchpowercorr.(celltypes{tt})(cc,:) = synchcoupling(tt).powercorr;
    end
end

%% Gamma-phase coupling in the best band(s) for each cell with each channel

gamma.fband = [50 90]; %REM best

% gamma = bz_Filter(lfp,'passband',gammaband);
clear ratepowercorr
clear spikephasemag
clear spikephaseangle
for cc = 1:length(ctxchan) 
    cc
    [~,synchcoupling,ratepowercorr(cc,:),...
        spikephasemag(cc,:),spikephaseangle(cc,:)] = ...
        bz_GenSpikeLFPCoupling(spikes.times,lfp.data(:,cc),'sf_LFP',lfp.samplingRate,...
        'int',SleepState.ints.(state),'DOWNSAMPLE',1,'frange',gamma.fband,'ncyc',5,...
        'subpop',CellClass.pE+2.*CellClass.pI,'nfreqs',1);
    for tt = 1:length(celltypes)
        synchgammaphasecoupling.(celltypes{tt})(cc) = synchcoupling(tt).phasemag;
        synchgammapowercorr.(celltypes{tt})(cc) = synchcoupling(tt).powercorr;
    end
end

%Pick the channel with the best coupling to the pE population
[~,gamma.pEchanIDX] = max(synchgammaphasecoupling.pE);
gamma.pEchan = lfp.channels(gamma.pEchanIDX);

%% Figure: pop-phase coupling by channel
figure
    for tt = 1:length(celltypes)
        subplot(2,3,tt)
        imagesc(log2(freqs),[1 length(ctxchan)],synchphasecoupling.(celltypes{tt})(spikegroupsort,:));
        hold on
        plot(log2([gamma.fband ; gamma.fband]),[1 length(ctxchan);1 length(ctxchan)]',...
            'w--')
        plot(bz_NormToRange(synchgammaphasecoupling.(celltypes{tt})(spikegroupsort),log2(gamma.fband)),...
            1:length(ctxchan),'color',cellcolor{tt})
        plot(log2(gamma.fband(2)),find(ismember(spikegroupordering,gamma.pEchan)),'+')
        title([celltypes{tt},' Pop. Rate'])
        LogScale('x',2)
        colorbar('northoutside')

        xlabel('freq (Hz)');
        if tt==1
        ylabel({'Phase Coupling','Channel'})
        end 
        %xlim(log2([16 192]))
        
        
        subplot(2,3,tt+3)
        imagesc(log2(freqs),[1 length(ctxchan)],synchpowercorr.(celltypes{tt})(spikegroupsort,:));
        hold on
        plot(log2([gamma.fband ; gamma.fband]),[1 length(ctxchan);1 length(ctxchan)]',...
            'w--')
        plot(bz_NormToRange(synchgammapowercorr.(celltypes{tt})(spikegroupsort),log2(gamma.fband)),...
            1:length(ctxchan),'color',cellcolor{tt})
        title([celltypes{tt},' Pop. Rate'])
        LogScale('x',2)
        colorbar('northoutside')
        xlabel('freq (Hz)');
        if tt==1
        ylabel({'Rate-Power Corr.','Channel'})
        end 

    end
    
    
    subplot(2,3,3)
    imagesc([1 spikes.numcells],[1 length(ctxchan)],spikephasemag(spikegroupsort,cellchanneltypesort));
    hold on
    idxvec = 0;
    for tt = 1:length(celltypes)
        idxvec=idxvec(end)+[1:sum(CellClass.(celltypes{tt}))];
         plot(idxvec,cellchannelIDX(CellClass.(celltypes{tt})),...
             '.','color',cellcolor{tt})
    end
    title('All Cells')
            colorbar('northoutside')

   % LogScale('x',2)
    xlabel('Cell');%ylabel('Channel')
    subplot(2,3,6)
    imagesc([1 spikes.numcells],[1 length(ctxchan)],ratepowercorr(spikegroupsort,cellchanneltypesort));
    hold on
    idxvec = 0;
    for tt = 1:length(celltypes)
        idxvec=idxvec(end)+[1:sum(CellClass.(celltypes{tt}))];
         plot(idxvec,cellchannelIDX(CellClass.(celltypes{tt})),...
             '.','color',cellcolor{tt})
    end
    title('All Cells')
            colorbar('northoutside')

   % LogScale('x',2)
    xlabel('Cell');%ylabel('Channel')
NiceSave(['GammaCoupling',state],figfolder,baseName,'includeDate',true)

%%
figure
plot(log10(ISIStats.summstats.(state).meanrate),spikephasemag,'.')
