function [ ] = Analysis20190202_Gamma(basePath,figfolder)
% Date 02/02/2019
%
%Question: How does population-gamma coupling vary by channel? Note: this
%would be better to do with a cortical (or hippocampal) dataset that has
%full depth information... Yuta's or William's
%This will be used to best select gamma band/channel for further ISI
%analysis
%
%Plots
%-Wavelet-pop synch coupling 
%-Channel-channel gamma coupling
%-Cell-gamma coupling to each channel
%-Two observed gamma bands: gamma and high gamma (spiking)
%-Calcualted best filter order (5 cycles, not much difference)
%-Gamma in the return map.... need power and pMRL
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/Dino_mPFC/Dino_061814';
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

state = states{3}; %was using REM?! oops....
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);


%% Spike coupling with wavelets
%Calculate pop-phase coupling for all channels
nfreqs = 100;
%synchphasecoupling = zeros(length(ctxchan),nfreqs);
%synchpowercorr = zeros(length(ctxchan),nfreqs);
for cc = 1:length(ctxchan)
    cc
[freqs,synchcoupling] = ...
    bz_GenSpikeLFPCoupling(spikes.times,lfp.data(:,cc),'sf_LFP',lfp.samplingRate,...
    'int',SleepState.ints.(state),'DOWNSAMPLE',1,'frange',[20 200],'ncyc',15,...
    'subpop',CellClass.pE+2.*CellClass.pI,'synchwin',0.002,'synchdt',0.002,...
    'nfreqs',nfreqs);
    for tt = 1:length(celltypes)
        synchphasecoupling.(celltypes{tt})(cc,:) = synchcoupling(tt).phasemag;
        synchpowercorr.(celltypes{tt})(cc,:) = synchcoupling(tt).powercorr;
    end
end


%% Gamma-phase coupling in the best band(s) for each cell with each channel

gamma.fband = [50 90];
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
NiceSave('GammaCoupling',figfolder,baseName)


%% High gamma. reflects spiking?
 hgamma.fband = [110 200];

 %PUT THIS LOOP INTO GenSpikeLFPCoupling... for multiple channels. No need
 %to calculate synchrony etc multiple times
for cc = 1:length(ctxchan) 
    cc
    [~,synchcoupling,hgamma.ratepowercorr(cc,:),...
        hgamma.spikephasemag(cc,:),hgamma.spikephaseangle(cc,:)] = ...
        bz_GenSpikeLFPCoupling(spikes.times,lfp.data(:,cc),'sf_LFP',lfp.samplingRate,...
        'int',SleepState.ints.(state),'DOWNSAMPLE',1,'frange',hgamma.fband,'ncyc',5,...
        'subpop',CellClass.pE+2.*CellClass.pI,'nfreqs',1);
    for tt = 1:length(celltypes)
        hgamma.synchphasecoupling.(celltypes{tt})(cc) = synchcoupling(tt).phasemag;
        hgamma.synchpowercorr.(celltypes{tt})(cc) = synchcoupling(tt).powercorr;
    end
end
%% Figure: pop-phase coupling by channel: HGamma
figure
    for tt = 1:length(celltypes)
        subplot(2,3,tt)
        imagesc(log2(freqs),[1 length(ctxchan)],synchphasecoupling.(celltypes{tt})(spikegroupsort,:));
        hold on
        plot(log2([hgamma.fband ; hgamma.fband]),[1 length(ctxchan);1 length(ctxchan)]',...
            'w--')
        plot(bz_NormToRange(hgamma.synchphasecoupling.(celltypes{tt})(spikegroupsort),log2(hgamma.fband)),...
            1:length(ctxchan),'color',cellcolor{tt})
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
        plot(log2([hgamma.fband ; hgamma.fband]),[1 length(ctxchan);1 length(ctxchan)]',...
            'w--')
        plot(bz_NormToRange(hgamma.synchpowercorr.(celltypes{tt})(spikegroupsort),log2(hgamma.fband)),...
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
    imagesc([1 spikes.numcells],[1 length(ctxchan)],hgamma.spikephasemag(spikegroupsort,cellchanneltypesort));
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
    imagesc([1 spikes.numcells],[1 length(ctxchan)],hgamma.ratepowercorr(spikegroupsort,cellchanneltypesort));
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
    
NiceSave('HGammaCoupling',figfolder,baseName)

%% Question: what filter order?
gammaordersweep.forders = 2:15;

gammaordersweep.fband = gamma.fband;
for oo = 1:length(gammaordersweep.forders) 
    oo
    [~,synchcoupling,gammaordersweep.ratepowercorr(oo,:),...
        gammaordersweep.spikephasemag(oo,:),gammaordersweep.spikephaseangle(oo,:)] = ...
        bz_GenSpikeLFPCoupling(spikes.times,lfp.data(:,gamma.pEchanIDX),'sf_LFP',lfp.samplingRate,...
        'int',SleepState.ints.(state),'DOWNSAMPLE',1,'frange',gammaband,'ncyc',gammaordersweep.forders(oo),...
        'subpop',CellClass.pE+2.*CellClass.pI,'nfreqs',1);
    for tt = 1:length(celltypes)
        gammaordersweep.synchphasecoupling.(celltypes{tt})(oo) = synchcoupling(tt).phasemag;
        gammaordersweep.synchpowercorr.(celltypes{tt})(oo) = synchcoupling(tt).powercorr;
    end
end

%%
figure
subplot(2,2,1)
hold on
for tt = 1:length(celltypes)
    plot(gammaordersweep.forders,gammaordersweep.synchphasecoupling.(celltypes{tt}),'color',cellcolor{tt})
end
xlabel('Filter Order (# cycles)');ylabel('Pop-Phase Coupling')
axis tight
box off 

subplot(2,2,2)
hold on
for tt = 1:length(celltypes)
    plot(gammaordersweep.forders,gammaordersweep.synchpowercorr.(celltypes{tt}),'color',cellcolor{tt})
end
xlabel('Filter Order (# cycles)');ylabel('Rate-Power Correlation')
axis tight
box off

subplot(2,2,3)
imagesc(gammaordersweep.forders,[1 length(ctxchan)],gammaordersweep.spikephasemag(:,cellchanneltypesort)')
xlabel('Filter Order (# cycles)');ylabel('Cell')
title('Phase-Coupling')

subplot(2,2,4)
imagesc(gammaordersweep.forders,[1 length(ctxchan)],gammaordersweep.ratepowercorr(:,cellchanneltypesort)')
xlabel('Filter Order (# cycles)');ylabel('Cell')
title('Rate-Power Corr')

NiceSave('GammaFilterOrder',figfolder,baseName)


%% Filter gamma

gammaLFP = bz_Filter(lfp,'passband',gamma.fband,'order',5,'filter','fir1',...
    'fast',true,'channels',gamma.pEchan);
gammaLFP.normamp = NormToInt(log10(gammaLFP.amp),'modZ',SleepState.ints.(state),gammaLFP.samplingRate);

gamma.powerphasemaps = bz_PowerPhaseRatemap(spikes,gammaLFP,'ints',SleepState.ints.(state));

%% Calculate gamma power/phase at every spike -
%show phase and power for all spikes in the return map and calculate weighted mrl 
%every point in the return map, using gaussian weighting.


ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end);nan],ISIStats.allspikes.ISIs,...
    'UniformOutput',false);

ISIStats.allspikes.gammapower = ...
    cellfun(@(Y) interp1(gammaLFP.timestamps,gammaLFP.normamp,Y,'nearest'),...
    ISIStats.allspikes.times,...
    'UniformOutput',false);
ISIStats.allspikes.gammahilb = ...
    cellfun(@(Y) interp1(gammaLFP.timestamps,gammaLFP.hilb,Y,'nearest'),...
    ISIStats.allspikes.times,...
    'UniformOutput',false);
ISIStats.allspikes.gammaphase = ...
    cellfun(@(Y) interp1(gammaLFP.timestamps,gammaLFP.phase,Y,'nearest'),...
    ISIStats.allspikes.times,...
    'UniformOutput',false);

%% Gamma in the return map
%Gamma power in the return map
ISIbounds = [-3 1.5];
[meanZ,countmap ] = cellfun(@(X,Y,Z,Q) ConditionalHist3(log10(X(Q)),...
    log10(Y(Q)),(Z(Q)),...
    'numXbins',80,'numYbins',80,'Xbounds',ISIbounds,'Ybounds',ISIbounds,...
    'minXY',25,'sigma',0.2),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,ISIStats.allspikes.gammapower,...
    ISIStats.allspikes.instate,'UniformOutput',false);

for tt = 1:length(celltypes)
    meangammapowerreturn.(celltypes{tt}) = nanmean(cat(3,meanZ{CellClass.(celltypes{tt})}),3);
    meanreturn.(celltypes{tt}) = nanmean(cat(3,countmap{CellClass.(celltypes{tt})}),3);

end
%%
%Gamma pMRL in the return map
ISIbounds = [-3 1.5];
[meanZ,countmap ] = cellfun(@(X,Y,Z,Q) ConditionalHist3(log10(X(Q)),...
    log10(Y(Q)),(Z(Q)),...
    'numXbins',80,'numYbins',80,'Xbounds',ISIbounds,'Ybounds',ISIbounds,...
    'minXY',10,'sigma',0.2),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,ISIStats.allspikes.gammahilb,...
    ISIStats.allspikes.instate,'UniformOutput',false);
%%
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
    scatter(log10(ISIStats.allspikes.ISIs{excell}(ISIStats.allspikes.instate{excell})),...
        log10(ISIStats.allspikes.ISInp1{excell}(ISIStats.allspikes.instate{excell})),1,...
        (ISIStats.allspikes.gammapower{excell}(ISIStats.allspikes.instate{excell})))
    colorbar
    
subplot(2,2,1)
    imagesc(ISIbounds,ISIbounds,abs(meanZ{excell})')
    axis xy
    LogScale('xy',10)
    
    
%% Figure: population return map
cmap = [1 1 1;colormap(gca)];

figure
colormap(cmap)
for tt = 1:length(celltypes)
    subplot(2,2,tt)
        s = imagesc(ISIbounds,ISIbounds,abs(meangammapowerreturn.(celltypes{tt}))');
        axis xy
        alpha(s,meanreturn.(celltypes{tt})'*25)
        axis xy
        LogScale('xy',10)
        colorbar
end

%% Conditional CV2 distribution on Gamma Power

[ CONDXY ] = ConditionalHist(ISIStats.allspikes.gammapower,...
    ISIStats.allspikes.CV2,'Xbounds',[-3 1.65],'numXbins',13,...
    'numYbins',20,'Ybounds',[0 2]);

%%

for tt = 1:length(celltypes)
    CV2distbyPOP.(celltypes{tt}) = nanmean(CONDXY.pYX(:,:,CellClass.(celltypes{tt})),3);
    meanCV2byPOP.(celltypes{tt}) = nanmean(CONDXY.meanYX(:,:,CellClass.(celltypes{tt})),3);
end
%%
figure
for tt = 1:length(celltypes)
subplot(2,2,tt)
imagesc(CONDXY.Xbins(:,:,1),CONDXY.Ybins(:,:,1),CV2distbyPOP.(celltypes{tt})')
hold on
plot(CONDXY.Xbins(:,:,excell),meanCV2byPOP.(celltypes{tt}),'w-')
plot(CONDXY.Xbins(:,[1 end],excell),[1 1],'w--')
title(celltypes{tt})
axis xy
xlabel('Gamma Power');ylabel('CV2')
end



for tt = 1:length(celltypes)
    subplot(2,2,tt+2)
    colormap(gca,cmap)
        s = imagesc(ISIbounds,ISIbounds,abs(meangammapowerreturn.(celltypes{tt}))');
        axis xy
        alpha(s,meanreturn.(celltypes{tt})'*25)
        axis xy
        LogScale('xy',10)
        colorbar
        title(celltypes{tt})

end

NiceSave('CV2byGammaPower',figfolder,baseName)
