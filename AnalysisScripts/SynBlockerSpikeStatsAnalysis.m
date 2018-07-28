function [ ] = SynBlockerSpikeStatsAnalysis(basePath,figfolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%

%% DEV
basePath = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/Datasets/MVData/pvch13_180719_143357';
figfolder = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/SynBlockerSpikeStatsAnalysis';
%figfolder = '/mnt/data1/Dropbox/research/Current Projects/FRHET_temp/SpikeStatsAnalysis';
%REMOVE THIS SOON, JUST TO GET STUFF TO DP
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
%CellClass = bz_LoadCellinfo(basePath,'CellClass');
[sessionInfo] = bz_getSessionInfo(basePath, 'noPrompts', false)

%%
% [spikes] = bz_LoadPhy('basepath',basePath,'kilosort_path',basePath,...
%     'getWaveforms',false);
%%
%spikes.times = cellfun(@(X) X./20000,spikes.ts,'UniformOutput',false);
%%
timeints.lightStim = [27 43].*60.*2/3;
timeints.dark = [0 27;43 50].*60.*2/3;
timeints.Eonly = [57 92;102 175].*60.*2/3;
timeints.Eonly = [102 175].*60.*2/3;
timeints.Eonly_light = [92 102].*60.*2/3;
timeints.EI = [180 230].*60.*2/3;
%%
%Waveform classification
%CellClass = bz_CellClassification(basePath);

%%
[ ISIstats ] = bz_ISIStats( spikes,'ints',timeints,...
    'figfolder',figfolder,'forceRedetect',true,...
    'savecellinfo',true,'basePath',basePath);
%%

ISIstats.summstats.Eonly.meanrate(isnan(ISIstats.summstats.Eonly.meanrate))=0;
ISIstats.summstats.dark.meanrate(isnan(ISIstats.summstats.dark.meanrate))=0;

ISIstats.summstats.lightchange = ISIstats.summstats.lightStim.meanrate./...
    ISIstats.summstats.dark.meanrate;
[~,ISIstats.sorts.lightchange] = sort(ISIstats.summstats.lightchange);

ISIstats.summstats.Esynchange = ISIstats.summstats.Eonly.meanrate./...
    ISIstats.summstats.dark.meanrate;
[~,ISIstats.sorts.Esynchange] = sort(ISIstats.summstats.Esynchange);

ISIstats.summstats.Esynchange_CV2 = ISIstats.summstats.Eonly.meanCV2-...
    ISIstats.summstats.dark.meanCV2;
%%
figure
subplot(2,2,1)
plot(log10(ISIstats.summstats.dark.meanrate),log10(ISIstats.summstats.lightStim.meanrate),'k.')
hold on
plot([-4 2],[-4 2],'k')
LogScale('xy',10)
xlabel('Spike Rate: Dark (Hz)');ylabel('Spike Rate: Light (Hz)')
title('Pre Treat')

subplot(2,2,2)
plot((ISIstats.summstats.dark.meanCV2),(ISIstats.summstats.lightStim.meanCV2),'k.')
hold on
plot([0 2],[0 2],'k')
xlabel('CV2: Dark');ylabel('CV2: Light')
title('Pre Treat')

subplot(2,2,3)
plot(log10(ISIstats.summstats.Eonly.meanrate),log10(ISIstats.summstats.Eonly_light.meanrate),'k.')
hold on
plot([-4 2],[-4 2],'k')
LogScale('xy',10)
xlabel('Spike Rate: Dark (Hz)');ylabel('Spike Rate: Light (Hz)')
title('E-Block')

subplot(2,2,4)
plot((ISIstats.summstats.Eonly.meanCV2),(ISIstats.summstats.Eonly_light.meanCV2),'k.')
hold on
plot([0 2],[0 2],'k')
xlabel('CV2: Dark');ylabel('CV2: Light')
title('E-Block')
%%
figure
subplot(2,2,1)
plot(log10(ISIstats.summstats.dark.meanrate),log10(ISIstats.summstats.Eonly.meanrate),'k.')
hold on
plot([-4 2],[-4 2],'k')
LogScale('xy',10)
xlabel('Rate: Dark (Hz)');ylabel('Rate: E Block (Hz)')


subplot(2,2,2)
plot((ISIstats.summstats.dark.meanCV2),(ISIstats.summstats.Eonly.meanCV2),'k.')
hold on
plot([0 2],[0 2],'k')
xlabel('CV2: Dark');ylabel('CV2: E Block')
%LogScale('xy',10)

%%
figure
subplot(2,2,1)
    plot(log10(ISIstats.summstats.dark.meanrate), log10(ISIstats.summstats.Esynchange),'k.')
    hold on
    plot(get(gca,'xlim'),[0 0],'k')
    LogScale('xy',10)
    xlabel('Mean Rate (Dark)');ylabel('Rate (E block)/Rate (Dark)')
subplot(2,2,2)
    plot(log10(ISIstats.summstats.Esynchange),ISIstats.summstats.Esynchange_CV2,'k.')
    hold on
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    LogScale('x',10)
    xlabel('Rate (E block)/Rate (Dark)');ylabel('CV2 (E block)-CV2 (Dark)')
subplot(2,2,3)
    plot(log10(ISIstats.summstats.Eonly.meanrate),ISIstats.summstats.Eonly.meanCV2,'k.')
    hold on
    plot(get(gca,'xlim'),[1 1],'k')
    LogScale('x',10)
    xlabel('Rate (Eblock)');ylabel('CV2 (Eblock)')
%%
spkwinsize = 20;
spkmat = bz_SpktToSpkmat(spikes,'binsize',spkwinsize,'overlap',4);
%%
CV2mat.winsize = spkwinsize;
CV2mat.timestamps = spkmat.timestamps;
CV2mat.binedges = bsxfun(@(X,Y) X+Y,spkmat.timestamps,[-0.5 0.5].*CV2mat.winsize);
    [~,CV2mat.cells] = ...
        cellfun(@(X,Y) BinDataTimes(X,Y,CV2mat.binedges),...
        ISIstats.allspikes.CV2,ISIstats.allspikes.times,...
        'UniformOutput',false);
%%
CV2mat.cells = cat(2,CV2mat.cells{:});
%%
CV2mat.mean = nanmean(CV2mat.cells,2);
CV2mat.std = nanstd(CV2mat.cells,[],2);

%%
spkmat.meanrate = mean(spkmat.data./spkwinsize,2);
spkmat.stdrate = std(spkmat.data./spkwinsize,[],2);
%%
excell = randsample(length(spikes.times),1);

excell = 7;

figure
subplot(3,1,1)
imagesc(spkmat.timestamps,[0 spikes.numcells],   log10(spkmat.data(:,ISIstats.sorts.Esynchange)./spkwinsize)')
hold on
StateScorePlot( timeints,{'r','k','b','r','g'} )
axis xy
xlabel('t (s)')
%xlim([0 12000])
ylabel('Cell (Sorted by Dark Rate)')
%colorbar('east')
LogScale('c',10)

subplot(3,1,2)
imagesc(CV2mat.timestamps,[0 1],CV2mat.cells(:,ISIstats.sorts.dark.rate)')
%xlim([0 12000])
axis xy

subplot(3,1,2)
plot(ISIstats.allspikes.times{excell},ISIstats.allspikes.CV2{excell},'.')
hold on
plot(CV2mat.timestamps,CV2mat.cells(:,excell),'o-')
plot(get(gca,'xlim'),[1 1],'k')
%xlim([0 12000])
xlabel('t (s)')
ylabel(['CV2: cell ',num2str(excell)])

subplot(3,1,3)
plot(CV2mat.timestamps,CV2mat.mean,'k','linewidth',2)
hold on
errorshade(CV2mat.timestamps,CV2mat.mean,CV2mat.std,CV2mat.std,'k','scalar')
ylim([0 2])
%xlim([0 12000])
ylabel('Mean CV2')
xlabel('t (s)')


%%
ratecolormap = makeColorMap([1 1 1],[0.5 0.5 0.5],[0.8 0 0]);
cv2colormap = makeColorMap([0 0 0.8],[0.5 0.5 0.5],[0.8 0 0]);
figure
subplot(3,1,1)
imagesc(spkmat.timestamps,[0 spikes.numcells],   log10(spkmat.data(:,ISIstats.sorts.Esynchange)./spkwinsize)')
hold on
StateScorePlot( timeints,{'r','k','b','r','g'} )
axis xy
xlabel('t (s)')
%xlim([0 12000])
ylabel('Cell (Sorted by Rate Change)')
colorbar('eastoutside')
colormap(gca,ratecolormap)
LogScale('c',10)

subplot(3,1,2)
h = imagesc(CV2mat.timestamps,[0 1],CV2mat.cells(:,ISIstats.sorts.Esynchange)');
set(h, 'AlphaData', ~isnan(CV2mat.cells(:,ISIstats.sorts.Esynchange)'));
colormap(gca,cv2colormap)
%xlim([0 12000])
caxis([0.6 1.4])
colorbar('eastoutside')
axis xy

%%
figure
subplot(4,1,1)
plot(spkmat.timestamps,spkmat.meanrate,'k')
hold on
errorshade(spkmat.timestamps,spkmat.meanrate,spkmat.stdrate,spkmat.stdrate,'k','scalar')
ylim([0 20])

%%
%figure
%plot(CV2mat.cells(:,6),spkmat.data(:,6),'k.')
%%
spikes.waveforms = cat(1,spikes.rawWaveform{:});

%% channel
[~,celldepth] = ismember(spikes.maxWaveformCh,sessionInfo.AnatGrps.Channels);
[~,depthsort] = sort(celldepth);

%%
figure
subplot(2,2,1)
plot(log10(ISIstats.summstats.dark.meanrate),-celldepth,'.')
LogScale('x',10)
xlabel('Rate (Dark, Hz)');ylabel('Depth')
subplot(2,2,2)
plot(log10(ISIstats.summstats.Esynchange),-celldepth,'.')
LogScale('x',10)
xlabel('Rate Change with E Block');ylabel('Depth')
%%
figure
subplot(2,2,1)
imagesc(spikes.waveforms)

subplot(2,2,2)
plot(spikes.rawWaveform{17})
xlabel('t');ylabel('mV?')


%%
figure
plot(ISIstats.summstats.dark.meanCV2,log10(ISIstats.summstats.Esynchange),'k.')
hold on
plot([1 1],get(gca,'ylim'),'k')
plot(get(gca,'xlim'),[0 0],'k')
xlabel('CV2 (Dark)');ylabel('Rate Change with E block (Ratio)')
LogScale('y',10)
