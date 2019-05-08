function [ ] = SpikeStatsbyCellActivityAnalysis(basePath,figfolder)

%% DEV
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/Dino_mPFC/Dino_061814';s
basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyCellActivityAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
SleepState.ints.ALL = [0 Inf];
% lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
%     'basepath',basePath);

%% Cell types and states
[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};
statenames = fieldnames(SleepState.ints);

%% Calculate spike count matrix

%IDEA: Instead of constant time bin.... constant spike count bin.
dt = 1;
binsize = 10;
spikemat = bz_SpktToSpkmat(spikes,'dt',1,'binsize',binsize);


for cc = 1:spikes.numcells
    spikemat.cellrate{cc} = spikemat.data(:,cc);
end


%% Calculate Cell rate for each spike
    ISIStats.allspikes.cellrate = cellfun(@(X,Y) interp1(spikemat.timestamps,X,Y,'nearest'),...
        spikemat.cellrate,ISIStats.allspikes.times,'UniformOutput',false);
    
    
%% Calculate stuff in each state
%for ss = 1:length(statenames)
state = statenames{4};

instatespiketimes = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);
instateratetimes = InIntervals(spikemat.timestamps,double(SleepState.ints.(state)));


%% Rate-CV2 histogram
numbins = 20;
cellhists.ratebins = 1:20;
cellhists.cv2bins = linspace(0,2,numbins);
cellhists.logISIbins = linspace(-2.3,1.3,numbins);

%Calculate the histogram of spikes in bins with rate X
cellhists.ratehist = cellfun(@(X) hist(X,cellhists.ratebins),ISIStats.allspikes.cellrate,...
    'UniformOutput',false);
cellhists.meanISIhist = cellfun(@(X) hist(log10(X),cellhists.logISIbins),ISIStats.allspikes.meanISI,...
    'UniformOutput',false);

cellhists.rateVcv2 = cellfun(@(X,Y,Z) hist3([X(Z),Y(Z)],{cellhists.ratebins,cellhists.cv2bins}),...
    ISIStats.allspikes.cellrate,ISIStats.allspikes.CV2,instatespiketimes,...
    'UniformOutput',false);
cellhists.rateVcv2 = cellfun(@(X,Y) bsxfun(@(x,y) x./y,X,Y'),...
    cellhists.rateVcv2,cellhists.ratehist,...
    'UniformOutput',false);

cellhists.meanisiVcv2 = cellfun(@(X,Y,Z) hist3([log10(X(Z)),Y(Z)],{cellhists.logISIbins,cellhists.cv2bins}),...
    ISIStats.allspikes.meanISI,ISIStats.allspikes.CV2,instatespiketimes,...
    'UniformOutput',false);
cellhists.meanisiVcv2 = cellfun(@(X,Y) bsxfun(@(x,y) x./y,X,Y'),...
    cellhists.meanisiVcv2,cellhists.meanISIhist,...
    'UniformOutput',false);




%%
for tt=1:length(celltypes)
    cellhists.meanhist.(celltypes{tt}) = nanmean(cat(3,cellhists.rateVcv2{CellClass.(celltypes{tt})}),3);
    cellhists.meanISIvCV2hist.(celltypes{tt}) = nanmean(cat(3,cellhists.meanisiVcv2{CellClass.(celltypes{tt})}),3);
end
%% Rate-CV2 Correlation
ISIStats.summstats.(state).rateCV2corr = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
    ISIStats.allspikes.cellrate,ISIStats.allspikes.CV2,instatespiketimes);
ISIStats.summstats.(state).ISICV2corr = cellfun(@(X,Y,Z) corr(log10(1./X(Z)),Y(Z),'type','spearman'),...
    ISIStats.allspikes.meanISI,ISIStats.allspikes.CV2,instatespiketimes);


%%
figure
subplot(2,2,1)
for tt = 1:length(celltypes)
    plot(log10(ISIStats.summstats.(state).meanrate(CellClass.(celltypes{tt}))),...
        ISIStats.summstats.(state).rateCV2corr(CellClass.(celltypes{tt})),...
        '.','color',cellcolor{tt})
    hold on
end
    plot(get(gca,'xlim'),[0 0],'k--')
    LogScale('x',10)
xlabel('Mean Rate (Hz)');ylabel('CV2-Rate Correlation')

for tt = 1:length(celltypes)
subplot(2,2,tt+2)
imagesc(cellhists.ratebins,cellhists.cv2bins,cellhists.meanhist.(celltypes{tt})')
axis xy
colorbar
caxis([0 0.08])
title(celltypes{tt})
xlabel(['Spike Count, (',num2str(binsize),'s bins']);
ylabel('CV2')
end

NiceSave('RateCV2correlation',figfolder,baseName)
%%
figure
for tt = 1:length(celltypes)
subplot(2,2,tt)
imagesc(cellhists.logISIbins,cellhists.cv2bins,cellhists.meanISIvCV2hist.(celltypes{tt})')
axis xy
colorbar
caxis([0 0.1])
xlabel('meanISI');ylabel('CV2')
LogScale('x',10)
title(celltypes{tt})
end

NiceSave('MeanISIvCV2',figfolder,baseName)


%%
excell = randsample(spikes.numcells,1);
figure
subplot(2,2,2)
plot(ISIStats.allspikes.cellrate{excell},log10(ISIStats.allspikes.meanISI{excell}),'.')
xlabel(['Spike Count, (',num2str(binsize),'s bins']);
ylabel('Mean ISI')

subplot(2,2,1)
plot(ISIStats.allspikes.cellrate{excell},(ISIStats.allspikes.CV2{excell}),'.')
xlabel(['Spike Count, (',num2str(binsize),'s bins']);
ylabel('CV2')

subplot(2,2,3)
plot(log10(ISIStats.allspikes.meanISI{excell}),ISIStats.allspikes.CV2{excell},'.')
xlabel('Mean ISI');ylabel('CV2')

subplot(2,2,4)
imagesc(cellhists.ratebins,cellhists.cv2bins,log10(cellhists.rateVcv2{excell})')
axis xy
%%
subplot(2,2,2)
for tt = 1:length(celltypes)
    plot((ISIStats.summstats.(state).meanCV2(CellClass.(celltypes{tt}))),...
        ISIStats.summstats.(state).rateCV2corr(CellClass.(celltypes{tt})),...
        '.','color',cellcolor{tt})
    hold on
end
%plot(

%% Constant Spike Count Bins

nspkintervals = 10;

cspkbinstats.meanISI = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.ISIs,'UniformOutput',false);
cspkbinstats.rate = cellfun(@(X) 1./X ,cspkbinstats.meanISI,'UniformOutput',false);
cspkbinstats.bindur = cellfun(@(X) movsum(X,nspkintervals),ISIStats.allspikes.ISIs,'UniformOutput',false);
cspkbinstats.stdISI = cellfun(@(X) movstd(X,nspkintervals),ISIStats.allspikes.ISIs,'UniformOutput',false);
cspkbinstats.CVISI = cellfun(@(X,Y) X./Y,cspkbinstats.stdISI,cspkbinstats.meanISI,'UniformOutput',false);
cspkbinstats.meanCV2 = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.CV2,'UniformOutput',false);
cspkbinstats.times = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.times,'UniformOutput',false);

cspkbinstats.summstats.rateCV2corr = cellfun(@(X,Y) corr(X,Y,'type','spearman'),...
    cspkbinstats.rate,cspkbinstats.meanCV2);
cspkbinstats.summstats.binderCV = cellfun(@(X) std((X))/mean((X)),...
    cspkbinstats.bindur);

instatebintimes = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    cspkbinstats.times,'UniformOutput',false);


cspkbinstats.cellhists.cv2bins = linspace(0,2,30);
cspkbinstats.cellhists.ratebins = linspace(-1.5,2.5,40);
numbinthresh = 25;

cspkbinstats.cellhists.ratehist = cellfun(@(X,Z) hist(log10(X(Z)),cspkbinstats.cellhists.ratebins),...
    cspkbinstats.rate,instatebintimes,...
    'UniformOutput',false);

cspkbinstats.cellhists.ratehist_norm = cellfun(@(X) X./sum(X),...
    cspkbinstats.cellhists.ratehist,'UniformOutput',false);
cspkbinstats.cellhists.ratehist_norm = cat(1,cspkbinstats.cellhists.ratehist_norm{:});

underthresh = cellfun(@(X) X<numbinthresh,cspkbinstats.cellhists.ratehist,...
    'UniformOutput',false);
for cc = 1:spikes.numcells
    cspkbinstats.cellhists.ratehist{cc}(underthresh{cc}) = nan;
end


cspkbinstats.cellhists.rateVcv2 = cellfun(@(X,Y,Z) hist3([log10(X(Z)),Y(Z)],...
    {cspkbinstats.cellhists.ratebins,cspkbinstats.cellhists.cv2bins}),...
    cspkbinstats.rate,cspkbinstats.meanCV2,instatebintimes,...
    'UniformOutput',false);

cspkbinstats.cellhists.rateVcv2 = cellfun(@(X,Y) bsxfun(@(x,y) x./y,X,Y'),...
    cspkbinstats.cellhists.rateVcv2,cspkbinstats.cellhists.ratehist,...
    'UniformOutput',false);

% cspkbinstats.cellhists.rateVcv2 = cellfun(@(X) X./sum(X(:)),...
%     cspkbinstats.cellhists.rateVcv2,...
%     'UniformOutput',false);

%Two normalizations needed : joint probability and conditional on rate
%probabiltiy....

for tt=1:length(celltypes)
    cspkbinstats.cellhists.meanhist.(celltypes{tt}) = nanmean(cat(3,cspkbinstats.cellhists.rateVcv2{CellClass.(celltypes{tt})}),3);
end

%%
figure
for tt = 1:length(celltypes)
    subplot(3,2,tt)
        imagesc(cspkbinstats.cellhists.ratebins,cspkbinstats.cellhists.cv2bins,cspkbinstats.cellhists.meanhist.(celltypes{tt})')
        hold on
        plot(get(gca,'xlim'),[1 1],'k')
        axis xy
        LogScale('x',10)
        colorbar
        caxis([0 0.18])
        title(celltypes{tt})
        xlabel(['Rate (',num2str(nspkintervals),'spk bins)']);
        ylabel('CV2')
end


subplot(2,2,3)
    imagesc(cspkbinstats.cellhists.ratebins,[0 spikes.numcells],...
        cspkbinstats.cellhists.ratehist_norm(ISIStats.sorts.ALL.ratebyclass,:))
    hold on
    plot(get(gca,'xlim'),sum(CellClass.pE).*[1 1],'r')
    LogScale('x',10)
    colorbar
    %caxis([0 0.18])
    xlabel(['Rate (',num2str(nspkintervals),'spk bins)']);
    ylabel('Cell (Sorted by Rate)')

NiceSave('MeanRateCV2',figfolder,baseName)


%%
figure

for tt = 1:length(celltypes)
subplot(3,3,1)
    plot(log10(ISIStats.summstats.ALL.meanrate(CellClass.(celltypes{tt}))),...
        cspkbinstats.summstats.rateCV2corr(CellClass.(celltypes{tt})),...
        '.','color',cellcolor{tt},'markersize',8)
    hold on
    box off
    axis tight
    plot(get(gca,'xlim'),[0 0],'k')
    LogScale('x',10)
    xlabel('Cell Rate (Hz)');ylabel('Rate-CV2 correlation')
subplot(3,3,2)
    plot(log10(ISIStats.summstats.ALL.meanrate(CellClass.(celltypes{tt}))),...
        cspkbinstats.summstats.binderCV(CellClass.(celltypes{tt})),...
        '.','color',cellcolor{tt})
    hold on
    xlabel('Cell Rate');ylabel('CV of bin duration')
    
subplot(3,3,3)
    plot((ISIStats.summstats.ALL.meanCV2(CellClass.(celltypes{tt}))),...
        cspkbinstats.summstats.rateCV2corr(CellClass.(celltypes{tt})),...
        '.','color',cellcolor{tt},'markersize',8)
    hold on
    box off 
    axis tight
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean CV2');ylabel('Rate-CV2 correlation')
end
NiceSave('BinRateCV2Correlation',figfolder,baseName)

%%
figure
subplot(2,2,1)
plot(log10(1./cspkbinstats.meanISI{excell}),cspkbinstats.CVISI{excell},'.')
xlabel('Mean ISI (s)')
ylabel('CV ISI')
LogScale('x',10)



%%
% twin = bz_RandomWindowInIntervals(SleepState.ints.NREMstate,60);
% excell = randsample(spikes.numcells,1);
% figure
% subplot(3,1,1)
% bz_MultiLFPPlot(lfp,'timewin',twin,'spikes',spikes,'plotcells',excell)
% % subplot(3,1,2)
% %     plot(ISIStats.allspikes.times{excell},log10(ISIStats.allspikes.ISIs{excell}),'.')
% %     xlim(twin)
% %     LogScale('y',10)
% subplot(3,1,2)
%     plot(cspkbinstats.times{excell},log10(cspkbinstats.rate{excell}),'o-')
%     axis tight
%     xlim(twin)
%     LogScale('y',10)
%     
% subplot(3,1,3)
%     plot(cspkbinstats.times{excell},(cspkbinstats.meanCV2{excell}),'o-')
%     axis tight
%     xlim(twin)
%     ylim([0 2])
    


%%
figure
subplot(2,2,1)
imagesc(cspkbinstats.cellhists.ratebins,cspkbinstats.cellhists.cv2bins,cspkbinstats.cellhists.rateVcv2{excell}')
hold on
plot(get(gca,'xlim'),[1 1],'k')
plot(log10(ISIStats.summstats.ALL.meanrate(excell)),ISIStats.summstats.ALL.meanCV2(excell),'r+')
axis xy
caxis([0 0.18])
LogScale('x',10)
xlabel('Rate (Hz)');ylabel('CV2')
title({['Cell: ',num2str(excell)],['Type: ',CellClass.label{excell}]})

subplot(2,2,2)
hist(log10(cspkbinstats.bindur{excell}))
LogScale('x',10)
xlabel('Bin Duration (s)')

subplot(2,2,3)
plot(log10(cspkbinstats.meanISI{excell}),cspkbinstats.meanCV2{excell},'.')
hold on
plot(get(gca,'xlim'),[1 1],'k')
xlabel('Mean ISI (s)')
ylabel('Mean CV2')
LogScale('x',10)
ylim([0 2])

subplot(2,2,4)
plot(log10(1./cspkbinstats.meanISI{excell}),cspkbinstats.meanCV2{excell},'.')
hold on
plot(get(gca,'xlim'),[1 1],'k')
plot(log10(ISIStats.summstats.ALL.meanrate(excell)),ISIStats.summstats.ALL.meanCV2(excell),'r+')
xlabel('Rate (Hz)')
ylabel('Mean CV2')
LogScale('x',10)
ylim([0 2])

NiceSave('ExCellRateCV2 ',figfolder,baseName)
end