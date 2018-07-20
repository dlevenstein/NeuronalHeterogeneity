function [ ] = SpikeStatsbyPopActivityAnalysis(basePath,figfolder)

%% DEV
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/Dino_mPFC/Dino_061814';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyPopActivityAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
SleepState.ints.ALL = [0 Inf];
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath);

%% Cell types and states
[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};
statenames = fieldnames(SleepState.ints);

%% Calculate spike count matrix
binsize = 0.25; %s
binsize = 0.5; %s
overlap = 10;
spikemat = bz_SpktToSpkmat(spikes,'binsize',binsize,'overlap',overlap);


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
caxis([0 0.15])
title(celltypes{tt})
xlabel(['Spike Count, (',num2str(binsize),'s bins']);
ylabel('CV2')
end

%%
figure
for tt = 1:length(celltypes)
subplot(2,2,tt)
imagesc(cellhists.logISIbins,cellhists.cv2bins,cellhists.meanISIvCV2hist.(celltypes{tt})')
axis xy
colorbar
caxis([0 0.15])
xlabel('meanISI');ylabel('CV2')
LogScale('x',10)
title(celltypes{tt})
end

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
end