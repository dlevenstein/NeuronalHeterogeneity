function [ ] = SynBlockerSpikeStatsAnalysis(basePath,figfolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%

%% DEV
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/MVData/pvch13_180719_143357';
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
%timeints.Eonly = [102 175].*60.*2/3;
timeints.Eonly_light = [92 102].*60.*2/3;
timeints.EI = [180 230].*60.*2/3;
%%
%Waveform classification
%CellClass = bz_CellClassification(basePath);

[~,depthID] = ismember(spikes.maxWaveformCh,sessionInfo.AnatGrps.Channels);
[~,depthsort] = sort(depthID);

depthlayer = cell(size(depthID));
depthlayer(depthsort(1:47)) = {'Sup'};
depthlayer(depthsort(48:end)) = {'Deep'};


%%
[ ISIstats ] = bz_ISIStats( spikes,'ints',timeints,...
    'figfolder',figfolder,'forceRedetect',true,...
    'savecellinfo',true,'basePath',basePath,'cellclass',depthlayer);

%%

%ISIstats.summstats.Eonly.meanrate(isnan(ISIstats.summstats.Eonly.meanrate))=0;
%ISIstats.summstats.dark.meanrate(isnan(ISIstats.summstats.dark.meanrate))=0;

ISIstats.summstats.lightchange = ISIstats.summstats.lightStim.meanrate./...
    ISIstats.summstats.dark.meanrate;
[~,ISIstats.sorts.lightchange] = sort(ISIstats.summstats.lightchange);

ISIstats.summstats.Esynchange = ISIstats.summstats.Eonly.meanrate./...
    ISIstats.summstats.dark.meanrate;
[~,ISIstats.sorts.Esynchange] = sort(ISIstats.summstats.Esynchange);

ISIstats.summstats.Esynchange_CV2 = ISIstats.summstats.Eonly.meanCV2-...
    ISIstats.summstats.dark.meanCV2;
[~,ISIstats.sorts.Esynchange_CV2] = sort(ISIstats.summstats.Esynchange_CV2);



%%
celltypes = unique(depthlayer);
cellcolors = {'k','r'};
for cc = 1:length(celltypes)
    layerclass.(celltypes{cc})= strcmp(depthlayer,celltypes{cc});

end

%% GET GS rate

ISIrate.dt = 0.002;
ISIrate.timestamps = [0:ISIrate.dt:max(cat(1,spikes.times{:}))]';

%Bug Fix two spikes same time
[~,ISIstats.allspikes.unique] = cellfun(@(X) unique(X),ISIstats.allspikes.times,'UniformOutput',false);

ISIrate.ISI = cellfun(@(X,Y,Z) interp1(X(Z),Y(Z),ISIrate.timestamps,'next'),...
    ISIstats.allspikes.times,ISIstats.allspikes.ISIs,ISIstats.allspikes.unique,'UniformOutput',false);
ISIrate.ISI = cat(2,ISIrate.ISI{:});

%%
states = fieldnames(timeints);
for ss = 1:length(states)
    state = states{ss};
%     ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
%         ISIStats.allspikes.times,'UniformOutput',false);
    ISIrate.instate = InIntervals(ISIrate.timestamps,timeints.(state));    
    

    %% Calculate occupancy statistics
    OccupancyStats.(state).median = nanmedian(ISIrate.ISI(ISIrate.instate,:));
    OccupancyStats.(state).GSrate = 1./OccupancyStats.(state).median;


    [~,OccupancyStats.sorts.(state).median] = sort(OccupancyStats.(state).median);
    
    
        for cl = 1:length(celltypes)
            inclasscells{cl} = layerclass.(celltypes{cl});
            
            numspikethresh = 50;
            enoughspikes = ISIstats.summstats.dark.numspikes>numspikethresh;
            ISIstats.meanISIhist.(celltypes{cl}).(state).log = nanmean(ISIstats.ISIhist.(state).log(inclasscells{cl}&enoughspikes,:),1);

            OccupancyStats.sorts.(state).(['median',celltypes{cl}]) = ...
                intersect(OccupancyStats.sorts.(state).median,...
                find(inclasscells{cl}),'stable');

                if cl==1
                    OccupancyStats.sorts.(state).(['median','byclass'])=[];
                   
                end
            OccupancyStats.sorts.(state).(['median','byclass']) = ...
                [OccupancyStats.sorts.(state).(['median','byclass']),...
               OccupancyStats. sorts.(state).(['median',celltypes{cl}])];

        end
    
end
%%
histcolors = flipud(gray);
figure
subplot(4,3,1)

imagesc(ISIstats.ISIhist.logbins,[1 length(OccupancyStats.sorts.dark.medianSup)],...
    ISIstats.ISIhist.dark.log(OccupancyStats.sorts.dark.medianSup,:))
hold on
plot(ISIstats.ISIhist.logbins,bz_NormToRange(ISIstats.meanISIhist.Sup.dark.log,0.3),'k','linewidth',1)
plot(log10(OccupancyStats.dark.median(OccupancyStats.sorts.dark.medianSup)),...
    [1:length(OccupancyStats.sorts.dark.medianSup)],'k.','markersize',2)
colormap(gca,histcolors)
caxis([0 0.15])
LogScale('x',10,'exp',true)
ylabel({'Sup Cells','Sort by GS Rate'})
axis xy
title('Baseline')

subplot(4,3,2)

imagesc(ISIstats.ISIhist.logbins,[1 length(OccupancyStats.sorts.dark.medianSup)],ISIstats.ISIhist.Eonly.log(OccupancyStats.sorts.dark.medianSup,:))
hold on
plot(ISIstats.ISIhist.logbins,bz_NormToRange(ISIstats.meanISIhist.Sup.Eonly.log,0.3),'k','linewidth',1)
colormap(gca,histcolors)
caxis([0 0.15])
LogScale('x',10,'exp',true)
axis xy
title('Synaptic Block')

subplot(4,3,4)

imagesc(ISIstats.ISIhist.logbins,[1 length(OccupancyStats.sorts.dark.medianDeep)],ISIstats.ISIhist.dark.log(OccupancyStats.sorts.dark.medianDeep,:))
hold on
plot(ISIstats.ISIhist.logbins,bz_NormToRange(ISIstats.meanISIhist.Deep.dark.log,0.3),'k','linewidth',1)

plot(log10(OccupancyStats.dark.median(OccupancyStats.sorts.dark.medianDeep)),...
    [1:length(OccupancyStats.sorts.dark.medianDeep)],'k.','markersize',2)
colormap(gca,histcolors)
LogScale('x',10,'exp',true)
ylabel({'Deep Cells','Sort by GS Rate'})
xlabel('ISI (s)')
caxis([0 0.15])
axis xy

subplot(4,3,5)

imagesc(ISIstats.ISIhist.logbins,[1 length(OccupancyStats.sorts.Eonly.medianDeep)],ISIstats.ISIhist.Eonly.log(OccupancyStats.sorts.dark.medianDeep,:))
colormap(gca,histcolors)
hold on
plot(ISIstats.ISIhist.logbins,bz_NormToRange(ISIstats.meanISIhist.Deep.Eonly.log,0.3),'k','linewidth',1)

caxis([0 0.15])
LogScale('x',10,'exp',true)
axis xy
xlabel('ISI (s)')


subplot(4,3,7)

imagesc(ISIstats.ISIhist.logbins,[1 length(ISIstats.sorts.dark.rateSup)],...
    ISIstats.ISIhist.dark.log(ISIstats.sorts.dark.rateSup,:))
hold on
plot(log10(1./ISIstats.summstats.dark.meanrate(ISIstats.sorts.dark.rateSup)),...
    [1:length(ISIstats.sorts.dark.rateSup)],'k.','markersize',2)
colormap(gca,histcolors)
caxis([0 0.15])
LogScale('x',10,'exp',true)
ylabel({'Sup Cells','Sort by Mean Rate'})
%axis xy

subplot(4,3,8)

imagesc(ISIstats.ISIhist.logbins,[1 length(ISIstats.sorts.dark.rateSup)],ISIstats.ISIhist.Eonly.log(ISIstats.sorts.dark.rateSup,:))
hold on

colormap(gca,histcolors)
caxis([0 0.15])

%axis xy

subplot(4,3,10)

imagesc(ISIstats.ISIhist.logbins,[1 length(ISIstats.sorts.dark.rateDeep)],...
    ISIstats.ISIhist.dark.log(ISIstats.sorts.dark.rateDeep,:))
hold on
plot(log10(1./ISIstats.summstats.dark.meanrate(ISIstats.sorts.dark.rateDeep)),...
    [1:length(ISIstats.sorts.dark.rateDeep)],'k.','markersize',2)
colormap(gca,histcolors)
caxis([0 0.15])
LogScale('x',10,'exp',true)
ylabel({'Deep Cells','Sort by Mean Rate'})
%axis xy

subplot(4,3,11)

imagesc(ISIstats.ISIhist.logbins,[1 length(ISIstats.sorts.dark.rateDeep)],...
    ISIstats.ISIhist.Eonly.log(ISIstats.sorts.dark.rateDeep,:))
colormap(gca,histcolors)
caxis([0 0.15])
%axis xy
NiceSave('ISIdists',figfolder,baseName)


%%
figure
subplot(2,2,1)
hold on
for cc = 1:length(celltypes)
plot(log10(ISIstats.summstats.dark.meanrate(layerclass.(celltypes{cc}))),...
    log10(ISIstats.summstats.lightStim.meanrate(layerclass.(celltypes{cc}))),'.','color',cellcolors{cc})
end
plot([-3 2],[-3 2],'k')
LogScale('xy',10)
xlabel('Spike Rate: Dark (Hz)');ylabel('Spike Rate: Light (Hz)')
title('Pre Treat')

subplot(2,2,2)
hold on
for cc = 1:length(celltypes)
plot((ISIstats.summstats.dark.meanCV2(layerclass.(celltypes{cc}))),...
    (ISIstats.summstats.lightStim.meanCV2(layerclass.(celltypes{cc}))),'.','color',cellcolors{cc})
end
plot([0 2],[0 2],'k')
xlabel('CV2: Dark');ylabel('CV2: Light')
title('Pre Treat')

subplot(2,2,3)
hold on
for cc = 1:length(celltypes)
plot(log10(ISIstats.summstats.Eonly.meanrate(layerclass.(celltypes{cc}))),...
    log10(ISIstats.summstats.Eonly_light.meanrate(layerclass.(celltypes{cc}))),'.','color',cellcolors{cc})
end
plot([-3 2],[-3 2],'k')
LogScale('xy',10)
xlabel('Spike Rate: Dark (Hz)');ylabel('Spike Rate: Light (Hz)')
title('E-Block')

subplot(2,2,4)
hold on
for cc = 1:length(celltypes)
plot((ISIstats.summstats.Eonly.meanCV2(layerclass.(celltypes{cc}))),...
    (ISIstats.summstats.Eonly_light.meanCV2(layerclass.(celltypes{cc}))),'.','color',cellcolors{cc})
end
plot([0 2],[0 2],'k')
xlabel('CV2: Dark');ylabel('CV2: Light')
title('E-Block')
%%
figure
subplot(3,3,1)
hold on
for cc = 1:length(celltypes)
plot(log10(ISIstats.summstats.dark.meanrate(layerclass.(celltypes{cc}))),...
    log10(ISIstats.summstats.Eonly.meanrate(layerclass.(celltypes{cc}))),'.','color',cellcolors{cc})
end
plot([-2 2],[-2 2],'k')
axis tight
LogScale('xy',10)
xlabel('Rate: Baseline (Hz)');ylabel('Rate: E Block (Hz)')


subplot(3,3,2)
hold on
for cc = 1:length(celltypes)
plot((ISIstats.summstats.dark.meanCV2(layerclass.(celltypes{cc}))),...
    (ISIstats.summstats.Eonly.meanCV2(layerclass.(celltypes{cc}))),'.','color',cellcolors{cc})
end
plot([0 2],[0 2],'k')
xlabel('CV2: Baseline');ylabel('CV2: E Block')
axis tight
%LogScale('xy',10)

subplot(3,3,3)
hold on
for cc = 1:length(celltypes)
plot(log10(ISIstats.summstats.dark.meanrate(layerclass.(celltypes{cc}))),...
    (ISIstats.summstats.Eonly.meanCV2(layerclass.(celltypes{cc}))),'.','color',cellcolors{cc})
end
plot([-1 2],[1 1],'k')
xlabel('Rate: Baseline');ylabel('CV2: E Block')
axis tight
LogScale('x',10)


subplot(3,3,4)
hold on
for cc = 1:length(celltypes)
plot(log10(ISIstats.summstats.dark.meanrate(layerclass.(celltypes{cc}))),...
    (ISIstats.summstats.dark.meanCV2(layerclass.(celltypes{cc}))),'.','color',cellcolors{cc})
end
axis tight
plot(get(gca,'xlim'),[1 1],'k')
xlabel('Rate: Baseline');ylabel('CV2: Baseline')

LogScale('x',10)
ylim([0 2])
  NiceSave('BaseStats',figfolder,baseName)


%%
figure
subplot(2,2,1)
hold on
for cc = 1:length(celltypes)
    plot(log10(ISIstats.summstats.dark.meanrate(layerclass.(celltypes{cc}))),...
        log10(ISIstats.summstats.Esynchange(layerclass.(celltypes{cc}))),'.','color',cellcolors{cc})
end
    plot(get(gca,'xlim'),[0 0],'k')
    LogScale('xy',10)
    xlabel('Mean Rate (Dark)');ylabel('Rate (E block)/Rate (Dark)')
subplot(3,3,3)
hold on
for cc = 1:length(celltypes)
    plot(log10(ISIstats.summstats.Esynchange(layerclass.(celltypes{cc}))),...
        ISIstats.summstats.Esynchange_CV2(layerclass.(celltypes{cc})),'.','color',cellcolors{cc})
end
    axis tight
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
   % LogScale('x',10)
    
    xlabel('Fold Rate Change');ylabel('CV2 Change')
    
subplot(2,2,3)
hold on
for cc = 1:length(celltypes)
    plot(log10(ISIstats.summstats.dark.meanrate(layerclass.(celltypes{cc}))),...
        ISIstats.summstats.Esynchange_CV2(layerclass.(celltypes{cc})),'.','color',cellcolors{cc})
end
    plot(get(gca,'xlim'),[0 0],'k')
    LogScale('x',10)
    xlabel('Rate (Dark)');ylabel('Change CV2')
    
subplot(2,2,4)
hold on
for cc = 1:length(celltypes)
    plot((ISIstats.summstats.dark.meanCV2(layerclass.(celltypes{cc}))),...
        ISIstats.summstats.Esynchange_CV2(layerclass.(celltypes{cc})),'.','color',cellcolors{cc})
end
    plot(get(gca,'xlim'),[0 0],'k')
    %LogScale('x',10)
    xlabel('CV2 (Dark)');ylabel('Change CV2')
    
  NiceSave('RelateCV2Rate',figfolder,baseName)

    
%% Everything by depth!
figure
subplot(3,4,1)
hold on
for cc = 1:length(celltypes)
    plot(log10(ISIstats.summstats.Esynchange(layerclass.(celltypes{cc}))),...
        -depthID(layerclass.(celltypes{cc})),'.','color',cellcolors{cc})
end
    plot([0 0],get(gca,'ylim'),'k')
    %LogScale('x',10)
    xlabel('Fold Rate Change')
    ylabel('Depth')
    axis tight
    box off
subplot(3,4,2)
hold on
for cc = 1:length(celltypes)
    plot((ISIstats.summstats.Esynchange_CV2(layerclass.(celltypes{cc}))),...
        -depthID(layerclass.(celltypes{cc})),'.','color',cellcolors{cc})
end
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('CV2 Change')
    ylabel('Depth')
    axis tight
    box off
 
    
subplot(3,4,3)
hold on
for cc = 1:length(celltypes)
    plot(log10(ISIstats.summstats.dark.meanrate(layerclass.(celltypes{cc}))),...
        -depthID(layerclass.(celltypes{cc})),'.','color',cellcolors{cc})
end
    %plot([0 0],get(gca,'ylim'),'k')
    xlabel('Rate')
    ylabel('Depth')
    axis tight
    box off
    
subplot(3,4,4)
hold on
for cc = 1:length(celltypes)
    plot((ISIstats.summstats.dark.meanCV2(layerclass.(celltypes{cc}))),...
        -depthID(layerclass.(celltypes{cc})),'.','color',cellcolors{cc})
end
    %plot([0 0],get(gca,'ylim'),'k')
    xlabel('CV2')
    ylabel('Depth')
    axis tight
    box off
    
    
subplot(3,4,5)
hold on
for cc = 1:length(celltypes)
    plot(log10(OccupancyStats.dark.GSrate(layerclass.(celltypes{cc}))),...
        log10(ISIstats.summstats.Esynchange(layerclass.(celltypes{cc}))),...
        '.','color',cellcolors{cc})
end
    plot(get(gca,'xlim'),[0 0],'k')
    %LogScale('x',10)
    ylabel('Fold Rate Change')
    xlabel('GS Rate')
    axis tight
    LogScale('x',10,'exp',true)
    box off
subplot(3,4,6)
hold on
for cc = 1:length(celltypes)
    plot(log10(OccupancyStats.dark.GSrate(layerclass.(celltypes{cc}))),...
        (ISIstats.summstats.Esynchange_CV2(layerclass.(celltypes{cc}))),...
        '.','color',cellcolors{cc})
end
    plot(get(gca,'xlim'),[0 0],'k')
    ylabel('CV2 Change')
    xlabel('GS Rate')
    
    axis tight
    LogScale('x',10,'exp',true)
    box off    
    
    
  NiceSave('ChangeByDepth',figfolder,baseName)

%%
spkwinsize = 25; %s
spkmat = bz_SpktToSpkmat(spikes,'binsize',spkwinsize,'dt',2.5);
%%
CV2mat.winsize = spkwinsize;
CV2mat.timestamps = spkmat.timestamps;
CV2mat.binedges = bsxfun(@(X,Y) X+Y,spkmat.timestamps,[-0.5 0.5].*CV2mat.winsize);
    [~,CV2mat.cells] = ...
        cellfun(@(X,Y) BinDataTimes(X,Y,CV2mat.binedges),...
        ISIstats.allspikes.CV2,ISIstats.allspikes.times,...
        'UniformOutput',false);
    
%% PCA for CV2 Patterns
%[COEFF, SCORE, LATENT] = pca(CV2mat.cells','Rows','pairwise');

%%
%%
CV2mat.cells = cat(2,CV2mat.cells{:});
%%

for cc = 1:length(celltypes)
CV2mat.(celltypes{cc}).mean = nanmean(CV2mat.cells(:,layerclass.(celltypes{cc})),2);
CV2mat.(celltypes{cc}).std = nanstd(CV2mat.cells(:,layerclass.(celltypes{cc})),[],2);
end

%%
spkmat.meanrate = mean(spkmat.data./spkwinsize,2);
spkmat.stdrate = std(spkmat.data./spkwinsize,[],2);
%%
excell = randsample(length(spikes.times),1);

excell = 7;

figure
subplot(3,1,1)
imagesc(spkmat.timestamps,[0 spikes.numcells],   log10(spkmat.data(:,depthsort)./spkwinsize)')
hold on
StateScorePlot( timeints,{'r','k','b','r','g'} )
%axis xy
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
hold on
for cc = 1:length(celltypes)
plot(CV2mat.timestamps,CV2mat.(celltypes{cc}).mean,cellcolors{cc},'linewidth',2)

errorshade(CV2mat.timestamps,CV2mat.(celltypes{cc}).mean,...
    CV2mat.(celltypes{cc}).std,CV2mat.(celltypes{cc}).std,cellcolors{cc},'scalar')
end
ylim([0 2])
%xlim([0 12000])
ylabel('Mean CV2')
xlabel('t (s)')


%% Entire Recording
ratecolormap = makeColorMap([1 1 1],[0.5 0.5 0.5],[0.8 0 0]);
cv2colormap = makeColorMap([0 0 0.8],[0.5 0.5 0.5],[0.8 0 0]);

whichsort = depthsort;
figure
subplot(3,1,1)
imagesc(spkmat.timestamps,[0 spikes.numcells],   log10(spkmat.data(:,whichsort)./spkwinsize)')
hold on
StateScorePlot( timeints,{'r','k','b','r','g'} )
%axis xy
xlabel('t (s)')
%xlim([0 12000])
ylabel('Cell (Sorted by Depth)')
colorbar('eastoutside')
colormap(gca,ratecolormap)
LogScale('c',10)

subplot(3,1,2)
h = imagesc(CV2mat.timestamps,[0 1],CV2mat.cells(:,whichsort)');
set(h, 'AlphaData', ~isnan(CV2mat.cells(:,whichsort)'));
colormap(gca,cv2colormap)
%xlim([0 12000])
caxis([0.5 1.5])
colorbar('eastoutside')
%axis xy
xlabel('t (s)')
ylabel('Cell (Sorted by Depth)')

NiceSave('CellsActivityByDepth',figfolder,baseName)


%% Onset
ratecolormap = makeColorMap([1 1 1],[0.5 0.5 0.5],[0.8 0 0]);
cv2colormap = makeColorMap([0 0 0.8],[0.5 0.5 0.5],[0.8 0 0]);

winshow = timeints.Eonly(1,1)+[-1000 3000];

whichsort = depthsort;
figure
subplot(3,1,1)
imagesc(spkmat.timestamps,[0 spikes.numcells],   log10(spkmat.data(:,whichsort)./spkwinsize)')
hold on
plot(timeints.Eonly(1,1).*[1 1],get(gca,'ylim'),'k--')
%StateScorePlot( timeints,{'r','k','b','r','g'} )
%axis xy
xlabel('t (s)')
box off
%xlim([0 12000])
ylabel('Cell (Sorted by Depth)')
%colorbar('eastoutside')
ColorbarWithAxis([-1.5 1.5],'Rate (Hz)')
colormap(gca,ratecolormap)
LogScale('c',10,'nohalf',true)
xlim(winshow);ylim([-3 spikes.numcells])
bz_ScaleBar('s')

subplot(3,1,2)
h = imagesc(CV2mat.timestamps,[0 spikes.numcells],CV2mat.cells(:,whichsort)');
set(h, 'AlphaData', ~isnan(CV2mat.cells(:,whichsort)'));
colormap(gca,cv2colormap)
hold on
plot(timeints.Eonly(1,1).*[1 1],get(gca,'ylim'),'k--')
%xlim([0 12000])
ColorbarWithAxis([0.4 1.6],'CV2')
%axis xy
box off
xlabel('t (s)')
ylabel('Cell (Sorted by Depth)')
xlim(winshow);ylim([-3 spikes.numcells])
bz_ScaleBar('s')

NiceSave('CellsActivityByDepth_Onset',figfolder,baseName)
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
