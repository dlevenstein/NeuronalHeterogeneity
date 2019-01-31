function [ PSShist,ratePSScorr,CV2PSScorr,...
    PSSpECV2hist,PSSpICV2hist,PSSpECVhist,PSSpICVhist,PSSEIhist,PSScellhist,...
    PSSpEpopratehist,PSSpIpopratehist,PSSpEsynchhist,PSSpIsynchhist] = PSSCorticalStateAnalysis( basePath,figfolder )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
repoRoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity'; %desktop
%repoRoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/NyuShare/Buzsakilabspace/Datasets/GrosmarkAD/Cicero/Cicero_09102014';
basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/FRHetAndDynamics/AnalysisScripts/AnalysisFigs';
%figfolder = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/PSSCorticalStateAnalysis';
figfolder = [repoRoot,'/AnalysisScripts/AnalysisFigs/PSSCorticalStateAnalysis'];
%%
baseName = bz_BasenameFromBasepath(basePath);


spikes = bz_GetSpikes('basePath',basePath);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath);

states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

%% Pick the "best" window for further analysis
dt = 0.5;
winsize = 8;
specslope = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',true,...
    'saveMat',basePath);


numbins = 40;
PSShist.bins = linspace(-2,0,numbins);
for ss = 1:length(states)
    specslope.timeidx.(states{ss}) = InIntervals(specslope.timestamps,SleepState.ints.(states{ss}));
    
    PSShist.(states{ss}) = hist(specslope.data(specslope.timeidx.(states{ss})),PSShist.bins);
end


%%
rescolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
exwinsize = 2000;
exwin = bz_RandomWindowInIntervals(specslope.timestamps([1 end])',exwinsize);
figure
subplot(4,1,1)
    imagesc(specslope.timestamps,log2(specslope.freqs),specslope.specgram)
    ylabel({'Specgram','f (Hz)'})
    axis xy
    StateScorePlot({SleepState.ints.NREMstate,SleepState.ints.REMstate,SleepState.ints.WAKEstate},...
        {'b','r','k'})
    LogScale('y',2)
    hold on
    yyaxis right
    plot(specslope.timestamps,specslope.data,'w','linewidth',0.5)
    ylabel('PSS')
    xlim(exwin)
    
subplot(4,1,2)
    imagesc(specslope.timestamps,log2(specslope.freqs),specslope.resid')
    %colorbar
    colormap(gca,rescolormap)
    caxis(1.*[-1 1])
    ylabel({'Resid.','f (Hz)'})
    axis xy
    yyaxis right
    plot(specslope.timestamps,specslope.rsq,'k')
    xlim(exwin)
    LogScale('y',2)
    ylabel('Rsq')

subplot(6,2,11)
    for ss = 1:length(states)
    plot(PSShist.bins,PSShist.(states{ss}),'color',statecolors{ss},'linewidth',2)
    hold on
    end
    xlabel('PSS')
    legend(states{:},'location','eastoutside')
    axis tight
    box off
    xlim([-1.75 -0.25])
    
NiceSave('PSSbyState',figfolder,baseName)
%% Relate PSS and Spiking

ISIStats.allspikes.PSS = cellfun(@(X) interp1(specslope.timestamps,specslope.data,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);
spkwinsize = winsize;
overlap = spkwinsize./dt;
spikemat = bz_SpktToSpkmat(spikes,'binsize',spkwinsize,'overlap',overlap);
spikemat.PSS = interp1(specslope.timestamps,specslope.data,spikemat.timestamps);

cellclasses = {'pE','pI'};
classcolors = {'k','r'};
for tt = 1:length(cellclasses)
    spikemat.poprate.(cellclasses{tt}) = ...
        sum(spikemat.data(:,CellClass.(cellclasses{tt})),2)./winsize./sum(CellClass.(cellclasses{tt}));
    spikemat.popCV.(cellclasses{tt}) = ...
        std(spikemat.data(:,CellClass.(cellclasses{tt})),[],2)./mean(spikemat.data(:,CellClass.(cellclasses{tt})),2);
    spikemat.popvar.(cellclasses{tt}) = ...
        std(spikemat.data(:,CellClass.(cellclasses{tt})),[],2);
end
spikemat.poprate.EI = (spikemat.poprate.pE-spikemat.poprate.pI)./...
    (spikemat.poprate.pE+spikemat.poprate.pI);
%%
ratePSScorr.ALL = corr(spikemat.PSS,spikemat.data,'type','spearman','rows','complete');
CV2PSScorr.ALL = cellfun(@(X,Y) corr(X,Y,'type','spearman','rows','complete'),...
    ISIStats.allspikes.PSS,ISIStats.allspikes.CV2);
%% Mean binned CV2...
clear CV2mat
CV2mat.winsize = spkwinsize;
CV2mat.timestamps = spikemat.timestamps;
CV2mat.binedges = bsxfun(@(X,Y) X+Y,spikemat.timestamps,[-0.5 0.5].*CV2mat.winsize);
for tt = 1:length(cellclasses)
    allspikes.CV2.(cellclasses{tt}) = cat(1,ISIStats.allspikes.CV2{CellClass.(cellclasses{tt})});
    allspikes.times.(cellclasses{tt}) = cat(1,ISIStats.allspikes.times{CellClass.(cellclasses{tt})});
    [CV2mat.timestamps,CV2mat.(cellclasses{tt})] = ...
        BinDataTimes(allspikes.CV2.(cellclasses{tt}),allspikes.times.(cellclasses{tt}),CV2mat.binedges);
    CV2mat.rate.(cellclasses{tt}) = interp1(spikemat.timestamps,spikemat.poprate.(cellclasses{tt}),CV2mat.timestamps);
end
CV2mat.PSS = interp1(specslope.timestamps,specslope.data,CV2mat.timestamps);

%%
figure
subplot(2,2,1)
plot(CV2mat.PSS,CV2mat.pI,'r.')
xlabel('PSS');ylabel('Pop CV2')

subplot(2,2,2)
plot(CV2mat.PSS,CV2mat.pE,'k.')
xlabel('PSS');ylabel('Pop CV2')

subplot(2,2,3)
plot(log10(CV2mat.rate.pI),CV2mat.pI,'r.')
xlabel('Pop Rate');ylabel('Pop CV2')

subplot(2,2,4)
plot(log10(CV2mat.rate.pE),(CV2mat.pE),'k.')
xlabel('Pop Rate');ylabel('Pop CV2')


%%
exwinsize = 100;
exwin = bz_RandomWindowInIntervals(specslope.timestamps([1 end])',exwinsize);
figure
for tt = 1:length(cellclasses)
    subplot(4,1,tt)
    plot(allspikes.times.(cellclasses{tt}),allspikes.CV2.(cellclasses{tt}),'.','color',classcolors{tt})
    hold on
    plot(CV2mat.timestamps,CV2mat.(cellclasses{tt}),'color',classcolors{tt},'linewidth',2)
    xlim(exwin)
end

%%
nspklim = 20;

PSScorrhist.bins = linspace(-1,1,40);
for ss = 1:length(states)
%ss = 1;
    spikemat.timeidx.(states{ss}) = InIntervals(spikemat.timestamps,SleepState.ints.(states{ss}));
    CV2mat.timeidx.(states{ss}) = InIntervals(CV2mat.timestamps,SleepState.ints.(states{ss}));
    instatespikes = cellfun(@(X) InIntervals(X,SleepState.ints.(states{ss})),ISIStats.allspikes.times,'UniformOutput',false);

    toofewspikes = cellfun(@sum,instatespikes) < nspklim;
    instatespikes(toofewspikes) = {1}; %To account for 0 spikes

    ratePSScorr.(states{ss}) = corr(spikemat.PSS(spikemat.timeidx.(states{ss})),...
        spikemat.data(spikemat.timeidx.(states{ss}),:),'type','spearman','rows','complete');
    CV2PSScorr.(states{ss}) = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman','rows','complete'),...
        ISIStats.allspikes.PSS,ISIStats.allspikes.CV2,instatespikes);
    
    CV2PSScorr.(states{ss})(toofewspikes) = nan;
    
%     for tt = 1:length(cellclasses)
%         PSScorrhist.rate.(states{ss}).(cellclasses{tt}) = ...
%             hist(ratePSScorr.(states{ss})(CellClass.(cellclasses{tt})),PSScorrhist.bins);
%         PSScorrhist.CV2.(states{ss}).(cellclasses{tt}) = ...
%             hist(CV2PSScorr.(states{ss})(CellClass.(cellclasses{tt})),PSScorrhist.bins);
%     end
 
end    
%% Figure: PSS and Spiking
    figure
for ss = 1:length(states)
    subplot(4,4,ss)
    for tt = 1:length(cellclasses)
    	plot(ISIStats.summstats.(states{ss}).meanCV2(CellClass.(cellclasses{tt})),ratePSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    axis tight
    plot(get(gca,'xlim'),[0 0],'k')
    title(states{ss})
    %ylim([-0.3 0.3])
    xlabel('<CV2>');ylabel('Rate-PSS Corr')
    

    subplot(4,4,ss+4)
    for tt = 1:length(cellclasses)
    	plot(ISIStats.summstats.(states{ss}).meanCV2(CellClass.(cellclasses{tt})),CV2PSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    axis tight
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('<CV2>');ylabel('CV2-PSS Corr')
    

    subplot(4,4,ss+8)
    for tt = 1:length(cellclasses)
    	plot(log10(ISIStats.summstats.(states{ss}).meanrate(CellClass.(cellclasses{tt}))),ratePSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    axis tight
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate');ylabel('Rate-PSS Corr')
    
    subplot(4,4,ss+12)
    for tt = 1:length(cellclasses)
    	plot(log10(ISIStats.summstats.(states{ss}).meanrate(CellClass.(cellclasses{tt}))),CV2PSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    axis tight
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate');ylabel('CV2-PSS Corr')
end

NiceSave('PSSandCells',figfolder,baseName)


%%
figure
subplot(4,4,1)
 for ss = 1:3 
%         ScatterWithLinFit(spikemat.PSS(spikemat.timeidx.(states{ss})),...
%             log2(spikemat.poprate.pE(spikemat.timeidx.(states{ss}))),...
%             'color',statecolors{ss},'markersize',1,'showtext',false)
        plot(spikemat.PSS(spikemat.timeidx.(states{ss})),...
            log2(spikemat.poprate.pE(spikemat.timeidx.(states{ss}))),...
            '.','color',statecolors{ss},'markersize',1)
        hold on
 end
 axis tight
 box off
        xlabel('PSS');ylabel('pE Rate (spk/cell/s)')
       xlim([-1.8 -0.2])
        LogScale('y',2)
        
subplot(4,4,2)
 for ss = 1:3
%         ScatterWithLinFit(spikemat.PSS(spikemat.timeidx.(states{ss})),...
%             log2(spikemat.poprate.pI(spikemat.timeidx.(states{ss}))),...
%             'color',statecolors{ss},'markersize',1,'showtext',false)
        plot(spikemat.PSS(spikemat.timeidx.(states{ss})),...
            log2(spikemat.poprate.pI(spikemat.timeidx.(states{ss}))),...
            '.','color',statecolors{ss},'markersize',1)        
        hold on
 end
        axis tight
        box off
        xlabel('PSS');ylabel('pI Rate (spk/cell/s)')
        xlim([-1.8 -0.2])
        LogScale('y',2)
        
subplot(4,4,3)
 for ss = 1:3  
%         ScatterWithLinFit(spikemat.PSS(spikemat.timeidx.(states{ss})),...
%             log10(spikemat.poprate.pE(spikemat.timeidx.(states{ss}))./spikemat.poprate.pI(spikemat.timeidx.(states{ss}))),...
%             'color',statecolors{ss},'markersize',1,'showtext',false)
        plot(spikemat.PSS(spikemat.timeidx.(states{ss})),...
            log2(spikemat.poprate.pE(spikemat.timeidx.(states{ss}))./spikemat.poprate.pI(spikemat.timeidx.(states{ss}))),...
            '.','color',statecolors{ss},'markersize',1)
        hold on
 end
 axis tight
 box off
        xlabel('PSS');ylabel('E/I SpkRatio')
        xlim([-1.8 -0.2])
        LogScale('y',2)
        
for ss = 1:4
    subplot(6,4,ss*4)
    for tt = 1:length(cellclasses)
%         plot(PSScorrhist.bins,PSScorrhist.rate.(states{ss}).(cellclasses{tt}),...
%             classcolors{tt},'linewidth',2)
            HistWithMean(ratePSScorr.(states{ss})(CellClass.(cellclasses{tt})),...
                'numbins',8,'color',classcolors{tt},'showtext',false)
        hold on
    end
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('PSS-Rate Corr');
    ylabel('# cells')
    title(states{ss})
    xlim([-0.4 0.4])
end

for ss = 1:4
    subplot(6,4,ss*4+5)
    for tt = 1:length(cellclasses)
%         plot(PSScorrhist.bins,PSScorrhist.CV2.(states{ss}).(cellclasses{tt}),...
%             classcolors{tt},'linewidth',2)
            HistWithMean(CV2PSScorr.(states{ss})(CellClass.(cellclasses{tt})),...
                'numbins',8,'color',classcolors{tt},'showtext',false)
        hold on
    end
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('PSS-CV2 Corr');
    ylabel('# cells')
    xlim([-0.25 0.25])
end

for tt = 1:length(cellclasses)
    subplot(4,4,tt+13)
    for ss = 1:3
    	plot(CV2mat.PSS(CV2mat.timeidx.(states{ss})),CV2mat.(cellclasses{tt})(CV2mat.timeidx.(states{ss})),...
            '.','color',statecolors{ss},'markersize',2)
        hold on
    end
    axis tight
    box off
    xlabel('PSS');ylabel(['<CV2>, ',cellclasses{tt}])

end

subplot(4,4,16)
for tt = 1:length(cellclasses)
    plot(log2(CV2mat.rate.(cellclasses{tt})),CV2mat.(cellclasses{tt}),...
        '.','color',classcolors{tt},'markersize',1)
    hold on
end
box off
axis tight
LogScale('x',2)
xlabel('Pop. Rate');ylabel('<CV2>')

NiceSave('PSSandSpiking',figfolder,baseName,'tiff')


%%
figure
plot(log10(ISIStats.summstats.ALL.meanrate),ratePSScorr.ALL,'.')
hold on
plot(xlim(gca),[0 0],'k')
%%
%note good window: 1500 2500
exwinsize = 800;
exwin = bz_RandomWindowInIntervals(specslope.timestamps([1 end])',exwinsize);
figure

subplot(6,1,1)
    imagesc(specslope.timestamps,log2(specslope.freqs),specslope.specgram)
    hold on
    StateScorePlot({SleepState.ints.NREMstate,SleepState.ints.REMstate,SleepState.ints.WAKEstate},...
        {'b','r','k'})
    axis xy
    xlim(exwin)
    ylabel({'Specgram','f (Hz)'})
    LogScale('y',2)
        bz_ScaleBar('s')

    yyaxis right
    plot(specslope.timestamps,specslope.data,'w','linewidth',1)
    ylabel('PSS')

    
    
subplot(6,1,[2 3])
    bz_MultiLFPPlot(lfp,'spikes',spikes,'timewin',exwin,...
        'cellgroups',{CellClass.pE,CellClass.pI},...
        'sortmetric',ISIStats.summstats.ALL.meanrate,...
        'scaleLFP',0.4,'scalespikes',0.1,'spikeside','bottom')
    xlabel('')
    box off
    bz_ScaleBar('s')

subplot(6,1,4)
for tt = 1:length(cellclasses)
    plot(spikemat.timestamps,spikemat.popCV.(cellclasses{tt}),classcolors{tt},'linewidth',1)
    hold on
end
    plot(exwin,[1 1],'k--')
    axis tight
    %ylim([0.65 1.35])
    xlim(exwin)
    box off
    ylabel('CV')
    
% subplot(6,1,5)
% for tt = 1:length(cellclasses)
%     plot(spikemat.timestamps,log2(spikemat.poprate.(cellclasses{tt})),classcolors{tt},'linewidth',1)
%     hold on
% end
%     axis tight
%     xlim(exwin)
%     box off
%     ylabel({'Pop. Rate', '(Spk/Cell/S)'})
%     legend(cellclasses{:})
%     LogScale('y',2)
%     
subplot(3,1,3)
for tt = 1:length(cellclasses)
    plot(CV2mat.timestamps,CV2mat.(cellclasses{tt}),classcolors{tt},'linewidth',1)
    hold on
end
    plot(exwin,[1 1],'k--')
    axis tight
    ylim([0.65 1.35])
    xlim(exwin)
    box off
    ylabel('<CV2>')
    
    NiceSave('PSSandSpikingExample',figfolder,baseName,'tiff')

%%
NREMex = randsample(Restrict(CV2mat.timestamps,SleepState.ints.NREMstate),1);
WAKEex = randsample(Restrict(CV2mat.timestamps,SleepState.ints.WAKEstate),1);

figure
subplot(3,3,4)
    bz_MultiLFPPlot(lfp,'spikes',spikes,'timewin',NREMex+specslope.detectionparms.winsize.*[-0.2 0.2],...
        'cellgroups',{CellClass.pE,CellClass.pI},...
        'sortmetric',ISIStats.summstats.ALL.meanrate,...
        'scaleLFP',0.4,'scalespikes',0.1,'spikeside','bottom')
    xlabel('')
    box off
    bz_ScaleBar('s')
    title({['PSS: ',num2str(round(CV2mat.PSS(CV2mat.timestamps==NREMex),3))],...
        ['CV2pE: ',num2str(round(CV2mat.pE(CV2mat.timestamps==NREMex),3))],...
        ['CV2pI: ',num2str(round(CV2mat.pI(CV2mat.timestamps==NREMex),3))]})
    
subplot(3,3,1)
    bz_MultiLFPPlot(lfp,'spikes',spikes,'timewin',WAKEex+specslope.detectionparms.winsize.*[-0.2 0.2],...
        'cellgroups',{CellClass.pE,CellClass.pI},...
        'sortmetric',ISIStats.summstats.ALL.meanrate,...
        'scaleLFP',0.4,'scalespikes',0.1,'spikeside','bottom')
    xlabel('')
    box off
    bz_ScaleBar('s')
    title({['PSS: ',num2str(round(CV2mat.PSS(CV2mat.timestamps==WAKEex),3))],...
        ['CV2pE: ',num2str(round(CV2mat.pE(CV2mat.timestamps==WAKEex),3))],...
        ['CV2pI: ',num2str(round(CV2mat.pI(CV2mat.timestamps==WAKEex),3))]})
    
    NiceSave('PSSandSpikingExample2',figfolder,baseName,'tiff')

%% Heatmaps
numXbins = 60;
numYbins = 150;

Xbounds = [-2 0];
Ybounds = [0 2];

minx = 100;
mincells = 4;

[ PSSpECV2hist ] = ConditionalHist(CV2mat.PSS,CV2mat.pE,...
    'numXbins',60,'numYbins',200,'Xbounds',[-2 0],'Ybounds',[0 2],...
    'minX',minx);
[ PSSpICV2hist ] = ConditionalHist(CV2mat.PSS,CV2mat.pI,...
    'numXbins',60,'numYbins',200,'Xbounds',[-2 0],'Ybounds',[0 2],...
    'minX',minx);

[ PSSpECVhist ] = ConditionalHist(spikemat.PSS,spikemat.popCV.pE,...
    'numXbins',60,'numYbins',150,'Xbounds',[-2 0],'Ybounds',[0 4],...
    'minX',minx);
[ PSSpICVhist ] = ConditionalHist(spikemat.PSS,spikemat.popCV.pI,...
    'numXbins',60,'numYbins',150,'Xbounds',[-2 0],'Ybounds',[0 4],...
    'minX',minx);
if sum(CellClass.pI)<mincells
    PSSpICVhist.pYX = nan(size(PSSpICVhist.pYX));
end

[ PSSpEpophist ] = ConditionalHist(spikemat.PSS,spikemat.poprate.pE,...
    'numXbins',60,'numYbins',150,'Xbounds',[-2 0],'Ybounds',[0 4],...
    'minX',minx);
[ PSSpIpophist ] = ConditionalHist(spikemat.PSS,spikemat.poprate.pI,...
    'numXbins',60,'numYbins',150,'Xbounds',[-2 0],'Ybounds',[0 25],...
    'minX',minx);

% [ PSSpEhist ] = ConditionalHist(spikemat.PSS,(spikemat.data(:,CellClass.pE)./winsize),...
%     'numXbins',60,'numYbins',50,'Xbounds',[-2 0],'Ybounds',[]);
% [ PSSpIhist ] = ConditionalHist(spikemat.PSS,(spikemat.data(:,CellClass.pI)./winsize),...
%     'numXbins',60,'numYbins',50,'Xbounds',[-2 0],'Ybounds',[]);

[ PSSEIhist ] = ConditionalHist(spikemat.PSS,spikemat.poprate.EI,...
    'numXbins',60,'numYbins',150,'Xbounds',[-2 0],'Ybounds',[-1 1],...
    'minX',minx);

[ PSScellhist ] = ConditionalHist(spikemat.PSS,spikemat.data./spikemat.binsize,...
    'numXbins',60,'numYbins',50,'Xbounds',[-2 0],'Ybounds',[],...
    'minX',minx);
%%
for tt = 1:length(cellclasses) 
    [PSScellhist.ratedist.(cellclasses{tt}),PSScellhist.ratebins.(cellclasses{tt})] = ...
        hist(squeeze(log10(PSScellhist.meanYX(:,:,CellClass.(cellclasses{tt}))))');
    PSScellhist.meanrate.(cellclasses{tt}) = ...
        mean(squeeze(log10(PSScellhist.meanYX(:,:,CellClass.(cellclasses{tt}))))',1);
    PSScellhist.stdrate.(cellclasses{tt}) = ...
        std(squeeze(log10(PSScellhist.meanYX(:,:,CellClass.(cellclasses{tt}))))',[],1);
end
%spikemat.poprate.EI

%%
figure
subplot(2,2,1)
imagesc(squeeze(log10(PSScellhist.meanYX(:,:,ISIStats.sorts.NREMstate.ratebyclass)))')
subplot(2,2,2)
imagesc(PSScellhist.Xbins(:,:,1),PSScellhist.ratebins.pE,PSScellhist.ratedist.pE)
hold on
plot(PSScellhist.Xbins(:,:,1),PSScellhist.meanrate.pE,'o-')
plot(PSScellhist.Xbins(:,:,1),PSScellhist.stdrate.pE,'o-')
LogScale('y',10)
axis xy
%% Figure

figure
cmap = [1 1 1;colormap(gca)];
colormap(cmap)

subplot(10,4,[1 5])
imagesc(PSSpECV2hist.Xbins,PSSpECV2hist.Ybins,PSSpECV2hist.pYX')
axis xy
hold on
plot(PSSpECV2hist.Xbins,PSSpECV2hist.meanYX,'-k')
xlim([-1.6 -0.3])
ylabel({'CV_2', 'pE Pop.'})
ylim([0.9 1.4])
plot(get(gca,'xlim'),[1 1],'k--')
box off

subplot(10,4,[9 13])
imagesc(PSSpICV2hist.Xbins,PSSpICV2hist.Ybins,PSSpICV2hist.pYX')
axis xy
hold on
plot(PSSpICV2hist.Xbins,PSSpICV2hist.meanYX,'-k')
xlim([-1.6 -0.3])
ylim([0.75 1.15])
ylabel({'CV_2',' pI Pop.'})
plot(get(gca,'xlim'),[1 1],'k--')
box off

subplot(8,4,17)
    for ss = 1:3
    plot(PSShist.bins,PSShist.(states{ss}),'color',statecolors{ss},'linewidth',2)
    hold on
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
    set(gca,'ytick',[])
    %legend(states{:},'location','eastoutside')
    axis tight
    box off
    xlim([-1.6 -0.3])
    
subplot(10,4,[3 7])
imagesc(PSSpECVhist.Xbins,PSSpECVhist.Ybins,PSSpECVhist.pYX')
axis xy
hold on
plot(PSSpECVhist.Xbins,PSSpECVhist.meanYX,'-k')
xlim([-1.6 -0.3])
box off
ylabel({'Rate CV', 'pE Pop.'})
ylim([0.5 3.5])
%plot(get(gca,'xlim'),[1 1],'w--')

subplot(10,4,[11 15])
imagesc(PSSpICVhist.Xbins,PSSpICVhist.Ybins,PSSpICVhist.pYX')
axis xy
hold on
plot(PSSpICVhist.Xbins,PSSpICVhist.meanYX,'-k')
box off
xlim([-1.6 -0.3])
ylim([0.5 1.5])
ylabel({'Rate CV',' pI Pop.'})
%plot(get(gca,'xlim'),[1 1],'w--')
    
 
subplot(8,4,19)
    for ss = 1:3
    plot(PSShist.bins,PSShist.(states{ss}),'color',statecolors{ss},'linewidth',2)
    hold on
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
    set(gca,'ytick',[])
    %legend(states{:},'location','eastoutside')
    axis tight
    box off
    xlim([-1.6 -0.3])
    
subplot(10,4,[2 6])
imagesc(PSSpEpophist.Xbins,PSSpEpophist.Ybins,PSSpEpophist.pYX')
axis xy
hold on
plot(PSSpEpophist.Xbins,PSSpEpophist.meanYX,'-k')
xlim([-1.6 -0.3])
ylabel('pE Pop Rate')
ylim([0 3])
box off

subplot(10,4,[10 14])
imagesc(PSSpIpophist.Xbins,PSSpIpophist.Ybins,PSSpIpophist.pYX')
axis xy
hold on
plot(PSSpIpophist.Xbins,PSSpIpophist.meanYX,'-k')
xlim([-1.6 -0.3])
ylim([0.5 25])
ylabel('pI Pop Rate')
box off

subplot(8,4,18)
    for ss = 1:3
    plot(PSShist.bins,PSShist.(states{ss}),'color',statecolors{ss},'linewidth',2)
    hold on
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
    set(gca,'ytick',[])
    %legend(states{:},'location','eastoutside')
    axis tight
    box off
    xlim([-1.6 -0.3])
    
% subplot(5,4,15)
% imagesc(PSSpEhist.Xbins,PSSpEhist.Ybins,log10(PSSpEhist.pYX)')
% axis xy
% hold on
% xlim([-1.6 -0.3])
% ylabel({'Cell FR', 'pE Pop.'})
% %ylim([0.5 1.6])
% plot(get(gca,'xlim'),[1 1],'w--')
% 
% subplot(5,4,19)
% imagesc(PSSpIhist.Xbins,PSSpIhist.Ybins,log10(PSSpIhist.pYX)')
% axis xy
% hold on
% xlim([-1.6 -0.3])
% %ylim([0.5 1.6])
% ylabel({'Cell FR',' pI Pop.'})
% plot(get(gca,'xlim'),[1 1],'w--')

subplot(5,4,20)
imagesc(PSSEIhist.Xbins,PSSEIhist.Ybins,PSSEIhist.pYX')
axis xy
hold on
xlim([-1.6 -0.3])
%ylim([0.5 1.6])
ylabel('E-I Ratio')
plot(get(gca,'xlim'),[1 1],'w--')

subplot(6,4,21)
    for tt = 1:length(cellclasses)
%         plot(PSScorrhist.bins,PSScorrhist.CV2.(states{ss}).(cellclasses{tt}),...
%             classcolors{tt},'linewidth',2)
            HistWithMean(CV2PSScorr.ALL(CellClass.(cellclasses{tt})),...
                'numbins',8,'color',classcolors{tt},'showtext',false)
        hold on
    end
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('PSS-CV2 Corr');
    ylabel('# cells')
    xlim([-0.25 0.25])
    %legend('pE cells','pI ce
    
    NiceSave('CV2byPSSstats',figfolder,baseName,'tiff')

    
    
%% Faster Time Scale Spiking 


spikemat_fast = bz_SpktToSpkmat(spikes,'binsize',0.08,'overlap',4);
spikemat_fast.PSS = interp1(specslope.timestamps,specslope.data,spikemat_fast.timestamps);

cellclasses = {'pE','pI'};
classcolors = {'k','r'};
for tt = 1:length(cellclasses)
    spikemat_fast.popspikes.(cellclasses{tt}) = ...
        sum(spikemat_fast.data(:,CellClass.(cellclasses{tt})),2);
    spikemat_fast.popsynch.(cellclasses{tt}) = ...
        sum(spikemat_fast.data(:,CellClass.(cellclasses{tt}))>0,2);

end

%% Fast Time Scale
minx = 500;
[ PSSpEpopratehist ] = ConditionalHist(spikemat_fast.PSS,spikemat_fast.popspikes.pE,...
    'numXbins',80,'numYbins',sum(CellClass.pE)+1,'Xbounds',[-2 0],'Ybounds',[-0.5 sum(CellClass.pE)+0.5],...
    'minX',minx);
[ PSSpIpopratehist ] = ConditionalHist(spikemat_fast.PSS,spikemat_fast.popspikes.pI,...
    'numXbins',80,'numYbins',4.*sum(CellClass.pI)+1,'Xbounds',[-2 0],'Ybounds',[-0.5 4.*sum(CellClass.pI)+0.5],...
    'minX',minx);
PSSpEpopratehist.Ybins = PSSpEpopratehist.Ybins./sum(CellClass.pE)./spikemat_fast.binsize;
PSSpIpopratehist.Ybins = PSSpIpopratehist.Ybins./sum(CellClass.pI)./spikemat_fast.binsize;
PSSpEpopratehist.meanYX = PSSpEpopratehist.meanYX./sum(CellClass.pE)./spikemat_fast.binsize;
PSSpIpopratehist.meanYX = PSSpIpopratehist.meanYX./sum(CellClass.pI)./spikemat_fast.binsize;


[ PSSpEsynchhist ] = ConditionalHist(spikemat_fast.PSS,spikemat_fast.popsynch.pE,...
    'numXbins',80,'numYbins',sum(CellClass.pE)+1,'Xbounds',[-2 0],'Ybounds',[-0.5 sum(CellClass.pE)+0.5],...
    'minX',minx);
[ PSSpIsynchhist ] = ConditionalHist(spikemat_fast.PSS,spikemat_fast.popsynch.pI,...
    'numXbins',80,'numYbins',sum(CellClass.pI)+1,'Xbounds',[-2 0],'Ybounds',[-0.5 sum(CellClass.pI)+0.5],...
    'minX',minx);
PSSpEsynchhist.Ybins = PSSpEsynchhist.Ybins./sum(CellClass.pE);
PSSpIsynchhist.Ybins = PSSpIsynchhist.Ybins./sum(CellClass.pI);
PSSpEsynchhist.meanYX = PSSpEsynchhist.meanYX./sum(CellClass.pE);
PSSpIsynchhist.meanYX = PSSpIsynchhist.meanYX./sum(CellClass.pI);

% [ PSSpEhist ] = ConditionalHist(spikemat_fast.PSS,(spikemat_fast.data(:,CellClass.pE)),...
%     'numXbins',100,'numYbins',16,'Xbounds',[-2 0],'Ybounds',[0 15]);
% [ PSSpIhist ] = ConditionalHist(spikemat_fast.PSS,(spikemat_fast.data(:,CellClass.pI)),...
%     'numXbins',100,'numYbins',16,'Xbounds',[-2 0],'Ybounds',[0 15]);
% PSSpEhist.Ybins = PSSpEhist.Ybins./spikemat_fast.binsize;
% PSSpIhist.Ybins = PSSpIhist.Ybins./spikemat_fast.binsize;



%% Figure: fast time scale rate

figure
cmap = [1 1 1;colormap(gca)];

colormap(cmap)
subplot(10,4,[1 5])
imagesc(PSSpEpopratehist.Xbins,PSSpEpopratehist.Ybins,PSSpEpopratehist.pYX')
axis xy
hold on
plot(PSSpEpopratehist.Xbins,PSSpEpopratehist.meanYX,'-k')
xlim([-1.6 -0.3])
ylabel({'pE Pop Rate', ['(',num2str(spikemat_fast.binsize*1000),'ms bins)']})
ylim([0 5])
box off
   set(gca,'xticklabel',[])


subplot(10,4,[9 13])
imagesc(PSSpIpopratehist.Xbins,PSSpIpopratehist.Ybins,PSSpIpopratehist.pYX')
axis xy
hold on
plot(PSSpIpopratehist.Xbins,PSSpIpopratehist.meanYX,'-k')
xlim([-1.6 -0.3])
ylim([0 40])
box off
ylabel({'pI Pop Rate', ['(',num2str(spikemat_fast.binsize*1000),'ms bins)']})

subplot(8,4,17)
    for ss = 1:3
    plot(PSShist.bins,PSShist.(states{ss}),'color',statecolors{ss},'linewidth',2)
    hold on
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
    set(gca,'ytick',[])
    %legend(states{:},'location','eastoutside')
    axis tight
    box off
    xlim([-1.6 -0.3])
    
    
    
subplot(10,4,[3 7])
imagesc(PSSpEsynchhist.Xbins,PSSpEsynchhist.Ybins,PSSpEsynchhist.pYX')
axis xy
hold on
plot(PSSpEsynchhist.Xbins,PSSpEsynchhist.meanYX,'-k')
xlim([-1.6 -0.3])
ylabel({'pE Synch', ['(',num2str(spikemat_fast.binsize*1000),'ms bins)']})
ylim([0 0.5])

subplot(10,4,[11 15])
imagesc(PSSpIsynchhist.Xbins,PSSpIsynchhist.Ybins,PSSpIsynchhist.pYX')
axis xy
hold on
plot(PSSpIsynchhist.Xbins,PSSpIsynchhist.meanYX,'-k')
xlim([-1.6 -0.3])
ylim([0 1])
ylabel({'pI Synch', ['(',num2str(spikemat_fast.binsize*1000),'ms bins)']})

subplot(8,4,19)
    for ss = 1:3
    plot(PSShist.bins,PSShist.(states{ss}),'color',statecolors{ss},'linewidth',2)
    hold on
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
    set(gca,'ytick',[])
    %legend(states{:},'location','eastoutside')
    axis tight
    box off
    xlim([-1.6 -0.3])
    
    
%     
%  subplot(5,4,3)
% imagesc(PSSpEhist.Xbins,PSSpEhist.Ybins,log10(PSSpEhist.pYX)')
% axis xy
% hold on
% xlim([-1.6 -0.3])
% ylabel({'pE Cell Rate', ['(',num2str(spikemat_fast.binsize*1000),'ms bins)']})
% %ylim([0.5 1.6])
% 
% subplot(5,4,7)
% imagesc(PSSpIhist.Xbins,PSSpIhist.Ybins,log10(PSSpIhist.pYX)')
% axis xy
% hold on
% xlim([-1.6 -0.3])
% %ylim([0.5 1.6])
% ylabel({'pI Cell Rate', ['(',num2str(spikemat_fast.binsize*1000),'ms bins)']})
% 
% subplot(10,4,19)
%     for ss = 1:3
%     plot(PSShist.bins,PSShist.(states{ss}),'color',statecolors{ss},'linewidth',2)
%     hold on
%     end
%     xlabel('PSS')
%     ylabel({'Time', 'Occupancy'})
%     set(gca,'ytick',[])
%     %legend(states{:},'location','eastoutside')
%     axis tight
%     box off
%     xlim([-1.6 -0.3])

    NiceSave('FastRatebyPSSstats',figfolder,baseName,'tiff')
    
%% ISI Stats
[ PSSISIhist_allspkE ] = ConditionalHist(cat(1,ISIStats.allspikes.PSS{CellClass.pE}),...
    log10(cat(1,ISIStats.allspikes.ISIs{CellClass.pE})),...
    'numXbins',100,'numYbins',50,'Xbounds',[-2 0],'Ybounds',[-3 2]);
[ PSSISIhist_allspkI ] = ConditionalHist(cat(1,ISIStats.allspikes.PSS{CellClass.pI}),...
    log10(cat(1,ISIStats.allspikes.ISIs{CellClass.pI})),...
    'numXbins',100,'numYbins',50,'Xbounds',[-2 0],'Ybounds',[-3 2]);


PSSISIhist  = ConditionalHist(ISIStats.allspikes.PSS,...
    cellfun(@(X) log10(X),ISIStats.allspikes.ISIs,'UniformOutput',false),...
    'numXbins',100,'numYbins',50,'Xbounds',[-2 0],'Ybounds',[-3 2],'minX',25);

for tt = 1:length(cellclasses)
     PSSISIhist.mean.(cellclasses{tt}) = ...
         nanmean(PSSISIhist.pYX(:,:,CellClass.(cellclasses{tt})),3);
end
%%
figure
subplot(5,4,1)
imagesc(PSSISIhist_allspkE.Xbins,(PSSISIhist_allspkE.Ybins),PSSISIhist_allspkE.pYX')
axis xy
hold on
xlim([-1.6 -0.3])
ylabel('all pE ISIs')
%ylim([0 5])
LogScale('y',10)


subplot(5,4,5)
imagesc(PSSISIhist_allspkI.Xbins,(PSSISIhist_allspkI.Ybins),PSSISIhist_allspkI.pYX')
axis xy
hold on
xlim([-1.6 -0.3])
ylabel('all pE ISIs')%ylim([0 5])
LogScale('y',10)



subplot(10,4,17)
    for ss = 1:3
    plot(PSShist.bins,PSShist.(states{ss}),'color',statecolors{ss},'linewidth',2)
    hold on
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
    set(gca,'ytick',[])
    %legend(states{:},'location','eastoutside')
    axis tight
    box off
    xlim([-1.6 -0.3])
    
    
subplot(5,4,2)
imagesc(PSSISIhist_allspkE.Xbins,(PSSISIhist_allspkE.Ybins),PSSISIhist.mean.pE')
axis xy
hold on
xlim([-1.6 -0.3])
ylabel('all pE ISIs')
%ylim([0 5])
LogScale('y',10)


subplot(5,4,6)
imagesc(PSSISIhist_allspkI.Xbins,(PSSISIhist_allspkI.Ybins),PSSISIhist.mean.pI')
axis xy
hold on
xlim([-1.6 -0.3])
ylabel('all pE ISIs')%ylim([0 5])
LogScale('y',10)



subplot(10,4,18)
    for ss = 1:3
    plot(PSShist.bins,PSShist.(states{ss}),'color',statecolors{ss},'linewidth',2)
    hold on
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
    set(gca,'ytick',[])
    %legend(states{:},'location','eastoutside')
    axis tight
    box off
    xlim([-1.6 -0.3])

% 
