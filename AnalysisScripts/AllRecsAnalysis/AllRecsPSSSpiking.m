reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/PSSCorticalStateAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC';
%regions = {'fCTX','CA1'};
%datasetPath = figfolder;

[PSSandSpikingAll,baseNames] = GetMatResults(figfolder,'PSSCorticalStateAnalysis','select',true);
PSSandSpikingAll = bz_CollapseStruct(PSSandSpikingAll);
thisregion = 'CA1';

CellClass = bz_LoadCellinfo(datasetPath.(thisregion),'CellClass','dataset',true,...
    'baseNames',baseNames,'catall',true);
ISIStats = bz_LoadCellinfo(datasetPath.(thisregion),'ISIStats','dataset',true,...
    'baseNames',baseNames,'catall',true);

%%
PSScellhist = bz_CollapseStruct(PSSandSpikingAll.PSScellhist,'match');
PSScellhist.meanYX(PSScellhist.meanYX==0) = nan;

%%
PSScellhist.allratebins.pE = linspace(-2.5,1,30);
PSScellhist.allratebins.pI = linspace(-1.5,2,30);

cellclasses = unique(CellClass.label);
for tt = 1:length(cellclasses) 
    [PSScellhist.allratedist.(cellclasses{tt})] = ...
        hist(squeeze(log10(PSScellhist.meanYX(:,:,CellClass.(cellclasses{tt}))))',...
        PSScellhist.allratebins.(cellclasses{tt}));
    PSScellhist.allratedist.(cellclasses{tt}) = bsxfun(@(X,Y) X./Y,...
        PSScellhist.allratedist.(cellclasses{tt}),sum(PSScellhist.allratedist.(cellclasses{tt}),1));
    PSScellhist.allmeanrate.(cellclasses{tt}) = ...
        nanmean(squeeze(log10(PSScellhist.meanYX(:,:,CellClass.(cellclasses{tt}))))',1);
    PSScellhist.allstdrate.(cellclasses{tt}) = ...
        nanstd(squeeze(log10(PSScellhist.meanYX(:,:,CellClass.(cellclasses{tt}))))',[],1);
end




%%
ratePSScorr = bz_CollapseStruct(PSSandSpikingAll.ratePSScorr,'match');

for tt = 1:length(cellclasses) 
    [ratePSScorr.popcorr.(cellclasses{tt}),ratePSScorr.popcorr.pval.(cellclasses{tt})]...
        = corr(log10(ISIStats.summstats.ALL.meanrate(CellClass.(cellclasses{tt})))',...
    ratePSScorr.ALL(CellClass.(cellclasses{tt}))','type','spearman');
end

CV2PSScorr = bz_CollapseStruct(PSSandSpikingAll.CV2PSScorr,'match');

for tt = 1:length(cellclasses) 
    [CV2PSScorr.popcorr.(cellclasses{tt}),CV2PSScorr.popcorr.pval.(cellclasses{tt})]...
        = corr(log10(ISIStats.summstats.ALL.meanrate(CellClass.(cellclasses{tt})))',...
    CV2PSScorr.ALL(CellClass.(cellclasses{tt}))','type','spearman');
end

%% PSS DIst
PSShist = bz_CollapseStruct(PSSandSpikingAll.PSShist,1);
PSShist.mean = bz_CollapseStruct(PSSandSpikingAll.PSShist,1,'mean');
PSShist.std = bz_CollapseStruct(PSSandSpikingAll.PSShist,1,'std');


states = {'WAKEstate','NREMstate','REMstate','ALL'};
statecolors = {'k','b','r',[0.6 0.6 0.6]};
%%
classcolors = {'k','r'};
cmap = [1 1 1;colormap(gca)];

figure
colormap(cmap)

%subplot(2,2,1)
% imagesc(squeeze(log10(PSScellhist.meanYX(:,:,ISIStats.sorts.NREMstate.ratebyclass)))')
subplot(3,3,1)
imagesc(PSScellhist.Xbins(:,:,1),PSScellhist.allratebins.pE,PSScellhist.allratedist.pE)
hold on
%plot(PSScellhist.Xbins(:,:,1),PSScellhist.allmeanrate.pE,'o-')
%plot(PSScellhist.Xbins(:,:,1),PSScellhist.allstdrate.pE,'o-')
LogScale('y',10)
axis xy

subplot(3,3,4)
imagesc(PSShist.mean.bins,[0 1],PSShist.ALL)


subplot(2,2,2)
hold on
for tt = 1:length(cellclasses) 
plot(log10(ISIStats.summstats.ALL.meanrate(CellClass.(cellclasses{tt}))),...
    ratePSScorr.ALL(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
end
axis tight
box off
xlim([-2.5 2])
plot(xlim(gca),[0 0],'k')
xlabel('sFR (Hz)');ylabel('Rate-PSS Corr')
LogScale('x',10)
title({['pE - rho:',num2str(round(ratePSScorr.popcorr.pE,2)),...
    '  p:',num2str(round(ratePSScorr.popcorr.pval.pE,2))],...
    ['pI - rho:',num2str(round(ratePSScorr.popcorr.pI,2)),...
    '  p:',num2str(round(ratePSScorr.popcorr.pval.pI,2))]})


%for ss = 1:4
subplot(6,4,24)
for tt = 1:length(cellclasses)
%         plot(PSScorrhist.bins,PSScorrhist.CV2.(states{ss}).(cellclasses{tt}),...
%             classcolors{tt},'linewidth',2)
        HistWithMean(CV2PSScorr.ALL(CellClass.(cellclasses{tt})),...
            'numbins',15,'color',classcolors{tt},'showtext',false)
    hold on
end
plot([0 0],get(gca,'ylim'),'k')
xlabel('PSS-CV2 Corr');
ylabel('# cells')
%xlim([-0.25 0.25])
%end





NiceSave(['RatebyPSS_',(thisregion)],figfolder,[])

%% CV/CV2 Stats
PSSpECV2hist = bz_CollapseStruct(PSSandSpikingAll.PSSpECV2hist,3,'mean');
PSSpICV2hist = bz_CollapseStruct(PSSandSpikingAll.PSSpICV2hist,3,'mean');
PSSpECVhist = bz_CollapseStruct(PSSandSpikingAll.PSSpECVhist,3,'mean');
PSSpICVhist = bz_CollapseStruct(PSSandSpikingAll.PSSpICVhist,3,'mean');
%%
figure
colormap(cmap)

subplot(10,4,[1 5])
imagesc(PSSpECV2hist.Xbins,PSSpECV2hist.Ybins,PSSpECV2hist.pYX')
axis xy
hold on
plot(PSSpECV2hist.Xbins,PSSpECV2hist.meanYX,'-k')
xlim([-1.6 -0.3])
ylabel({'CV_2', 'pE Pop.'})
ylim([0.6 1.5])
plot(get(gca,'xlim'),[1 1],'k--')
box off
   set(gca,'xticklabel',[])


subplot(10,4,[9 13])
imagesc(PSSpICV2hist.Xbins,PSSpICV2hist.Ybins,PSSpICV2hist.pYX')
axis xy
hold on
plot(PSSpICV2hist.Xbins,PSSpICV2hist.meanYX,'-k')
xlim([-1.6 -0.3])
ylim([0.5 1.25])
ylabel({'CV_2',' pI Pop.'})
plot(get(gca,'xlim'),[1 1],'k--')
box off



subplot(8,4,17)
    hold on
for ss = 1:3
   errorshade(PSShist.mean.bins,PSShist.mean.(states{ss}),...
       PSShist.std.(states{ss}),PSShist.std.(states{ss}),statecolors{ss},'scalar')
end
    for ss = 1:3
    plot(PSShist.mean.bins,PSShist.mean.(states{ss}),'color',statecolors{ss},'linewidth',2)
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
   set(gca,'ytick',[])
    %legend(states{:},'location','eastoutside')
    axis tight
    y = ylim(gca);
    ylim([0 y(2)])
    box off
    xlim([-1.6 -0.3])
    

subplot(10,4,[3 7])
imagesc(PSSpECVhist.Xbins,PSSpECVhist.Ybins,PSSpECVhist.pYX')
axis xy
hold on
plot(PSSpECVhist.Xbins,PSSpECVhist.meanYX','w')
%xlim([-1.6 -0.3])
%colorbar
box off
caxis([0 0.05])
ylabel({'Rate CV', 'pE Pop.'})
ylim([0.5 3.8])
%plot(get(gca,'xlim'),[1 1],'w--')
   set(gca,'xticklabel',[])


subplot(10,4,[11 15])
imagesc(PSSpICVhist.Xbins,PSSpICVhist.Ybins,PSSpICVhist.pYX')
axis xy
hold on
plot(PSSpICVhist.Xbins,PSSpICVhist.meanYX','w')
%xlim([-1.6 -0.3])
%ylim([0.5 1.6])
box off
ylim([0.25 2])
ylabel({'Rate CV',' pI Pop.'})
%plot(get(gca,'xlim'),[1 1],'w--')


subplot(8,4,19)
    hold on
for ss = 1:3
   errorshade(PSShist.mean.bins,PSShist.mean.(states{ss}),...
       PSShist.std.(states{ss}),PSShist.std.(states{ss}),statecolors{ss},'scalar')
end
    for ss = 1:3
    plot(PSShist.mean.bins,PSShist.mean.(states{ss}),'color',statecolors{ss},'linewidth',2)
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
    set(gca,'ytick',[])
    %legend(states{:},'location','eastoutside')
    axis tight
    y = ylim(gca);
    ylim([0 y(2)])
    box off
    xlim([-1.6 -0.3])

subplot(6,4,21)
for tt = 1:length(cellclasses)
%         plot(PSScorrhist.bins,PSScorrhist.CV2.(states{ss}).(cellclasses{tt}),...
%             classcolors{tt},'linewidth',2)
        HistWithMean(CV2PSScorr.ALL(CellClass.(cellclasses{tt})),...
            'numbins',15,'color',classcolors{tt},'showtext',false)
    hold on
end
plot([0 0],get(gca,'ylim'),'k')
xlabel('PSS-CV2 Corr');
ylabel('# cells')
%xlim([-0.25 0.25])
    
NiceSave(['CVbyPSS_',thisregion],figfolder,[])

%% ISI Stats by PSS

PSSISIhist  =  bz_CollapseStruct(PSSandSpikingAll.PSSISIhist,3);

for tt = 1:length(cellclasses)
     PSSISIhist.allcellmean.(cellclasses{tt}) = ...
         nanmean(PSSISIhist.pYX(:,:,CellClass.(cellclasses{tt})),3);
end

%%
figure

for tt = 1:length(cellclasses)
subplot(4,3,tt*3)
imagesc(PSSISIhist.Xbins(1,:,1),(PSSISIhist.Ybins(1,:,1)),PSSISIhist.allcellmean.(cellclasses{tt})')
axis xy
hold on
xlim([-1.5 -0.4])
ylabel([cellclasses{tt},' ISI (s)'])
%ylim([0 5])
LogScale('y',10)
if tt==1
caxis([0 0.05])
elseif tt==2
    caxis([0 0.08])
end
end



subplot(8,3,15)
    hold on
for ss = 1:3
   errorshade(PSShist.mean.bins,PSShist.mean.(states{ss}),...
       PSShist.std.(states{ss}),PSShist.std.(states{ss}),statecolors{ss},'scalar')
end
    for ss = 1:3
    plot(PSShist.mean.bins,PSShist.mean.(states{ss}),'color',statecolors{ss},'linewidth',2)
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
    set(gca,'ytick',[])
    
    %legend(states{:},'location','eastoutside')
    axis tight
    y = ylim(gca);
    ylim([0 y(2)])
    box off
    xlim([-1.5 -0.4])
    
    
    NiceSave(['ISIbyPSS_',thisregion],figfolder,[])
    