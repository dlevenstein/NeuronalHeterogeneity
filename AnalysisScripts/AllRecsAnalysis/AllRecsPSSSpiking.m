reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/PSSCorticalStateAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC';
%regions = {'fCTX','CA1'};
%datasetPath = figfolder;

[PSSandSpikingAll,baseNames] = GetMatResults(figfolder,'PSSCorticalStateAnalysis','select',true);
PSSandSpikingAll = bz_CollapseStruct(PSSandSpikingAll);
%%
CellClass = bz_LoadCellinfo(datasetPath.CA1,'CellClass','dataset',true,...
    'baseNames',baseNames,'catall',true);
ISIStats = bz_LoadCellinfo(datasetPath.CA1,'ISIStats','dataset',true,...
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
        PSScellhist.allratedist.(cellclasses{tt}),sum(PSScellhist.allratedist.(cellclasses{tt}),1))
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
%%
classcolors = {'k','r'};
figure
%subplot(2,2,1)
% imagesc(squeeze(log10(PSScellhist.meanYX(:,:,ISIStats.sorts.NREMstate.ratebyclass)))')
subplot(2,2,1)
imagesc(PSScellhist.Xbins(:,:,1),PSScellhist.allratebins.pE,PSScellhist.allratedist.pE)
hold on
%plot(PSScellhist.Xbins(:,:,1),PSScellhist.allmeanrate.pE,'o-')
%plot(PSScellhist.Xbins(:,:,1),PSScellhist.allstdrate.pE,'o-')
LogScale('y',10)
axis xy


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

NiceSave('RatebyPSS_HPC',figfolder,[])
