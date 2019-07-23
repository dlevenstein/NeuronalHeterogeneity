%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/OptimizeSleepScoreParms'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC';
%regions = {'fCTX','CA1'};
%datasetPath = figfolder;

[OptimizeSleepScoreParmsAll,baseNames] = GetMatResults(figfolder,'OptimizeSleepScoreParms','select',false);
OptimizeSleepScoreParmsAll = bz_CollapseStruct(OptimizeSleepScoreParmsAll);
thisregion = 'vCTX';

%%
dipmap = bz_CollapseStruct(OptimizeSleepScoreParmsAll.dipmap,3,'median',true);
dipmap_all = bz_CollapseStruct(OptimizeSleepScoreParmsAll.dipmap,1,'justcat',true);

%%
dipmap.wins = 1:12;
dipmap.swins = 1:25;

%%
figure
subplot(3,2,1)
imagesc(dipmap_all.SWhist)
caxis([0 0.3])
colorbar
subplot(3,2,3)
imagesc(dipmap_all.SWhist_IRASA)
colorbar
caxis([0 0.3])
subplot(3,2,5)
plot(dipmap.bins,dipmap.SWhist)
hold on
plot(dipmap.bins,dipmap.SWhist_IRASA)

subplot(3,2,2)
imagesc(dipmap_all.THhist_used)
caxis([0 0.25])
colorbar
subplot(3,2,4)
imagesc(dipmap_all.THhist_IRASA_used)
colorbar
caxis([0 0.25])
subplot(3,2,6)
plot(dipmap.bins,dipmap.THhist_used)
hold on
plot(dipmap.bins,dipmap.THhist_IRASA_used)
%%

examples = [10 10; 2 15];

figure

subplot(4,4,1)
    imagesc(dipmap.wins,dipmap.swins,dipmap.SW')
    hold all
    plot(dipmap.wins(10),dipmap.swins(10),'k+')
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.06])
    title('Bimodality: Linear Fit')

subplot(4,4,2)
    imagesc(dipmap.wins,dipmap.swins,dipmap.SW_IRASA')
    hold all
        plot(dipmap.wins(2),dipmap.swins(15),'r+')
    
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.06])
    title('Bimodality: IRASA') 
    
subplot(4,4,3)
    imagesc(dipmap.wins,dipmap.swins,dipmap.TH')
    hold all
        plot(dipmap.wins(10),dipmap.swins(10),'k+')
    
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.0075])
    title('Bimodality: Theta')
    
subplot(4,4,4)
    imagesc(dipmap.wins,dipmap.swins,dipmap.TH_IRASA')
    hold all
        plot(dipmap.wins(2),dipmap.swins(15),'r+')
    
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.0075])
    title('Bimodality: Theta IRASA')


    
    
 NiceSave('WinIRASAOptimization',figfolder,[])