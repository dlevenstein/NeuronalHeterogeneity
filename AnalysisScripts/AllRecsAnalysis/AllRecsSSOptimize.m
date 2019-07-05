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
dipmap = bz_CollapseStruct(OptimizeSleepScoreParmsAll.dipmap,3,'mean',true);

%%
dipmap.wins = 1:12;
dipmap.swins = 1:25;

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