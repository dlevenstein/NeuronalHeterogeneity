reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/BehaviorTransitionsAnalysis'];
datasetPath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox';
[BehAnalysis,baseNames] = bz_LoadAnalysisResults(datasetPath,'BehaviorTransitionsAnalysis','dataset',true);
BehAnalysis = bz_CollapseStruct(BehAnalysis,3,'justcat',true);
%%

figure
subplot(2,2,1)
imagesc(mean(BehAnalysis.transprob,3))
set(gca,'ytick',[1 2 3 4]);set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabels',{'W','N','R','M'});
set(gca,'yticklabels',{'W','N','R','M'});
crameri bilbao
colorbar

subplot(2,2,2)
imagesc(mean(BehAnalysis.normprob,3))
set(gca,'ytick',[1 2 3 4]);set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabels',{'W','N','R','M'});
set(gca,'yticklabels',{'W','N','R','M'});
crameri('vik','pivot',1)
colorbar

subplot(2,2,3)
hold on
for ss = 1:3
    plot(durhist.bins,mean(BehAnalysis.durhist.(statenames{ss}),3),statecolors{ss})

end
        hold on
        axis tight
        plot(log10(MAthresh).*[1 1],ylim(gca),'r--')
        title(statenames{ss})
        LogScale('x',10)
        
        NiceSave('BehaviorTransition',figfolder,[]);
