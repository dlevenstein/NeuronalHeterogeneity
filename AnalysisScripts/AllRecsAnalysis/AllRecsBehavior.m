reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/BehaviorTransitionsAnalysis'];
datasetPath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox';
[BehAnalysis,baseNames] = bz_LoadAnalysisResults(datasetPath,'BehaviorTransitionsAnalysis','dataset',true);
BehAnalysis = bz_CollapseStruct(BehAnalysis,3,'justcat',true);
%%
statenames = {'WAKEstate','NREMstate','REMstate','MA','LWAKE'};
statecolors = {'k','b','r','k'};

for ss = 1:5
    tottime.(statenames{ss}) = mean(BehAnalysis.tottime.(statenames{ss}),3);
end
% tottime.MA = mean(BehAnalysis.tottime.MA,3);
% tottime.LWAKE = mean(BehAnalysis.tottime.LWAKE,3);

%%
alltime = struct2cell(BehAnalysis.tottime);
alltime = cellfun(@squeeze,alltime,'UniformOutput',false);
%%
durhist.MAthresh = 180;

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

subplot(4,2,5)
hold on
for ss = 1:3
    plot(mean(BehAnalysis.durhist.bins,3),mean(BehAnalysis.durhist.(statenames{ss}),3),statecolors{ss},'linewidth',2)

end
        hold on
        axis tight
        plot(log10(durhist.MAthresh).*[1 1],ylim(gca),'r--','linewidth',2)
        LogScale('x',10,'nohalf',true)
        xlabel('Duration (s)')
    
colors = [0 0 0;0 0 1;1 0 0;0 0 0];
subplot(4,2,6)     
	BoxAndScatterPlot(alltime([end,3:end-1]),'colors',colors,'labels',{'WAKE','NREM','REM','MA'})   
    ylim([0 1])
    ylabel('P[time]')
        NiceSave('BehaviorTransition',figfolder,[]);

%%