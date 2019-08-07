projectfolder = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
projectfolder = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';

basePath = [projectfolder,'Datasets/onProbox/AG_HPC/Cicero_09012014'];
figfolder = [projectfolder,'AnalysisScripts/AnalysisFigs/Illustrations'];
baseName = bz_BasenameFromBasepath(basePath);
%%
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');

%%


%%
cc =2;

shortthresh = 0.05;
longthresh = 2;
[shorthist,shortbins] = hist(ISIStats.allspikes.ISIs{cc}(ISIStats.allspikes.ISIs{cc}<=shortthresh),100);
shorthist = shorthist./length(ISIStats.allspikes.ISIs{cc});
[longhist,longbins] = hist(ISIStats.allspikes.ISIs{cc}(ISIStats.allspikes.ISIs{cc}<longthresh),200);
timehist = longhist.*longbins;
timehist = timehist./sum(ISIStats.allspikes.ISIs{cc});
longhist = longhist./length(ISIStats.allspikes.ISIs{cc});

figure
subplot(3,2,1)
bar(shortbins,shorthist)
box off
axis tight
xlabel('ISI (s)')
ylabel('P(ISI)')

subplot(3,4,3)
bar([sum(ISIStats.allspikes.ISIs{cc}<=shortthresh),sum(ISIStats.allspikes.ISIs{cc}>shortthresh)]./length(ISIStats.allspikes.ISIs{cc}))
box off
set(gca,'XTickLabels',{['< ',num2str(shortthresh.*1000),'ms'],['> ',num2str(shortthresh.*1000),'ms']})
ylabel('P(ISI)')

subplot(3,2,3)
bar(longbins,longhist)
box off
axis tight
xlabel('ISI (s)')
ylabel('P(ISI)')

subplot(3,4,7)
bar([sum(ISIStats.allspikes.ISIs{cc}<=longthresh),sum(ISIStats.allspikes.ISIs{cc}>longthresh)]./length(ISIStats.allspikes.ISIs{cc}))
box off
set(gca,'XTickLabels',{['< ',num2str(longthresh),'s'],['> ',num2str(longthresh),'s']})
ylabel('P(ISI)')

subplot(3,2,5)
bar(longbins,timehist)
box off
axis tight
xlabel('ISI (s)')
ylabel('P_t(ISI)')


subplot(3,2,6)
hist(log10(ISIStats.allspikes.ISIs{cc}),100)
xlabel('Log_1_0(ISI)')
box off
axis tight

NiceSave('IllustrateLogISI',figfolder,baseName);


%%

rate1ISI = exprnd(1,2500,1);
rate2ISI = exprnd(50,5000,1);

%%
linbins = linspace(0,300,100);
logbins = linspace(-2.5,3,100);
figure
subplot(3,2,1)
hist(rate1ISI,100)
%xlim(linbins([1 end]))
box off
title('2500 Poisson Spikes, Rate: 1')

subplot(3,2,3)
hist(rate2ISI,linbins)
xlim(linbins([1 end]))
box off
title('5000 Poisson Spikes, Rate: 50')



subplot(3,2,2)
hist(log10(rate1ISI),logbins)
xlim(logbins([1 end]))
box off
title('Poisson, Rate: 1')



subplot(3,2,4)
hist(log10(rate2ISI),logbins)
xlim(logbins([1 end]))
box off
title('Poisson, Rate: 50')



subplot(3,2,5)
hist([rate1ISI;rate2ISI],linbins)
xlim(linbins([1 end]))
box off


subplot(3,2,6)
hist(log10([rate1ISI;rate2ISI]),logbins)
xlim(logbins([1 end]))
box off

NiceSave('SimulateLogISI',figfolder,[]);
