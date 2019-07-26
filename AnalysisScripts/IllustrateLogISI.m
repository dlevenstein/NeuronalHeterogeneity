basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC/Cicero_09012014';
figfolder = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/Illustrations';
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