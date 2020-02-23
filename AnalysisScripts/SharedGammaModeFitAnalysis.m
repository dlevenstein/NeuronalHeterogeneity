function [GammaFit] = SharedGammaModeFitAnalysis(basePath,figfolder)
% Date XX/XX/20XX
%
%Question:
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = pwd;
%basePath = '/Users/dlevenstein/Dropbox/research/Datasets/Cicero_09102014';
%basePath = '/Users/dlevenstein/Dropbox/research/Datasets/20140526_277um';
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
%ISIStats = bz_LoadCellinfo(basePath,'ISIStats');

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};

%%
statenames = {'NREMstate','WAKEstate','REMstate'};
%%

[ ISIstats ] = bz_ISIStats( spikes,'ints',SleepState.ints,...
    'cellclass',CellClass.label,'shuffleCV2',false,...
    'savecellinfo',false,'basePath',basePath,'forceRedetect',true,...
    'numISIbins',150,'logISIbounds',[0.0005 300]);


%%
numAS.NREMstate = 4;
numAS.WAKEstate = 4;
numAS.REMstate = 4;

%%

for ss = 1:3
   % numspks = sum(ISIstats.allspikes.instate.NREMstate.(statenames{ss}));
    logtimebins = ISIstats.ISIhist.logbins;
    logISIhist = ISIstats.ISIhist.(statenames{ss}).log(CellClass.pE,:)';
    logISIhist = logISIhist./mode(diff(logtimebins));
    GammaFit.(statenames{ss}) = bz_FitISISharedGammaModes(logISIhist,logtimebins,...
        'numAS',numAS.(statenames{ss}),...
        'figfolder',figfolder,'basePath',basePath,...
        'AScost',0.2,'ASguess',true);
end



GScolor = [0.6 0.4 0];

%% Example cell: 3 states
numex=1;
excell = randi(GammaFit.NREMstate.numcells,numex);
figure
for ss = 1:3
    %excell = excells(ee);
subplot(6,3,ss+3)
plot(GammaFit.(statenames{ss}).logtimebins,...
    GammaFit.(statenames{ss}).ISIdists(:,excell),...
    'color',[0.5 0.5 0.5],'linewidth',2)
hold on
plot(GammaFit.(statenames{ss}).logtimebins,...
    GSASmodel(GammaFit.(statenames{ss}).singlecell(excell),...
    GammaFit.(statenames{ss}).taubins),...
    'k','linewidth',2)
hold on
plot(GammaFit.(statenames{ss}).logtimebins,...
    LogGamma(GammaFit.(statenames{ss}).singlecell(excell).GSlogrates,...
    GammaFit.(statenames{ss}).singlecell(excell).GSCVs,...
    GammaFit.(statenames{ss}).singlecell(excell).GSweights',...
    GammaFit.(statenames{ss}).taubins'),'color',GScolor);
for aa = 1:numAS.(statenames{ss})
    plot(GammaFit.(statenames{ss}).logtimebins,...
        LogGamma(GammaFit.(statenames{ss}).singlecell(excell).ASlogrates(aa),...
        GammaFit.(statenames{ss}).singlecell(excell).ASCVs(aa),...
        GammaFit.(statenames{ss}).singlecell(excell).ASweights(aa)',...
        GammaFit.(statenames{ss}).taubins'),'k');
end
box off
axis tight
title(statenames{ss})

subplot(3,3,3+ss)
scatter(-GammaFit.(statenames{ss}).singlecell(excell).ASlogrates(:),...
    log10(GammaFit.(statenames{ss}).singlecell(excell).ASCVs(:)),...
    100*GammaFit.(statenames{ss}).singlecell(excell).ASweights(:)+0.00001,'k','filled')
hold on
scatter(-GammaFit.(statenames{ss}).singlecell(excell).GSlogrates,...
    log10(GammaFit.(statenames{ss}).singlecell(excell).GSCVs),...
    100*GammaFit.(statenames{ss}).singlecell(excell).GSweights+0.00001,GScolor,'filled')
ylabel('CV');xlabel('mean ISI')
xlim(logtimebins([1 end]))
LogScale('x',10)
box on

end
if figfolder
    NiceSave('CellExample_states',figfolder,baseName);
end
%% Mean and all points (single cell and group)
%Mean dist with group AS. use mean weight.

%% Mean rate and GS rate
% figure
% plot(sharedfit.GSlogrates,log10(ISIStats.summstats.WAKEstate.meanrate(ISIStats.sorts.WAKEstate.ratepE)),'.')
% hold on
% %UnityLine
% xlabel('GS Rate');
% ylabel('Mean Rate')









%%
%cc = 1
% Nmodes = 5;
% maxNmodes = 12;
% numcells = length(ISIStats.summstats.WAKEstate.meanrate);
% clear ISIfits
% for ss = 1:3
%     for cc = 1:numcells
%      
%     bz_Counter(cc,numcells,'Cell')
%     fitISIs = InIntervals(ISIStats.allspikes.times{cc},SleepState.ints.(statenames{ss}));
%     fitISIs = ISIStats.allspikes.ISIs{cc}(fitISIs);
%     [ISIfits.(statenames{ss}).lambdas(:,cc),ISIfits.(statenames{ss}).ks(:,cc),...
%         ISIfits.(statenames{ss}).weights(:,cc),ISIfits.(statenames{ss}).fiterror(cc,:),...
%         ISIfits.(statenames{ss}).Nmodes(cc)] = ...
%         bz_FitISIGammaModes(fitISIs,...
%         'showfig',false,'returnNmodes',Nmodes,'maxNmodes',maxNmodes,...
%         'sequentialreduce',true,'Nestimatemethod','descending');
%     end
% 
%     ISIfits.(statenames{ss}).rates = 1./(ISIfits.(statenames{ss}).ks./ISIfits.(statenames{ss}).lambdas);
%     ISIfits.(statenames{ss}).CVs = 1./ISIfits.(statenames{ss}).ks;
%     ISIfits.(statenames{ss}).weights(ISIfits.(statenames{ss}).weights<0.01) = nan;
%     ISIfits.(statenames{ss}).rates(isnan(ISIfits.(statenames{ss}).weights)) = nan;
%     
% end
% 
% 
% %% Example cell
% 
% cc = randi(numcells);
% for ss =1:3
% %
% fitISIs = InIntervals(ISIStats.allspikes.times{cc},SleepState.ints.(statenames{ss}));
% fitISIs = ISIStats.allspikes.ISIs{cc}(fitISIs);
% [~] = ...
%     bz_FitISIGammaModes(fitISIs,...
%     'showfig',true,'sequentialreduce',true,...
%     'maxNmodes',maxNmodes,'returnNmodes',Nmodes,'autoNmodes','LargeInflection',...
%     'Nestimatemethod','descending');
% 
%     NiceSave(['ISImodefits_ExCell_',num2str(cc),'_',(statenames{ss})],figfolder,baseName)
% end
% %%
% for ss = 1:3
% figure
% subplot(2,2,1)
% hold on
% for cc = 1:2
%     for mm = 1:Nmodes
%     scatter(log10(ISIfits.(statenames{ss}).rates(mm,CellClass.(celltypes{cc}))),...
%         ISIfits.(statenames{ss}).CVs(mm,CellClass.(celltypes{cc})),...
%         10.*ISIfits.(statenames{ss}).weights(mm,CellClass.(celltypes{cc})),cellcolor{cc})
%     end
% %     scatter(log10(1./(ks(:,CellClass.(celltypes{cc}))./rates(:,CellClass.(celltypes{cc})))),...
% %         (1./ks(:,CellClass.(celltypes{cc}))),weights(:,CellClass.(celltypes{cc})),...
% %         repmat(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc})),4,1))
%     LogScale('x',10)
%     xlabel('Rate');ylabel('CV')
% end
% title((statenames{ss}))
% 
% subplot(2,2,2)
% %imagesc(log10(fiterror))
% hold on
% for cc = 1:2
%     plot(1:maxNmodes,nanmean(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),1),...
%         '-o','linewidth',2,'color',cellcolor{cc})
%     errorshade(1:maxNmodes,nanmean(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),1),...
%         nanstd(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),[],1),...
%         nanstd(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),[],1),cellcolor{cc},'scalar');
% end
% LogScale('y',10)
% xlabel('N Modes');ylabel('Total Squared Error')
% 
% subplot(4,2,5)
% hold on
% for cc = 1:2
%     plot(log10(ISIfits.(statenames{ss}).rates(:,CellClass.(celltypes{cc}))),...
%         log10(ISIfits.(statenames{ss}).weights(:,CellClass.(celltypes{cc}))),'.','color',cellcolor{cc})
% 
% end
%     LogScale('xy',10)
%     xlabel('Rate');ylabel('Weight')
% 
% subplot(4,2,7)
% hold on
% for cc = 1:2
%     for mm = 1:Nmodes
%     scatter(log10(ISIfits.(statenames{ss}).rates(mm,CellClass.(celltypes{cc}))),...
%         log10(ISIStats.summstats.(statenames{ss}).meanrate(CellClass.(celltypes{cc}))),...
%         10.*ISIfits.(statenames{ss}).weights(mm,CellClass.(celltypes{cc})),cellcolor{cc})
%     end
%     
% %    plot(log10(ISIfits.(statenames{ss}).rates(:,CellClass.(celltypes{cc}))),...
% %        repmat(log10(ISIStats.summstats.(statenames{ss}).meanrate(CellClass.(celltypes{cc}))),Nmodes,1),'.','color',cellcolor{cc})
% 
% end
%     LogScale('xy',10)
%     UnityLine
%     ylim(log10([min(ISIStats.summstats.(statenames{ss}).meanrate) max(ISIStats.summstats.(statenames{ss}).meanrate)]))
%     xlabel('Mode Rate');ylabel(' Cell Rate')
% 
%     for cc = 1:2
% subplot(4,2,6+(cc-1)*2)
%     hist(ISIfits.(statenames{ss}).Nmodes(CellClass.(celltypes{cc})))
%     end
%     
%     NiceSave(['ISImodefits',(statenames{ss})],figfolder,baseName)
% 
% end
%%
% figure
%
% subplot(2,2,1)
% hold on
% for cc = 1
%     plot(log10(rates(:,CellClass.(celltypes{cc}))),...
%         repmat(log10(ISIStats.summstats.WAKEstate.meanrate(CellClass.(celltypes{cc}))),Nmodes,1),'.','color',cellcolor{cc})
% %     scatter(log10(1./(ks(:,CellClass.(celltypes{cc}))./rates(:,CellClass.(celltypes{cc})))),...
% %         (1./ks(:,CellClass.(celltypes{cc}))),weights(:,CellClass.(celltypes{cc})),...
% %         repmat(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc})),4,1))
% end
%     LogScale('xy',10)
%     UnityLine
%     %ylim(log10([min(ISIStats.summstats.NREMstate.meanrate) max(ISIStats.summstats.NREMstate.meanrate)]))
%     xlabel('Mode Rate');ylabel(' Cell Rate')
%
% title('CA1 - NREM')
%
% subplot(2,2,2)
% hold on
% for cc = 1
%     for mm = 1:Nmodes
% %     plot(log10(weights(:,CellClass.(celltypes{cc}))),...
% %         repmat(log10(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc}))),3,1),'.','color',cellcolor{cc})
%     scatter(log10(rates(mm,CellClass.(celltypes{cc}))),...
%         log10(ISIStats.summstats.WAKEstate.meanrate(CellClass.(celltypes{cc}))),...
%         5,20*weights(mm,CellClass.(celltypes{cc})))
%     end
% end
%     LogScale('xy',10)
%     %UnityLine
%     xlabel('Weight');ylabel('mean rate')

end
