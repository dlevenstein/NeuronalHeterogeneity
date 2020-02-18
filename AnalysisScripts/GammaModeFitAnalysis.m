function [ISIfits] = GammaModeFitAnalysis(basePath,figfolder)
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
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
%spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};

%%
statenames = {'WAKEstate','NREMstate','REMstate'};
%%
%cc = 1
Nmodes = 4;
maxNmodes = 10;
logbins = ISIStats.ISIhist.logbins;
numcells = length(ISIStats.summstats.WAKEstate.meanrate);
clear lambdas ks weights fiterror
for ss = 1:3
    for cc = 1:numcells
    bz_Counter(cc,numcells,'Cell')
    fitISIs = InIntervals(ISIStats.allspikes.times{cc},SleepState.ints.(statenames{ss}));
    fitISIs = ISIStats.allspikes.ISIs{cc}(fitISIs);
    [ISIfits.(statenames{ss}).lambdas(:,cc),ISIfits.(statenames{ss}).ks(:,cc),...
        ISIfits.(statenames{ss}).weights(:,cc),ISIfits.(statenames{ss}).fiterror(cc,:)] = ...
        bz_FitISIGammaModes(fitISIs,...
        'showfig',false,'returnNmodes',Nmodes,'maxNmodes',maxNmodes);
    end
    
    ISIfits.(statenames{ss}).rates = 1./(ISIfits.(statenames{ss}).ks./ISIfits.(statenames{ss}).lambdas);
    ISIfits.(statenames{ss}).CVs = 1./ISIfits.(statenames{ss}).ks;
    ISIfits.(statenames{ss}).rates(ISIfits.(statenames{ss}).weights == 0 | isnan(ISIfits.(statenames{ss}).weights)) = nan;
    ISIfits.(statenames{ss}).weights(ISIfits.(statenames{ss}).weights<0.01) = nan;
end

%% Example cell
cc = randi(numcells);
cc=6
fitISIs = InIntervals(ISIStats.allspikes.times{cc},SleepState.ints.WAKEstate);
fitISIs = ISIStats.allspikes.ISIs{cc}(fitISIs);
[~] = ...
    bz_FitISIGammaModes(fitISIs,...
    'showfig',true,'returnNmodes',Nmodes);
    NiceSave(['ISImodefits_ExCell_',(statenames{ss})],figfolder,baseName)

%%
for ss = 1:3
figure
subplot(2,2,1)
hold on
for cc = 1:2
    plot(log10(ISIfits.(statenames{ss}).rates(:,CellClass.(celltypes{cc}))),...
        ISIfits.(statenames{ss}).CVs(:,CellClass.(celltypes{cc})),'.','color',cellcolor{cc})
%     scatter(log10(1./(ks(:,CellClass.(celltypes{cc}))./rates(:,CellClass.(celltypes{cc})))),...
%         (1./ks(:,CellClass.(celltypes{cc}))),weights(:,CellClass.(celltypes{cc})),...
%         repmat(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc})),4,1))
    LogScale('x',10)
    xlabel('Rate');ylabel('CV')
end
title((statenames{ss}))

subplot(2,2,2)
%imagesc(log10(fiterror))
hold on
for cc = 1:2
    plot(1:maxNmodes,nanmean(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),1),...
        '-o','linewidth',2,'color',cellcolor{cc})
    errorshade(1:maxNmodes,nanmean(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),1),...
        nanstd(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),[],1),...
        nanstd(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),[],1),cellcolor{cc},'scalar');
end
LogScale('y',10)
xlabel('N Modes');ylabel('MSE')

subplot(2,2,3)
hold on
for cc = 1:2
    plot(log10(ISIfits.(statenames{ss}).rates(:,CellClass.(celltypes{cc}))),...
        log10(ISIfits.(statenames{ss}).weights(:,CellClass.(celltypes{cc}))),'.','color',cellcolor{cc})

end
    LogScale('xy',10)
    xlabel('Rate');ylabel('Weight')

subplot(4,2,6)
hold on
for cc = 1:2
    plot(log10(ISIfits.(statenames{ss}).rates(:,CellClass.(celltypes{cc}))),...
        repmat(log10(ISIStats.summstats.(statenames{ss}).meanrate(CellClass.(celltypes{cc}))),Nmodes,1),'.','color',cellcolor{cc})
%     scatter(log10(1./(ks(:,CellClass.(celltypes{cc}))./rates(:,CellClass.(celltypes{cc})))),...
%         (1./ks(:,CellClass.(celltypes{cc}))),weights(:,CellClass.(celltypes{cc})),...
%         repmat(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc})),4,1))
end
    LogScale('xy',10)
    UnityLine
    ylim(log10([min(ISIStats.summstats.(statenames{ss}).meanrate) max(ISIStats.summstats.(statenames{ss}).meanrate)]))
    xlabel('Mode Rate');ylabel(' Cell Rate')

    NiceSave(['ISImodefits',(statenames{ss})],figfolder,baseName)

end
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