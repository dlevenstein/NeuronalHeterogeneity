function [ ] = AnalysisXXXXXXXX(basePath,figfolder)
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
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = pwd;
basePath = '/Users/dlevenstein/Dropbox/research/Datasets/Cicero_09102014';
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
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
%cc = 1
Nmodes = 4;
maxNmodes = 10;
logbins = ISIStats.ISIhist.logbins;
numcells = length(ISIStats.summstats.WAKEstate.meanrate);
clear lambdas ks weights fiterror
for cc = 1:numcells
bz_Counter(cc,numcells,'Cell')
logISIhist = ISIStats.ISIhist.WAKEstate.log(cc,:);
[lambdas(:,cc),ks(:,cc),weights(:,cc),fiterror(cc,:)] = ...
    bz_FitISIGammaModes(logbins,logISIhist,...
    'showfig',false,'returnNmodes',Nmodes,'maxNmodes',maxNmodes);
end

%% Example cell
cc = 17;
logISIhist = ISIStats.ISIhist.WAKEstate.log(cc,:);
[~] = ...
    bz_FitISIGammaModes(logbins,logISIhist,...
    'showfig',true,'returnNmodes',1);
%%

rates = 1./(ks./lambdas);
CVs = 1./ks;
rates(weights == 0 | isnan(weights)) = nan;
weights(weights == 0) = nan;

figure
subplot(2,2,1)
hold on
for cc = 1:2
    plot(log10(rates(:,CellClass.(celltypes{cc}))),...
        CVs(:,CellClass.(celltypes{cc})),'.','color',cellcolor{cc})
%     scatter(log10(1./(ks(:,CellClass.(celltypes{cc}))./rates(:,CellClass.(celltypes{cc})))),...
%         (1./ks(:,CellClass.(celltypes{cc}))),weights(:,CellClass.(celltypes{cc})),...
%         repmat(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc})),4,1))
    LogScale('x',10)
    xlabel('Rate');ylabel('CV')
end
title('CA1 - NREM')

subplot(2,2,2)
%imagesc(log10(fiterror))
hold on
for cc = 1:2
    plot(1:maxNmodes,mean(log10(fiterror(CellClass.(celltypes{cc}),:)),1),...
        '-o','linewidth',2,'color',cellcolor{cc})
    errorshade(1:maxNmodes,mean(log10(fiterror(CellClass.(celltypes{cc}),:)),1),...
        std(log10(fiterror(CellClass.(celltypes{cc}),:)),[],1),...
        std(log10(fiterror(CellClass.(celltypes{cc}),:)),[],1),cellcolor{cc},'scalar');
end
LogScale('y',10)
xlabel('N Modes');ylabel('MSE')

subplot(2,2,3)
hold on
for cc = 1:2
    plot(log10(rates(:,CellClass.(celltypes{cc}))),...
        log10(weights(:,CellClass.(celltypes{cc}))),'.','color',cellcolor{cc})
    LogScale('x',10)
    xlabel('Rate');ylabel('Weight')
end

subplot(4,2,6)
hold on
for cc = 1:2
    plot(log10(rates(:,CellClass.(celltypes{cc}))),...
        repmat(log10(ISIStats.summstats.WAKEstate.meanrate(CellClass.(celltypes{cc}))),Nmodes,1),'.','color',cellcolor{cc})
%     scatter(log10(1./(ks(:,CellClass.(celltypes{cc}))./rates(:,CellClass.(celltypes{cc})))),...
%         (1./ks(:,CellClass.(celltypes{cc}))),weights(:,CellClass.(celltypes{cc})),...
%         repmat(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc})),4,1))
end
    LogScale('xy',10)
    UnityLine
    ylim(log10([min(ISIStats.summstats.WAKEstate.meanrate) max(ISIStats.summstats.WAKEstate.meanrate)]))
    xlabel('Mode Rate');ylabel(' Cell Rate')


%%
figure

subplot(2,2,1)
hold on
for cc = 1
    plot(log10(rates(:,CellClass.(celltypes{cc}))),...
        repmat(log10(ISIStats.summstats.WAKEstate.meanrate(CellClass.(celltypes{cc}))),Nmodes,1),'.','color',cellcolor{cc})
%     scatter(log10(1./(ks(:,CellClass.(celltypes{cc}))./rates(:,CellClass.(celltypes{cc})))),...
%         (1./ks(:,CellClass.(celltypes{cc}))),weights(:,CellClass.(celltypes{cc})),...
%         repmat(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc})),4,1))
end
    LogScale('xy',10)
    UnityLine
    %ylim(log10([min(ISIStats.summstats.NREMstate.meanrate) max(ISIStats.summstats.NREMstate.meanrate)]))
    xlabel('Mode Rate');ylabel(' Cell Rate')

title('CA1 - NREM')

subplot(2,2,2)
hold on
for cc = 1
    for mm = 1:Nmodes
%     plot(log10(weights(:,CellClass.(celltypes{cc}))),...
%         repmat(log10(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc}))),3,1),'.','color',cellcolor{cc})
    scatter(log10(rates(mm,CellClass.(celltypes{cc}))),...
        log10(ISIStats.summstats.WAKEstate.meanrate(CellClass.(celltypes{cc}))),...
        5,20*weights(mm,CellClass.(celltypes{cc})))
    end
end
    LogScale('xy',10)
    %UnityLine
    xlabel('Weight');ylabel('mean rate')

    
%% Try fitting the ISIs directly...
excell = 1;
ISIs2fit = InIntervals(ISIStats.allspikes.times{excell},SleepState.ints.NREMstate);
ISIs2fit = ISIStats.allspikes.ISIs{excell}(ISIs2fit);

FitISIGammaModes2(ISIs2fit)

end