function [  ] = AutoCorrelationAnalysis(basePath,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/Dino_mPFC/Dino_061814';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyPopActivityAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
%SleepState.ints.ALL = [0 Inf];
statenames = fieldnames(SleepState.ints);
numstates = length(statenames);
statecolors = {'k','b','r'};
classnames = unique(CellClass.label);
numclasses = length(classnames);
classcolors = {'k','r'};
% lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
%     'basepath',basePath);

%%
sessionInfo = bz_getSessionInfo(basePath);

%% Separate out E and I population spikes
popspikes.Espikes = sort(cat(1,spikes.times{CellClass.pE}));
popspikes.Ispikes = sort(cat(1,spikes.times{CellClass.pI}));


%% Calculate the CCGs (and thus, ACGs)
ccg.binsize = 0.001;
[ccg.hists,ccg.t] = CCG([spikes.times,popspikes.Espikes,popspikes.Ispikes],...
    [],'binSize',ccg.binsize,'duration',10,'norm','rate');

%% ACGs
ccg.popccgs = ccg.hists(:,[end-1:end],[end-1:end]);
ccg.popccgs(:,1,1) = ccg.popccgs(:,1,1)./sum(CellClass.pE);
ccg.popccgs(:,2,2) = ccg.popccgs(:,2,2)./sum(CellClass.pI);

ccg.popccgs(:,1,2) = ccg.popccgs(:,1,2)./sum(CellClass.pI);
ccg.popccgs(:,2,1) = ccg.popccgs(:,2,1)./sum(CellClass.pE);

acgs = zeros(length(ccg.t),spikes.numcells);
for cc = 1:spikes.numcells
    acgs(:,cc) = ccg.hists(:,cc,cc);
end
%% Plot the ACGs/CCGs
figure
subplot(2,2,1)
    imagesc(ccg.t,[1 spikes.numcells],log10(acgs(:,ISIStats.sorts.ALL.ratebyclass))')
    axis xy
    xlim([-0.2 0.2])
    
subplot(4,2,2)
    plot(ccg.t,ccg.popccgs(:,1,1),'k')
    hold on
    plot(ccg.t,ccg.popccgs(:,2,2),'r')
    xlim([-0.2 0.2])
    
subplot(4,2,4)
    plot(ccg.t,ccg.popccgs(:,2,1),'k')
    hold on
    plot(ccg.t,ccg.popccgs(:,1,2),'r')
    xlim([-0.2 0.2])
    
    
    %% Log Spaced ACG/CCG
    

    [ logacg ] = bz_LogACG( spikes.times,...
        'ints',SleepState.ints);
    %%
    [ poplogacg ] = bz_LogACG( {popspikes.Espikes,popspikes.Ispikes},...
        'ints',SleepState.ints,'numbins',60);
    for ss = 1:numstates
        for cc = 1:numclasses
        poplogacg.(statenames{ss}).(classnames{cc}) = ...
            poplogacg.(statenames{ss}).acg(:,cc)./sum(CellClass.(classnames{cc}));
        end
    end
%end
%%
figure
for ss = 1:numstates
subplot(3,3,ss*3)
imagesc(logacg.t,[1 spikes.numcells],log10(logacg.(statenames{ss}).acg(:,ISIStats.sorts.(statenames{ss}).ratepE)'))
colorbar
caxis([-1.5 1.5])
title((statenames{ss}))
LogScale('x',10)
ylabel({'Unit:','Sort by Rate/Type'})
end

for cc = 1:numclasses
subplot(2,2,cc*2-1)
    for ss = 1:numstates
        plot(poplogacg.t,(poplogacg.(statenames{ss}).(classnames{cc})),...
            statecolors{ss},'linewidth',2)
        hold on
    end
    title(['Population ACG: ',classnames{cc},' cells'])
    xlabel('t lag (s)');ylabel('Rate (Hz/Cell)')
    axis tight
    box off
    LogScale('x',10)
end

%LogScale('y',10)
%plot(poplogacg.t,log10(poplogacg.(statenames{ss}).pI),'r')
%ylim([0 1.2])
end