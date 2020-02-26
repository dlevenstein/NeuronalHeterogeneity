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
SAVECELLINFO = true;

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
statecolor = {'b','k','r'};

%%
statenames = {'NREMstate','WAKEstate','REMstate'};
%%

[ ISIstats ] = bz_ISIStats( spikes,'ints',SleepState.ints,...
    'cellclass',CellClass.label,'shuffleCV2',false,...
    'savecellinfo',false,'basePath',basePath,'forceRedetect',true,...
    'numISIbins',150,'logISIbounds',[0.0001 500]);

cellinfofilename = fullfile(basePath,[baseName,'.GammaFit.cellinfo.mat']);
%%
numAS.NREMstate = 5;
numAS.WAKEstate = 5;
numAS.REMstate = 5;

%%
spkthresh = 250;
for ss = 1:3
    numspks = cellfun(@sum,ISIstats.allspikes.instate.(statenames{ss}));
    logtimebins = ISIstats.ISIhist.logbins;
    logISIhist = ISIstats.ISIhist.(statenames{ss}).log;
    usecells{ss} = find(CellClass.pE & numspks>spkthresh);
    logISIhist = logISIhist(usecells{ss},:)';
    logISIhist = logISIhist./mode(diff(logtimebins));
    GammaFit.(statenames{ss}) = bz_FitISISharedGammaModes(logISIhist,logtimebins,...
        'numAS',numAS.(statenames{ss}),...
        'figfolder',figfolder,'basePath',basePath,...
        'AScost_lambda',0.13,'AScost_p',1/2,'ASguess',true,'MScost',2.2,'figname',(statenames{ss}));
    
    
    GammaFit.(statenames{ss}).cellstats.meanrate = ...
        ISIstats.summstats.(statenames{ss}).meanrate(usecells{ss});
    GammaFit.(statenames{ss}).cellstats.UID = spikes.UID(usecells{ss});
    if isfield(spikes,'region')
        GammaFit.(statenames{ss}).cellstats.region = spikes.region(usecells{ss});
    end
    if length(GammaFit.(statenames{ss}).cellstats.meanrate) ~= ...
            length(GammaFit.(statenames{ss}).singlecell)
        error('bad number of cells')
    end
    


end

%Wake/NREM joint indexing
GammaFit.(statenames{1}).cellstats.NW = false(size(GammaFit.(statenames{1}).cellstats.meanrate));
GammaFit.(statenames{2}).cellstats.NW = false(size(GammaFit.(statenames{2}).cellstats.meanrate));
[~,NW1,NW2]=intersect(usecells{1},usecells{2});
GammaFit.(statenames{1}).cellstats.NW(NW1)=true;
GammaFit.(statenames{2}).cellstats.NW(NW2)=true;


if SAVECELLINFO
    save(cellinfofilename,'GammaFit')
end

GScolor = [0.6 0.4 0];

%%


%%
diffrate = log10(GammaFit.WAKEstate.cellstats.meanrate(GammaFit.WAKEstate.cellstats.NW))-...
    log10(GammaFit.NREMstate.cellstats.meanrate(GammaFit.NREMstate.cellstats.NW));
diffGS = [GammaFit.WAKEstate.singlecell(GammaFit.WAKEstate.cellstats.NW).GSlogrates]-...
     [GammaFit.NREMstate.singlecell(GammaFit.NREMstate.cellstats.NW).GSlogrates];
 diffAR = (1-[GammaFit.WAKEstate.singlecell(GammaFit.WAKEstate.cellstats.NW).GSweights]) - ...
      (1-[GammaFit.NREMstate.singlecell(GammaFit.NREMstate.cellstats.NW).GSweights]);
%%
figure
subplot(2,3,1)
plot(log10(GammaFit.WAKEstate.cellstats.meanrate(GammaFit.WAKEstate.cellstats.NW)),...
    log10(GammaFit.NREMstate.cellstats.meanrate(GammaFit.NREMstate.cellstats.NW)),'.');
hold on
UnityLine
xlabel('WAKE');ylabel('NREM')
subplot(2,3,2)
plot([GammaFit.WAKEstate.singlecell(GammaFit.WAKEstate.cellstats.NW).GSlogrates],...
    [GammaFit.NREMstate.singlecell(GammaFit.NREMstate.cellstats.NW).GSlogrates],'.');
hold on
UnityLine
xlabel('WAKE ');ylabel('NREM')
subplot(2,3,3)
plot(1-[GammaFit.WAKEstate.singlecell(GammaFit.WAKEstate.cellstats.NW).GSweights],...
    1-[GammaFit.NREMstate.singlecell(GammaFit.NREMstate.cellstats.NW).GSweights],'.');
hold on
UnityLine
xlabel('WAKE');ylabel('NREM')
subplot(2,3,4)
    plot(diffrate,diffGS,'.')
    hold on
    UnityLine
    plot([0 0],ylim(gca),'k')
    plot(xlim(gca),[0 0],'k')
    xlabel('NW Rate Diff');ylabel('NW GSRate Diff')
subplot(2,3,5)
    plot(diffrate,diffAR,'.')
    hold on
    UnityLine
    plot([0 0],ylim(gca),'k')
    plot(xlim(gca),[0 0],'k')
    xlabel('NW Rate Diff');ylabel('NW ASWeight Diff')
subplot(2,3,6)
    plot(diffGS,diffAR,'.')
    hold on
    UnityLine
    plot([0 0],ylim(gca),'k')
    plot(xlim(gca),[0 0],'k')
    xlabel('NW GSRate Diff');ylabel('NW ASWeight Diff')
        NiceSave('NWDiff',figfolder,baseName);

%% Example cell: 3 states
for ff = 1:3
numex=2;
excells = randi(GammaFit.(statenames{ss}).numcells,numex);
figure
for ee = 1:2
for ss = 1:3
    excell = excells(ee);
subplot(6,3,ss+(ee-1)*9)
    plot(GammaFit.(statenames{ss}).logtimebins,...
        GammaFit.(statenames{ss}).ISIdists(:,excell),...
        'color',[0.5 0.5 0.5],'linewidth',2)
    hold on
    plot(GammaFit.(statenames{ss}).logtimebins,...
        GSASmodel(GammaFit.(statenames{ss}).singlecell(excell),...
        GammaFit.(statenames{ss}).taubins),...
        statecolor{ss},'linewidth',2)
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
    if ee == 1
        title(statenames{ss})
    end
    if ss == 1
        ylabel(['UID: ',num2str(GammaFit.(statenames{ss}).cellstats.UID(excell))])
    end
    xlim([-3 2])

subplot(6,3,[3 6]+ss+(ee-1)*9)
    scatter(-GammaFit.(statenames{ss}).singlecell(excell).ASlogrates(:),...
        log10(GammaFit.(statenames{ss}).singlecell(excell).ASCVs(:)),...
        100*GammaFit.(statenames{ss}).singlecell(excell).ASweights(:)+0.00001,'k','filled')
    hold on
    scatter(-GammaFit.(statenames{ss}).singlecell(excell).GSlogrates,...
        log10(GammaFit.(statenames{ss}).singlecell(excell).GSCVs),...
        100*GammaFit.(statenames{ss}).singlecell(excell).GSweights+0.00001,GScolor,'filled')
    plot(GammaFit.(statenames{ss}).logtimebins([1 end]),[0 0],'k--')
    ylabel('CV');xlabel('mean ISI (s)')
    xlim([-3 1.9])
    ylim([-2 0.75])
    LogScale('x',10,'exp',true)
    LogScale('y',10)
    box on

end
end
%if figfolder
    NiceSave(['CellExample_states',num2str(ff)],figfolder,baseName);
%end
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
