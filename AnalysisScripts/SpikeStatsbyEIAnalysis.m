function [ ] = SpikeStatsbyEIAnalysis(basePath,figfolder)

%% DEV
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/Dino_mPFC/Dino_061814';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyEIAnalysis'];
fitfolder = [reporoot,'/AnalysisScripts/modelfits'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
SleepState.ints.ALL = [0 Inf];
SlowWaves = bz_LoadEvents(basePath,'SlowWaves');
lfp = bz_GetLFP(SlowWaves.detectorinfo.detectionchannel,...
    'basepath',basePath);

%% Cell types and states
[celltypes,~,typeidx] = unique(CellClass.label);

cellcolor = {'k','r'};
statenames = fieldnames(SleepState.ints);

%%
thisstate = 'NREMstate';
stateints = SleepState.ints.(thisstate);
binsize = 0.1; %s
overlap = 10;
[ EIMUARateMap,EIBalRateMap,ISIStats.allspikes.poprate,spikemat ] = bz_EIMUARateMap( ISIStats.allspikes,CellClass,...
    'intervals',stateints,'metric',ISIStats.allspikes.CV2,'binsize',binsize,'overlap',overlap);
%celltypes = [unique(CellClass.label),'ALL'];
%CellClass.ALL = ones(size(CellClass.(celltypes{1})));


%% Correlate CV2, rate with E/I rate
instatespiketimes = cellfun(@(X) InIntervals(X,stateints),...
    ISIStats.allspikes.times,'UniformOutput',false);
nanEI = cellfun(@(X) isnan(X),...
    ISIStats.allspikes.poprate.EIratio,'UniformOutput',false);

for tt = 1:length(celltypes)
    CV2popcorr.(celltypes{tt}) = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
        ISIStats.allspikes.poprate.(celltypes{tt}),ISIStats.allspikes.CV2,instatespiketimes);
    ratepopcorr.(celltypes{tt}) = cellfun(@(X,Y,Z) corr(X(Z),1./Y(Z),'type','spearman'),...
        ISIStats.allspikes.poprate.(celltypes{tt}),ISIStats.allspikes.ISIs,instatespiketimes);
end
    CV2popcorr.ALL = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
        ISIStats.allspikes.poprate.ALL,ISIStats.allspikes.CV2,instatespiketimes);
    ratepopcorr.ALL = cellfun(@(X,Y,Z) corr(X(Z),1./Y(Z),'type','spearman'),...
        ISIStats.allspikes.poprate.ALL,ISIStats.allspikes.ISIs,instatespiketimes);
    
    CV2popcorr.EI = cellfun(@(X,Y,Z,W) corr(X(Z & ~W),Y(Z & ~W),'type','spearman'),...
        ISIStats.allspikes.poprate.EIratio,ISIStats.allspikes.CV2,instatespiketimes,nanEI);
    ratepopcorr.EI = cellfun(@(X,Y,Z,W) corr(X(Z & ~W),1./Y(Z & ~W),'type','spearman'),...
        ISIStats.allspikes.poprate.EIratio,ISIStats.allspikes.ISIs,instatespiketimes,nanEI);
    


%%
cv2color = [makeColorMap([0.6 0.6 1],[0 0 0.8],[0.5 0.5 0.5]);makeColorMap([0.5 0.5 0.5],[0.8 0 0],[1 0.6 0.6])];

figure

    subplot(4,4,1)
    for tt = 1:length(celltypes)
        plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),CV2popcorr.pE(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-eMUA Corr.')
    LogScale('x',10)
    %title(state)
    
    subplot(4,4,2)
    for tt = 1:length(celltypes)
        plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),CV2popcorr.pI(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-iMUA Corr.')
    LogScale('x',10)
    
    
    subplot(4,4,8)
    for tt = 1:length(celltypes)
        plot(CV2popcorr.pE(CellClass.(celltypes{tt})),CV2popcorr.pI(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight
    box off
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('CV2-eMUA Corr.');ylabel('CV2-iMUA Corr.')
    
    subplot(4,4,3)
    for tt = 1:length(celltypes)
        plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),CV2popcorr.ALL(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-MUA Corr.')
    LogScale('x',10)
    
    subplot(4,4,4)
    for tt = 1:length(celltypes)
        plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),CV2popcorr.EI(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-EI Corr.')
    LogScale('x',10)
    
    


subplot(4,4,9)
colormap(gca,cv2color)
h = imagesc(EIMUARateMap.bins{1}./sum(CellClass.pE)./binsize,...
    EIMUARateMap.bins{2}./(sum(CellClass.pI)-1)./binsize,...
    (EIMUARateMap.meanmetric.pI)');
set(h,'AlphaData',~isnan(EIMUARateMap.meanmetric.pI'));
axis xy
colorbar
%LogScale('c',10)
xlabel('pE Rate (Hz/cell)');ylabel('pI Rate (Hz/cell)')
caxis([0.8 1.2])
title('CV2 - pI cells')

subplot(4,4,10)
colormap(gca,cv2color)
h = imagesc(EIMUARateMap.bins{1}./(sum(CellClass.pE)-1)./binsize,...
    EIMUARateMap.bins{2}./(sum(CellClass.pI))./binsize,...
    (EIMUARateMap.meanmetric.pE)');
set(h,'AlphaData',~isnan(EIMUARateMap.meanmetric.pE'));
axis xy
colorbar
xlabel('pE Rate (Hz/cell)');ylabel('pI Rate (Hz/cell)')
title('CV2 - pE cells')
caxis([0.8 1.2])
%LogScale('c',10)

subplot(4,4,13)
colormap(gca,cv2color)
h = imagesc(EIBalRateMap.bins{1},...
    EIBalRateMap.bins{2}./(spikes.numcells-1)./binsize,...
    (EIBalRateMap.meanmetric.pI)');
set(h,'AlphaData',~isnan(EIBalRateMap.meanmetric.pI'));
axis xy
colorbar
%LogScale('c',10)
xlabel('E/I Ratio');ylabel('MUA (Hz/cell)')
caxis([0.8 1.2])
title('CV2 - pI cells')

subplot(4,4,14)
colormap(gca,cv2color)
h = imagesc(EIBalRateMap.bins{1},...
    EIBalRateMap.bins{2}./(spikes.numcells-1)./binsize,...
    (EIBalRateMap.meanmetric.pE)');
set(h,'AlphaData',~isnan(EIBalRateMap.meanmetric.pE'));
axis xy
colorbar
xlabel('E/I Ratio');ylabel('MUA (Hz/cell)')
title('CV2 - pE cells')
caxis([0.8 1.2])
%LogScale('c',10)

NiceSave('CV2byPopRate',figfolder,baseName)

%%
figure
    subplot(4,4,1)
    for tt = 1:length(celltypes)
        plot((ISIStats.summstats.(thisstate).meanCV2(CellClass.(celltypes{tt}))),CV2popcorr.pE(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    plot([1 1],get(gca,'ylim'),'k')
    xlabel('Mean CV2');ylabel('CV2-eMUA Corr.')
    %LogScale('x',10)
    %title(state)
    
    subplot(4,4,2)
    for tt = 1:length(celltypes)
        plot((ISIStats.summstats.(thisstate).meanCV2(CellClass.(celltypes{tt}))),CV2popcorr.pI(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    plot([1 1],get(gca,'ylim'),'k')
    xlabel('Mean CV2');ylabel('CV2-iMUA Corr.')
    %LogScale('x',10)
    
    subplot(4,4,4)
    for tt = 1:length(celltypes)
        plot((ISIStats.summstats.(thisstate).meanCV2(CellClass.(celltypes{tt}))),CV2popcorr.ALL(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    plot([1 1],get(gca,'ylim'),'k')
    xlabel('Mean CV2');ylabel('CV2-MUA Corr.')
    %LogScale('x',10)    
    
    subplot(4,4,5)
    for tt = 1:length(celltypes)
        plot((ISIStats.summstats.(thisstate).meanCV2(CellClass.(celltypes{tt}))),ratepopcorr.pE(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    plot([1 1],get(gca,'ylim'),'k')
    xlabel('Mean CV2');ylabel('Rate-eMUA Corr.')
    %LogScale('x',10)
    %title(state)
    
    subplot(4,4,6)
    for tt = 1:length(celltypes)
        plot((ISIStats.summstats.(thisstate).meanCV2(CellClass.(celltypes{tt}))),ratepopcorr.pI(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    plot([1 1],get(gca,'ylim'),'k')
    xlabel('Mean CV2');ylabel('Rate-iMUA Corr.')
    %LogScale('x',10)
    
    subplot(4,4,8)
    for tt = 1:length(celltypes)
        plot((ISIStats.summstats.(thisstate).meanCV2(CellClass.(celltypes{tt}))),ratepopcorr.ALL(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    plot([1 1],get(gca,'ylim'),'k')
    xlabel('Mean CV2');ylabel('Rate-MUA Corr.')
    %LogScale('x',10)
    
    
    subplot(4,4,11)
        for tt = 1:length(celltypes)
            plot(ratepopcorr.ALL(CellClass.(celltypes{tt})),ratepopcorr.EI(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        plot([0 0],get(gca,'ylim'),'k')
        ylabel('Rate-EI Corr.');xlabel('Rate-MUA Corr.')
    
    subplot(4,4,12)
    hold on
        for tt = 1:length(celltypes)
            plot(ratepopcorr.ALL(CellClass.(celltypes{tt})),CV2popcorr.ALL(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            
        end
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        plot([0 0],get(gca,'ylim'),'k')
        ylabel('CV2-MUA Corr.');xlabel('Rate-MUA Corr.')
        
    subplot(4,4,16)
        for tt = 1:length(celltypes)
            plot(ratepopcorr.EI(CellClass.(celltypes{tt})),CV2popcorr.EI(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        plot([0 0],get(gca,'ylim'),'k')
        ylabel('CV2-EI Corr.');xlabel('Rate-EI Corr.')

NiceSave('CV2SummStats',figfolder,baseName)

%% 

figure
    subplot(4,4,1)
        for tt = 1:length(celltypes)
            plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),ratepopcorr.pE(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Cell Rate-pE Rate Corr.')
        LogScale('x',10)
        %title(state)

    subplot(4,4,2)
        for tt = 1:length(celltypes)
            plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),ratepopcorr.pI(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Cell Rate-pI Rate Corr.')
        LogScale('x',10)

    subplot(4,4,8)
        for tt = 1:length(celltypes)
            plot(ratepopcorr.pE(CellClass.(celltypes{tt})),ratepopcorr.pI(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        xlabel('Rate-eMUA Corr');ylabel('Rate-iMUA Corr')
        plot(get(gca,'xlim'),[0 0],'k')
        plot([0 0],get(gca,'ylim'),'k')
        
    subplot(4,4,3)
        for tt = 1:length(celltypes)
            plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),ratepopcorr.ALL(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Rate-MUA Corr.')
        LogScale('x',10)
        
        
    subplot(4,4,4)
        for tt = 1:length(celltypes)
            plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),ratepopcorr.EI(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Rate-EI Corr.')
        LogScale('x',10)
        
        
    subplot(4,4,9)
        h = imagesc(EIMUARateMap.bins{1}./sum(CellClass.pE)./binsize,...
            EIMUARateMap.bins{2}./(sum(CellClass.pI)-1)./binsize,...
            log10(EIMUARateMap.meanrate.pI)');
        set(h,'AlphaData',~isnan(EIMUARateMap.meanrate.pI'));
        xlabel('pE Rate (Hz/cell)');ylabel('pI Rate (Hz/cell)')
        title('Rate - pI cells')
        axis xy
        colorbar
        caxis([-0.5 0.5])
        LogScale('c',10)
        %

    subplot(4,4,10)
        h = imagesc(EIMUARateMap.bins{1}./(sum(CellClass.pE)-1)./binsize,...
            EIMUARateMap.bins{2}./(sum(CellClass.pI))./binsize,...
            log10(EIMUARateMap.meanrate.pE)');
        set(h,'AlphaData',~isnan(EIMUARateMap.meanrate.pE'));
        xlabel('pE Rate (Hz/cell)');ylabel('pI Rate (Hz/cell)')
        title('Rate - pE cells')
        axis xy
        colorbar
        caxis([-1.5 0])
        LogScale('c',10)
        
    subplot(4,4,13)
        h = imagesc(EIBalRateMap.bins{1},...
            EIBalRateMap.bins{2}./(spikes.numcells)./binsize,...
            log10(EIBalRateMap.meanrate.pI)');
        set(h,'AlphaData',~isnan(EIBalRateMap.meanrate.pI'));
        xlabel('E/I Ratio');ylabel('MUA (Hz/cell)')
        title('Rate - pI cells')
        axis xy
        colorbar
        caxis([-0.5 0.5])
        LogScale('c',10)
        %

    subplot(4,4,14)
        h = imagesc(EIBalRateMap.bins{1},...
            EIBalRateMap.bins{2}./(spikes.numcells)./binsize,...
            log10(EIBalRateMap.meanrate.pE)');
        set(h,'AlphaData',~isnan(EIBalRateMap.meanrate.pE'));
        xlabel('E/I Ratio');ylabel('MUA (Hz/cell)')
        title('Rate - pE cells')
        axis xy
        colorbar
        caxis([-1.5 0])
        LogScale('c',10)

        
NiceSave('RatebyPopRate',figfolder,baseName)




%% Constant Spike Count Bins

nspkintervals = 10;

%Get the stuff in sliding windows
cspkbinstats.meanISI = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.ISIs,'UniformOutput',false);
cspkbinstats.rate = cellfun(@(X) 1./X ,cspkbinstats.meanISI,'UniformOutput',false);
cspkbinstats.bindur = cellfun(@(X) movsum(X,nspkintervals),ISIStats.allspikes.ISIs,'UniformOutput',false);
cspkbinstats.stdISI = cellfun(@(X) movstd(X,nspkintervals),ISIStats.allspikes.ISIs,'UniformOutput',false);
cspkbinstats.CVISI = cellfun(@(X,Y) X./Y,cspkbinstats.stdISI,cspkbinstats.meanISI,'UniformOutput',false);
cspkbinstats.meanCV2 = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.CV2,'UniformOutput',false);
cspkbinstats.times = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.times,'UniformOutput',false);
cspkbinstats.MUA = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.poprate.ALL,'UniformOutput',false);
cspkbinstats.eMUA = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.poprate.pE,'UniformOutput',false);
cspkbinstats.iMUA = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.poprate.pI,'UniformOutput',false);
cspkbinstats.eibal = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.poprate.EIratio,'UniformOutput',false);

%Only the times in the state interval
instatebintimes = cellfun(@(X) InIntervals(X,stateints),...
    cspkbinstats.times,'UniformOutput',false);

%Calculate some summary statistics
cspkbinstats.summstats.MUACV2corr = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
    cspkbinstats.MUA,cspkbinstats.meanCV2,instatebintimes);
cspkbinstats.summstats.eMUACV2corr = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
    cspkbinstats.eMUA,cspkbinstats.meanCV2,instatebintimes);
cspkbinstats.summstats.iMUACV2corr = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
    cspkbinstats.iMUA,cspkbinstats.meanCV2,instatebintimes);
cspkbinstats.summstats.MUAratecorr = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
    cspkbinstats.MUA,cspkbinstats.rate,instatebintimes);
cspkbinstats.summstats.EICV2corr = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
    cspkbinstats.eibal,cspkbinstats.meanCV2,instatebintimes);
cspkbinstats.summstats.EIratecorr = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
    cspkbinstats.eibal,cspkbinstats.rate,instatebintimes);


%% Conditional heatmap on MUA
%Bins
cspkbinstats.cellhists.cv2bins = linspace(0,2,40);
cspkbinstats.cellhists.ratebins = linspace(-1.5,2.5,50);
cspkbinstats.cellhists.MUAbins = linspace(0,100,50);
numbinthresh = 10;

%MUA hist for normalization
cspkbinstats.cellhists.MUAhist = cellfun(@(X,Z) hist((X(Z)),cspkbinstats.cellhists.MUAbins),...
    cspkbinstats.MUA,instatebintimes,...
    'UniformOutput',false);

%Threshold for not enough bins.
underthresh = cellfun(@(X) X<numbinthresh,cspkbinstats.cellhists.MUAhist,...
    'UniformOutput',false);
for cc = 1:spikes.numcells
    cspkbinstats.cellhists.MUAhist{cc}(underthresh{cc}) = nan;
end

%CV2 by MUA
cspkbinstats.cellhists.muaVcv2 = cellfun(@(X,Y,Z) hist3([(X(Z)),Y(Z)],...
    {cspkbinstats.cellhists.MUAbins,cspkbinstats.cellhists.cv2bins}),...
    cspkbinstats.MUA,cspkbinstats.meanCV2,instatebintimes,...
    'UniformOutput',false);

cspkbinstats.cellhists.muaVcv2 = cellfun(@(X,Y) bsxfun(@(x,y) x./y,X,Y'),...
    cspkbinstats.cellhists.muaVcv2,cspkbinstats.cellhists.MUAhist,...
    'UniformOutput',false);

for tt=1:length(celltypes)
    cspkbinstats.cellhists.meanhist.(celltypes{tt}) = nanmean(cat(3,cspkbinstats.cellhists.muaVcv2{CellClass.(celltypes{tt})}),3);
end

%Rate by MUA
cspkbinstats.cellhists.muaVrate = cellfun(@(X,Y,Z) hist3([(X(Z)),log10(Y(Z))],...
    {cspkbinstats.cellhists.MUAbins,cspkbinstats.cellhists.ratebins}),...
    cspkbinstats.MUA,cspkbinstats.rate,instatebintimes,...
    'UniformOutput',false);

cspkbinstats.cellhists.muaVrate = cellfun(@(X,Y) bsxfun(@(x,y) x./y,X,Y'),...
    cspkbinstats.cellhists.muaVrate,cspkbinstats.cellhists.MUAhist,...
    'UniformOutput',false);

for tt=1:length(celltypes)
    cspkbinstats.cellhists.meanratehist.(celltypes{tt}) = nanmean(cat(3,cspkbinstats.cellhists.muaVrate{CellClass.(celltypes{tt})}),3);
end

%% Conditional on EI

cspkbinstats.cellhists.EIbins = linspace(-1,1,40);
numbinthresh = 25;

%EI hist for normalization
cspkbinstats.cellhists.EIhist = cellfun(@(X,Z) hist((X(Z)),cspkbinstats.cellhists.EIbins),...
    cspkbinstats.eibal,instatebintimes,...
    'UniformOutput',false);

%Threshold for not enough bins.
underthresh = cellfun(@(X) X<numbinthresh,cspkbinstats.cellhists.EIhist,...
    'UniformOutput',false);
for cc = 1:spikes.numcells
    cspkbinstats.cellhists.EIhist{cc}(underthresh{cc}) = nan;
end

%CV2 by EI
cspkbinstats.cellhists.eiVcv2 = cellfun(@(X,Y,Z) hist3([(X(Z)),Y(Z)],...
    {cspkbinstats.cellhists.EIbins,cspkbinstats.cellhists.cv2bins}),...
    cspkbinstats.eibal,cspkbinstats.meanCV2,instatebintimes,...
    'UniformOutput',false);

cspkbinstats.cellhists.eiVcv2 = cellfun(@(X,Y) bsxfun(@(x,y) x./y,X,Y'),...
    cspkbinstats.cellhists.eiVcv2,cspkbinstats.cellhists.EIhist,...
    'UniformOutput',false);

for tt=1:length(celltypes)
    cspkbinstats.cellhists.meaneihist.(celltypes{tt}) = nanmean(cat(3,cspkbinstats.cellhists.eiVcv2{CellClass.(celltypes{tt})}),3);
end

%Rate by EI
cspkbinstats.cellhists.eiVrate = cellfun(@(X,Y,Z) hist3([(X(Z)),log10(Y(Z))],...
    {cspkbinstats.cellhists.EIbins,cspkbinstats.cellhists.ratebins}),...
    cspkbinstats.eibal,cspkbinstats.rate,instatebintimes,...
    'UniformOutput',false);

cspkbinstats.cellhists.eiVrate = cellfun(@(X,Y) bsxfun(@(x,y) x./y,X,Y'),...
    cspkbinstats.cellhists.eiVrate,cspkbinstats.cellhists.EIhist,...
    'UniformOutput',false);

for tt=1:length(celltypes)
    cspkbinstats.cellhists.meaneiratehist.(celltypes{tt}) = nanmean(cat(3,cspkbinstats.cellhists.eiVrate{CellClass.(celltypes{tt})}),3);
end

%%
figure
imagesc(cat(1,cspkbinstats.cellhists.MUAhist{ISIStats.sorts.(thisstate).CV2}))

%% COnditional probabilities
figure
for tt = 1:length(celltypes)
    subplot(3,2,tt)
        imagesc(cspkbinstats.cellhists.MUAbins./spikes.numcells./binsize,cspkbinstats.cellhists.cv2bins,cspkbinstats.cellhists.meanhist.(celltypes{tt})')
        hold on
        plot(get(gca,'xlim'),[1 1],'k')
        axis xy
        colorbar
        caxis([0 0.15])
        title(celltypes{tt})
        xlabel(['MUA (Hz)']);
        ylabel('CV2')
end

for tt = 1:length(celltypes)
    subplot(3,2,tt+2)
        imagesc(cspkbinstats.cellhists.MUAbins./spikes.numcells./binsize,cspkbinstats.cellhists.ratebins,cspkbinstats.cellhists.meanratehist.(celltypes{tt})')
        hold on
        LogScale('y',10)
        axis xy
        colorbar
        caxis([0 0.15])
        title(celltypes{tt})
        xlabel(['MUA (Hz)']);
        ylabel('Rate')
end

NiceSave('ConditionalCVRatebyMUA',figfolder,baseName)


%% COnditional probabilities
figure
for tt = 1:length(celltypes)
    subplot(3,2,tt)
        imagesc(cspkbinstats.cellhists.EIbins,cspkbinstats.cellhists.cv2bins,cspkbinstats.cellhists.meaneihist.(celltypes{tt})')
        hold on
        plot(get(gca,'xlim'),[1 1],'k')
        axis xy
        colorbar
        caxis([0 0.15])
        title(celltypes{tt})
        xlabel(['EI (',num2str(nspkintervals),'spk bins)']);
        ylabel('CV2')
end

for tt = 1:length(celltypes)
    subplot(3,2,tt+2)
        imagesc(cspkbinstats.cellhists.EIbins,cspkbinstats.cellhists.ratebins,cspkbinstats.cellhists.meaneiratehist.(celltypes{tt})')
        hold on
        LogScale('y',10)
        axis xy
        colorbar
        caxis([0 0.15])
        title(celltypes{tt})
        xlabel(['EI (',num2str(nspkintervals),'spk bins)']);
        ylabel('Rate')
end

NiceSave('ConditionalCVRatebyEI',figfolder,baseName)

%% Example cells
excells = randsample(spikes.numcells,2);

figure
for ee = 1:length(excells)
    subplot(3,2,ee)
        imagesc(cspkbinstats.cellhists.MUAbins,cspkbinstats.cellhists.cv2bins,cspkbinstats.cellhists.muaVcv2{excells(ee)}')
        hold on
        plot(get(gca,'xlim'),[1 1],'k')
        axis xy
        colorbar
        caxis([0 0.18])
        %title(celltypes{tt})
        xlabel(['MUA (',num2str(nspkintervals),'spk bins)']);
        ylabel('CV2')
end

for ee = 1:length(excells)
    subplot(3,2,ee+2)
        imagesc(cspkbinstats.cellhists.MUAbins,cspkbinstats.cellhists.ratebins,cspkbinstats.cellhists.muaVrate{excells(ee)}')
        hold on
        LogScale('y',10)
        axis xy
        colorbar
        caxis([0 0.18])
        %title(celltypes{tt})
        xlabel(['MUA (',num2str(nspkintervals),'spk bins)']);
        ylabel('Rate')
end


figure
for ee = 1:length(excells)
    subplot(3,2,ee)
        imagesc(cspkbinstats.cellhists.EIbins,cspkbinstats.cellhists.cv2bins,cspkbinstats.cellhists.eiVcv2{excells(ee)}')
        hold on
        plot(get(gca,'xlim'),[1 1],'k')
        axis xy
        colorbar
        caxis([0 0.18])
        %title(celltypes{tt})
        xlabel(['EI']);
        ylabel('CV2')
end

for ee = 1:length(excells)
    subplot(3,2,ee+2)
        imagesc(cspkbinstats.cellhists.EIbins,cspkbinstats.cellhists.ratebins,cspkbinstats.cellhists.eiVrate{excells(ee)}')
        hold on
        LogScale('y',10)
        axis xy
        colorbar
        caxis([0 0.18])
        %title(celltypes{tt})
        xlabel(['EI']);
        ylabel('Rate')
end

%% Figure: Summary Statistics
figure
    subplot(3,3,1)
    plot(cspkbinstats.summstats.MUACV2corr,CV2popcorr.ALL,'.')
    hold on
    %plot([-0.3 0.3],[-0.3 0.3],'k')
    xlabel('10-spk bin Calculation');ylabel('1-spk bin calculation')

subplot(3,3,2)
    hold on
    for tt = 1:length(celltypes)
    plot(cspkbinstats.summstats.eMUACV2corr(CellClass.(celltypes{tt})),...
        cspkbinstats.summstats.iMUACV2corr(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
    end
    axis tight;box off
    %plot([-0.3 0.3],[-0.3 0.3],'k')
    xlabel('eMUA-CV2 correlation');ylabel('iMUA-CV2 correlation')
    
subplot(3,3,4)
    hold on
    for tt = 1:length(celltypes)
    plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),...
        cspkbinstats.summstats.MUACV2corr(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
    end
    plot(get(gca,'xlim'),[0 0],'k')
    LogScale('x',10)
    xlabel('sFR (Hz)');ylabel('MUA-CV2 Correlation')


    
subplot(3,3,5)
    hold on
    for tt = 1:length(celltypes)
    plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),...
        cspkbinstats.summstats.MUAratecorr(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
    end
    plot(get(gca,'xlim'),[0 0],'k')
    LogScale('x',10)
    xlabel('sFR (Hz)');ylabel('MUA-Rate Correlation')
    
subplot(3,3,6)
    hold on
    for tt = 1:length(celltypes)
    plot(cspkbinstats.summstats.MUAratecorr(CellClass.(celltypes{tt})),...
        cspkbinstats.summstats.MUACV2corr(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
    end
    axis tight
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    %LogScale('x',10)
    xlabel('MUA-Rate Correlation');ylabel('MUA-CV2 Correlation')
      
    
subplot(3,3,7)
    hold on
    for tt = 1:length(celltypes)
    plot((ISIStats.summstats.(thisstate).meanCV2(CellClass.(celltypes{tt}))),...
        cspkbinstats.summstats.MUACV2corr(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
    end
    axis tight
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('CV2 (Hz)');ylabel('MUA-CV2 Correlation')


    
subplot(3,3,8)
    hold on
    for tt = 1:length(celltypes)
    plot((ISIStats.summstats.(thisstate).meanCV2(CellClass.(celltypes{tt}))),...
        cspkbinstats.summstats.MUAratecorr(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
    end
    axis tight
    plot(get(gca,'xlim'),[0 0],'k')
%    LogScale('x',10)
    xlabel('CV2 (Hz)');ylabel('MUA-Rate Correlation')
    
    
%% Figure: Pop Rate

%Example cell... show example window: high pE, medium pI, high pI medium
%pE, low pE pI

excell = randsample(spikes.numcells,1);

bigwin = bz_RandomWindowInIntervals( stateints,10);

figure
subplot(4,1,1)
bz_MultiLFPPlot(lfp,'timewin',bigwin,'spikes',spikes,'sortmetric',ISIStats.summstats.NREMstate.meanrate)
set(gca,'XTickLabels',[])
xlabel('')

subplot(6,1,3)
plot(spikemat.timestamps,spikemat.poprate.pE,'k')
hold on
plot(spikemat.timestamps,spikemat.poprate.pI,'r')
plot(spikemat.timestamps,spikemat.poprate.ALL,'k','linewidth',2)
xlim(bigwin)
ylabel('MUA')
legend('eMUA','iMUA','MUA')
set(gca,'XTickLabels',[])
box off


subplot(6,1,4)
plot(spikemat.timestamps,spikemat.poprate.EIratio,'k')
xlim(bigwin)
ylabel('E/I Ratio')
set(gca,'XTickLabels',[])
box off



subplot(6,1,5)
plot(cspkbinstats.times{excell},cspkbinstats.meanCV2{excell},'o-')
%plot(cat(1,ISIStats.allspikes.times{:}),cat(1,ISIStats.allspikes.CV2{:}),'.')
hold on
plot(bigwin,[1 1],'k--')
xlim(bigwin)

NiceSave('ExWin',figfolder,baseName)


%%
excell = 20;
figure
plot(cspkbinstats.MUA{excell},log10(cspkbinstats.rate{excell}),'.')
LogScale('y',10)
%% Delta for GLM 
deLFP = bz_Filter(lfp,'passband',[0.5 6],'order',1,'filter','fir1');
predLFP = deLFP;
predLFP.data = predLFP.hilb;
predLFP.freqs = 4;

%%

% [ GLMmodelfit ] = GLMEI_LFP(spikes,predLFP,CellClass,'intervals',stateints,...
%     'refcell',9,'smoothwin',binsize);
%% GLM
%for cc = 1:spikes.numcells
clear GLMmodelfit
for cc = 1:spikes.numcells
    cc
% [ GLMmodelfit(cc) ] = GLMLFP_param(spikes.times(cc),predLFP,...
%     'intervals',SleepState.ints.NREMstate );
%     [ GLMmodelfit(cc) ] = GLMEI(spikes,CellClass,'intervals',stateints,...
%         'refcell',cc,'smoothwin',binsize);
    [ GLMmodelfit(cc) ] = GLMEI_LFP(spikes,predLFP,CellClass,'intervals',stateints,...
        'refcell',cc,'smoothwin',binsize);
end
%save([fitfolder,'EIMUA'],GLMmodelfit)

%% Take a look at some parameters
figure
subplot(4,4,1)
    plot(log10(ISIStats.summstats.(thisstate).meanrate),log10([GLMmodelfit.R0]),'.')
subplot(4,4,2)
    scatter([GLMmodelfit.RE],[GLMmodelfit.RI],2,log10(ISIStats.summstats.(thisstate).meanrate))
    xlabel('E weight');ylabel('I weight')
    hold on
    axis tight;box off
    plot([0 0],get(gca,'ylim'),'k');plot(get(gca,'xlim'),[0 0],'k');
subplot(4,4,3)
hold on
for tt = 1:length(celltypes)
    plot([GLMmodelfit(CellClass.(celltypes{tt})).RE],[GLMmodelfit(CellClass.(celltypes{tt})).RI],'.','color',cellcolor{tt})
end
    xlabel('E weight');ylabel('I weight')
    hold on;axis tight;box off
    plot([0 0],get(gca,'ylim'),'k');plot(get(gca,'xlim'),[0 0],'k');
    
subplot(4,4,4)
hold on
for tt = 1:length(celltypes)
    plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),...
        ([GLMmodelfit(CellClass.(celltypes{tt})).RE]),'.','color',cellcolor{tt})
end
    hold on;axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k');
    xlabel('Rate');ylabel('E weight')
    LogScale('x',10)
subplot(4,4,5)
hold on
for tt = 1:length(celltypes)
    plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),...
        ([GLMmodelfit(CellClass.(celltypes{tt})).RI]),'.','color',cellcolor{tt})
end
    hold on
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k');
    LogScale('x',10)
    xlabel('Rate');ylabel('I weight')
    
subplot(4,4,6)
hold on
for tt = 1:length(celltypes)
    plot(ratepopcorr.pI(CellClass.(celltypes{tt})),...
        ([GLMmodelfit(CellClass.(celltypes{tt})).RI]),'.','color',cellcolor{tt})
end
    hold on
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k');
    plot([0 0],get(gca,'ylim'),'k');

    %LogScale('x',10)
    xlabel('iMUA-Rate Corr');ylabel('I weight')
subplot(4,4,7)
hold on
for tt = 1:length(celltypes)
    plot(ratepopcorr.pE(CellClass.(celltypes{tt})),...
        ([GLMmodelfit(CellClass.(celltypes{tt})).RE]),'.','color',cellcolor{tt})
end
    hold on
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k');
    plot([0 0],get(gca,'ylim'),'k');

    %LogScale('x',10)
    xlabel('eMUA-Rate Corr');ylabel('E weight')

NiceSave('GLMParms',figfolder,baseName)

%% Simulate spikes from GLM
clear Poissmodelfit
clear simspikes
clear simspikes_poiss
for cc = 1:spikes.numcells
    cc
    GLMspkmat = rand(size(GLMmodelfit(cc).predRate))<=GLMmodelfit(cc).predRate;
    Poissspkmat  = rand(size(GLMmodelfit(cc).predRate))<=ISIStats.summstats.(thisstate).meanrate(cc).*GLMmodelfit(cc).dt;
%     for tt = 1:length(GLMmodelfit(cc).timestamps)
%         GLMspkmat(tt) = rand(1)<=GLMmodelfit(cc).predRate(tt);
%         Poissspkmat(tt) = rand(1)<=ISIStats.summstats.NREMstate.meanrate(cc).*GLMmodelfit(cc).dt;
%     end
    simspikes.times{cc} = GLMmodelfit(cc).timestamps(GLMspkmat);
    simspikes_poiss.times{cc} = GLMmodelfit(cc).timestamps(Poissspkmat);
end
simspikes.UID = spikes.UID;
simspikes_poiss.UID = spikes.UID;
%% Calculate ISI stats for simulated spikes
[ ISIstats_sim ] = bz_ISIStats( simspikes,'ints',SleepState.ints,...
    'cellclass',CellClass.label);%,'figfolder',figfolder);

[ ISIstats_poiss ] = bz_ISIStats( simspikes_poiss,'ints',SleepState.ints,...
    'cellclass',CellClass.label);%,'figfolder',figfolder);

%Mean ISI histograms
for tt= 1:length(celltypes)
    ISIstats_sim.ISIhist.NREMstate.popmean.(celltypes{tt}) = ...
        mean(ISIstats_sim.ISIhist.NREMstate.log(CellClass.(celltypes{tt}),:),1);
    ISIstats_poiss.ISIhist.NREMstate.popmean.(celltypes{tt}) = ...
        mean(ISIstats_poiss.ISIhist.NREMstate.log(CellClass.(celltypes{tt}),:),1);
    ISIStats.ISIhist.NREMstate.popmean.(celltypes{tt}) = ...
        mean(ISIStats.ISIhist.NREMstate.log(CellClass.(celltypes{tt}),:),1);
end
%% Calculate MUA stats for the simulated spikes
[ EIMUARateMap_sim,EIBalRateMap_sim,ISIstats_sim.allspikes.poprate,spikemat_sim ] = bz_EIMUARateMap( ISIstats_sim.allspikes,CellClass,...
    'intervals',stateints,'metric',ISIstats_sim.allspikes.CV2,'binsize',binsize,'overlap',overlap,...
    'SHOWFIG',true,'spikemat',spikemat);

%Correlate CV2, rate with E/I rate
instatespiketimes = cellfun(@(X) InIntervals(X,stateints),...
    ISIstats_sim.allspikes.times,'UniformOutput',false);
nanEI = cellfun(@(X) isnan(X),...
    ISIstats_sim.allspikes.poprate.EIratio,'UniformOutput',false);

for tt = 1:length(celltypes)
    CV2popcorr_sim.(celltypes{tt}) = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
        ISIstats_sim.allspikes.poprate.(celltypes{tt}),ISIstats_sim.allspikes.CV2,instatespiketimes);
    ratepopcorr_sim.(celltypes{tt}) = cellfun(@(X,Y,Z) corr(X(Z),1./Y(Z),'type','spearman'),...
        ISIstats_sim.allspikes.poprate.(celltypes{tt}),ISIstats_sim.allspikes.ISIs,instatespiketimes);
end
    CV2popcorr_sim.ALL = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
        ISIstats_sim.allspikes.poprate.ALL,ISIstats_sim.allspikes.CV2,instatespiketimes);
    ratepopcorr_sim.ALL = cellfun(@(X,Y,Z) corr(X(Z),1./Y(Z),'type','spearman'),...
        ISIstats_sim.allspikes.poprate.ALL,ISIstats_sim.allspikes.ISIs,instatespiketimes);
    
    CV2popcorr_sim.EI = cellfun(@(X,Y,Z,W) corr(X(Z & ~W),Y(Z & ~W),'type','spearman'),...
        ISIstats_sim.allspikes.poprate.EIratio,ISIstats_sim.allspikes.CV2,instatespiketimes,nanEI);
    ratepopcorr_sim.EI = cellfun(@(X,Y,Z,W) corr(X(Z & ~W),1./Y(Z & ~W),'type','spearman'),...
        ISIstats_sim.allspikes.poprate.EIratio,ISIstats_sim.allspikes.ISIs,instatespiketimes,nanEI);
    

%%
figure
subplot(3,3,1)
    plot(log10(ISIStats.summstats.(thisstate).meanrate),log10(ISIstats_sim.summstats.(thisstate).meanrate),'.')
subplot(3,3,2)
    plot(ISIStats.summstats.(thisstate).meanCV2,ISIstats_sim.summstats.(thisstate).meanCV2,'.')
    hold on;box off
        plot([0.8 1.4],[0.8 1.4],'k')
        plot([0.8 1.4],[1 1],'k:')
        xlabel('CV2 - Observed');ylabel('CV2 - Simulated')

subplot(3,3,3)
plot(log10(ISIStats.summstats.(thisstate).meanrate),ISIstats_sim.summstats.(thisstate).meanCV2,'.')
hold on;box off
LogScale('x',10)
    plot(get(gca,'xlim'),[1 1],'k:')
    xlabel('FR');ylabel('CV2')
% subplot(2,2,4)
% plot(log10(ISIStats.summstats.(thisstate).meanrate),ISIstats_poiss.summstats.(thisstate).meanCV2,'.')
% hold on
% LogScale('x',10)
%     plot(get(gca,'xlim'),[1 1],'k:')


subplot(3,3,4)
plot(log10(ISIStats.summstats.(thisstate).meanrate),...
    (ISIstats_sim.summstats.(thisstate).meanCV2-ISIStats.summstats.(thisstate).meanCV2)./...
    (ISIStats.summstats.(thisstate).meanCV2-1),'.')
hold on
    plot(get(gca,'xlim'),[0 0],'k-')
    plot(get(gca,'xlim'),[1 1],'k:')
    plot(get(gca,'xlim'),[-1 -1],'k:')

LogScale('x',10)
    xlabel('FR');ylabel('FracCV2')


for tt = 1:length(celltypes)    
subplot(4,3,9+tt)
    plot(ISIstats_sim.ISIhist.logbins,...
        ISIstats_sim.ISIhist.NREMstate.popmean.(celltypes{tt}),'--','color',cellcolor{tt})
    hold on
    plot(ISIstats_poiss.ISIhist.logbins,...
        ISIstats_poiss.ISIhist.NREMstate.popmean.(celltypes{tt}),':','color',cellcolor{tt})
    plot(ISIStats.ISIhist.logbins,...
        ISIStats.ISIhist.NREMstate.popmean.(celltypes{tt}),'-','linewidth',2,'color',cellcolor{tt})
    xlabel('ISI (s)');
    axis tight
    box off
    LogScale('x',10)
    
end

NiceSave('GLMSimCV2',figfolder,baseName)

%% MUA stats compare sim/data
figure
    subplot(4,4,9)
        h = imagesc(EIBalRateMap_sim.bins{1},...
            EIBalRateMap_sim.bins{2}./(spikes.numcells)./binsize,...
            log10(EIBalRateMap_sim.meanrate.pI)');
        set(h,'AlphaData',~isnan(EIBalRateMap_sim.meanrate.pI'));
        xlabel('E/I Ratio');ylabel('MUA (Hz/cell)')
        title('Rate - pI cells')
        axis xy
        colorbar
        caxis([-0.5 0.5])
        LogScale('c',10)
        %

    subplot(4,4,10)
        h = imagesc(EIBalRateMap_sim.bins{1},...
            EIBalRateMap_sim.bins{2}./(spikes.numcells)./binsize,...
            log10(EIBalRateMap_sim.meanrate.pE)');
        set(h,'AlphaData',~isnan(EIBalRateMap_sim.meanrate.pE'));
        xlabel('E/I Ratio');ylabel('MUA (Hz/cell)')
        title('Rate - pE cells')
        axis xy
        colorbar
        caxis([-1.5 0])
        LogScale('c',10)
        
        
    subplot(4,4,13)
        colormap(gca,cv2color)
        h = imagesc(EIBalRateMap_sim.bins{1},...
            EIBalRateMap_sim.bins{2}./(spikes.numcells-1)./binsize,...
            (EIBalRateMap_sim.meanmetric.pI)');
        set(h,'AlphaData',~isnan(EIBalRateMap_sim.meanmetric.pI'));
        axis xy
        colorbar
        %LogScale('c',10)
        xlabel('E/I Ratio');ylabel('MUA (Hz/cell)')
        caxis([0.8 1.2])
        title('CV2 - pI cells')

    subplot(4,4,14)
        colormap(gca,cv2color)
        h = imagesc(EIBalRateMap_sim.bins{1},...
            EIBalRateMap_sim.bins{2}./(spikes.numcells-1)./binsize,...
            (EIBalRateMap_sim.meanmetric.pE)');
        set(h,'AlphaData',~isnan(EIBalRateMap_sim.meanmetric.pE'));
        axis xy
        colorbar
        xlabel('E/I Ratio');ylabel('MUA (Hz/cell)')
        title('CV2 - pE cells')
        caxis([0.8 1.2])
        %LogScale('c',10)

    subplot(4,4,3)
    for tt = 1:length(celltypes)
        plot(log10(ISIstats_sim.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),...
            CV2popcorr_sim.ALL(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-MUA Corr.')
    LogScale('x',10)
    
    subplot(4,4,4)
    for tt = 1:length(celltypes)
        plot(log10(ISIstats_sim.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),...
            CV2popcorr_sim.EI(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-EI Corr.')
    LogScale('x',10)
    
    
    subplot(4,4,7)
        for tt = 1:length(celltypes)
            plot(log10(ISIstats_sim.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),...
                ratepopcorr_sim.ALL(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Rate-MUA Corr.')
        LogScale('x',10)
        
        
    subplot(4,4,8)
        for tt = 1:length(celltypes)
            plot(log10(ISIstats_sim.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),...
                ratepopcorr_sim.EI(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Rate-EI Corr.')
        LogScale('x',10)
        
        
        
    subplot(4,4,1)
    for tt = 1:length(celltypes)
        plot(CV2popcorr.ALL(CellClass.(celltypes{tt})),...
            CV2popcorr_sim.ALL(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    UnityLine
    plot(get(gca,'xlim'),[0 0],'k');plot([0 0],get(gca,'ylim'),'k')
    xlabel('CV2-MUA Corr.');ylabel('CV2-MUA Corr. (sim)')
    
    subplot(4,4,2)
    for tt = 1:length(celltypes)
        plot(CV2popcorr.EI(CellClass.(celltypes{tt})),...
            CV2popcorr_sim.EI(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    UnityLine
    plot(get(gca,'xlim'),[0 0],'k');plot([0 0],get(gca,'ylim'),'k')
    xlabel('CV2-EI Corr.');ylabel('CV2-EI Corr. (sim)')
    
    
    subplot(4,4,5)
        for tt = 1:length(celltypes)
            plot(ratepopcorr.ALL(CellClass.(celltypes{tt})),...
                ratepopcorr_sim.ALL(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
    UnityLine
    plot(get(gca,'xlim'),[0 0],'k');plot([0 0],get(gca,'ylim'),'k')
        xlabel('Rate-MUA Corr.');ylabel('Rate-MUA Corr. (sim)')
        
        
    subplot(4,4,6)
        for tt = 1:length(celltypes)
            plot(ratepopcorr.EI(CellClass.(celltypes{tt})),...
                ratepopcorr_sim.EI(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
    UnityLine
    plot(get(gca,'xlim'),[0 0],'k');plot([0 0],get(gca,'ylim'),'k')
        xlabel('Rate-EI Corr.');ylabel('Rate-EI Corr. (sim)')
        
NiceSave('GLMPopCoupling',figfolder,baseName)
    
    
%% Plot the simulated data
%% Figure
viewwin = bz_RandomWindowInIntervals(SleepState.ints.(thisstate),5);
figure
subplot(5,1,4)
    hold on
    for cc = 1:spikes.numcells
        whichcell = ISIStats.sorts.(thisstate).rate(cc);
        plot(simspikes.times{whichcell},cc.*ones(size(simspikes.times{whichcell})),'k.')
    end
    xlim(viewwin);ylim([0 spikes.numcells])
    ylabel('Simulated Spikes');
    box off
    set(gca,'xticklabels',[])
    
subplot(5,1,5)
    hold on
    for cc = 1:spikes.numcells
        whichcell = ISIStats.sorts.(thisstate).rate(cc);
        plot(simspikes_poiss.times{whichcell},cc.*ones(size(simspikes_poiss.times{whichcell})),'k.')
    end
    xlim(viewwin);ylim([0 spikes.numcells])
    ylabel('Poisson Spikes');
    box off
    set(gca,'xticklabels',[])
%  

subplot(5,1,3)
plot([GLMmodelfit(:).timestamps],log([GLMmodelfit(:).predRate]))
xlim(viewwin)

subplot(5,1,2)
    plot(lfp.timestamps,lfp.data,'k')
    hold on
    %plot(deLFP.timestamps,deLFP.data,'b')
    xlim(viewwin)
    ylabel('LFP');
    box off
    set(gca,'xticklabels',[]);set(gca,'yticklabels',[])
    
subplot(5,1,1)
    hold on
    for cc = 1:spikes.numcells
        whichcell = ISIStats.sorts.(thisstate).rate(cc);
        plot(spikes.times{whichcell},cc.*ones(size(spikes.times{whichcell})),'k.')
    end
    for tt = 1:length(celltypes)
        plot(spikemat.timestamps,spikemat.poprate.(celltypes{tt}),cellcolor{tt})
    end
    xlim(viewwin);ylim([0 spikes.numcells])
    ylabel('Observed Spikes');
    box off
    set(gca,'xticklabels',[])
NiceSave('SimPopEI',figfolder,baseName)


%% GLM: timescale
timescales = logspace(-2.5,0.5,15);
cc = 20;
%clear GLMmodelfit
for tt = 1:length(timescales)
    tt
%     [ GLMmodelfit_ts(tt) ] = GLMEI(spikes,CellClass,'intervals',stateints,...
%         'refcell',cc,'smoothwin',timescales(tt));
    [ GLMmodelfit_ts(tt) ] = GLMEI_LFP(spikes,predLFP,CellClass,'intervals',stateints,...
        'refcell',cc,'smoothwin',timescales(tt));
end
%%
figure
plot(log10(timescales),[GLMmodelfit_ts.nlogL],'o-')
LogScale('x',10)
xlabel('MUA bin size');ylabel('nLogL')
NiceSave('BinSize',figfolder,baseName)
title(['Cell', num2str(cc)])
NiceSave(['BinSizeCell', num2str(cc)],figfolder,baseName)

end