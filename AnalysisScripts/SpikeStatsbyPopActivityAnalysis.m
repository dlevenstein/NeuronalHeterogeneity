function [popratehist_joint,popratehist,ISIbySynch,SynchbyISI,CV2popcorr,ratepopcorr,Ncells ] = SpikeStatsbyPopActivityAnalysis(basePath,figfolder)

%% DEV
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/AGData/Cicero/Cicero_09012014';
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyPopActivityAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
SleepState.ints.ALL = [0 Inf];
% lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
%     'basepath',basePath);

%%

ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end);nan],ISIStats.allspikes.ISIs,...
    'UniformOutput',false);

%% Cell types and states
% try
% celltypes = CellClass.celltypes;
% catch
%     celltypes = unique(CellClass.label);
% end
celltypes = {'pE','pI'};
cellcolor = {'k','r'};
statenames = fieldnames(SleepState.ints);

%% Calculate spike count matrix
binsize = 0.1; %s
dt = 0.005;
spikemat = bz_SpktToSpkmat(spikes,'binsize',binsize,'dt',dt,'bintype','gaussian','units','rate');

%% For each cell, calculate E and I pop rates of all OTHER cells
for tt = 1:length(celltypes)
    Ncells.(celltypes{tt}) = sum(CellClass.(celltypes{tt}));
    %Mean Rate of active cells
    spikemat.poprate.(celltypes{tt}) = mean(spikemat.data(:,CellClass.(celltypes{tt})),2);%./...
            %sum(CellClass.(celltypes{tt}))./binsize;
    %Percentage of cells active (more than half a gaussian kernal)
    spikemat.cellsync.(celltypes{tt}) = mean(spikemat.data(:,CellClass.(celltypes{tt}))>0.5,2);%./...
            %sum(CellClass.(celltypes{tt}));
end
spikemat.poprate.ALL = mean(spikemat.data,2);

for cc = 1:spikes.numcells
    cc
    thiscell = false(size(CellClass.pE));
    thiscell(cc) = true;
    spikemat.cellrate{cc} = spikemat.data(:,cc);
    for tt = 1:length(celltypes)
        popratehist.(celltypes{tt})(cc) = sum(CellClass.(celltypes{tt}) & ~thiscell);
        spikemat.bycellpoprate.(celltypes{tt}){cc} = mean(spikemat.data(:,CellClass.(celltypes{tt}) & ~thiscell),2);%./...
            %sum(CellClass.(celltypes{tt}) & ~thiscell)./binsize;
    end
    spikemat.bycellpoprate.ALL{cc} = mean(spikemat.data(:,~thiscell),2);%./...
            %sum(~thiscell)./binsize;
end


%% Calculate E and I pop rate (of other cells) for each spike
synchtypes = [celltypes,'ALL'];
for tt = 1:length(synchtypes)
    ISIStats.allspikes.poprate.(synchtypes{tt}) = ...
        cellfun(@(X,Y) interp1(spikemat.timestamps,X,Y,'nearest'),...
        spikemat.bycellpoprate.(synchtypes{tt}),ISIStats.allspikes.times,...
        'UniformOutput',false);
end

%% Calculate Cell rate for each spike
    ISIStats.allspikes.cellrate = cellfun(@(X,Y) interp1(spikemat.timestamps,X,Y,'nearest'),...
        spikemat.cellrate,ISIStats.allspikes.times,'UniformOutput',false);
    
%%
for ss = 1:length(statenames)
    spikemat.instate.(statenames{ss}) = InIntervals(spikemat.timestamps,double(SleepState.ints.(statenames{ss})));
    ISIStats.allspikes.instate.(statenames{ss}) = cellfun(@(X) InIntervals(X,double(SleepState.ints.(statenames{ss}))),...
        ISIStats.allspikes.times,'UniformOutput',false);
end
%% Calculate Population Rate Histograms
%     maxrate.pE = 10;
%     maxrate.pI = 80;
%     maxrate.ALL = 25;
    
nbins = 100;
% popratehist.bins.pE = linspace(0,maxrate.pE,nbins+1);
% popratehist.bins.pI = linspace(0,maxrate.pI,nbins+1);
% popratehist.bins.ALL = linspace(0,maxrate.ALL,nbins+1);
% 
popratehist.bins.pE = linspace(-1.25,1,nbins+1);
popratehist.bins.pI = linspace(-0.25,1.75,nbins+1);
popratehist.bins.ALL = linspace(-0.75,1.25,nbins+1);


popratehist_joint.Ebins = popratehist.bins.pE; 
popratehist_joint.Ibins = popratehist.bins.pI;

for ss = 1:3      
    %Pop rate distribution
    [popratehist_joint.(statenames{ss}).alltime]...
        = histcounts2(log10(spikemat.poprate.pE(spikemat.instate.(statenames{ss}))),log10(spikemat.poprate.pI(spikemat.instate.(statenames{ss}))),...
        popratehist_joint.Ebins,popratehist_joint.Ibins);
    popratehist_joint.(statenames{ss}).alltime = popratehist_joint.(statenames{ss}).alltime./sum(popratehist_joint.(statenames{ss}).alltime(:));
    for st = 1:length(synchtypes)
        popratehist.(statenames{ss}).(synchtypes{st}) = hist(log10(spikemat.poprate.(synchtypes{st})(spikemat.instate.(statenames{ss}))),popratehist.bins.(synchtypes{st}));
        popratehist.(statenames{ss}).(synchtypes{st}) = popratehist.(statenames{ss}).(synchtypes{st})./sum(popratehist.(statenames{ss}).(synchtypes{st}));
    end
    
end

%%
figure
for ss =1:3
    subplot(3,3,ss)
        h = imagesc(popratehist_joint.Ebins,popratehist_joint.Ibins,popratehist_joint.(statenames{ss}).alltime');
        axis xy
        set(h,'AlphaData',~(popratehist_joint.(statenames{ss}).alltime'==0));

        title(statenames{ss})
        
    subplot(3,3,ss+3)
    hold on
    for st = 1:length(synchtypes)
        plot(popratehist.bins.(synchtypes{st}),popratehist.(statenames{ss}).(synchtypes{st}))
    end
end
    %% Spiking statistics wrt pop rate
clear countmap
for ss = 1:3  
    %Mean CV2 by pop rate
    [popratehist_joint.(statenames{ss}).cellCV2s,countmap.(statenames{ss}).numspikes ] = cellfun(@(Erate,Irate,CV2,instate)...
        ConditionalHist3(log10(Erate(instate)),log10(Irate(instate)),CV2(instate),...
        'numXbins',25,'numYbins',25,'Xbounds',popratehist_joint.Ebins([1 end]),'Ybounds',popratehist_joint.Ibins([1 end]),...
        'minXY',25,'countnorm',false),...
        ISIStats.allspikes.poprate.pE,ISIStats.allspikes.poprate.pI,ISIStats.allspikes.CV2,...
        ISIStats.allspikes.instate.(statenames{ss}),'UniformOutput',false);
    popratehist_joint.(statenames{ss}).cellCV2s = cat(3,popratehist_joint.(statenames{ss}).cellCV2s{:});
%     
%     %Mean ISI by pop rate
%     [popratehist.(statenames{ss}).geomeanISIs ] = cellfun(@(X,Y,Z,Q,W) ConditionalHist3([X(Q);X(Q)],...
%         [Y(Q);Y(Q)],log10([Z(Q);W(Q)]),...
%         'numXbins',25,'numYbins',25,'Xbounds',[0 maxrate.pE],'Ybounds',[0 maxrate.pI],...
%         'minXY',25),...
%         ISIStats.allspikes.poprate.pE,ISIStats.allspikes.poprate.pI,ISIStats.allspikes.ISIs,...
%         ISIStats.allspikes.instate.(statenames{ss}),ISIStats.allspikes.ISInp1,'UniformOutput',false);
%     popratehist.(statenames{ss}).geomeanISIs = cat(3,popratehist.(statenames{ss}).geomeanISIs{:});
%     
    %Per-cell pop rate distribution
    [popratehist_joint.(statenames{ss}).cellsync,countmap.(statenames{ss}).numtimebins ] = cellfun(@(X,Y) ConditionalHist3(...
        log10(X(spikemat.instate.(statenames{ss}))),log10(Y(spikemat.instate.(statenames{ss}))),...
        spikemat.cellsync.pE(spikemat.instate.(statenames{ss})),...
        'numXbins',25,'numYbins',25,'Xbounds',popratehist_joint.Ebins([1 end]),'Ybounds',popratehist_joint.Ibins([1 end]),...
        'minXY',100,'countnorm',false),...
        spikemat.bycellpoprate.pE,spikemat.bycellpoprate.pI,...
        'UniformOutput',false);
    popratehist_joint.(statenames{ss}).cellsync = cat(3,popratehist_joint.(statenames{ss}).cellsync{:});
    
    %P(spike|poprate)
    popratehist_joint.(statenames{ss}).pSpk = cellfun(@(X,Y) X./(Y.*dt),...
        countmap.(statenames{ss}).numspikes,countmap.(statenames{ss}).numtimebins,...
        'UniformOutput',false);
    popratehist_joint.(statenames{ss}).pSpk = cat(3,popratehist_joint.(statenames{ss}).pSpk{:});

    for tt = 1:length(celltypes)
        popratehist_joint.(statenames{ss}).(celltypes{tt}).cellCV2s = nanmean(popratehist_joint.(statenames{ss}).cellCV2s(:,:,CellClass.(celltypes{tt})),3);
        popratehist_joint.(statenames{ss}).(celltypes{tt}).pSpk = nanmean(popratehist_joint.(statenames{ss}).pSpk(:,:,CellClass.(celltypes{tt})),3);
       % popratehist.(statenames{ss}).(celltypes{tt}).geomeanISIs = nanmean(popratehist.(statenames{ss}).geomeanISIs(:,:,CellClass.(celltypes{tt})),3);

    end
end


%%
figure
for ss = 1:3
    subplot(3,3,ss)
        h = imagesc(popratehist_joint.Ebins,popratehist_joint.Ibins,popratehist_joint.(statenames{ss}).alltime');
        axis xy
        set(h,'AlphaData',~(popratehist_joint.(statenames{ss}).alltime'==0));

        title(statenames{ss})
    for tt = 1:length(celltypes)
    subplot(6,6,(ss-1)*2+12+tt)
        h = imagesc(popratehist_joint.Ebins,popratehist_joint.Ibins,log10(popratehist_joint.(statenames{ss}).(celltypes{tt}).pSpk)');
        axis xy
        set(h,'AlphaData',~isnan(popratehist_joint.(statenames{ss}).(celltypes{tt}).pSpk'));
        colorbar
        %crameri lajolla
        caxis([-0.5 1.75])
        
%     subplot(6,6,(ss-1)*2+18+tt)
%             h = imagesc(popratehist.Ebins,popratehist.Ibins,1./(popratehist.(statenames{ss}).(celltypes{tt}).geomeanISIs)');
%         axis xy
%         set(h,'AlphaData',~isnan(popratehist.(statenames{ss}).(celltypes{tt}).geomeanISIs'));
%         colorbar
%         %crameri lajolla
%         %caxis([-1 1])
        
    subplot(6,6,(ss-1)*2+24+tt)
        h = imagesc(popratehist_joint.Ebins,popratehist_joint.Ibins,popratehist_joint.(statenames{ss}).(celltypes{tt}).cellCV2s');
        axis xy
        set(h,'AlphaData',~isnan(popratehist_joint.(statenames{ss}).(celltypes{tt}).cellCV2s'));
        colorbar
        crameri berlin
        caxis([0.7 1.3])
    end
end

NiceSave('PopRateHists',figfolder,baseName)



    
%% Calculate Conditional distributions on synchrony (pop rate) in each state
for ss = 1:3
    %statenames{ss} = statenames{ss};

    %% Calculate conditional ISI/Synch distributions


    for st = 1:length(synchtypes)
        [ ISIbySynch.(synchtypes{st}).(statenames{ss}) ] = cellfun(@(X,Y,Z,W) ConditionalHist( log10([Z(W);Z(W)]),log10([X(W);Y(W)]),...
            'Xbounds',popratehist.bins.(synchtypes{st})([1 end]),'numXbins',25,'Ybounds',[-3 2],'numYbins',125,'minX',100),...
            ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
            ISIStats.allspikes.poprate.(synchtypes{st}),ISIStats.allspikes.instate.(statenames{ss}),...
            'UniformOutput',false);
        ISIbySynch.(synchtypes{st}).(statenames{ss}) = cat(1,ISIbySynch.(synchtypes{st}).(statenames{ss}){:});
        ISIbySynch.(synchtypes{st}).(statenames{ss}) = bz_CollapseStruct( ISIbySynch.(synchtypes{st}).(statenames{ss}),3);

        [ SynchbyISI.(synchtypes{st}).(statenames{ss}) ] = cellfun(@(X,Y,Z,W) ConditionalHist( log10([X(W);Y(W)]),log10([Z(W);Z(W)]),...
            'Ybounds',popratehist.bins.(synchtypes{st})([1 end]),'numYbins',25,'Xbounds',[-3 2],'numXbins',125,'minX',100),...
            ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
            ISIStats.allspikes.poprate.(synchtypes{st}),ISIStats.allspikes.instate.(statenames{ss}),...
            'UniformOutput',false);
        SynchbyISI.(synchtypes{st}).(statenames{ss}) = cat(1,SynchbyISI.(synchtypes{st}).(statenames{ss}){:});
        SynchbyISI.(synchtypes{st}).(statenames{ss}) = bz_CollapseStruct( SynchbyISI.(synchtypes{st}).(statenames{ss}),3);

        for cc = 1:length(celltypes)
            ISIbySynch.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) = nanmean(ISIbySynch.(synchtypes{st}).(statenames{ss}).pYX(:,:,CellClass.(celltypes{cc})),3);
            SynchbyISI.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) = nanmean(SynchbyISI.(synchtypes{st}).(statenames{ss}).pYX(:,:,CellClass.(celltypes{cc})),3);

            ISIbySynch.(synchtypes{st}).(statenames{ss}).celltypeidx.(celltypes{cc}) = CellClass.(celltypes{cc});
            SynchbyISI.(synchtypes{st}).(statenames{ss}).celltypeidx.(celltypes{cc}) = CellClass.(celltypes{cc});
        end
    end
end
%%
figure
for ss = 1:3
for tt = 1:length(celltypes)
for st = 1:3
subplot(6,3,(ss-1)+(tt-1)*3+(st-1)*6+1)
    imagesc(ISIbySynch.(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),ISIbySynch.(synchtypes{st}).(statenames{ss}).Ybins(1,:,1), ISIbySynch.(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    ylabel([(celltypes{tt}),' ISI (log(s))']);xlabel([(synchtypes{st}),' Synch'])
    if tt==1 & st == 1
        title(statenames{ss})
    end
    %title((celltypes{tt}))
   % colorbar
    if tt ==1 
        caxis([0 0.02])
    elseif tt==2
         caxis([0 0.03])
    end
    xlim(ISIbySynch.(synchtypes{st}).(statenames{ss}).Xbins(1,[1 end],1))
  
end
end
end
   
NiceSave('ISIbySynch',figfolder,baseName)
%%
figure
for ss = 1:3
for tt = 1:length(celltypes)
for st = 1:3
subplot(6,3,(ss-1)+(tt-1)*3+(st-1)*6+1)
    
    imagesc(SynchbyISI.(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),SynchbyISI.(synchtypes{st}).(statenames{ss}).Ybins(1,:,1), SynchbyISI.(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    xlabel([(celltypes{tt}),' ISI (log(s))']);ylabel([(synchtypes{st}),' Synch'])
    %title((celltypes{tt}))
    if tt==1 & st == 1
        title(statenames{ss})
    end
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    xlim(SynchbyISI.(synchtypes{tt}).(statenames{ss}).Xbins(1,[1 end],1))
end 
end
end
NiceSave('SynchbyISI',figfolder,baseName)




%% Correlate CV2, rate with E/I rate
for ss = 1:length(statenames)
    %statenames{ss} = statenames{ss};


    for tt = 1:length(celltypes)
        CV2popcorr.(statenames{ss}).(celltypes{tt}) = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
            ISIStats.allspikes.poprate.(celltypes{tt}),ISIStats.allspikes.CV2,ISIStats.allspikes.instate.(statenames{ss}));
        ratepopcorr.(statenames{ss}).(celltypes{tt}) = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman'),...
            ISIStats.allspikes.poprate.(celltypes{tt}),ISIStats.allspikes.cellrate,ISIStats.allspikes.instate.(statenames{ss}));
    end
end

%%
for ss = 1:3
figure
    subplot(3,3,1)
        for tt = 1:length(celltypes)
            plot(log10(ISIStats.summstats.(statenames{ss}).meanrate(CellClass.(celltypes{tt}))),ratepopcorr.(statenames{ss}).pE(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Cell Rate-pE Rate Corr.')
        LogScale('x',10)
        title((statenames{ss}))

    subplot(3,3,2)
        for tt = 1:length(celltypes)
            plot(log10(ISIStats.summstats.(statenames{ss}).meanrate(CellClass.(celltypes{tt}))),ratepopcorr.(statenames{ss}).pI(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Cell Rate-pI Rate Corr.')
        LogScale('x',10)
        
        
    subplot(3,3,3)
        for tt = 1:length(celltypes)
            plot(ratepopcorr.(statenames{ss}).pE(CellClass.(celltypes{tt})),ratepopcorr.(statenames{ss}).pI(CellClass.(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        xlabel('Rate-pE Corr');ylabel('Rate-pI Corr')
        axis tight
        plot(get(gca,'xlim'),[0 0],'k')
        plot([0 0],get(gca,'ylim'),'k')
        
        
    subplot(3,3,4)
    for tt = 1:length(celltypes)
        plot(log10(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{tt}))),CV2popcorr.(statenames{ss}).pE(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-pE Rate Corr.')
    LogScale('x',10)
    %title(statenames{ss})
    
    subplot(3,3,5)
    for tt = 1:length(celltypes)
        plot(log10(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{tt}))),CV2popcorr.(statenames{ss}).pI(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-pI Rate Corr.')
    LogScale('x',10)
    
    
    subplot(3,3,6)
    for tt = 1:length(celltypes)
        plot(CV2popcorr.(statenames{ss}).pE(CellClass.(celltypes{tt})),CV2popcorr.(statenames{ss}).pI(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('CV2-pE Corr.');ylabel('CV2-pI Corr.')
    
    subplot(3,3,7)
    for tt = 1:length(celltypes)
        plot(ratepopcorr.(statenames{ss}).pE(CellClass.(celltypes{tt})),CV2popcorr.(statenames{ss}).pE(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('Rate-pE Corr.');ylabel('CV2-pE Corr.')
    
    subplot(3,3,8)
    for tt = 1:length(celltypes)
        plot(ratepopcorr.(statenames{ss}).pI(CellClass.(celltypes{tt})),CV2popcorr.(statenames{ss}).pI(CellClass.(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('Rate-pI Corr.');ylabel('CV2-pI Corr.')
    
    NiceSave(['RateCV2PopCorr_',(statenames{ss})],figfolder,baseName)
end

%%
% %% Average Synchrony (Pop Rate) during spikes in the ISI return map
% excell = randsample(spikes.numcells,1);
% 
% ISIbounds = [-3 1.5];
% [meanZ,countmap ] = cellfun(@(X,Y,Z,Q) ConditionalHist3(log10(X(Q)),...
%     log10(Y(Q)),Z(Q),...
%     'numXbins',80,'numYbins',80,'Xbounds',ISIbounds,'Ybounds',ISIbounds,...
%     'minXY',30,'sigma',0.1),...
%     ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,ISIStats.allspikes.poprate.pI,...
%     ISIStats.allspikes.instate.(statenames{ss}),'UniformOutput',false);
% meansynchall = cellfun(@mean,spikemat.bycellpoprate.pI,'UniformOutput',false);
% meanZ_norm = cellfun(@(X,Y) X./Y,meanZ,meansynchall,'UniformOutput',false);
% 
% %%
% for tt = 1:length(celltypes)
%     popratereturn.(celltypes{tt}) = nanmean(cat(3,meanZ_norm{CellClass.(celltypes{tt})}),3);
%     meanreturn.(celltypes{tt}) = nanmean(cat(3,countmap{CellClass.(celltypes{tt})}),3);
% end
% 
% %%
% posnegcolormap = makeColorMap([0 0 0.8],[0 0 0],[0.8 0 0]);
% histcolormap = makeColorMap([1 1 1],[0 0 0]);
% 
% figure
% subplot(2,2,1)
% colormap(gca,histcolormap)
% 
% s = imagesc(ISIbounds,ISIbounds,meanreturn.pE');
% axis xy
% LogScale('xy',10)
% colorbar
% %caxis([0.5 1.5])
% 
% subplot(2,2,2)
% colormap(gca,histcolormap)
% 
% s = imagesc(ISIbounds,ISIbounds,meanreturn.pI')
% axis xy
% LogScale('xy',10)
% colorbar
% %caxis([0.5 1.5])
% 
% subplot(2,2,3)
% colormap(gca,posnegcolormap)
% 
% s = imagesc(ISIbounds,ISIbounds,popratereturn.pE');
% alpha(s,meanreturn.pE'*25)
% axis xy
% caxis([0.5 1.5])
% LogScale('xy',10)
% colorbar
% 
% subplot(2,2,4)
% colormap(gca,posnegcolormap)
% s = imagesc(ISIbounds,ISIbounds,popratereturn.pI')
% alpha(s,meanreturn.pI'*25)
% axis xy
% caxis([0.5 1.5])
% LogScale('xy',10)
% colorbar
% 
% NiceSave('ReturnMapPopRate',figfolder,baseName)
% 
% %% Conditional CV2 distribution given pop synchrony
% [ CONDXY ] = ConditionalHist(ISIStats.allspikes.poprate.pI,...
%     ISIStats.allspikes.CV2,'Xbounds',[0 12],'numXbins',13,...
%     'numYbins',20,'Ybounds',[0 2]);
% %%
% for tt = 1:length(celltypes)
%     CV2distbyPOP.(celltypes{tt}) = nanmean(CONDXY.pYX(:,:,CellClass.(celltypes{tt})),3);
%     meanCV2byPOP.(celltypes{tt}) = nanmean(CONDXY.meanYX(:,:,CellClass.(celltypes{tt})),3);
%     %PopsynchhistbyPOP = .(celltypes{tt}) = 
%     
% end
% %%
% figure
% subplot(2,2,1)
% imagesc(CONDXY.Xbins(:,:,1),CONDXY.Ybins(:,:,1),CV2distbyPOP.pE')
% hold on
% plot(CONDXY.Xbins(:,:,excell),meanCV2byPOP.pE,'w-')
% plot(CONDXY.Xbins(:,[1 end],excell),[1 1],'w--')
% 
% axis xy
% xlabel('Pop Synchrony');ylabel('CV2')
% %%
% excell = randsample(spikes.numcells,1);
% figure
% imagesc(CONDXY.Xbins(:,:,excell),CONDXY.Ybins(:,:,excell),CONDXY.pYX(:,:,excell)')
% hold on
% plot(CONDXY.Xbins(:,:,excell),CONDXY.meanYX(:,:,excell),'w-')
% plot(CONDXY.Xbins(:,:,excell),CONDXY.meanYX(:,:,excell),'w-')
% axis xy
% %%
% excell = randsample(spikes.numcells,1);
% 
% figure
% subplot(2,2,1)
% plot(ISIStats.allspikes.cellrate{excell},ISIStats.allspikes.CV2{excell},'.')
% subplot(2,2,3)
% scatter(log10(ISIStats.allspikes.ISIs{excell}(ISIStats.allspikes.instate{excell})),...
%     log10(ISIStats.allspikes.ISInp1{excell}(ISIStats.allspikes.instate{excell})),1,...
%     log10(ISIStats.allspikes.poprate.pE{excell}(ISIStats.allspikes.instate{excell})))
% colorbar
% 
% subplot(2,2,4)
% scatter(log10(ISIStats.allspikes.ISIs{excell}(ISIStats.allspikes.instate{excell})),...
%     log10(ISIStats.allspikes.ISInp1{excell}(ISIStats.allspikes.instate{excell})),1,...
%     log10(ISIStats.allspikes.poprate.pI{excell}(ISIStats.allspikes.instate{excell})))
% colorbar
% 
% subplot(2,2,2)
% imagesc((meanZ_norm{excell})')
% axis xy
% colorbar
% caxis([0 2])
% %% Figure: Pop Rate
% 
% %Example cell... show example window: high pE, medium pI, high pI medium
% %pE, low pE pI
% 
% excell = randsample(spikes.numcells,1);
% 
% bigwin = bz_RandomWindowInIntervals( SleepState.ints.(state),10 );
% 
% figure
% subplot(3,1,1)
% bz_MultiLFPPlot(lfp,'timewin',bigwin,'spikes',spikes,'cellgroups',{CellClass.pE,CellClass.pI})
% 
% subplot(6,1,3)
% plot(spikemat.timestamps,spikemat.bycellpoprate.pE{excell},'k')
% hold on
% plot(spikemat.timestamps,spikemat.bycellpoprate.pI{excell},'r')
% xlim(bigwin)
% 
% subplot(6,1,4)
% plot(ISIStats.allspikes.times{excell},ISIStats.allspikes.CV2{excell},'o-')
% hold on
% plot(bigwin,[1 1],'k--')
% xlim(bigwin)
% 
% 
% 
% 
% 
% 
% %%
% figure
% subplot(3,3,1)
% imagesc(popratehist.bins{1},popratehist.bins{2},log10(popratehist.all)')
% axis xy
% xlabel('pE Spikes');ylabel('pI Spikes')
% colorbar
% LogScale('c',10)
% 
% subplot(3,3,2)
% imagesc(popratehist.bins{1},popratehist.bins{2},log10(popratehist.popstats.meanrate.pI)')
% axis xy
% colorbar
% LogScale('c',10)
% %caxis([0.5 1.5])
% 
% subplot(3,3,3)
% imagesc(popratehist.bins{1},popratehist.bins{2},log10(popratehist.popstats.meanrate.pE)')
% axis xy
% colorbar
% LogScale('c',10)
% %caxis([0.5 1.5])
% 
% subplot(3,3,5)
% imagesc(popratehist.bins{1},popratehist.bins{2},(popratehist.popstats.meanCV2.pI)')
% axis xy
% colorbar
% %LogScale('c',10)
% %caxis([0.5 1.5])
% 
% subplot(3,3,6)
% imagesc(popratehist.bins{1},popratehist.bins{2},(popratehist.popstats.meanCV2.pE)')
% axis xy
% colorbar
% %LogScale('c',10)
% 
% subplot(3,3,7)
% imagesc(popratehist.bins{1},popratehist.bins{2},log10(1./popratehist.meanISI_bycell{excell})')
% axis xy
% colorbar
% caxis([-1 1])
% LogScale('c',10)
% 
% subplot(3,3,8)
% imagesc(popratehist.bins{1},popratehist.bins{2},popratehist.meanCV2_bycell{excell}')
% axis xy
% colorbar
% caxis([0.8 1.4])


end