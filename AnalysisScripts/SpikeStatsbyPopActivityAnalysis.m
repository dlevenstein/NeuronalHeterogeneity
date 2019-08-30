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
    spikemat.poprate.(celltypes{tt}) = mean(spikemat.data(:,CellClass.(celltypes{tt})),2);
    %Percentage of cells active (more than half a gaussian kernal)
    spikemat.cellsync.(celltypes{tt}) = mean(spikemat.data(:,CellClass.(celltypes{tt}))>0.5,2);%./...
            %sum(CellClass.(celltypes{tt}));
end
spikemat.poprate.ALL = mean(spikemat.data(:,CellClass.pI|CellClass.pE),2);
%%
for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell');
    thiscell = false(size(CellClass.pE));
    thiscell(cc) = true;
    spikemat.cellrate{cc} = spikemat.data(:,cc);
    for tt = 1:length(celltypes)
        popratehist.(celltypes{tt})(cc) = sum(CellClass.(celltypes{tt}) & ~thiscell);
        spikemat.bycellpoprate.(celltypes{tt}){cc} = mean(spikemat.data(:,CellClass.(celltypes{tt}) & ~thiscell),2);

    end
    spikemat.bycellpoprate.ALL{cc} = mean(spikemat.data(:,(CellClass.pI|CellClass.pE) & ~thiscell),2);
            %sum(~thiscell)./binsize;
end

%% Normalizations
normtypes = {'lin','log','lognorm','norm'};
synchtypes = [celltypes,'ALL'];

for tt = 1:length(synchtypes)
    spikemat.poprate.lognorm.(synchtypes{tt}) = log10(spikemat.poprate.(synchtypes{tt})./mean(spikemat.poprate.(synchtypes{tt})));
    spikemat.poprate.norm.(synchtypes{tt}) = (spikemat.poprate.(synchtypes{tt})./mean(spikemat.poprate.(synchtypes{tt})));

    spikemat.poprate.log.(synchtypes{tt}) = log10(spikemat.poprate.(synchtypes{tt}));
    spikemat.poprate.lin.(synchtypes{tt}) = spikemat.poprate.(synchtypes{tt});

    spikemat.bycellpoprate.lognorm.(synchtypes{tt}) = cellfun(@(X) ...
        log10(X./mean(X)),spikemat.bycellpoprate.(synchtypes{tt}),'UniformOutput',false);
    spikemat.bycellpoprate.norm.(synchtypes{tt}) = cellfun(@(X) ...
        (X./mean(X)),spikemat.bycellpoprate.(synchtypes{tt}),'UniformOutput',false);
    spikemat.bycellpoprate.log.(synchtypes{tt}) = cellfun(@(X) ...
        log10(X),spikemat.bycellpoprate.(synchtypes{tt}),'UniformOutput',false);
    spikemat.bycellpoprate.lin.(synchtypes{tt}) = spikemat.bycellpoprate.(synchtypes{tt});
end

%% Calculate E and I pop rate (of other cells) for each spike

for tt = 1:length(synchtypes)
    for nn = 1:length(normtypes)
    ISIStats.allspikes.poprate.(normtypes{nn}).(synchtypes{tt}) = ...
        cellfun(@(X,Y) interp1(spikemat.timestamps,X,Y,'nearest'),...
        spikemat.bycellpoprate.(normtypes{nn}).(synchtypes{tt}),ISIStats.allspikes.times,...
        'UniformOutput',false);
    end
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

nbins = 100;
popratehist.log.bins.pE = linspace(-1.5,1.5,nbins+1);
popratehist.log.bins.pI = linspace(-0.5,2,nbins+1);
popratehist.log.bins.ALL = linspace(-1,1.5,nbins+1);

popratehist.lognorm.bins.pE = linspace(-1.25,1,nbins+1);
popratehist.lognorm.bins.pI = popratehist.lognorm.bins.pE;
popratehist.lognorm.bins.ALL = popratehist.lognorm.bins.pE;

popratehist.norm.bins.pE = linspace(0,5,nbins+1);
popratehist.norm.bins.pI = popratehist.norm.bins.pE;
popratehist.norm.bins.ALL = popratehist.norm.bins.pE;

popratehist.lin.bins.pE = linspace(0,20,nbins+1);
popratehist.lin.bins.pI = linspace(-0,50,nbins+1);
popratehist.lin.bins.ALL = linspace(0,20,nbins+1);


for nn = 1:length(normtypes)
    for st = 1:length(synchtypes)
       popratehist_joint.(normtypes{nn}).bins.(synchtypes{st}) =   popratehist.(normtypes{nn}).bins.(synchtypes{st});
    end
for ss = 1:3      
    %Pop rate distribution
    [popratehist_joint.(normtypes{nn}).(statenames{ss}).alltime]...
        = histcounts2((spikemat.poprate.(normtypes{nn}).pE(spikemat.instate.(statenames{ss}))),...
        (spikemat.poprate.(normtypes{nn}).pI(spikemat.instate.(statenames{ss}))),...
        popratehist_joint.(normtypes{nn}).bins.pE,popratehist_joint.(normtypes{nn}).bins.pI);
    
    
    popratehist_joint.(normtypes{nn}).(statenames{ss}).alltime = ...
        popratehist_joint.(normtypes{nn}).(statenames{ss}).alltime./...
        sum(popratehist_joint.(normtypes{nn}).(statenames{ss}).alltime(:));
    for st = 1:length(synchtypes)
        popratehist.(normtypes{nn}).(statenames{ss}).(synchtypes{st}) = ...
            hist((spikemat.poprate.(normtypes{nn}).(synchtypes{st})(spikemat.instate.(statenames{ss}))),...
            popratehist.(normtypes{nn}).bins.(synchtypes{st}));
        popratehist.(normtypes{nn}).(statenames{ss}).(synchtypes{st}) = ...
            popratehist.(normtypes{nn}).(statenames{ss}).(synchtypes{st})./...
            sum(popratehist.(normtypes{nn}).(statenames{ss}).(synchtypes{st}));
        
        bycellpopratehist.(normtypes{nn}).bins.(synchtypes{st}) = popratehist.(normtypes{nn}).bins.(synchtypes{st});
        bycellpopratehist.(normtypes{nn}).(statenames{ss}).(synchtypes{st}) = ...
            cellfun(@(BCrate) hist(BCrate(spikemat.instate.(statenames{ss})),...
            bycellpopratehist.(normtypes{nn}).bins.(synchtypes{st})),...
            spikemat.bycellpoprate.(normtypes{nn}).(synchtypes{st}),'UniformOutput',false);
        bycellpopratehist.(normtypes{nn}).(statenames{ss}).(synchtypes{st}) = ...
            cellfun(@(BCrate) BCrate./sum(BCrate),...
            bycellpopratehist.(normtypes{nn}).(statenames{ss}).(synchtypes{st}),'UniformOutput',false);
        bycellpopratehist.(normtypes{nn}).(statenames{ss}).(synchtypes{st}) = ...
            cat(1,bycellpopratehist.(normtypes{nn}).(statenames{ss}).(synchtypes{st}){:});
    end
    
end
end
%%
figure
for nn = 1:length(normtypes)
for ss =1:3
        
    subplot(4,3,ss+(nn-1)*3)
    hold on
    for st = 1:length(synchtypes)
        plot(popratehist.(normtypes{nn}).bins.(synchtypes{st}),popratehist.(normtypes{nn}).(statenames{ss}).(synchtypes{st}))
    end
    xlabel('Pop Rate (Hz)')
end
end
    %% Spiking statistics wrt pop rate
clear countmap
for nn = 1:length(normtypes)
for ss = 1:3  
    %Mean CV2 by pop rate
    [popratehist_joint.(normtypes{nn}).(statenames{ss}).cellCV2s,countmap.(normtypes{nn}).(statenames{ss}).numspikes ] = cellfun(@(Erate,Irate,CV2,instate)...
        ConditionalHist3((Erate(instate)),(Irate(instate)),CV2(instate),...
        'numXbins',25,'numYbins',25,'Xbounds',popratehist_joint.(normtypes{nn}).bins.pE([1 end]),'Ybounds',popratehist_joint.(normtypes{nn}).bins.pI([1 end]),...
        'minXY',25,'countnorm',false),...
        ISIStats.allspikes.poprate.(normtypes{nn}).pE,ISIStats.allspikes.poprate.(normtypes{nn}).pI,ISIStats.allspikes.CV2,...
        ISIStats.allspikes.instate.(statenames{ss}),'UniformOutput',false);
    popratehist_joint.(normtypes{nn}).(statenames{ss}).cellCV2s = cat(3,popratehist_joint.(normtypes{nn}).(statenames{ss}).cellCV2s{:});
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
    [popratehist_joint.(normtypes{nn}).(statenames{ss}).cellsync,countmap.(normtypes{nn}).(statenames{ss}).numtimebins ] = cellfun(@(X,Y) ConditionalHist3(...
        (X(spikemat.instate.(statenames{ss}))),(Y(spikemat.instate.(statenames{ss}))),...
        spikemat.cellsync.pE(spikemat.instate.(statenames{ss})),...
        'numXbins',25,'numYbins',25,'Xbounds',popratehist_joint.(normtypes{nn}).bins.pE([1 end]),...
        'Ybounds',popratehist_joint.(normtypes{nn}).bins.pI([1 end]),...
        'minXY',100,'countnorm',false),...
        spikemat.bycellpoprate.(normtypes{nn}).pE,spikemat.bycellpoprate.(normtypes{nn}).pI,...
        'UniformOutput',false);
    popratehist_joint.(normtypes{nn}).(statenames{ss}).cellsync = cat(3,popratehist_joint.(normtypes{nn}).(statenames{ss}).cellsync{:});
    
    %P(spike|poprate)
    popratehist_joint.(normtypes{nn}).(statenames{ss}).pSpk = cellfun(@(X,Y) X./(Y.*dt),...
        countmap.(normtypes{nn}).(statenames{ss}).numspikes,countmap.(normtypes{nn}).(statenames{ss}).numtimebins,...
        'UniformOutput',false);
    popratehist_joint.(normtypes{nn}).(statenames{ss}).pSpk = cat(3,popratehist_joint.(normtypes{nn}).(statenames{ss}).pSpk{:});

    for tt = 1:length(celltypes)
        popratehist_joint.(normtypes{nn}).(statenames{ss}).(celltypes{tt}).cellCV2s = nanmean(popratehist_joint.(normtypes{nn}).(statenames{ss}).cellCV2s(:,:,CellClass.(celltypes{tt})),3);
        popratehist_joint.(normtypes{nn}).(statenames{ss}).(celltypes{tt}).pSpk = nanmean(popratehist_joint.(normtypes{nn}).(statenames{ss}).pSpk(:,:,CellClass.(celltypes{tt})),3);
       % popratehist.(statenames{ss}).(celltypes{tt}).geomeanISIs = nanmean(popratehist.(statenames{ss}).geomeanISIs(:,:,CellClass.(celltypes{tt})),3);

    end
end
end

%%
for nn = 1:length(normtypes)
figure
for ss = 1:3
    subplot(3,3,ss)
        h = imagesc(popratehist_joint.(normtypes{nn}).bins.pE,popratehist_joint.(normtypes{nn}).bins.pI,...
            popratehist_joint.(normtypes{nn}).(statenames{ss}).alltime');
        axis xy
        set(h,'AlphaData',~(popratehist_joint.(normtypes{nn}).(statenames{ss}).alltime'==0));

        title(statenames{ss})
    for tt = 1:length(celltypes)
    subplot(6,6,(ss-1)*2+12+tt)
        h = imagesc(popratehist_joint.(normtypes{nn}).bins.pE,popratehist_joint.(normtypes{nn}).bins.pI,...
            log10(popratehist_joint.(normtypes{nn}).(statenames{ss}).(celltypes{tt}).pSpk)');
        axis xy
        set(h,'AlphaData',~isnan(popratehist_joint.(normtypes{nn}).(statenames{ss}).(celltypes{tt}).pSpk'));
        colorbar
        %crameri lajolla
        caxis([-0.5 1.75])
        
%     subplot(6,6,(ss-1)*2+18+tt)
%             h = imagesc(popratehist.bins.pE,popratehist.bins.pI,1./(popratehist.(statenames{ss}).(celltypes{tt}).geomeanISIs)');
%         axis xy
%         set(h,'AlphaData',~isnan(popratehist.(statenames{ss}).(celltypes{tt}).geomeanISIs'));
%         colorbar
%         %crameri lajolla
%         %caxis([-1 1])
        
    subplot(6,6,(ss-1)*2+24+tt)
        h = imagesc(popratehist_joint.(normtypes{nn}).bins.pE,popratehist_joint.(normtypes{nn}).bins.pI,...
            popratehist_joint.(normtypes{nn}).(statenames{ss}).(celltypes{tt}).cellCV2s');
        axis xy
        set(h,'AlphaData',~isnan(popratehist_joint.(normtypes{nn}).(statenames{ss}).(celltypes{tt}).cellCV2s'));
        colorbar
        crameri berlin
        caxis([0.7 1.3])
    end
end

NiceSave(['PopRateHists_',(normtypes{nn})],figfolder,baseName)
end


    
%% Calculate Conditional distributions on synchrony (pop rate) in each state
for nn = 1:length(normtypes)
for ss = 1:3
    %statenames{ss} = statenames{ss};

    %% Calculate conditional ISI/Synch distributions


    for st = 1:length(synchtypes)
        [ ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) ] = cellfun(@(X,Y,Z,W) ConditionalHist( ([Z(W);Z(W)]),log10([X(W);Y(W)]),...
            'Xbounds',popratehist.(normtypes{nn}).bins.(synchtypes{st})([1 end]),'numXbins',25,'Ybounds',[-3 2],'numYbins',125,'minX',100),...
            ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
            ISIStats.allspikes.poprate.(normtypes{nn}).(synchtypes{st}),ISIStats.allspikes.instate.(statenames{ss}),...
            'UniformOutput',false);
        ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = cat(1,ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}){:});
        ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = bz_CollapseStruct( ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}),3);

        [ SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) ] = cellfun(@(X,Y,Z,W) ConditionalHist( log10([X(W);Y(W)]),([Z(W);Z(W)]),...
            'Ybounds',popratehist.(normtypes{nn}).bins.(synchtypes{st})([1 end]),'numYbins',50,'Xbounds',[-3 2],'numXbins',125,'minX',50),...
            ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
            ISIStats.allspikes.poprate.(normtypes{nn}).(synchtypes{st}),ISIStats.allspikes.instate.(statenames{ss}),...
            'UniformOutput',false);
        SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = cat(1,SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}){:});
        SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = bz_CollapseStruct( SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}),3);

        for cc = 1:length(celltypes)
            ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) = nanmean(ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pYX(:,:,CellClass.(celltypes{cc})),3);
            SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) = nanmean(SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pYX(:,:,CellClass.(celltypes{cc})),3);

            ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).celltypeidx.(celltypes{cc}) = CellClass.(celltypes{cc});
            SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).celltypeidx.(celltypes{cc}) = CellClass.(celltypes{cc});
        end
    end
end
end
%%
for nn = 1:length(normtypes)
figure
for ss = 1:3
for tt = 1:length(celltypes)
for st = 1:3
subplot(6,3,(ss-1)+(tt-1)*3+(st-1)*6+1)
    imagesc(ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Ybins(1,:,1), ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
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
    xlim(ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,[1 end],1))
  
end
end
end

NiceSave(['ISIbySynch_',(normtypes{nn})],figfolder,baseName)

end
%%
for nn = 1:length(normtypes)
figure
for ss = 1:3
for tt = 1:length(celltypes)
for st = 1:3
subplot(6,3,(ss-1)+(tt-1)*3+(st-1)*6+1)
    
    imagesc(SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Ybins(1,:,1), SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
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
    xlim(SynchbyISI.(normtypes{nn}).(synchtypes{tt}).(statenames{ss}).Xbins(1,[1 end],1))
end 
end
end

NiceSave(['SynchbyISI_',(normtypes{nn})],figfolder,baseName)
end



%% Correlate CV2, rate with E/I rate
for ss = 1:length(statenames)
    %statenames{ss} = statenames{ss};


    for tt = 1:length(synchtypes)
        for cc = 1:spikes.numcells
            if sum(ISIStats.allspikes.instate.(statenames{ss}){cc})==0
                CV2popcorr.(statenames{ss}).(synchtypes{tt})(cc) = nan;
                ratepopcorr.(statenames{ss}).(synchtypes{tt})(cc) = nan;
                continue
            end
            
            
            CV2popcorr.(statenames{ss}).(synchtypes{tt})(cc) = ...
                corr(ISIStats.allspikes.poprate.lin.(synchtypes{tt}){cc}(ISIStats.allspikes.instate.(statenames{ss}){cc}),...
                ISIStats.allspikes.CV2{cc}(ISIStats.allspikes.instate.(statenames{ss}){cc}),'type','spearman');
            
            ratepopcorr.(statenames{ss}).(synchtypes{tt})(cc) = ...
                corr(ISIStats.allspikes.poprate.lin.(synchtypes{tt}){cc}(ISIStats.allspikes.instate.(statenames{ss}){cc}),...
                ISIStats.allspikes.cellrate{cc}(ISIStats.allspikes.instate.(statenames{ss}){cc}),'type','spearman');
            
        end
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