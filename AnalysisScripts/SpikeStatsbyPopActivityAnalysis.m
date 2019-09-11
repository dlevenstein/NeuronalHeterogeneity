function [popratehist_joint,popratehist,bycellpopratehist,ISIbySynch,...
    normISIbySynch,CV2popcorr,ratepopcorr,cellinfo,Ncells,...
    PopRatebyPSS,PopRatebyTheta] = SpikeStatsbyPopActivityAnalysis(basePath,figfolder)

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
OccupancyStats = bz_LoadCellinfo(basePath,'OccupancyStats');
SleepState = bz_LoadStates(basePath,'SleepState');
SleepState.ints.ALL = [0 Inf];
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath);
%%
cellinfo.CellClass = CellClass;
cellinfo.ISIStats = ISIStats.summstats;
cellinfo.OccupancyStats = OccupancyStats;
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
    spikemat.totpoprate.(celltypes{tt}) = sum(spikemat.data(:,CellClass.(celltypes{tt})),2);
    spikemat.poprate.(celltypes{tt}) = spikemat.totpoprate.(celltypes{tt})./Ncells.(celltypes{tt});
    
    spikemat.cellsync.(celltypes{tt}) = mean(spikemat.data(:,CellClass.(celltypes{tt}))>0.5,2);
end
Ncells.ALL = Ncells.pE + Ncells.pI;
spikemat.totpoprate.ALL = sum(spikemat.data(:,(CellClass.pI|CellClass.pE)),2);
spikemat.poprate.ALL = spikemat.totpoprate.ALL./Ncells.ALL;

for cc = 1:spikes.numcells %weird roundabout way to calculate is much faster
    bz_Counter(cc,spikes.numcells,'Cell');
    thiscell = false(size(CellClass.pE));
    thiscell(cc) = true;
    spikemat.cellrate{cc} = spikemat.data(:,cc);
    for tt = 1:length(celltypes)
        popratehist.(celltypes{tt})(cc) = sum(CellClass.(celltypes{tt}) & ~thiscell);
        if CellClass.(celltypes{tt})(cc) %if it's in theclass, subtract off the current cell
            spikemat.bycellpoprate.(celltypes{tt}){cc} = ...
                (spikemat.totpoprate.(celltypes{tt})-spikemat.cellrate{cc})./...
               popratehist.(celltypes{tt})(cc);
        else
            spikemat.bycellpoprate.(celltypes{tt}){cc} = ...
                spikemat.totpoprate.(celltypes{tt})./popratehist.(celltypes{tt})(cc);
        end
    end
    
    if CellClass.pI(cc)||CellClass.pE(cc)
        spikemat.bycellpoprate.ALL{cc} = (spikemat.totpoprate.ALL-spikemat.cellrate{cc})./...
            (Ncells.ALL-1);
    else
        spikemat.bycellpoprate.ALL{cc} = spikemat.totpoprate.ALL./Ncells.ALL;
    end
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
    
    cellinfo.instatespikes.(statenames{ss}) = cellfun(@(X) sum(X),...
        ISIStats.allspikes.instate.(statenames{ss}));
end
%% Calculate Population Rate Histograms

nbins = 100;
popratehist.log.bins.pE = linspace(-1.5,1.5,nbins+1);
popratehist.log.bins.pI = linspace(-0.5,2,nbins+1);
popratehist.log.bins.ALL = linspace(-1,1.5,nbins+1);

popratehist.lognorm.bins.pE = linspace(-1,0.75,nbins+1);
popratehist.lognorm.bins.pI = linspace(-1,0.75,nbins+1);
popratehist.lognorm.bins.ALL = linspace(-0.75,0.75,nbins+1);

popratehist.norm.bins.pE = linspace(0,5,nbins+1);
popratehist.norm.bins.pI = linspace(0,4,nbins+1);
popratehist.norm.bins.ALL = linspace(0,5,nbins+1);

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
        
        bycellpopratehist.(normtypes{nn}).bins.(synchtypes{st}) = ...
            popratehist.(normtypes{nn}).bins.(synchtypes{st});
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

for nn = 3:4 %only do normalized rates
for ss = 1:3
    %statenames{ss} = statenames{ss};

    %% Calculate conditional ISI/Synch distributions


    for st = 1:length(synchtypes)
        [ ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) ] = cellfun(@(X,Y,Z,W) ConditionalHist( ([Z(W);Z(W)]),log10([X(W);Y(W)]),...
            'Xbounds',popratehist.(normtypes{nn}).bins.(synchtypes{st})([1 end]),'numXbins',30,'Ybounds',[-3 2],'numYbins',125,'minX',100),...
            ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
            ISIStats.allspikes.poprate.(normtypes{nn}).(synchtypes{st}),ISIStats.allspikes.instate.(statenames{ss}),...
            'UniformOutput',false);
        ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = cat(1,ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}){:});
        ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = bz_CollapseStruct( ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}),3);

%         [ SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) ] = cellfun(@(X,Y,Z,W) ConditionalHist( log10([X(W);Y(W)]),([Z(W);Z(W)]),...
%             'Ybounds',popratehist.(normtypes{nn}).bins.(synchtypes{st})([1 end]),'numYbins',50,'Xbounds',[-3 2],'numXbins',125,'minX',50),...
%             ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
%             ISIStats.allspikes.poprate.(normtypes{nn}).(synchtypes{st}),ISIStats.allspikes.instate.(statenames{ss}),...
%             'UniformOutput',false);
%         SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = cat(1,SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}){:});
%         SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = bz_CollapseStruct( SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}),3);

        for cc = 1:length(celltypes)
            ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) = nanmean(ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pYX(:,:,CellClass.(celltypes{cc})),3);
%            SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) = nanmean(SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pYX(:,:,CellClass.(celltypes{cc})),3);

            ISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).celltypeidx.(celltypes{cc}) = CellClass.(celltypes{cc});
%            SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).celltypeidx.(celltypes{cc}) = CellClass.(celltypes{cc});
        end
        
        
        [ normISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) ] = cellfun(@(X,Y,Z,W,MTO) ConditionalHist( ([Z(W);Z(W)]),log10([X(W);Y(W)]./MTO),...
            'Xbounds',popratehist.(normtypes{nn}).bins.(synchtypes{st})([1 end]),'numXbins',30,'Ybounds',[-4 1],'numYbins',100,'minX',100),...
            ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
            ISIStats.allspikes.poprate.(normtypes{nn}).(synchtypes{st}),...
            ISIStats.allspikes.instate.(statenames{ss}),num2cell(OccupancyStats.(statenames{ss}).median),...
            'UniformOutput',false);
        normISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = cat(1,normISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}){:});
        normISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = bz_CollapseStruct( normISIbySynch.(normtypes{nn}).(synchtypes{st}).(statenames{ss}),3);
        
%         [ SynchbynormISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) ] = cellfun(@(X,Y,Z,W,MTO) ConditionalHist( log10([X(W);Y(W)]./MTO),([Z(W);Z(W)]),...
%             'Ybounds',popratehist.(normtypes{nn}).bins.(synchtypes{st})([1 end]),'numYbins',50,'Xbounds',[-4 1],'numXbins',100,'minX',50),...
%             ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
%             ISIStats.allspikes.poprate.(normtypes{nn}).(synchtypes{st}),...
%             ISIStats.allspikes.instate.(statenames{ss}),num2cell(OccupancyStats.(statenames{ss}).median),...
%             'UniformOutput',false);
%         SynchbynormISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = cat(1,SynchbynormISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}){:});
%         SynchbynormISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}) = bz_CollapseStruct( SynchbynormISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}),3);
    end
end
end
%%
for nn = 3:4
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
%% Synch by ISI figure
% for nn = 3:4
% figure
% for ss = 1:3
% for tt = 1:length(celltypes)
% for st = 1:3
% subplot(6,3,(ss-1)+(tt-1)*3+(st-1)*6+1)
%     
%     imagesc(SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Ybins(1,:,1), SynchbyISI.(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
%     %hold on
%     %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
%     axis xy
%     %LogScale('y',10)
%     xlabel([(celltypes{tt}),' ISI (log(s))']);ylabel([(synchtypes{st}),' Synch'])
%     %title((celltypes{tt}))
%     if tt==1 & st == 1
%         title(statenames{ss})
%     end
% %     if tt ==1 
% %         caxis([0 0.02])
% %     elseif tt==2
% %          caxis([0 0.03])
% %     end
%     xlim(SynchbyISI.(normtypes{nn}).(synchtypes{tt}).(statenames{ss}).Xbins(1,[1 end],1))
% end 
% end
% end
% 
% NiceSave(['SynchbyISI_',(normtypes{nn})],figfolder,baseName)
% end



%% Correlate CV2, rate with E/I rate
for ss = 1:length(statenames)
    %statenames{ss} = statenames{ss};


    for tt = 1:length(synchtypes)
        for cc = 1:spikes.numcells
            ratepopcorr.(statenames{ss}).cellrate = ISIStats.summstats.(statenames{ss}).meanrate;
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

%% Pop Synch wrt state variables

spikemat.BSmetrics.PSS = interp1(...
    SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus,...
    SleepState.detectorinfo.detectionparms.SleepScoreMetrics.broadbandSlowWave,...
    spikemat.timestamps,'nearest');
spikemat.BSmetrics.thratio = interp1(...
    SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus,...
    SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio,...
    spikemat.timestamps,'nearest');

%%
for nn = 1:4
for st = 1:length(synchtypes)
PopRatebyPSS.(normtypes{nn}).(synchtypes{st})  = ConditionalHist(...
    spikemat.BSmetrics.PSS,spikemat.poprate.(normtypes{nn}).(synchtypes{st}),...
    'Xbounds',[0 1],'numXbins',25,'numYbins',100,'minX',2000,...
    'Ybounds',popratehist.(normtypes{nn}).bins.(synchtypes{st})([1 end]));

PopRatebyTheta.notNREM.(normtypes{nn}).(synchtypes{st})  = ConditionalHist(...
    spikemat.BSmetrics.thratio(spikemat.instate.WAKEstate|spikemat.instate.REMstate),...
    spikemat.poprate.(normtypes{nn}).(synchtypes{st})(spikemat.instate.WAKEstate|spikemat.instate.REMstate),...
    'Xbounds',[0 1],'numXbins',25,'numYbins',100,'minX',2000,...
    'Ybounds',popratehist.(normtypes{nn}).bins.(synchtypes{st})([1 end]));

PopRatebyTheta.WAKE.(normtypes{nn}).(synchtypes{st})  = ConditionalHist(...
    spikemat.BSmetrics.thratio(spikemat.instate.WAKEstate),...
    spikemat.poprate.(normtypes{nn}).(synchtypes{st})(spikemat.instate.WAKEstate),...
    'Xbounds',[0 1],'numXbins',25,'numYbins',100,'minX',2000,...
    'Ybounds',popratehist.(normtypes{nn}).bins.(synchtypes{st})([1 end]));

PopRatebyTheta.REM.(normtypes{nn}).(synchtypes{st})  = ConditionalHist(...
    spikemat.BSmetrics.thratio(spikemat.instate.REMstate),...
    spikemat.poprate.(normtypes{nn}).(synchtypes{st})(spikemat.instate.REMstate),...
    'Xbounds',[0 1],'numXbins',25,'numYbins',100,'minX',2000,...
    'Ybounds',popratehist.(normtypes{nn}).bins.(synchtypes{st})([1 end]));

%[-1.5 1.75]
% [ PopRatebyTheta.(synchtypes{st}) ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(~W);Z(~W)],log10([X(~W);Y(~W)]),...
%     'Xbounds',[0 1],'numXbins',20,'Ybounds',[-3 2],'numYbins',125,'minX',20),...
%     ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
%     ISIStats.allspikes.thetarat,ISIStats.allspikes.instate.NREMstate,...
%     'UniformOutput',false);
% ISIbytheta = cat(1,ISIbytheta{:});
% ISIbytheta = CollapseStruct( ISIbytheta,3);
% 
% [ ISIbytheta.WAKE ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(W);Z(W)],log10([X(W);Y(W)]),...
%     'Xbounds',[0 1],'numXbins',20,'Ybounds',[-3 2],'numYbins',125,'minX',20),...
%     ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
%     ISIStats.allspikes.thetarat,ISIStats.allspikes.instate.WAKEstate,...
%     'UniformOutput',false);
% ISIbytheta.WAKE = cat(1,ISIbytheta.WAKE{:});
% ISIbytheta.WAKE = CollapseStruct( ISIbytheta.WAKE,3);
% 
% [ ISIbytheta.REM ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(W);Z(W)],log10([X(W);Y(W)]),...
%     'Xbounds',[0 1],'numXbins',20,'Ybounds',[-3 2],'numYbins',125,'minX',20),...
%     ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
%     ISIStats.allspikes.thetarat,ISIStats.allspikes.instate.REMstate,...
%     'UniformOutput',false);
% ISIbytheta.REM = cat(1,ISIbytheta.REM{:});
% ISIbytheta.REM = CollapseStruct( ISIbytheta.REM,3);


end
end
%%
for nn = [1 3]
figure
subplot(2,2,1)
imagesc(PopRatebyPSS.(normtypes{nn}).(synchtypes{st}).Xbins,...
    PopRatebyPSS.(normtypes{nn}).(synchtypes{st}).Ybins,...
    PopRatebyPSS.(normtypes{nn}).(synchtypes{st}).pYX')
hold on
plot(PopRatebyPSS.(normtypes{nn}).(synchtypes{st}).Xbins,...
    bz_NormToRange(PopRatebyPSS.(normtypes{nn}).(synchtypes{st}).pX),'k')
axis xy
xlabel('PSS');ylabel('Pop Rate (norm)')

subplot(2,2,2)
imagesc(PopRatebyTheta.notNREM.(normtypes{nn}).(synchtypes{st}).Xbins,...
    PopRatebyTheta.notNREM.(normtypes{nn}).(synchtypes{st}).Ybins,...
    PopRatebyTheta.notNREM.(normtypes{nn}).(synchtypes{st}).pYX')
hold on
plot(PopRatebyTheta.notNREM.(normtypes{nn}).(synchtypes{st}).Xbins,...
    bz_NormToRange(PopRatebyTheta.notNREM.(normtypes{nn}).(synchtypes{st}).pX),'k')
axis xy
xlabel('Theta');ylabel('Pop Rate (norm)')


subplot(4,2,6)
imagesc(PopRatebyTheta.WAKE.(normtypes{nn}).(synchtypes{st}).Xbins,...
    PopRatebyTheta.WAKE.(normtypes{nn}).(synchtypes{st}).Ybins,...
    PopRatebyTheta.WAKE.(normtypes{nn}).(synchtypes{st}).pYX')
hold on
plot(PopRatebyTheta.WAKE.(normtypes{nn}).(synchtypes{st}).Xbins,...
    bz_NormToRange(PopRatebyTheta.WAKE.(normtypes{nn}).(synchtypes{st}).pX),'k')
axis xy
ylabel({'WAKE','Pop Rate'})

subplot(4,2,8)
imagesc(PopRatebyTheta.REM.(normtypes{nn}).(synchtypes{st}).Xbins,...
    PopRatebyTheta.REM.(normtypes{nn}).(synchtypes{st}).Ybins,...
    PopRatebyTheta.REM.(normtypes{nn}).(synchtypes{st}).pYX')
hold on
plot(PopRatebyTheta.REM.(normtypes{nn}).(synchtypes{st}).Xbins,...
    bz_NormToRange(PopRatebyTheta.REM.(normtypes{nn}).(synchtypes{st}).pX),'k')
axis xy
xlabel('Theta');ylabel({'REM','Pop Rate'})

    NiceSave(['PopRatebyState_',(normtypes{nn})],figfolder,baseName)

end
%% Examples
%excell = randsample(spikes.numcells,1);
try
for ss = 1:3
exwins.(statenames{ss}) = bz_RandomWindowInIntervals( SleepState.ints.(statenames{ss}),8 );
end

figure
for ss = 1:3
subplot(3,1,ss)
[ ywinrange ] = bz_MultiLFPPlot(lfp,'timewin',exwins.(statenames{ss}),'spikes',spikes,...
    'cellgroups',{CellClass.pE,CellClass.pI},'spikeside','bottom','scaleLFP',0.7);
bz_ScaleBar('s')
hold on
plot(spikemat.timestamps,bz_NormToRange(spikemat.poprate.pE,ywinrange(1).*[1 0.5]),'k')
plot(spikemat.timestamps,bz_NormToRange(spikemat.poprate.pI,ywinrange(1).*[1 0.5]),'r')
xlim(exwins.(statenames{ss}))
ylabel((statenames{ss}))
end
NiceSave('Examples',figfolder,baseName,'figtype','tiff')
catch
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