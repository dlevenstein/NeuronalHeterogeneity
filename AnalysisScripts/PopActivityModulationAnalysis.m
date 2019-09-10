function [PopMod,PopMod_MTO,cellinfo,PopCellCorr,MutInfo] =...
    PopActivityModulationAnalysis(basePath,figfolder)

%% DEV
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/AGData/Cicero/Cicero_09012014';
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/PopActivityModulationAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
OccupancyStats = bz_LoadCellinfo(basePath,'OccupancyStats');
SleepState = bz_LoadStates(basePath,'SleepState');
%SleepState.ints.ALL = [0 Inf];
% lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
%     'basepath',basePath);

%%
cellinfo.CellClass = CellClass;
cellinfo.ISIStats = ISIStats.summstats;
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
cellcolor = {[0 0 0],[1 0 0]};
statenames = fieldnames(SleepState.ints);
synchtypes = [celltypes,'ALL'];


for ss = 1:length(statenames)
    ISIStats.allspikes.instate.(statenames{ss}) = cellfun(@(X) InIntervals(X,double(SleepState.ints.(statenames{ss}))),...
        ISIStats.allspikes.times,'UniformOutput',false);
end

%%
alltime = [0 ...
    max([SleepState.ints.NREMstate(:,2);SleepState.ints.WAKEstate(:,2);SleepState.ints.REMstate(:,2)])];
%% Loop binsizes
numbins = 25;
binrange = [0.002 20];
mindt = 0.002; %minimum acceptable dt
% numbins = 10;
% binrange = [0.005 10];
% mindt = 0.005; %minimum acceptable dt
binsizes = logspace(log10(binrange(1)),log10(binrange(2)),numbins);
dt = max(binsizes/20,mindt);
%%
for bb = 1:numbins
    bb
%% Calculate spike count matrix
%binsize = 0.1; %s
clear spikemat
spikemat = bz_SpktToSpkmat(spikes,'binsize',binsizes(bb),'dt',dt(bb),...
    'bintype','gaussian','units','rate','win',alltime);

%% For each cell, calculate E and I pop rates of all OTHER cells

for tt = 1:length(celltypes)
    cellinfo.Ncells.(celltypes{tt}) = sum(CellClass.(celltypes{tt}));
    spikemat.totpoprate.(celltypes{tt}) = sum(spikemat.data(:,CellClass.(celltypes{tt})),2);
end
spikemat.totpoprate.ALL = sum(spikemat.data(:,(CellClass.pI|CellClass.pE)),2);
cellinfo.Ncells.ALL = cellinfo.Ncells.pE + cellinfo.Ncells.pI;

for cc = 1:spikes.numcells %weird roundabout way to calculate is much faster
    bz_Counter(cc,spikes.numcells,'Cell');
    thiscell = false(size(CellClass.pE));
    thiscell(cc) = true;
    spikemat.cellrate{cc} = spikemat.data(:,cc);
    for tt = 1:length(celltypes)
        PopMod.Ncells.(celltypes{tt})(cc) = sum(CellClass.(celltypes{tt}) & ~thiscell);
        if CellClass.(celltypes{tt})(cc) %if it's in theclass, subtract off the current cell
            spikemat.bycellpoprate.(celltypes{tt}){cc} = ...
                (spikemat.totpoprate.(celltypes{tt})-spikemat.cellrate{cc})./...
                PopMod.Ncells.(celltypes{tt})(cc);
        else
            spikemat.bycellpoprate.(celltypes{tt}){cc} = ...
                spikemat.totpoprate.(celltypes{tt})./PopMod.Ncells.(celltypes{tt})(cc);
        end
    end
    
    if CellClass.pI(cc)||CellClass.pE(cc)
        spikemat.bycellpoprate.ALL{cc} = (spikemat.totpoprate.ALL-spikemat.cellrate{cc})./...
            (cellinfo.Ncells.ALL-1);
    else
        spikemat.bycellpoprate.ALL{cc} = spikemat.totpoprate.ALL./cellinfo.Ncells.ALL;
    end
end
%% Normalizations

for tt = 1:length(synchtypes)
    spikemat.bycellpoprate.(synchtypes{tt}) = cellfun(@(X) ...
        (X./mean(X)),spikemat.bycellpoprate.(synchtypes{tt}),'UniformOutput',false);
end

%% Calculate E and I pop rate (of other cells) for each spike

for tt = 1:length(synchtypes)
    ISIStats.allspikes.poprate.(synchtypes{tt}) = ...
        cellfun(@(X,Y) interp1(spikemat.timestamps,X,Y,'nearest'),...
        spikemat.bycellpoprate.(synchtypes{tt}),ISIStats.allspikes.times,...
        'UniformOutput',false);
end

%%
for ss = 1:length(statenames)
    spikemat.instate.(statenames{ss}) = InIntervals(spikemat.timestamps,double(SleepState.ints.(statenames{ss})));
end 
%% Calculate Conditional distributions on synchrony (pop rate) in each state
for ss = 1:3

    for st = 1:length(synchtypes)

        [ SynchbyISI.(synchtypes{st}).(statenames{ss}) ] = cellfun(@(X,Y,Z,W) ...
            ConditionalHist( log10([X(W);Y(W)]),([Z(W);Z(W)]),...
            'Ybounds',[-1 1],'numYbins',50,'Xbounds',[-3 2],'numXbins',100,'minX',50),...
            ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
            ISIStats.allspikes.poprate.(synchtypes{st}),...
            ISIStats.allspikes.instate.(statenames{ss}),...
            'UniformOutput',false);
        SynchbyISI.(synchtypes{st}).(statenames{ss}) = ...
            cat(1,SynchbyISI.(synchtypes{st}).(statenames{ss}){:});
        SynchbyISI.(synchtypes{st}).(statenames{ss}) = ...
            bz_CollapseStruct( SynchbyISI.(synchtypes{st}).(statenames{ss}),3);

        [ SynchbynormISI.(synchtypes{st}).(statenames{ss}) ] = cellfun(@(X,Y,Z,W,MTO) ...
            ConditionalHist( log10([X(W);Y(W)]./MTO),([Z(W);Z(W)]),...
            'Ybounds',[-1 1],'numYbins',50,'Xbounds',[-4 1],'numXbins',100,'minX',50),...
            ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
            ISIStats.allspikes.poprate.(synchtypes{st}),...
            ISIStats.allspikes.instate.(statenames{ss}),...
            num2cell(OccupancyStats.(statenames{ss}).median),...
            'UniformOutput',false);
        SynchbynormISI.(synchtypes{st}).(statenames{ss}) = ...
            cat(1,SynchbynormISI.(synchtypes{st}).(statenames{ss}){:});
        SynchbynormISI.(synchtypes{st}).(statenames{ss}) = ...
            bz_CollapseStruct( SynchbynormISI.(synchtypes{st}).(statenames{ss}),3);
        
%         for cc = 1:length(celltypes)
%             SynchbyISI.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) = nanmean(SynchbyISI.(synchtypes{st}).(statenames{ss}).pYX(:,:,CellClass.(celltypes{cc})),3);
%             SynchbyISI.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) = nanmean(SynchbyISI.(synchtypes{st}).(statenames{ss}).pYX(:,:,CellClass.(celltypes{cc})),3);
%             SynchbyISI.(synchtypes{st}).(statenames{ss}).celltypeidx.(celltypes{cc}) = CellClass.(celltypes{cc});
%         end
        
        PopCellCorr.(synchtypes{st}).(statenames{ss}).allcells(bb,:) = cellfun(@(CR,PR) ...
            corr(CR(spikemat.instate.(statenames{ss})),PR(spikemat.instate.(statenames{ss})),'type','spearman'),...
            spikemat.cellrate,spikemat.bycellpoprate.(synchtypes{st}));
        
        PopMod.(synchtypes{st}).(statenames{ss}).allcells(bb,:,:) = SynchbyISI.(synchtypes{st}).(statenames{ss}).meanYX;
        PopMod_MTO.(synchtypes{st}).(statenames{ss}).allcells(bb,:,:) = SynchbynormISI.(synchtypes{st}).(statenames{ss}).meanYX;
        PopMod.(synchtypes{st}).(statenames{ss}).pISI(bb,:,:) = SynchbyISI.(synchtypes{st}).(statenames{ss}).pX;
        PopMod_MTO.(synchtypes{st}).(statenames{ss}).pISI(bb,:,:) = SynchbynormISI.(synchtypes{st}).(statenames{ss}).pX;
        PopMod.(synchtypes{st}).(statenames{ss}).nISI(bb,:,:) = SynchbyISI.(synchtypes{st}).(statenames{ss}).Xhist;
        PopMod_MTO.(synchtypes{st}).(statenames{ss}).nISI(bb,:,:) = SynchbynormISI.(synchtypes{st}).(statenames{ss}).Xhist;
        
        %Mutual info
        MutInfo.(synchtypes{st}).(statenames{ss}).allcells(bb,:) = ...
            cellfun(@(ISIs,ISInp1,PopRate,instate) ...
            mutualinfo(log10([ISIs(instate);ISInp1(instate)]),([PopRate(instate);PopRate(instate)])),...
            ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
            ISIStats.allspikes.poprate.(synchtypes{st}),...
            ISIStats.allspikes.instate.(statenames{ss}));
    end

end
end
%%
PopMod.bins.ISIbins = SynchbyISI.(synchtypes{st}).(statenames{ss}).Xbins(1,:,1);
PopMod.bins.BinSizeBins = binsizes;

PopMod_MTO.bins.ISIbins = SynchbynormISI.(synchtypes{st}).(statenames{ss}).Xbins(1,:,1);
PopMod_MTO.bins.BinSizeBins = binsizes;

PopCellCorr.bins.BinSizeBins = PopMod.bins.BinSizeBins;
MutInfo.bins.BinSizeBins = PopMod.bins.BinSizeBins;
%% Population average modulation
for ss = 1:3
    for st = 1:length(synchtypes)
        for cc = 1:length(celltypes)
                PopMod.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
                    nanmean(PopMod.(synchtypes{st}).(statenames{ss}).allcells(:,:,CellClass.(celltypes{cc})),3); 
                PopMod_MTO.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
                    nanmean(PopMod_MTO.(synchtypes{st}).(statenames{ss}).allcells(:,:,CellClass.(celltypes{cc})),3); 
                PopCellCorr.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
                    nanmean(PopCellCorr.(synchtypes{st}).(statenames{ss}).allcells(:,CellClass.(celltypes{cc})),2); 
                MutInfo.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) = ...
                    nanmean(MutInfo.(synchtypes{st}).(statenames{ss}).allcells(:,CellClass.(celltypes{cc})),2);
        end
    end
end

 
%%
figure
for ss = 1:3
for st = 1:2
for cc = 1:length(celltypes)
    subplot(4,3,ss+(st-1)*3+(cc-1)*6)
        imagesc((PopMod.bins.ISIbins),log10(PopMod.bins.BinSizeBins),...
            PopMod.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}))
        %colorbar
        
        caxis([0.9 1.8])
        %crameri('vik','pivot',1)
      
        %LogScale('y',2)
         LogScale('xy',10,'exp',true);
        axis xy
        
        if cc==1 &st==1
           title(statenames{ss}) 
        end
        if ss == 1
            ylabel({[(synchtypes{st}),' Modulation'],'Bin Size'})
        end
        if st == 2
            xlabel([(celltypes{cc}),' ISI (s)'])
        end
        
end
end 
end

NiceSave('PopRateModulation',figfolder,baseName)
%%
figure
for ss = 1:3
for st = 1:2
for cc = 1:length(celltypes)
    subplot(4,3,ss+(st-1)*3+(cc-1)*6)
        imagesc((PopMod_MTO.bins.ISIbins),log10(PopMod_MTO.bins.BinSizeBins),...
            PopMod_MTO.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}))
        %colorbar
        
        caxis([0.9 1.8])
        %crameri('vik','pivot',1)
      
        %LogScale('y',2)
         LogScale('xy',10,'exp',true);
        axis xy
        
        if cc==1 &st==1
           title(statenames{ss}) 
        end
        if ss == 1
            ylabel({[(synchtypes{st}),' Modulation'],'Bin Size'})
        end
        if st == 2
            xlabel([(celltypes{cc}),' ISI (MTOnorm)'])
        end
        
end
end 
end

NiceSave('PopRateModulation_MTO',figfolder,baseName)

%%
cc =1
figure
for ss = 1:3
for st = 1:3
%for cc = 1:length(celltypes)
    subplot(3,3,ss+(st-1)*3)
    hold on
    for cc = 1:length(celltypes)
plot(log10(PopCellCorr.bins.BinSizeBins),PopCellCorr.(synchtypes{st}).(statenames{ss}).allcells(:,CellClass.(celltypes{cc}))',...
    'color',min(1,cellcolor{cc}+0.5),'linewidth',0.1)

plot(log10(PopCellCorr.bins.BinSizeBins),PopCellCorr.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}),...
    'color',cellcolor{cc},'linewidth',2)
    end
    plot(log10(PopCellCorr.bins.BinSizeBins([1 end])),[0 0],'k')
%caxis([-1 1])
%ylim([-1 1])
xlim(log10(PopCellCorr.bins.BinSizeBins([1 end])))
xlabel('TimeScale');

        if ss == 1
            ylabel([(synchtypes{st}),' Corr'])
        end
        if cc==1 &st==1
           title(statenames{ss}) 
        end
%crameri('vik','pivot',0)
end
%end
end
NiceSave('PopRateCorr',figfolder,baseName)

%%
cc =1
figure
for ss = 1:3
for st = 1:3

    subplot(3,3,ss+(st-1)*3)
    hold on
    for cc = 1:length(celltypes)
plot(log10(MutInfo.bins.BinSizeBins),MutInfo.(synchtypes{st}).(statenames{ss}).allcells(:,CellClass.(celltypes{cc}))',...
    'color',min(1,cellcolor{cc}+0.5),'linewidth',0.1)

plot(log10(MutInfo.bins.BinSizeBins),MutInfo.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}),...
    'color',cellcolor{cc},'linewidth',2)
    end
%caxis([-1 1])
%ylim([-1 1])
xlim(log10(MutInfo.bins.BinSizeBins([1 end])))
xlabel('TimeScale');
LogScale('x',10)

        if ss == 1
            ylabel([(synchtypes{st}),' MI'])
        end
        if cc==1 &st==1
           title(statenames{ss}) 
        end
%crameri('vik','pivot',0)
end
%end
end
NiceSave('ISIMITimeScale',figfolder,baseName)
%%
exbin = 3;
figure
for ss = 1:3
for st = 1:3
subplot(3,3,ss+(st-1)*3)
h = imagesc(squeeze(PopMod.(synchtypes{st}).(statenames{ss}).allcells(exbin,:,ISIStats.sorts.(statenames{ss}).ratebyclass))');
 set(h,'AlphaData',50*(squeeze(PopMod.(synchtypes{st}).(statenames{ss}).pISI(exbin,:,ISIStats.sorts.(statenames{ss}).ratebyclass))'));
% PopMod.(synchtypes{st}).(statenames{ss}).pISI(bb,:,:)
%axis xy
colorbar
caxis([0.8 2])
end
end
%%
% pickbinsize = 5;
% figure
% imagesc(squeeze(PopMod.pE.NREMstate.allcells(pickbinsize,:,:))')
%%
% figure
% for ss = 1:3
% for tt = 1:length(celltypes)
% for st = 1:3
% subplot(6,3,(ss-1)+(tt-1)*3+(st-1)*6+1)
%     
%     imagesc(SynchbyISI.(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),SynchbyISI.(synchtypes{st}).(statenames{ss}).Ybins(1,:,1), SynchbyISI.(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
%     hold on
%     plot(SynchbyISI.(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),SynchbyISI.(synchtypes{st}).(statenames{ss}).meanYX,'w')
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
%     xlim(SynchbyISI.(synchtypes{tt}).(statenames{ss}).Xbins(1,[1 end],1))
% end 
% end
% end







end