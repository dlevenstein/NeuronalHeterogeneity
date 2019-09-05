function [PopMod,cellinfo,PopCellCorr] = PopActivityModulationAnalysis(basePath,figfolder)

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
SleepState.ints.ALL = [0 Inf];
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
cellcolor = {'k','r'};
statenames = fieldnames(SleepState.ints);
synchtypes = [celltypes,'ALL'];


for ss = 1:length(statenames)
    ISIStats.allspikes.instate.(statenames{ss}) = cellfun(@(X) InIntervals(X,double(SleepState.ints.(statenames{ss}))),...
        ISIStats.allspikes.times,'UniformOutput',false);
end
%% Loop binsizes
numbins = 20;
binrange = [0.003 30];
dt = 0.003;
binsizes = logspace(log10(binrange(1)),log10(binrange(2)),numbins);
%%
for bb = 1:numbins
    bb
%% Calculate spike count matrix
%binsize = 0.1; %s
clear spikemat
spikemat = bz_SpktToSpkmat(spikes,'binsize',binsizes(bb),'dt',dt,'bintype','gaussian','units','rate');

%% For each cell, calculate E and I pop rates of all OTHER cells

for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell');
    thiscell = false(size(CellClass.pE));
    thiscell(cc) = true;
    spikemat.cellrate{cc} = spikemat.data(:,cc);
    for tt = 1:length(celltypes)
        cellinfo.Ncells.(celltypes{tt}) = sum(CellClass.(celltypes{tt}));
        PopMod.Ncells.(celltypes{tt})(cc) = sum(CellClass.(celltypes{tt}) & ~thiscell);
        spikemat.bycellpoprate.(celltypes{tt}){cc} = mean(spikemat.data(:,CellClass.(celltypes{tt}) & ~thiscell),2);
    end
    spikemat.bycellpoprate.ALL{cc} = mean(spikemat.data(:,(CellClass.pI|CellClass.pE) & ~thiscell),2);
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
            'Ybounds',[-1 1],'numYbins',50,'Xbounds',[-3 2],'numXbins',125,'minX',100),...
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
            'Ybounds',[-1 1],'numYbins',50,'Xbounds',[-4 1],'numXbins',100,'minX',100),...
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
        
        PopCellCorr.(synchtypes{st}).(statenames{ss})(bb,:) = cellfun(@(CR,PR) ...
            corr(CR(spikemat.instate.(statenames{ss})),PR(spikemat.instate.(statenames{ss})),'type','spearman'),...
            spikemat.cellrate,spikemat.bycellpoprate.(synchtypes{st}));
        
        PopMod.(synchtypes{st}).(statenames{ss}).allcells(bb,:,:) = SynchbyISI.(synchtypes{st}).(statenames{ss}).meanYX;
        PopMod_MTO.(synchtypes{st}).(statenames{ss}).allcells(bb,:,:) = SynchbynormISI.(synchtypes{st}).(statenames{ss}).meanYX;
    end

end
end
%%
PopMod.bins.ISIbins = SynchbyISI.(synchtypes{st}).(statenames{ss}).Xbins(1,:,1);
PopMod.bins.BinSizeBins = binsizes;

PopMod_MTO.bins.ISIbins = SynchbynormISI.(synchtypes{st}).(statenames{ss}).Xbins(1,:,1);
PopMod_MTO.bins.BinSizeBins = binsizes;

PopCellCorr.bins.BinSizeBins = PopMod.bins.BinSizeBins;
%% Population average modulation
for ss = 1:3
    for st = 1:length(synchtypes)
        for cc = 1:length(celltypes)
                PopMod.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
                    nanmean(PopMod.(synchtypes{st}).(statenames{ss}).allcells(:,:,CellClass.(celltypes{cc})),3); 
                PopMod_MTO.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
                    nanmean(PopMod_MTO.(synchtypes{st}).(statenames{ss}).allcells(:,:,CellClass.(celltypes{cc})),3); 
%                 PopCellCorr.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
%                     nanmean(PopCellCorr.(synchtypes{st}).(statenames{ss}).allcells(:,CellClass.(celltypes{cc})),3); 
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
    subplot(3,3,ss+(st-1)*3+(cc-1)*6)
plot(log10(PopCellCorr.bins.BinSizeBins),PopCellCorr.(synchtypes{st}).(statenames{ss})')
hold on
plot(log10(PopCellCorr.bins.BinSizeBins([1 end])),[0 0],'k')
%caxis([-1 1])
ylim([-1 1])
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