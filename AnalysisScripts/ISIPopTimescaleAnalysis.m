function [] = ...
    ISIPopTimescaleAnalysis(basePath,figfolder)

%% DEV
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%reporoot = '/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/AGData/Cicero/Cicero_09012014';
basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse24-131213';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISIPopTimescaleAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
%ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
%OccupancyStats = bz_LoadCellinfo(basePath,'OccupancyStats');
SleepState = bz_LoadStates(basePath,'SleepState');
SleepState.ints.ALL = [0 Inf];
% lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
%     'basepath',basePath);
%%
cellinfo.CellClass = CellClass;
%cellinfo.ISIStats = ISIStats.summstats;
%cellinfo.OccupancyStats = OccupancyStats;
%% Load The Gamma stuff
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');


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

timescales = logspace(-2.5,1,15)
binsize = 0.1; %s
dt = 0.005;
spikemat = bz_SpktToSpkmat(spikes,'binsize',binsize,'dt',dt,'bintype','gaussian','units','rate');

for ss = 1:3
    spikemat.instate.(statenames{ss}) = InIntervals(spikemat.timestamps,SleepState.ints.(statenames{ss}));
end
%% For each cell, calculate E and I pop rates of all OTHER cells
ncellthresh = 5;

for tt = 1:length(celltypes)
    Ncells.(celltypes{tt}) = sum(CellClass.(celltypes{tt}));
    spikemat.totpoprate.(celltypes{tt}) = sum(spikemat.data(:,CellClass.(celltypes{tt})),2);
    spikemat.poprate.(celltypes{tt}) = spikemat.totpoprate.(celltypes{tt})./Ncells.(celltypes{tt});
    spikemat.totpopsynch.(celltypes{tt}) = sum(spikemat.data(:,CellClass.(celltypes{tt}))>0.5,2);
    spikemat.popsynch.(celltypes{tt}) = spikemat.totpopsynch.(celltypes{tt})./Ncells.(celltypes{tt});
    
%     if Ncells.(celltypes{tt})==0
%         celltypes(tt) = [];
%     end
end
Ncells.ALL = Ncells.pE + Ncells.pI;
spikemat.totpoprate.ALL = sum(spikemat.data(:,(CellClass.pI|CellClass.pE)),2);
spikemat.poprate.ALL = spikemat.totpoprate.ALL./Ncells.ALL;

for cc = 1:spikes.numcells %weird roundabout way to calculate is much faster
    bz_Counter(cc,spikes.numcells,'Cell');
    thiscell = false(size(CellClass.pE));
    thiscell(cc) = true;
    spikemat.cellrate{cc} = spikemat.data(:,cc);
    spikemat.cellspike{cc} = double(spikemat.data(:,cc)>0.5);
    for tt = 1:length(celltypes)
        MutInf.cellcount.(celltypes{tt})(cc) = sum(CellClass.(celltypes{tt}) & ~thiscell);
        if CellClass.(celltypes{tt})(cc) %if it's in theclass, subtract off the current cell
            spikemat.bycellpoprate.(celltypes{tt}){cc} = ...
                (spikemat.totpoprate.(celltypes{tt})-spikemat.cellrate{cc})./...
               MutInf.cellcount.(celltypes{tt})(cc);
           
            spikemat.bycellpopsynch.(celltypes{tt}){cc} = ...
                (spikemat.totpopsynch.(celltypes{tt})-spikemat.cellspike{cc})./...
               MutInf.cellcount.(celltypes{tt})(cc);
        else
            spikemat.bycellpoprate.(celltypes{tt}){cc} = ...
                spikemat.totpoprate.(celltypes{tt})./MutInf.cellcount.(celltypes{tt})(cc);
            
            spikemat.bycellpopsynch.(celltypes{tt}){cc} = ...
                spikemat.totpopsynch.(celltypes{tt})./MutInf.cellcount.(celltypes{tt})(cc);
        end
    end
    
    if CellClass.pI(cc)||CellClass.pE(cc)
        spikemat.bycellpoprate.ALL{cc} = (spikemat.totpoprate.ALL-spikemat.cellrate{cc})./...
            (Ncells.ALL-1);
    else
        spikemat.bycellpoprate.ALL{cc} = spikemat.totpoprate.ALL./Ncells.ALL;
    end
    
    

end

%%
% figure
% plot(spikemat.poprate.(celltypes{1}),spikemat.popsynch.(celltypes{1}),'.')
%% Correlation and Gamma parms

for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell');
    for ss = 1:3
        %Calculate Correlation
        for tt = 1:length(celltypes)
            MutInf.cellrate.(statenames{ss}).rate.(celltypes{tt})(cc) = ...
                mutualinfo(spikemat.cellrate{cc}(spikemat.instate.(statenames{ss})),...
                spikemat.bycellpoprate.(celltypes{tt}){cc}(spikemat.instate.(statenames{ss})));
            MutInf.cellspike.(statenames{ss}).synch.(celltypes{tt})(cc) = ...
                mutualinfo(spikemat.cellspike{cc}(spikemat.instate.(statenames{ss})),...
                spikemat.bycellpopsynch.(celltypes{tt}){cc}(spikemat.instate.(statenames{ss})));
            
            MutInf.cellrate.(statenames{ss}).synch.(celltypes{tt})(cc) = ...
                mutualinfo(spikemat.cellrate{cc}(spikemat.instate.(statenames{ss})),...
                spikemat.bycellpopsynch.(celltypes{tt}){cc}(spikemat.instate.(statenames{ss})));
            MutInf.cellspike.(statenames{ss}).rate.(celltypes{tt})(cc) = ...
                mutualinfo(spikemat.cellspike{cc}(spikemat.instate.(statenames{ss})),...
                spikemat.bycellpoprate.(celltypes{tt}){cc}(spikemat.instate.(statenames{ss})));
        end
    
        %Get the mean rate
        MutInf.(statenames{ss}).meanRate(cc) = ISIStats.summstats.(statenames{ss}).meanrate(cc);
        %Get the GS rate
        cellUID(cc) = spikes.UID(cc);
        GFIDX = find(GammaFit.(statenames{ss}).cellstats.UID==cellUID(cc));
        if isempty(GFIDX)
            MutInf.(statenames{ss}).GSrate(cc) = nan;
            MutInf.(statenames{ss}).GSCV(cc) = nan;
            MutInf.(statenames{ss}).GSweight(cc) = nan;
            continue
        end
        cellGamma = GammaFit.(statenames{ss}).singlecell(GFIDX);
        MutInf.(statenames{ss}).GSrate(cc) = cellGamma.GSlogrates;
        MutInf.(statenames{ss}).GSCV(cc) = cellGamma.GSCVs;
        MutInf.(statenames{ss}).GSweight(cc) = cellGamma.GSweights;
    end
end
MutInf.CellClass = CellClass;

%%
synchrate = {'rate','synch'};
%%
figure
for tt = 1:length(celltypes)
    for ss = 1:3
        for sr = 1:2
        subplot(4,3,(tt-1)*3+(sr-1).*6+ss)
        plot(MutInf.(statenames{ss}).GSrate,MutInf.(statenames{ss}).(synchrate{sr}).(celltypes{tt}),'k.')
        hold on
        plot(log10(MutInf.(statenames{ss}).meanRate(CellClass.pI)),...
            MutInf.(statenames{ss}).(synchrate{sr}).(celltypes{tt})(CellClass.pI),'r.')
        axis tight
        box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS/Mean Rate (Hz)');
        LogScale('x',10,'exp',true,'nohalf',true)
        if ss ==1
            ylabel([(celltypes{tt}),' ',(synchrate{sr}),' Corr'])
        end
        if tt == 1
            title(statenames{ss})
        end
        end
%         subplot(4,3,(tt-1)*3+ss+6)
%         plot(MutInf.(statenames{ss}).GSCV,MutInf.(statenames{ss}).rate.(celltypes{tt}),'.')
%         hold on
%         axis tight
%         box off
%         plot(xlim(gca),[0 0],'k--')
%         xlabel('GS CV');ylabel([(celltypes{tt}),' Corr'])
%         %LogScale('x',10,'exp',false,'nohalf',true)
        
    end 
    
end

NiceSave('MUACorrandGSRate',figfolder,baseName)

%% Conditional ISI distrobution 

%Should implement number cells threshold...
clear MUAConditionalISIDist_all
for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell');    
    for tt = 1:length(celltypes)
%         if Ncells.(celltypes{tt})==0
%             continue
%         end
        MUA.timestamps = spikemat.timestamps;
        MUA.data = spikemat.bycellpoprate.(celltypes{tt}){cc};
        for ss = 1:3
            [MUAConditionalISIDist_all.(statenames{ss}).rate.(celltypes{tt})(cc)] = ...
                bz_ConditionalISI(spikes.times{cc},MUA,...
                'ints',SleepState.ints.(statenames{ss}),...
                'showfig',false,'GammaFit',false,'minX',0,'numISIbins',100);
        end
        
        MUA.data = spikemat.bycellpopsynch.(celltypes{tt}){cc};
        for ss = 1:3
            [MUAConditionalISIDist_all.(statenames{ss}).synch.(celltypes{tt})(cc)] = ...
                bz_ConditionalISI(spikes.times{cc},MUA,...
                'ints',SleepState.ints.(statenames{ss}),...
                'showfig',false,'GammaFit',false,'minX',50,'numISIbins',100,...
                'normtype','none','Xwin',[0 1],'numXbins',20);
        end
        
    end
end

%%
for ss = 1:3
    for tt = 1:length(celltypes)
    for sr = 1:2
        MUAConditionalISIDist.(statenames{ss}).(synchrate{sr}).(celltypes{tt}) = ...
            bz_CollapseStruct(MUAConditionalISIDist_all.(statenames{ss}).(synchrate{sr}).(celltypes{tt}),3,'justcat',true);

        for tt2 = 1:length(celltypes)
            if sum(CellClass.(celltypes{tt2}))==0
                continue
            end
            MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}) = ...
                bz_CollapseStruct(MUAConditionalISIDist_all.(statenames{ss}).(synchrate{sr}).(celltypes{tt})(CellClass.(celltypes{tt2})),...
                3,'mean',true);
        end
    end
    end
end

%%
excell = 1;
figure
%imagesc(MUAConditionalISIDist.(statenames{2}).(synchrate{2}).(celltypes{1}).Dist.pYX(:,:,excell))
plot(MUAConditionalISIDist.(statenames{2}).(synchrate{2}).(celltypes{1}).Dist.Xhist(:,:,excell))
%% 
for sr = 1:2 
figure
for ss = 1:3
for tt = 1:length(celltypes) %Pop
    for tt2 = 1:length(celltypes) %Ref

        try
    subplot(4,3,(tt2-1)*6+(tt-1)*3+ss)
        imagesc(MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Ybins,...
            MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX')
        hold on
        plot(MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            -log10(MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.SpikeRate),...
            'r','LineWidth',2)
           
        xlabel([(celltypes{tt}),' ',(synchrate{sr})]);ylabel([(celltypes{tt2}),' ISI (s)'])
        LogScale('y',10,'exp',true,'nohalf',true)
        bz_AddRightRateAxis
        if tt ==1 & tt2 == 1
            title(statenames{ss})
        end
        catch
        end
    end 
end
end


NiceSave(['ISIbyMUA_',(synchrate{sr})],figfolder,baseName)

end
end