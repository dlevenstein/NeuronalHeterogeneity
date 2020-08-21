function [PopCorr,MUAConditionalISIDist,popratehist] = ...
    ISIModesbyPopActivityAnalysis(basePath,figfolder)

%% DEV
%reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%reporoot = '/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/AGData/Cicero/Cicero_09012014';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse24-131213';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISIModesbyPopActivityAnalysis'];
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
binsize = 0.1; %s
dt = 0.01;
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
    spikemat.cellsync.(celltypes{tt}) = mean(spikemat.data(:,CellClass.(celltypes{tt}))>0.5,2);
    
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
    for tt = 1:length(celltypes)
        PopCorr.cellcount.(celltypes{tt})(cc) = sum(CellClass.(celltypes{tt}) & ~thiscell);
        if CellClass.(celltypes{tt})(cc) %if it's in theclass, subtract off the current cell
            spikemat.bycellpoprate.(celltypes{tt}){cc} = ...
                (spikemat.totpoprate.(celltypes{tt})-spikemat.cellrate{cc})./...
               PopCorr.cellcount.(celltypes{tt})(cc);
        else
            spikemat.bycellpoprate.(celltypes{tt}){cc} = ...
                spikemat.totpoprate.(celltypes{tt})./PopCorr.cellcount.(celltypes{tt})(cc);
        end
    end
    
    if CellClass.pI(cc)||CellClass.pE(cc)
        spikemat.bycellpoprate.ALL{cc} = (spikemat.totpoprate.ALL-spikemat.cellrate{cc})./...
            (Ncells.ALL-1);
    else
        spikemat.bycellpoprate.ALL{cc} = spikemat.totpoprate.ALL./Ncells.ALL;
    end
    
    

end

%% Correlation and Gamma parms

for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell');
    for ss = 1:3
        %Calculate Correlation
        for tt = 1:length(celltypes)
            PopCorr.(statenames{ss}).(celltypes{tt})(cc) = ...
                corr(spikemat.cellrate{cc}(spikemat.instate.(statenames{ss})),...
                spikemat.bycellpoprate.(celltypes{tt}){cc}(spikemat.instate.(statenames{ss})),...
                'type','spearman');
        end

        %Get the GS rate
        cellUID(cc) = spikes.UID(cc);
        GFIDX = find(GammaFit.(statenames{ss}).cellstats.UID==cellUID(cc));
        if isempty(GFIDX)
            PopCorr.(statenames{ss}).GSrate(cc) = nan;
            PopCorr.(statenames{ss}).GSCV(cc) = nan;
            PopCorr.(statenames{ss}).GSweight(cc) = nan;
            continue
        end
        cellGamma = GammaFit.(statenames{ss}).singlecell(GFIDX);
        PopCorr.(statenames{ss}).GSrate(cc) = cellGamma.GSlogrates;
        PopCorr.(statenames{ss}).GSCV(cc) = cellGamma.GSCVs;
        PopCorr.(statenames{ss}).GSweight(cc) = cellGamma.GSweights;
    end
end
PopCorr.CellClass = CellClass;
%%
figure
for tt = 1:length(celltypes)
    for ss = 1:3
        subplot(4,4,(tt-1)*4+ss)
        plot(PopCorr.(statenames{ss}).GSrate,PopCorr.(statenames{ss}).(celltypes{tt}),'.')
        hold on
        axis tight
        box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS Rate (Hz)');ylabel([(celltypes{tt}),' Corr'])
        LogScale('x',10,'exp',true,'nohalf',true)
        if tt == 1
            title(statenames{ss})
        end
        
        subplot(4,3,(tt-1)*3+ss+6)
        plot(PopCorr.(statenames{ss}).GSCV,PopCorr.(statenames{ss}).(celltypes{tt}),'.')
        hold on
        axis tight
        box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS CV');ylabel([(celltypes{tt}),' Corr'])
        %LogScale('x',10,'exp',false,'nohalf',true)
        
    end 
    
end

NiceSave('MUACorrandGSRate',figfolder,baseName)

%% Conditional ISI distrobution 

%Should implement number cells threshold...

for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell');    
    for tt = 1:length(celltypes)
%         if Ncells.(celltypes{tt})==0
%             continue
%         end
        MUA.timestamps = spikemat.timestamps;
        MUA.data = spikemat.bycellpoprate.(celltypes{tt}){cc};
        for ss = 1:3
            [MUAConditionalISIDist_all.(statenames{ss}).(celltypes{tt})(cc)] = ...
                bz_ConditionalISI(spikes.times{cc},MUA,...
                'ints',SleepState.ints.(statenames{ss}),...
                'showfig',false,'GammaFit',false,'minX',0,'numISIbins',100);
        end
        
    end
end

%%
for ss = 1:3
    for tt = 1:length(celltypes)

        MUAConditionalISIDist.(statenames{ss}).(celltypes{tt}) = ...
            bz_CollapseStruct(MUAConditionalISIDist_all.(statenames{ss}).(celltypes{tt}),3,'justcat',true);

        for tt2 = 1:length(celltypes)
            if sum(CellClass.(celltypes{tt2}))==0
                continue
            end
            MeanCondISI.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}) = ...
                bz_CollapseStruct(MUAConditionalISIDist_all.(statenames{ss}).(celltypes{tt})(CellClass.(celltypes{tt2})),...
                3,'mean',true);
        end
        
    end
end


%% 

figure
for ss = 1:3
for tt = 1:length(celltypes) %Pop
    for tt2 = 1:length(celltypes) %Ref

        try
    subplot(4,3,(tt2-1)*6+(tt-1)*3+ss)
        imagesc(MeanCondISI.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            MeanCondISI.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.Ybins,...
            MeanCondISI.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX')
        hold on
        plot(MeanCondISI.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            -log10(MeanCondISI.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.SpikeRate),...
            'r','LineWidth',2)
           
        xlabel([(celltypes{tt}),' Rate']);ylabel([(celltypes{tt2}),' ISI (s)'])
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


NiceSave('ISIbyMUA',figfolder,baseName)


end