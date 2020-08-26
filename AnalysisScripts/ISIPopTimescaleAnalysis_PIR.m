function [MutInf] = ...
    ISIPopTimescaleAnalysis_PIR(basePath,figfolder)

%% DEV
%reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
reporoot = '/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/AGData/Cicero/Cicero_09012014';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse24-131213';
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISIPopTimescaleAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
%ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
%OccupancyStats = bz_LoadCellinfo(basePath,'OccupancyStats');
SleepState = bz_LoadStates(basePath,'SleepState');
SleepState.ints.ALL = [0 Inf];
% lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
%     'basepath',basePath);
%% Getting the right cells hack

LFPMapFolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISILFPMap'];
[ISILFPMap] = GetMatResults(LFPMapFolder,'ISILFPMap','baseNames',baseName);
region = 'pir';
lfpchannel = ISILFPMap.MIMap.selectedchans.(region).channel;
usecells = ISILFPMap.MIMap.(ISILFPMap.MIMap.selectedchans.(region).regname).UIDs;
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true,'UID',usecells);


CellClass.keep = ismember(CellClass.UID,usecells);
CellClass.pE = CellClass.pE(CellClass.keep);
CellClass.pI = CellClass.pI(CellClass.keep);
CellClass.UID = CellClass.UID(CellClass.keep);
CellClass.label = CellClass.label(CellClass.keep);
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
clear MutInf
timescales = logspace(-2.5,1.5,17);
timescales = logspace(-2,0.5,6);
for bb = 1:length(timescales)
    bb
binsize = timescales(bb); %s
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
        %Calculate MI - synch/rate and spike/rate
        for tt = 1:length(celltypes)
            MutInf.cellrate.(statenames{ss}).rate.(celltypes{tt})(cc,bb) = ...
                mutualinfo(spikemat.cellrate{cc}(spikemat.instate.(statenames{ss})),...
                spikemat.bycellpoprate.(celltypes{tt}){cc}(spikemat.instate.(statenames{ss})));
            MutInf.cellspike.(statenames{ss}).synch.(celltypes{tt})(cc,bb) = ...
                mutualinfo(spikemat.cellspike{cc}(spikemat.instate.(statenames{ss})),...
                spikemat.bycellpopsynch.(celltypes{tt}){cc}(spikemat.instate.(statenames{ss})));
            
            MutInf.cellrate.(statenames{ss}).synch.(celltypes{tt})(cc,bb) = ...
                mutualinfo(spikemat.cellrate{cc}(spikemat.instate.(statenames{ss})),...
                spikemat.bycellpopsynch.(celltypes{tt}){cc}(spikemat.instate.(statenames{ss})));
            MutInf.cellspike.(statenames{ss}).rate.(celltypes{tt})(cc,bb) = ...
                mutualinfo(spikemat.cellspike{cc}(spikemat.instate.(statenames{ss})),...
                spikemat.bycellpoprate.(celltypes{tt}){cc}(spikemat.instate.(statenames{ss})));
            
            
            MUA.timestamps = spikemat.timestamps;
            MUA.data = spikemat.bycellpoprate.(celltypes{tt}){cc};
            [temp] = ...
                bz_ConditionalISI(spikes.times{cc},MUA,...
                'ints',SleepState.ints.(statenames{ss}),...
                'showfig',false,'GammaFit',false,'ISIDist',false);
            MutInf.cellISI.(statenames{ss}).rate.(celltypes{tt})(cc,bb) = temp.MutInf;

            MUA.data = spikemat.bycellpopsynch.(celltypes{tt}){cc};
                [temp] = ...
                    bz_ConditionalISI(spikes.times{cc},MUA,...
                    'ints',SleepState.ints.(statenames{ss}),...
                    'showfig',false,'GammaFit',false,'ISIDist',false);
            MutInf.cellISI.(statenames{ss}).synch.(celltypes{tt})(cc,bb) = temp.MutInf;
            
        end
    
        %Get the mean rate
        %MutInf.(statenames{ss}).meanRate(cc) = ISIStats.summstats.(statenames{ss}).meanrate(cc);
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
end

%%
synchrate = {'rate','synch'};
spikerate = {'cellrate','cellspike','cellISI'};
%%
for cellsr = 1:3
    for ss = 1:3
        for popsr = 1:2
            for poptt = 1:2
                for celltt = 1:2
meanMI.(spikerate{cellsr}).(statenames{ss}).(synchrate{popsr}).(celltypes{poptt}).(celltypes{celltt}) = ...
    median(MutInf.(spikerate{cellsr}).(statenames{ss}).(synchrate{popsr}).(celltypes{poptt})(MutInf.CellClass.(celltypes{celltt}),:),1);
                end
            end
        end
    end
end
%% Mean all cells (e and I) for all metrics
%%
figure
% for sr = 1:2
% subplot(2,2,sr)
% plot(MutInf.(spikerate{sr}).(statenames{2}).synch.(celltypes{1}),...
%     MutInf.(spikerate{sr}).(statenames{2}).rate.(celltypes{1}),'.')
% axis tight
% box off
% hold on
% legend
% UnityLine
% xlabel('Pop Synch');ylabel('Pop Rate');title(spikerate{sr})
% end

for popsr = 1:2
for cellsr = 1:3
    for celltt = 1:2
subplot(4,3,(popsr-1)*6+cellsr+(celltt-1)*3)
hold on
%plot(log10(timescales),(MutInf.(spikerate{sr}).(statenames{1}).synch.(celltypes{1})(MutInf.CellClass.(celltypes{tt}),:)),'k')
%plot(log10(timescales),(MutInf.(spikerate{sr}).(statenames{1}).synch.(celltypes{2})(MutInf.CellClass.(celltypes{tt}),:)),'r')

for poptt = 1:2
plot(log10(timescales),(meanMI.(spikerate{cellsr}).(statenames{ss}).(synchrate{popsr}).(celltypes{poptt}).(celltypes{celltt})),...
    'color',cellcolor{poptt})
end

LogScale('x',10)
xlabel('Time Window (s)');
if cellsr == 1
    ylabel({celltypes{celltt},'MI'})
end
if celltt==1
title([ spikerate{cellsr}, ' by pop ',(synchrate{popsr})])
end
    end
end
end
NiceSave('CellModulationByPop',figfolder,baseName)

end