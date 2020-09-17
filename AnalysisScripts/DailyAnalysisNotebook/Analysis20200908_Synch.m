function [ ] = AnalysisXXXXXXXX(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%% DEV
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%reporoot = '/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/AGData/Cicero/Cicero_09012014';
basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse24-131213';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/YMV12_171211';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
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



%% Cell types and states
% try
% celltypes = CellClass.celltypes;
% catch
%     celltypes = unique(CellClass.label);
% end
celltypes = {'pE','pI'};
cellcolor = {'k','r'};
statenames = fieldnames(SleepState.ints);
statecolors = {[0 0 0],[0 0 1],[1 0 0]};

%% Calculate spike count matrix
bins = logspace(-2,1,21);
%bins = 0.08;
clear synchhist
for bb = 1:length(bins)
    bz_Counter(bb,length(bins),'Bin')
binsize = bins(bb); %s
dt = 0.01;
spikemat = bz_SpktToSpkmat(spikes,'binsize',binsize,'dt',dt,'bintype','gaussian','units','count');
spikemat.isspike = spikemat.data>0.5;
%spikemat.isspike(spikemat.isspike>1) = 1;
for ss = 1:3
    spikemat.instate.(statenames{ss}) = InIntervals(spikemat.timestamps,SleepState.ints.(statenames{ss}));
end

%% For each cell, calculate E and I pop rates of all OTHER cells
ncellthresh = 5;

for tt = 1:length(celltypes)
    Ncells.(celltypes{tt}) = sum(CellClass.(celltypes{tt}));
    spikemat.totpoprate.(celltypes{tt}) = sum(spikemat.data(:,CellClass.(celltypes{tt})),2);
    spikemat.poprate.(celltypes{tt}) = spikemat.totpoprate.(celltypes{tt})./Ncells.(celltypes{tt});
    spikemat.totpopsynch.(celltypes{tt}) = sum(spikemat.isspike(:,CellClass.(celltypes{tt})),2);
    spikemat.popsynch.(celltypes{tt}) = spikemat.totpopsynch.(celltypes{tt})./Ncells.(celltypes{tt});
    
%     if Ncells.(celltypes{tt})==0
%         celltypes(tt) = [];
%     end
end
%%
for tt = 1:length(celltypes)
    synchbins.(celltypes{tt}) = linspace(0,1,Ncells.(celltypes{tt})+1);
    for ss = 1:3

    synchhist.(statenames{ss}).(celltypes{tt})(bb,:) = hist(spikemat.popsynch.(celltypes{tt})(spikemat.instate.(statenames{ss})),synchbins.(celltypes{tt}));
    synchhist.(statenames{ss}).(celltypes{tt})(bb,:) = synchhist.(statenames{ss}).(celltypes{tt})(bb,:)./sum(synchhist.(statenames{ss}).(celltypes{tt})(bb,:));
    end
end
end
%%
figure
for tt = 1:length(celltypes)
    for ss = 1:3
    subplot(2,3,ss+(tt-1)*3)
    %hold on
    
        
        imagesc(log10(bins),synchbins.(celltypes{tt}),synchhist.(statenames{ss}).(celltypes{tt})')
        xlabel('Bin Size (s)')
        ylabel('% Active Neurons')
        LogScale('x',10)
        axis xy
    end
end
NiceSave('SynchTimescale',figfolder,baseName)

%%
%whichbin.pE = 
%whichbin.pI = 7;
figure
for tt = 1:length(celltypes)
    subplot(3,3,tt)
    hold on
    for ss = 1:3
        
        plot(synchbins.(celltypes{tt}),synchhist.(statenames{ss}).(celltypes{tt})(1,:),'color',statecolors{ss})
    end
end
NiceSave('80msbin',figfolder,baseName)
