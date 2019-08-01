function [ ISIoccupancy,OccupancyStats ] = GroundStateAnalysis( basePath,figfolder )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = [reporoot,'Datasets/onDesktop/AG_HPC/Achilles_10252013'];
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};


%% Load the LFP if needed

% lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID;
% downsamplefactor = 2;
% lfp = bz_GetLFP(lfpchan,...
%     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%Noralize the LFP
%lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

%% Interpolate containing ISI at each ms

ISIrate.dt = 0.002;
ISIrate.timestamps = [0:ISIrate.dt:max(cat(1,spikes.times{:}))]';

%Bug Fix two spikes same time
[~,ISIStats.allspikes.unique] = cellfun(@(X) unique(X),ISIStats.allspikes.times,'UniformOutput',false);

ISIrate.ISI = cellfun(@(X,Y,Z) interp1(X(Z),Y(Z),ISIrate.timestamps,'next'),...
    ISIStats.allspikes.times,ISIStats.allspikes.ISIs,ISIStats.allspikes.unique,'UniformOutput',false);
ISIrate.ISI = cat(2,ISIrate.ISI{:});


%% ISI occupancy
% 
ISIoccupancy.bins = linspace(0,20,100);
ISIoccupancy.logbins = linspace(-2.5,2.5,100);
for ss = 1:3
    state = states{ss};
%     ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
%         ISIStats.allspikes.times,'UniformOutput',false);
    ISIrate.instate = InIntervals(ISIrate.timestamps,SleepState.ints.(state));    
    
    if sum(ISIrate.instate)==0
        ISIoccupancy.(state).hist = nan(length(ISIoccupancy.bins),spikes.numcells);
        ISIoccupancy.(state).loghist = nan(length(ISIoccupancy.logbins),spikes.numcells);
        ISIoccupancy.(state).normhist = nan(length(ISIoccupancy.logbins),spikes.numcells);
       continue 
    end
    %% Calculate occupancy histogram


    ISIoccupancy.(state).hist = hist(ISIrate.ISI(ISIrate.instate,:),ISIoccupancy.bins);
    ISIoccupancy.(state).hist = ISIoccupancy.(state).hist./length(ISIrate.timestamps(ISIrate.instate));
    %ISIoccupancy.(state).hist(ISIoccupancy.(state).hist==0)=nan;

    ISIoccupancy.(state).loghist = hist(log10(ISIrate.ISI(ISIrate.instate,:)),ISIoccupancy.logbins);
    ISIoccupancy.(state).loghist = ISIoccupancy.(state).loghist./length(ISIrate.timestamps(ISIrate.instate));
    %ISIoccupancy.(state).loghist(ISIoccupancy.(state).loghist==0)=nan;
    
    ISIoccupancy.(state).normhist = hist(log10(ISIrate.ISI(ISIrate.instate,:)./ISIStats.summstats.(state).meanISI),...
        ISIoccupancy.logbins);
    ISIoccupancy.(state).normhist = ISIoccupancy.(state).normhist./length(ISIrate.timestamps(ISIrate.instate));
    
    %% Calculate occupancy statistics
    OccupancyStats.(state).median = nanmedian(ISIrate.ISI(ISIrate.instate,:));
    OccupancyStats.(state).std = std(ISIrate.ISI(ISIrate.instate,:));
    OccupancyStats.(state).mean = nanmean(log10(ISIrate.ISI(ISIrate.instate,:)));
    %OccupancyStats.(state).stdlog = std(log10(ISIrate.ISI(ISIrate.instate,:)));
    
    [~,OccupancyStats.sorts.(state).mean] = sort(OccupancyStats.(state).mean);
    [~,OccupancyStats.sorts.(state).median] = sort(OccupancyStats.(state).median);
    
    %% MedianOccupancy Normalization
    ISIoccupancy.(state).mednormhist = hist(log10(ISIrate.ISI(ISIrate.instate,:)./OccupancyStats.(state).median),...
        ISIoccupancy.logbins);
    ISIoccupancy.(state).mednormhist = ISIoccupancy.(state).mednormhist./length(ISIrate.timestamps(ISIrate.instate));
end

%%
% figure
% subplot(2,2,1)
% plot(OccupancyStats.NREMstate.mean,OccupancyStats.WAKEstate.mean,'.')
% 
% subplot(2,2,2)
% plot(OccupancyStats.NREMstate.mean,OccupancyStats.NREMstate.std,'.')
%%
%cmap = [1 1 1;colormap(parula)];

figure
%colormap(cmap)
for ss = 1:3
        state = states{ss};

subplot(3,2,ss*2-1)
    s = imagesc(ISIoccupancy.logbins,[1 spikes.numcells],...
        (ISIoccupancy.(state).loghist(:,ISIStats.sorts.(state).ratebyclass))');
    alpha(s,single(ISIoccupancy.(state).loghist(:,ISIStats.sorts.(state).ratebyclass)'~=0))

    hold on
    plot(log10(1./ISIStats.summstats.(state).meanrate(ISIStats.sorts.(state).ratebyclass)),...
        [1:length(ISIStats.sorts.(state).ratebyclass)],'.')
    LogScale('x',10)
    ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('ISI')
    ylabel(state)
% subplot(2,1,2)
%     imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
%         ISIoccupancy.(state).hist(:,ISIStats.sorts.(state).ratebyclass)')


subplot(3,2,ss*2)
    s = imagesc(ISIoccupancy.logbins,[1 spikes.numcells],...
        (ISIoccupancy.(state).normhist(:,ISIStats.sorts.(state).ratebyclass))');
    alpha(s,single(ISIoccupancy.(state).normhist(:,ISIStats.sorts.(state).ratebyclass)'~=0))

    hold on
    plot(log10(1),...
        [1:length(ISIStats.sorts.(state).ratebyclass)],'.')
    LogScale('x',10)
    ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('ISI')
    ylabel(state)
% subplot(2,1,2)
%     imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
%         ISIoccupancy.(state).hist(:,ISIStats.sorts.(state).ratebyclass)')
end
NiceSave(['ISIoccupancy'],figfolder,baseName)
