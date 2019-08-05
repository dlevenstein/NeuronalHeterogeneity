function [ ISIoccupancy,OccupancyStats,normISIhist ] = GroundStateAnalysis( basePath,figfolder )
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
ISIoccupancy.logbins = linspace(-3,2.5,100);
for ss = 1:3
    state = states{ss};
%     ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
%         ISIStats.allspikes.times,'UniformOutput',false);
    ISIrate.instate = InIntervals(ISIrate.timestamps,SleepState.ints.(state));    
    ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,SleepState.ints.(state)),...
        ISIStats.allspikes.times,'UniformOutput',false);    
    
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
    
    normISIs = cellfun(@(X,Y,Z) X(Y)./Z,ISIStats.allspikes.ISIs,ISIStats.allspikes.instate,...
        num2cell(OccupancyStats.(state).median),'UniformOutput',false);
    normISIhist.bins = linspace(-4,1,100);
    normISIhist.(state).mednorm = cellfun(@(X) hist(log10(X),normISIhist.bins),normISIs,'UniformOutput',false);
    normISIhist.(state).mednorm = cellfun(@(X) X./sum(X),normISIhist.(state).mednorm,'UniformOutput',false);
    normISIhist.(state).mednorm = cat(1,normISIhist.(state).mednorm{:});
    
    CV2s = cellfun(@(X,Y,Z) X(Y),ISIStats.allspikes.CV2,ISIStats.allspikes.instate,'UniformOutput',false);
    
    normISIhist.CV2bins = linspace(0,2,50+1);
    normISIhist.CV2bins = normISIhist.CV2bins(1:end-1)+0.5.*diff(normISIhist.CV2bins([1 2]));
    normISIhist.(state).jointCV2 = cellfun(@(X,Y) hist3([log10(X),Y],...
        {normISIhist.bins,normISIhist.CV2bins}),normISIs,CV2s,'UniformOutput',false);
    normISIhist.(state).jointCV2 = cellfun(@(X) X./sum(X(:)),normISIhist.(state).jointCV2,...
        'UniformOutput',false);
    normISIhist.(state).jointCV2 = cat(3,normISIhist.(state).jointCV2{:});
    normISIhist.(state).jointCV2 = shiftdim(normISIhist.(state).jointCV2,2);
 
    for tt = 1:length(celltypes)
        inclasscells{tt} = strcmp(celltypes{tt},CellClass.label);

        %Mean distributions
        meandists.(state).(celltypes{tt}).Jointdist = squeeze(nanmean(normISIhist.(state).jointCV2(inclasscells{tt},:,:),1));
    end
    
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
        (ISIoccupancy.(state).mednormhist(:,ISIStats.sorts.(state).ratebyclass))');
    alpha(s,single(ISIoccupancy.(state).normhist(:,ISIStats.sorts.(state).ratebyclass)'~=0))

    hold on
    plot(log10(1),...
        [1:length(ISIStats.sorts.(state).ratebyclass)],'.')
    LogScale('x',10)
    ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('norm ISI (medOcc)')
    ylabel(state)
% subplot(2,1,2)
%     imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
%         ISIoccupancy.(state).hist(:,ISIStats.sorts.(state).ratebyclass)')
end
NiceSave(['ISIoccupancy'],figfolder,baseName)


%%
histcolors = flipud(gray);
NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);
REMhistcolors = makeColorMap([1 1 1],[0.8 0 0]);
statecolormap = {histcolors,NREMhistcolors,REMhistcolors};


figure
%colormap(cmap)
for ss = 1:3
        state = states{ss};

subplot(3,4,(ss-1)*4+1)
colormap(gca,statecolormap{ss})
    s = imagesc(ISIStats.ISIhist.logbins(1,:),[1 length(OccupancyStats.sorts.(state).median)],...
        (ISIStats.ISIhist.(states{ss}).log(OccupancyStats.sorts.(state).median,:)));
    %alpha(s,single(ISIoccupancy.(state).loghist(:,OccupancyStats.sorts.(state).median)'~=0))

    hold on
    plot(log10(OccupancyStats.(states{ss}).median(OccupancyStats.sorts.(state).median)),...
        [1:length(OccupancyStats.sorts.(state).median)],'k.','markersize',4)
    LogScale('x',10)
    %caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('ISI')

            set(gca,'yticklabel',[])
%     if ss==1
%         title(regions{rr})
%     end
    caxis([0 0.1])
    LogScale('x',10,'exp',true)
% subplot(2,1,2)
%     imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
%         ISIoccupancy.(state).hist(:,ISIStats.sorts.(state).ratebyclass)')

subplot(3,4,(ss-1)*4+2)
colormap(gca,statecolormap{ss})
    s = imagesc(normISIhist.bins,[1 length(OccupancyStats.sorts.(state).median)],...
        (normISIhist.(state).mednorm(OccupancyStats.sorts.(state).median,:)));
    %alpha(s,single(ISIoccupancy.(state).loghist(:,OccupancyStats.sorts.(state).median)'~=0))

    hold on
    plot(0*log10(OccupancyStats.(states{ss}).median(OccupancyStats.sorts.(state).median)),...
        [1:length(OccupancyStats.sorts.(state).median)],'k.','markersize',4)
    LogScale('x',10)
    %caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('norm ISI (medOcc)')

            set(gca,'yticklabel',[])
%     if ss==1
%         title(regions{rr})
%     end
    caxis([0 0.1])
    LogScale('x',10,'exp',true)

for tt = 1:length(celltypes)
subplot(3,4,(ss-1)*4+2+tt)
imagesc(normISIhist.bins,normISIhist.CV2bins,meandists.(state).(celltypes{tt}).Jointdist')
axis xy
end

end
NiceSave(['ISIdists'],figfolder,baseName)
