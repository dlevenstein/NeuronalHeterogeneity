function [ ] = Analysis20190215(basePath,figfolder)
% Date 02/15/2019
%
%
%Question: Where do neurons spend their most time? What are the statistics
%of the ground state and how do we show it's poisson-like 
%
%Statement to support and expand on:
%"We find that forebrain neurons spend the majority of time in a low rate
%mode with poisson-like spiking statistics. This "ground state" determines
%the mean firing rate of the neuron and shows low amplitude fluctuations in
%rate the follow low-dimensional population-wide dynamics.
%
%Plots
%-ISI time occupancy
%-Conditional power distribution on ISI?
%-CV2 of low-rate state?... mean CV2 (or CV2 dist) as function of longest
%interval?
%
%% Load Header
%Initiate Paths
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
basePath = [reporoot,'Datasets/onDesktop/AG_HPC/Achilles_10252013'];
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
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
ISIrate.ISI = cellfun(@(X,Y) interp1(X,Y,ISIrate.timestamps,'next'),...
    ISIStats.allspikes.times,ISIStats.allspikes.ISIs,'UniformOutput',false);
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

    %% Calculate occupancy histogram


    ISIoccupancy.(state).hist = hist(ISIrate.ISI(ISIrate.instate,:),ISIoccupancy.bins);
    ISIoccupancy.(state).hist = ISIoccupancy.(state).hist./length(ISIrate.timestamps(ISIrate.instate));
    %ISIoccupancy.(state).hist(ISIoccupancy.(state).hist==0)=nan;

    ISIoccupancy.(state).loghist = hist(log10(ISIrate.ISI(ISIrate.instate,:)),ISIoccupancy.logbins);
    ISIoccupancy.(state).loghist = ISIoccupancy.(state).loghist./length(ISIrate.timestamps(ISIrate.instate));
    %ISIoccupancy.(state).loghist(ISIoccupancy.(state).loghist==0)=nan;
end
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
        [1:spikes.numcells],'.')
    LogScale('x',10)
    ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('ISI')
    ylabel(state)
% subplot(2,1,2)
%     imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
%         ISIoccupancy.(state).hist(:,ISIStats.sorts.(state).ratebyclass)')
end
NiceSave(['ISIoccupancy'],figfolder,baseName,'includeDate',true)
%%
excell = randsample(spikes.numcells,1);
viewwin = bz_RandomWindowInIntervals(SleepState.ints.(state),30);

figure
    subplot(4,1,1)
        plot(ISIrate.timestamps,log10(ISIrate.ISI{excell}),'k')
        axis tight
        hold on
        StateScorePlot(SleepState.ints,statecolors)
        box off 
        LogScale('y',10)
        %hold on
        %plot(ISIStats.allspikes.times{excell},zeros(size(ISIStats.allspikes.times{excell})),'+')
        %xlim(viewwin)
    subplot(4,1,2)
        plot(ISIrate.timestamps,log10(1./ISIrate.ISI{excell}),'k')
        %hold on
        axis tight
        box off
        LogScale('y',10)
        %plot(ISIStats.allspikes.times{excell},zeros(size(ISIStats.allspikes.times{excell})),'+')
        %xlim(viewwin)
    subplot(2,2,3)
        plot(ISIrate.timestamps,ISIrate.ISI{excell},'k')
        hold on
        plot(ISIStats.allspikes.times{excell},zeros(size(ISIStats.allspikes.times{excell})),'+')
        xlim(viewwin)
    subplot(4,2,6)
        bar(ISIoccupancy.bins,ISIoccupancy.(state).hist(excell,:))
        axis tight
        box off
    subplot(4,2,8)
        bar(ISIoccupancy.logbins,ISIoccupancy.(state).loghist(excell,:))
        axis tight
        box off
        LogScale('x',10)
