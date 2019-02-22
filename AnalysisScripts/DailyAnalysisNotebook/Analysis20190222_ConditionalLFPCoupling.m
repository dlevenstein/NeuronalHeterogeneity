function [ ] = Analysis20190222(basePath,figfolder)
% Date 02/20/2019
%
%Goal: Develop function for conditional LFP coupling (i.e. power and phase
%coupling as function of some condition (ISI, for example)
%Mean normalize after the fact...
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC/Achilles_10252013'
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

lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID;
downsamplefactor = 2;
lfp = bz_GetLFP(lfpchan,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%Noralize the LFP
%lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

%% Restrict to state

state = states{1};
%ints = SleepState.ints.(state);

%Take only subset of time (random intervals) so wavelets doesn't break
%computer (total 625s)
usetime = 2000;%2500
winsize = 25;
if sum(diff(SleepState.ints.(state),1,2))>usetime
    nwin = round(usetime./winsize);
    %winsize = 30; %s
    windows = bz_RandomWindowInIntervals( SleepState.ints.(state),winsize,nwin );
else
    windows = SleepState.ints.(state);
end

%%
% figure
% plot(lfp.timestamps,lfp.data,'k')
% hold on
% plot(windows',ones(size(windows')),'r','linewidth',4)
%% Get complex valued wavelet transform at each timestamp
wavespec = bz_WaveSpec(lfp,'intervals',windows,'showprogress',true,'ncyc',15,...
    'nfreqs',100,'frange',[1 312]);  %150 freqs, 1 312
%Mean-Normalize power within the interval
wavespec.meanpower = mean(abs(wavespec.data),1);

%%
ISIStats.allspikes.logISIs = cellfun(@(X) log10(X),ISIStats.allspikes.ISIs,'UniformOutput',false);
ISIStats.allspikes.logISIs_next = cellfun(@(X) log10([X(2:end);nan]),ISIStats.allspikes.ISIs,'UniformOutput',false);

doubleISIs.times = cellfun(@(X) [X;X],ISIStats.allspikes.times,'UniformOutput',false);
doubleISIs.ISIs = cellfun(@(X,Y) [X;Y],ISIStats.allspikes.logISIs,ISIStats.allspikes.logISIs_next,'UniformOutput',false);
%%
bz_ConditionalLFPCoupling( doubleISIs,doubleISIs.ISIs,wavespec,...
    'Xbounds',[-2.6 1],'intervals',windows,'showFig',false,... %true
'minX',25,'CellClass',CellClass,...
'saveFig',figfolder,'figName',['ISIConditionedLFP',state],'baseName',baseName);

%% next: coupling conditioned on CV2
%Add Power-ISI mutual information

%Get complex-valued filtered LFP at each spike time
for cc = 1:spikes.numcells
    cc
    ISIStats.allspikes.LFP{cc} = interp1(wavespec.timestamps,wavespec.data,ISIStats.allspikes.times{cc},'nearest');
end

%%
for cc = 1:spikes.numcells
    ISIStats.allspikes.LFP{cc} = bsxfun(@(X,Y) X./Y,ISIStats.allspikes.LFP{cc},wavespec.meanpower);
end
%%
ISIbins = linspace(-2.5,1,50);
powerbins = linspace(0,2,40);
excell=randsample(spikes.numcells,1);

joint = hist3([ISIStats.allspikes.logISIs{excell} abs(ISIStats.allspikes.LFP{excell})],{ISIbins,powerbins});


