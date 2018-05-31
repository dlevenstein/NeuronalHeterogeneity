%% DEV
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
figfolder = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs';
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
% for reformatting SleepState
%SleepState = SleepScoreMaster(basePath,'noPrompts',true);

%LFP
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath,'noPrompts',true);

ISIstats = bz_LoadCellinfo(basePath,'ISIstats');

%%
downsamplefactor = 5;
 [ lfpdown ] = bz_DownsampleLFP( lfp, downsamplefactor );
[wavespec] = bz_WaveSpec(lfpdown,'showprogress',true);

%% Interpolate power at each spike
wavespec.power = log10(abs(wavespec.data));

%% 
ISIstats.spikepower = cellfun(@(X) interp(wavespec.timestamps,wavespec.power,X,'nearest'),...
    ISIstats.allspikes.times,'UniformOutput',true)


