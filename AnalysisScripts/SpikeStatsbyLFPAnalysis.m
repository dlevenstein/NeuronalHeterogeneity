function [ jitterCV2,ISIstats ] = SpikeStatsbyLFPAnalysis(basePath,figfolder)

%% DEV
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyLFPAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
% for reformatting SleepState
%SleepState = SleepScoreMaster(basePath,'noPrompts',true);

%LFP
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath,'noPrompts',true);

ISIStats = bz_LoadCellinfo(basePath,'ISIStats');

%%
downsamplefactor = 5;
 [ lfpdown ] = bz_DownsampleLFP( lfp, downsamplefactor );
[wavespec] = bz_WaveSpec(lfpdown,'showprogress',true,'ncyc',4);

%% Interpolate power at each spike
wavespec.power = log10(abs(wavespec.data));


%%
statenames = fieldnames(SleepState.ints);
%% 
%for ss = 1:length(statenames)
state = statenames{2};
    normpower = NormToInt(wavespec.power,'modZ',SleepState.ints.(state),wavespec.samplingRate);

    %%
    instatespiketimes = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
        ISIStats.allspikes.times,'UniformOutput',false);
    ISIStats.spikepower = cellfun(@(X,Y) interp1(wavespec.timestamps,normpower,...
        X,'nearest'),...
        ISIStats.allspikes.times,'UniformOutput',false);
    
    %%
    %<CV2> by power....
    %CV2-power correlation in each f
    CV2freqcorr = zeros(spikes.numcells,length(wavespec.freqs));
    for cc = 1:spikes.numcells
        cc
        for ff = 1:length(wavespec.freqs)
        CV2freqcorr(cc,ff) = corr(ISIStats.spikepower{cc}(instatespiketimes{cc},ff),...
            ISIStats.allspikes.CV2{cc}(instatespiketimes{cc}),...
            'type','spearman');
        end
    end
    %%
    figure
    imagesc(log10(wavespec.freqs),[1 spikes.numcells],CV2freqcorr)
    colorbar
    LogScale('x',10)
%%
cellnum = 9;
freqidx = 30;
figure
subplot(2,2,1)
    plot(ISIStats.spikepower{cellnum}(instatespiketimes{cellnum},freqidx),...
        ISIStats.allspikes.CV2{cellnum}(instatespiketimes{cellnum}),'.')
    title(num2str(CV2freqcorr(cellnum,freqidx)))
%subplot(2,2,2)
%    plot(ISIStats.allspikes.CV2{cellnum}(instatespiketimes{cellnum})


%%

end