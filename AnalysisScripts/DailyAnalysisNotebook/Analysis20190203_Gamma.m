function [ ] = Analysis20190203(basePath,figfolder)
% Date 02/03/2019
%
%Question: Does gamma activity modulate CV2?
%
%Plots
%-CV2-power correlation with wavelets
%-
%
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/Dino_mPFC/Dino_061814';
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
downsamplefactor = 1;
lfp = bz_GetLFP(lfpchan,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%Noralize the LFP
lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

%% Restrict to state
state = states{2}; %REM... least time
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);

%% Get the wavelets (use ncyc=12, per 20190202
[wavespec] = bz_WaveSpec(lfp,'showprogress',true,'ncyc',15,'frange',[20 200],...
    'intervals',SleepState.ints.(state));
%wavespec.power = log10(abs(wavespec.data));
%wavespec.phase = angle(wavespec.data);

%% Interpolate power at each spike
ISIStats.allspikes.spikepower = cellfun(@(X,Y) interp1(wavespec.timestamps,log10(abs(wavespec.data)),...
    X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);


%%
CV2freqcorr = zeros(spikes.numcells,length(wavespec.freqs));
CV2freqcorr_pval = zeros(spikes.numcells,length(wavespec.freqs));
npCV2freqcorr = zeros(spikes.numcells,length(wavespec.freqs));
npCV2freqcorr_pval = zeros(spikes.numcells,length(wavespec.freqs));

for ff = 1:length(wavespec.freqs)
    ff
    [CV2freqcorr(:,ff),CV2freqcorr_pval(:,ff)] = cellfun(@(X,Y,Z) corr(X(Y,ff),Z(Y),...
        'type','spearman'),ISIStats.allspikes.spikepower,ISIStats.allspikes.instate,ISIStats.allspikes.CV2);
    [npCV2freqcorr(:,ff),npCV2freqcorr_pval(:,ff)] = cellfun(@(X,Y,Z) corr(X(Y,ff),abs(Z(Y)-1),...
        'type','spearman'),ISIStats.allspikes.spikepower,ISIStats.allspikes.instate,ISIStats.allspikes.CV2);

end

%%
    rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
    figure
    colormap(rwbcolormap)
    subplot(2,2,1)
        h = imagesc(log2(wavespec.freqs),[1 spikes.numcells],CV2freqcorr((ISIStats.sorts.NREMstate.ratebyclass),:));
        %set(h,'AlphaData',(CV2freqcorr_pval<0.05))
        colorbar('northoutside')
        LogScale('x',2)
        caxis([-0.15 0.15])
        ylabel('Cell - sorted by rate')
    subplot(2,2,2)
        h = imagesc(log2(wavespec.freqs),[1 spikes.numcells],npCV2freqcorr((ISIStats.sorts.NREMstate.ratebyclass),:));
        %set(h,'AlphaData',(npCV2freqcorr_pval<0.05))
        colorbar('northoutside')
        LogScale('x',2)
        caxis([-0.15 0.15])
        
