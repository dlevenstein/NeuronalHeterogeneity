function [ jitterCV2,ISIstats ] = SpikeStatsbyLFPAnalysis(basePath,figfolder)

%% DEV
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/Dino_mPFC/Dino_061814';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyLFPAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
% for reformatting SleepState
%SleepState = SleepScoreMaster(basePath,'noPrompts',true);
SlowWaves = bz_LoadEvents(basePath,'SlowWaves');
%LFP
%%
%Pick channel with most cells close to it...... for now;
usechannel = mode(spikes.maxWaveformCh);

lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath,'noPrompts',true);

ISIStats = bz_LoadCellinfo(basePath,'ISIStats');



%%
downsamplefactor = 5;
 [ lfpdown ] = bz_DownsampleLFP( lfp, downsamplefactor );
[wavespec] = bz_WaveSpec(lfpdown,'showprogress',true,'ncyc',5);

%% General Spike-LFP Coupling

%% Interpolate power at each spike
wavespec.power = log10(abs(wavespec.data));
wavespec.phase = angle(wavespec.data);

%%
[celltypes,~,typeidx] = unique(CellClass.label);
%%
statenames = fieldnames(SleepState.ints);
%% 
%for ss = 1:length(statenames)
state = statenames{2};
    normpower = NormToInt(wavespec.power,'modZ',SleepState.ints.(state),wavespec.samplingRate);
%% Spike-LFP Coupling.

%NOTE: ISSUE WITH WHICH LFP CHANNEL SELECTING!!! UGH.
%     [freqs,synchcoupling,ratepowercorr,...
%     spikephasemag,spikephaseangle,popcellind,cellpopidx] =...
%     GenSpikeLFPCoupling(spikes.times,lfpdown.data,'sf_LFP',lfpdown.samplingRate,'int',SleepState.ints.(state),...
%     'subpop',typeidx,'jittersig',false);
    %%
    instatespiketimes = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
        ISIStats.allspikes.times,'UniformOutput',false);
    ISIStats.spikepower = cellfun(@(X,Y) interp1(wavespec.timestamps,normpower,...
        X,'nearest'),...
        ISIStats.allspikes.times,'UniformOutput',false);
    %%
    %E and I population spikes (all)
    
    for tt = 1:length(celltypes)
        popspikes.(celltypes{tt}).power = cellfun(@(X,Y) X(Y,:),...
            ISIStats.spikepower(CellClass.(celltypes{tt})),...
            instatespiketimes(CellClass.(celltypes{tt})),'UniformOutput',false);
        popspikes.(celltypes{tt}).power = cat(1,popspikes.(celltypes{tt}).power{:});
  %     popspikes.(celltypes{tt}).power = [ISIStats.spikepower{CellClass.(celltypes{tt})}(instatespiketimes{CellClass.(celltypes{tt})})];
  %      popspikes.(celltypes{tt}).CV2 = [ISIStats.allspikes.CV2{CellClass.(celltypes{tt})}];
        popspikes.(celltypes{tt}).CV2 = cellfun(@(X,Y) X(Y,:),...
            ISIStats.allspikes.CV2(CellClass.(celltypes{tt})),...
            instatespiketimes(CellClass.(celltypes{tt})),'UniformOutput',false);
        popspikes.(celltypes{tt}).CV2 = cat(1,popspikes.(celltypes{tt}).CV2{:});
    end
    
    
    %%
    %<CV2> by power....
    %CV2-power correlation in each f
    CV2freqcorr = zeros(spikes.numcells,length(wavespec.freqs));
    CV2freqcorr_pval = zeros(spikes.numcells,length(wavespec.freqs));
    npCV2freqcorr = zeros(spikes.numcells,length(wavespec.freqs));
    npCV2freqcorr_pval = zeros(spikes.numcells,length(wavespec.freqs));
%     for cc = 1:spikes.numcells
%         cc
    for ff = 1:length(wavespec.freqs)
        ff
        [CV2freqcorr(:,ff),CV2freqcorr_pval(:,ff)] = cellfun(@(X,Y,Z) corr(X(Y,ff),Z(Y),...
            'type','spearman'),ISIStats.spikepower,instatespiketimes,ISIStats.allspikes.CV2);
        [npCV2freqcorr(:,ff),npCV2freqcorr_pval(:,ff)] = cellfun(@(X,Y,Z) corr(X(Y,ff),abs(Z(Y)-1),...
            'type','spearman'),ISIStats.spikepower,instatespiketimes,ISIStats.allspikes.CV2);
       
    end
    %%
        for tt = 1:length(celltypes)
           [popCV2freqcorr.(celltypes{tt}) popCV2freqcorr_pval.(celltypes{tt})]= ...
               corr(popspikes.(celltypes{tt}).power,popspikes.(celltypes{tt}).CV2,...
               'type','spearman');
           [nppopCV2freqcorr.(celltypes{tt}) nppopCV2freqcorr_pval.(celltypes{tt})]= ...
               corr(popspikes.(celltypes{tt}).power,abs(1-popspikes.(celltypes{tt}).CV2),...
               'type','spearman');
           
        end
        %%
        cellcolor = {'k','r'};
        figure

    %%
    rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
    figure
    colormap(rwbcolormap)
    subplot(2,2,1)
        h = imagesc(log10(wavespec.freqs),[1 spikes.numcells],CV2freqcorr((ISIStats.sorts.NREMstate.ratebyclass),:));
        %set(h,'AlphaData',(CV2freqcorr_pval<0.05))
        colorbar
        LogScale('x',10)
        caxis([-0.15 0.15])
    subplot(2,2,2)
        h = imagesc(log10(wavespec.freqs),[1 spikes.numcells],npCV2freqcorr((ISIStats.sorts.NREMstate.ratebyclass),:));
        %set(h,'AlphaData',(npCV2freqcorr_pval<0.05))
        colorbar
        LogScale('x',10)
        caxis([-0.15 0.15])
        
        subplot(4,2,5)
        plot(log2(wavespec.freqs([1 end])),[0 0],'k--')
        hold on
        for tt = 1:length(celltypes)
            plot(log2(wavespec.freqs),popCV2freqcorr.(celltypes{tt}),cellcolor{tt})
            plot(log2(wavespec.freqs(popCV2freqcorr_pval.(celltypes{tt})<0.05)),...
                popCV2freqcorr.(celltypes{tt})(popCV2freqcorr_pval.(celltypes{tt})<0.05),...
                '.','color',cellcolor{tt})
            
            LogScale('x',2)
        end
        subplot(4,2,6)
        plot(log2(wavespec.freqs([1 end])),[0 0],'k--')
        hold on
        for tt = 1:length(celltypes)
            plot(log2(wavespec.freqs),nppopCV2freqcorr.(celltypes{tt}),cellcolor{tt})
            plot(log2(wavespec.freqs(nppopCV2freqcorr_pval.(celltypes{tt})<0.05)),...
                nppopCV2freqcorr.(celltypes{tt})(nppopCV2freqcorr_pval.(celltypes{tt})<0.05),...
                '.','color',cellcolor{tt})
            LogScale('x',2)
        end

        NiceSave('LFPandCV2',figfolder,baseName)
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