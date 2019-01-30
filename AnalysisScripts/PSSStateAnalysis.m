function [ output_args ] = PSSStateAnalysis( basePath,figfolder )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% DEV
basePath = '/mnt/NyuShare/Buzsakilabspace/Datasets/GrosmarkAD/Achilles/Achilles_11012013';
repoRoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity'; %desktop

basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
repoRoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %laptop

figfolder = [repoRoot,'/AnalysisScripts/AnalysisFigs/PSSStateAnalysis'];
%%
baseName = bz_BasenameFromBasepath(basePath);
SleepState = bz_LoadStates(basePath,'SleepState');
sessionInfo = bz_getSessionInfo(basePath);
%% Calculate PSS and its bimodality for different bounds (SlowWave channel)
downsamplefactor = 5;

regions = {'CTX','HPC'};
for rr = 1:length(regions)
    regchans{rr}  = sessionInfo.channels(ismember(sessionInfo.region,regions{rr}));
end
usechans = [regchans{:}];
for rr = 1:length(regions)
    regchanIDX{rr} = ismember(usechans,regchans{rr});
end
regcolor = {[0 0 0.4],[0.5 0 0.2]};
numregchans = cellfun(@length,regchans);
repchan = [SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    randsample(regchans{2},1)];


%%
allLFP = bz_GetLFP(usechans,...
    'basepath',basePath,'downsample',downsamplefactor);

%% Calculate PSS correlation between channels
bounds = [3 120];
winsize = 10; %Why not use what your'e usingv for sleep scoring?....
dt = 5;
specslope = bz_PowerSpectrumSlope(allLFP,winsize,dt,'frange',bounds);


%%
PSScorr = corr(specslope.data,'type','spearman');

%Pick a HPC example channel with best correlation to CTX
%Best HPC chan
repchan(2) = max(PSScorr(ismember(specslope.channels,repchan(rr)),regchanIDX{2}));
repchan(2) = regchans{2}(repchan(2));
%%
spikegroupordering = [sessionInfo.spikeGroups.groups{:}];
[spikegroupordering,~,spikegroupsort] = intersect(spikegroupordering,usechans,'stable');


%% Compare Time Scales

numhistbins = 20;
histbins = linspace(0,1,numhistbins);

dt = 1;
winsizes = logspace(-0.33,1.7,15);
for ww = 1:length(winsizes)
    specslope_win = bz_PowerSpectrumSlope(allLFP,winsizes(ww),dt,...
        'frange',bounds,'channels',repchan);
    PSScorr_win(ww) = corr(specslope_win.data(:,1),specslope_win.data(:,2),'type','spearman');
    
        % Histogram and diptest of Slow Wave Power
        for rr = 1:length(regions)
            broadbandSlowWave = specslope_win.data(:,rr);
    %         broadbandSlowWave = smooth(broadbandSlowWave,smoothfact.*specslope_freqs.samplingRate);
            broadbandSlowWave = ...
                (broadbandSlowWave-min(broadbandSlowWave))./...
                max(broadbandSlowWave-min(broadbandSlowWave));
            
            [swhist.temp]= hist(broadbandSlowWave,histbins);

            swhists_win.(regions{rr})(:,ww) = swhist.temp;
            dipSW_win.(regions{rr})(ww) = bz_hartigansdiptest(broadbandSlowWave);
        end
end


%% Compare Frequency Bands
window = 10;
dt = 5; %Updated to speed up (don't need to sample at fine time resolution for channel selection)
%smoothfact = 10;

nbounds = 15;
upperbounds = logspace(1.66,2.33,nbounds);
lowerbounds = logspace(0,1.33,nbounds);


clear swhists
clear dipSW
for uu = 1:nbounds
    uu
    for ll = 1:nbounds
        [specslope_freqs] = bz_PowerSpectrumSlope(allLFP,window,dt,...
            'frange',[lowerbounds(ll) upperbounds(uu)],'channels',repchan);
        
        

        
        % Histogram and diptest of Slow Wave Power
        for rr = 1:length(regions)
            broadbandSlowWave = specslope_freqs.data(:,rr);
    %         broadbandSlowWave = smooth(broadbandSlowWave,smoothfact.*specslope_freqs.samplingRate);
            broadbandSlowWave = ...
                (broadbandSlowWave-min(broadbandSlowWave))./...
                max(broadbandSlowWave-min(broadbandSlowWave));
            
            [swhist.temp]= hist(broadbandSlowWave,histbins);

            swhists.(regions{rr})(:,ll,uu) = swhist.temp;
            dipSW.(regions{rr})(ll,uu) = bz_hartigansdiptest(broadbandSlowWave);
        end
        
        regioncorr(ll,uu) = corr(specslope_freqs.data(:,1),specslope_freqs.data(:,2),'type','spearman');

    end
end

%%
[~,bestboundsIdx] = max(dipSW(:));
[bestboundsIdx(1),bestboundsIdx(2)] = ind2sub(size(dipSW),bestboundsIdx);

bestbounds(2) = upperbounds(bestboundsIdx(2));
bestbounds(1) = lowerbounds(bestboundsIdx(1));

examples = [4,5];
%%

%% correlation with EMG
EMGfromLFP = bz_EMGFromLFP(basePath);
%%
specslope.EMG = interp1(EMGfromLFP.timestamps,EMGfromLFP.data,specslope.timestamps);
PSSEMGcorr = corr(specslope.data,specslope.EMG,'type','spearman');

%%
hist(PSSEMGcorr)

%% HPC-CTX heatmap
[multiregionhist.counts,multiregionhist.bins] = hist3([specslope.data(:,ismember(specslope.channels,repchan(1))),...
        specslope.data(:,ismember(specslope.channels,repchan(2)))],[50 50]);


%%

%add example windows
figure
subplot(5,4,13)
imagesc(PSScorr(spikegroupsort,spikegroupsort))
colorbar
ColorbarWithAxis([0 1],'Corr(PSS)')
set(gca,'xtick',cumsum(numregchans)-0.5.*numregchans)
set(gca,'xticklabel',regions)
set(gca,'ytick',cumsum(numregchans)-0.5.*numregchans)
set(gca,'yticklabel',regions)

subplot(5,4,17)
    imagesc(log2(lowerbounds),log2(upperbounds),regioncorr')
    hold on
    plot(log2(bounds(1)),log2(bounds(2)),'r+')

    %plot(log2(bestbounds(1)),log2(bestbounds(2)),'r+')
    %plot(log2(lowerbounds(examples(1))),log2(upperbounds(examples(2))),'k+')
    axis xy
    LogScale('xy',2)
    xlabel('L-Bound (Hz)');ylabel('U-Bound (Hz)')
    ColorbarWithAxis([0.4 0.9],'HPC-CTX Corr')

subplot(5,4,20)
%hist(PSSEMGcorr)
BoxAndScatterPlot({PSSEMGcorr(regchanIDX{1}),PSSEMGcorr(regchanIDX{2})},...
    'labels',regions,'colors',cat(1,regcolor{:}))
ylim([0 1])
ylabel('EMG-PSS corr')

subplot(5,4,16)
plot(log10(winsizes),PSScorr_win)
hold on
    plot(log10(winsize),0,'r+')

xlabel('Window Duration (s)');ylabel({'HPC-CTX', 'PSS corrrelation'})
xlim(log10(winsizes([1 end])))
ylim([0 1])
LogScale('x',10)

subplot(5,4,14)
    plot(specslope.data(:,ismember(specslope.channels,repchan(1))),...
        specslope.data(:,ismember(specslope.channels,repchan(2))),...
        '.')
    xlabel([regions{1},' PSS']);ylabel([regions{2},' PSS'])
    
subplot(5,4,18)
imagesc(multiregionhist.bins{1},multiregionhist.bins{2}, multiregionhist.counts')
axis xy
    xlabel([regions{1},' PSS']);ylabel([regions{2},' PSS'])

    
exwin = bz_RandomWindowInIntervals(specslope.timestamps([1 end]),2000);

for rr = 1:length(regions)
    subplot(6,1,rr)
        imagesc(specslope.timestamps,log2(specslope.freqs),...
            specslope.specgram(:,:,ismember(specslope.channels,repchan(rr))))
        hold on
        axis xy

        LogScale('y',2)
        ylabel({regions{rr},'f (Hz)'})

        yyaxis right
        plot(specslope.timestamps,...
            specslope.data(:,ismember(specslope.channels,repchan(rr))),...
            'color',regcolor{rr},'linewidth',1)
        ylabel('PSS')
        xlim(exwin)
          bz_ScaleBar('s')

end
    

subplot(6,1,3)
    plot(EMGfromLFP.timestamps,EMGfromLFP.data)
        xlim(exwin)
        box off
      bz_ScaleBar('s')
      ylabel('EMG')
      
NiceSave('PSSasGlobalState',figfolder,baseName)



%% Relate to firing rates, firing rate distribution width!
%% Calculate PSS and its bimodality for different bounds (SlowWave channel)


%%
figure
for rr = 1:length(regions)
subplot(5,4,4+rr)
    imagesc(log2(lowerbounds),log2(upperbounds),dipSW.(regions{rr})')
    hold on
    plot(log2(bounds(1)),log2(bounds(2)),'r+')
    %plot(log2(bestbounds(1)),log2(bestbounds(2)),'r+')
    %plot(log2(lowerbounds(examples(1))),log2(upperbounds(examples(2))),'k+')
    axis xy
    LogScale('xy',2)
    xlabel('Lower Bound (Hz)');ylabel('Upper Bound (Hz)')
    ColorbarWithAxis([0 max(dipSW.(regions{rr})(:))],'Dip')
    title(regions{rr})
end

subplot(5,4,2)
hold on
for rr = 1:length(regions)
    plot(log10(winsizes),dipSW_win.(regions{rr}),'color',regcolor{rr})
end
    plot(log10(winsize),0,'r+')

xlabel('Window Duration (s)');ylabel('Dip Test')
xlim(log10(winsizes([1 end])))
%ylim([0 1])
LogScale('x',10)

 
% 
% subplot(4,2,2)
%     plot(histbins,swhists(:,bestboundsIdx(1),bestboundsIdx(2)),'r')
%     hold on
%     plot(histbins,swhists(:,examples(1),examples(2)),'k')
%     
%histogram: good, OK, bad dips

NiceSave('PSSBimodality',figfolder,baseName)



%% PSS and UP/DOWN

SlowWaves = bz_LoadEvents(basePath,'SlowWaves');

updown = {'DOWN','UP'};
UDcolor = {'b','r'};
for ss = 1:2
    SlowWaves.dur.(updown{ss}) = diff(SlowWaves.ints.(updown{ss}),1,2);
    SlowWaves.midpoint.(updown{ss}) = mean(SlowWaves.ints.(updown{ss}),2);
    SlowWaves.PSS.(updown{ss}) = interp1(specslope.timestamps,specslope.data(:,ismember(specslope.channels,repchan(1))),SlowWaves.midpoint.(updown{ss}));
end

%%
figure
subplot(3,2,1)
for ss = 1:2
    plot(SlowWaves.PSS.(updown{ss}),log10(SlowWaves.dur.(updown{ss})),'.','color',UDcolor{ss},'markersize',3)
    hold on
end
xlabel('PSS');ylabel('Dur (s)')
axis tight
box off
LogScale('y',10)
legend(updown{:},'location','eastoutside')
NiceSave('PSSandUPDOWN',figfolder,baseName)
end

