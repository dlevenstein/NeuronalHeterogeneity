function [ output_args ] = PSSRangeBimodality( basePath,figfolder )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% DEV
basePath = '/mnt/NyuShare/Buzsakilabspace/Datasets/GrosmarkAD/Achilles/Achilles_11012013';
repoRoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity'; %desktop
figfolder = [repoRoot,'/AnalysisScripts/AnalysisFigs/PSSRangeBimodality'];
%%
baseName = bz_BasenameFromBasepath(basePath);
SleepState = bz_LoadStates(basePath,'SleepState');
%% Calculate PSS and its bimodality for different bounds
downsamplefactor = 10;
allLFP = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath,'downsample',downsamplefactor);

%%
%Parameters from SleepScoreMaster (optimize these too?)
window = 10;
noverlap = 5; %Updated to speed up (don't need to sample at fine time resolution for channel selection)
smoothfact = 10;

nbounds = 10;
upperbounds = logspace(1.66,2.33,nbounds);
lowerbounds = logspace(0,1.33,nbounds);

numhistbins = 21;
histbins = linspace(0,1,numhistbins);

for uu = 1:nbounds
    uu
    for ll = 1:nbounds
        [specslope,~] = bz_PowerSpectrumSlope(allLFP,window,window-noverlap,...
            'frange',[lowerbounds(ll) upperbounds(uu)]);
        
        
        broadbandSlowWave = specslope.data;
        broadbandSlowWave = smooth(broadbandSlowWave,smoothfact.*specslope.samplingRate);
        broadbandSlowWave = ...
            (broadbandSlowWave-min(broadbandSlowWave))./...
            max(broadbandSlowWave-min(broadbandSlowWave));
        
        % Histogram and diptest of Slow Wave Power
        [swhist]= hist(broadbandSlowWave,histbins);

        swhists(:,uu,ll) = swhist;
        dipSW(uu,ll) = bz_hartigansdiptest((broadbandSlowWave));
        
    end
end

%%
%upper and lower are backwards. silly
[~,bestboundsIdx] = max(dipSW(:));
[bestboundsIdx(1),bestboundsIdx(2)] = ind2sub(size(dipSW),bestboundsIdx);

bestbounds(1) = upperbounds(bestboundsIdx(1));
bestbounds(2) = lowerbounds(bestboundsIdx(2));

examples = [4,5];
%%
figure
subplot(2,2,1)
    imagesc(log2(lowerbounds),log2(upperbounds),dipSW)
    hold on
    plot(log2(bestbounds(2)),log2(bestbounds(1)),'r+')
    plot(log2(lowerbounds(examples(2))),log2(upperbounds(examples(1))),'k+')
    axis xy
    LogScale('xy',2)
    xlabel('Lower Bound (Hz)');ylabel('Upper Bound (Hz)')
    colorbar
subplot(4,2,2)
    plot(histbins,swhists(:,bestboundsIdx(1),bestboundsIdx(2)),'r')
    hold on
    plot(histbins,swhists(:,examples(1),examples(2)),'k')
    
%histogram: good, OK, bad dips

NiceSave('PSSBounds',figfolder,baseName)
end

