function [ ] = AnalysisXXXXXXXX(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC/Achilles_10252013';
basePath = '/mnt/proraidDL/Database/YSData/YutaData/YMV08_170922';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%%
load(fullfile(basePath,[baseName,'.SleepScoreLFP.LFP.mat']))
load(fullfile(basePath,[baseName,'.EMGFromLFP.LFP.mat']))
%%
sessionInfo = bz_getSessionInfo(basePath);
%%
lfp = bz_GetLFP(SleepScoreLFP.SWchanID,...
     'basepath',basePath,'noPrompts',true);
%% Load the SleepScoreLFP and EMG



%%
winsize = 10;
dt = 5;
[specslope_irasa,spec_irasa] = bz_PowerSpectrumSlope(lfp,winsize,dt,'IRASA',true,'showfig',true)

[specslope,spec] = bz_PowerSpectrumSlope(lfp,winsize,dt,'IRASA',false,'showfig',true)

%%
figure
plot(specslope_irasa.data,specslope.data,'.')

%% Compare IRASA with
lfp_down = bz_DownsampleLFP(lfp,5);
[specslope_down,spec_down] = bz_PowerSpectrumSlope(lfp_down,2,1,'IRASA',true,'showfig',true)

[specslope_high,spec_high] = bz_PowerSpectrumSlope(lfp,2,1,'IRASA',true,'showfig',true)

%%
figure
plot(specslope_high.data,specslope_down.data,'.')
%% Compare Bimodality IRASA vs no

%%
wins = 1:12;
swins = 1:25;


%% Get the metrics (PSS, theta)
for ww = 1:length(wins)
    bz_Counter(ww,length(wins),'Window')
    parfor ss = 1:length(swins)
        
[SleepScoreMetrics_IRASA(ww,ss),StatePlotMaterials_IRASA(ww,ss)] = ClusterStates_GetMetrics(...
                                           basePath,SleepScoreLFP,EMGFromLFP,true,...
                                           'window',wins(ww),'smoothfact',swins(ss),...
                                           'IRASA',true);
                                       
[SleepScoreMetrics(ww,ss),StatePlotMaterials(ww,ss)] = ClusterStates_GetMetrics(...
                                           basePath,SleepScoreLFP,EMGFromLFP,true,...
                                           'window',wins(ww),'smoothfact',swins(ss));
    end
end
       
%TRY MEDIAN SMOOTHING!
%% Test for bimodality
for ww = 1:length(wins)
    for ss = 1:length(swins)
        try
        dipSW(ww,ss) = bz_hartigansdiptest(sort(SleepScoreMetrics(ww,ss).broadbandSlowWave));
        dipTH(ww,ss) = bz_hartigansdiptest(sort(SleepScoreMetrics(ww,ss).thratio));
        catch
            dipSW(ww,ss) = nan;
            dipTH(ww,ss) = nan;
            
        end
        
        try
        dipSW_IRASA(ww,ss) = bz_hartigansdiptest(sort(SleepScoreMetrics_IRASA(ww,ss).broadbandSlowWave));
        dipTH_IRASA(ww,ss) = bz_hartigansdiptest(sort(SleepScoreMetrics_IRASA(ww,ss).thratio));
        catch
            dipSW_IRASA(ww,ss) = nan;
            dipTH_IRASA(ww,ss) = nan;
            
        end
    end
end

%%

examples = [10 10; 2 15];

figure

subplot(3,3,1)
    imagesc(wins,swins,dipSW')
    hold all
    for ee = 1:size(examples,1)
        plot(wins(examples(ee,1)),swins(examples(ee,2)),'+')
    end
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.065])
    title('Bimodality: Linear Fit')


    subplot(3,3,3)
    hold all
    for ee = 1:size(examples,1)
        plot(SleepScoreMetrics(examples(ee,1),examples(ee,2)).histsandthreshs.swhistbins,...
            SleepScoreMetrics(examples(ee,1),examples(ee,2)).histsandthreshs.swhist)
        plot(SleepScoreMetrics_IRASA(examples(ee,1),examples(ee,2)).histsandthreshs.swhistbins,...
            SleepScoreMetrics_IRASA(examples(ee,1),examples(ee,2)).histsandthreshs.swhist)
    end
    
    
subplot(3,3,2)
    imagesc(wins,swins,dipSW_IRASA')
    hold all
    for ee = 1:size(examples,1)
        plot(wins(examples(ee,1)),swins(examples(ee,2)),'+')
    end
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.065])
    title('Bimodality: IRASA')



subplot(4,1,3)
imagesc(SleepScoreMetrics_IRASA(examples(ee,1),examples(ee,2)).t_clus,...
    StatePlotMaterials_IRASA(examples(ee,1),examples(ee,2)).swFFTfreqs,...
    log10(StatePlotMaterials_IRASA(examples(ee,1),examples(ee,2)).swFFTspec))
axis xy
xlim(xwin)
SpecColorRange(log10(StatePlotMaterials(examples(ee,1),examples(ee,2)).swFFTspec))

    subplot(4,1,4)
    hold all
    for ee = 1:size(examples,1)
    plot(SleepScoreMetrics(examples(ee,1),examples(ee,2)).t_clus,...
        SleepScoreMetrics(examples(ee,1),examples(ee,2)).broadbandSlowWave)
    plot(SleepScoreMetrics_IRASA(examples(ee,1),examples(ee,2)).t_clus,...
        SleepScoreMetrics_IRASA(examples(ee,1),examples(ee,2)).broadbandSlowWave)
    end
    xlim(xwin)

    
%%
figure
subplot(4,1,1)
imagesc(SleepScoreMetrics(examples(ee,1),examples(ee,2)).t_clus,...
    StatePlotMaterials(examples(ee,1),examples(ee,2)).swFFTfreqs,...
    StatePlotMaterials(examples(ee,1),examples(ee,2)).swFFTspec)
axis xy
xlim(xwin)


    subplot(4,1,3)
    hold all
    for ee = 1:size(examples,1)
    plot(SleepScoreMetrics(examples(ee,1),examples(ee,2)).t_clus,...
        SleepScoreMetrics(examples(ee,1),examples(ee,2)).broadbandSlowWave)
    end
    xlim(xwin)
    subplot(4,1,4)
    hold all
    for ee = 1:size(examples,1)
    plot(SleepScoreMetrics(examples(ee,1),examples(ee,2)).t_clus,...
        SleepScoreMetrics(examples(ee,1),examples(ee,2)).thratio)
    end

    xlim(xwin)
