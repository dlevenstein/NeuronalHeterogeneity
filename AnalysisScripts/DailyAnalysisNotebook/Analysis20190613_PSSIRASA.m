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
 xwin = bz_RandomWindowInIntervals(SleepScoreMetrics(1,1).t_clus([1 end]),2000);


figure

subplot(4,3,1)
    imagesc(wins,swins,dipSW')
    hold all
    plot(wins(10),swins(10),'k+')
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.065])
    title('Bimodality: Linear Fit')


    subplot(4,2,3)
    hold all
        plot(SleepScoreMetrics(10,10).histsandthreshs.swhistbins,...
            SleepScoreMetrics(10,10).histsandthreshs.swhist,'k')
        plot(SleepScoreMetrics_IRASA(2,15).histsandthreshs.swhistbins,...
            SleepScoreMetrics_IRASA(2,15).histsandthreshs.swhist,'r')
        xlabel('Slow Wave (PSS)')
    
    subplot(4,2,4)
    hold all
        plot(SleepScoreMetrics(10,10).histsandthreshs.THhistbins,...
            SleepScoreMetrics(10,10).histsandthreshs.THhist,'k')
        plot(SleepScoreMetrics_IRASA(2,15).histsandthreshs.THhistbins,...
            SleepScoreMetrics_IRASA(2,15).histsandthreshs.THhist,'r')
        xlabel('Theta Ratio')
    
subplot(4,3,3)
    imagesc(wins,swins,dipTH')
    hold all
        plot(wins(2),swins(15),'r+')
        plot(wins(10),swins(10),'k+')
    
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.065])
    title('Bimodality: Theta')

subplot(4,3,2)
    imagesc(wins,swins,dipSW_IRASA')
    hold all
        plot(wins(2),swins(15),'r+')
    
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.065])
    title('Bimodality: IRASA')

subplot(4,1,3)
    imagesc(SleepScoreMetrics_IRASA(examples(ee,1),examples(ee,2)).t_clus,...
        log2(StatePlotMaterials_IRASA(examples(ee,1),examples(ee,2)).swFFTfreqs),...
        log10(StatePlotMaterials_IRASA(examples(ee,1),examples(ee,2)).swFFTspec))
    axis xy
        hold all
        plot(SleepScoreMetrics(10,10).t_clus,...
            bz_NormToRange(SleepScoreMetrics(10,10).broadbandSlowWave),'k')
        plot(SleepScoreMetrics_IRASA(2,15).t_clus,...
            bz_NormToRange(SleepScoreMetrics_IRASA(2,15).broadbandSlowWave),'r')

    xlim(xwin)
    SpecColorRange(log10(StatePlotMaterials(examples(ee,1),examples(ee,2)).swFFTspec),[1.5 1.5])
    LogScale('y',2)
    bz_ScaleBar('s')
    ylabel('PSS (f, Hz)')



subplot(4,1,4)
    imagesc(SleepScoreMetrics_IRASA(examples(ee,1),examples(ee,2)).t_clus,...
        log2(StatePlotMaterials_IRASA(examples(ee,1),examples(ee,2)).thFFTfreqs),...
        log10(StatePlotMaterials_IRASA(examples(ee,1),examples(ee,2)).thFFTspec))
    axis xy
        hold all
        plot(SleepScoreMetrics(10,10).t_clus,...
            bz_NormToRange(SleepScoreMetrics(10,10).thratio),'k')
        plot(SleepScoreMetrics_IRASA(2,15).t_clus,...
            bz_NormToRange(SleepScoreMetrics_IRASA(2,15).thratio),'r')

    xlim(xwin)
    SpecColorRange(log10(StatePlotMaterials(examples(ee,1),examples(ee,2)).thFFTspec),[1.5 1.5])
    LogScale('y',2)
    bz_ScaleBar('s')
    ylabel('Theta (f, Hz)')

 NiceSave('WinIRASAOptimization',figfolder,baseName,'includeDate',true)

