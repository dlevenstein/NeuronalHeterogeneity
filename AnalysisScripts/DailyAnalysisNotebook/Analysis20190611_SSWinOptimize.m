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

%% Load the SleepScoreLFP and EMG
load('YMV08_170922.SleepScoreLFP.LFP.mat')
load('YMV08_170922.EMGFromLFP.LFP.mat')

%%
wins = 2:15;
swins = 1:15;


%% Get the metrics (PSS, theta)
for ww = 1:length(wins)
    bz_Counter(ww,length(wins),'Window')
    parfor ss = 1:length(swins)
        
[SleepScoreMetrics(ww,ss)] = ClusterStates_GetMetrics(...
                                           basePath,SleepScoreLFP,EMGFromLFP,true,...
                                           'window',wins(ww),'smoothfact',swins(ss));
    end
end
       

%% Test for bimodality
for ww = 1:length(wins)
    for ss = 1:length(swins)
        try
        dipSW(ww,ss) = bz_hartigansdiptest(sort(SleepScoreMetrics(ww,ss).broadbandSlowWave));
        dipTH(ww,ss) = bz_hartigansdiptest(sort(SleepScoreMetrics(ww,ss).thratio));
        catch
            dipSW(ww,ss) = nan;
            dipTH(ww,ss) = nan;
            continue
        end
    end
end
   

%%
figure
subplot(3,3,1)
    imagesc(wins,swins,dipSW')
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    title('Slow Wave')
subplot(3,3,2)
    imagesc(wins,swins,dipTH')
    axis xy
    colorbar
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    title('Theta Ratio')

    subplot(3,3,3)
    plot(SleepScoreMetrics(9,10).histsandthreshs.swhistbins,SleepScoreMetrics(9,10).histsandthreshs.swhist,'k')
    hold on
    plot(SleepScoreMetrics(1,10).histsandthreshs.swhistbins,SleepScoreMetrics(1,10).histsandthreshs.swhist,'r')
    plot(SleepScoreMetrics(9,1).histsandthreshs.swhistbins,SleepScoreMetrics(9,1).histsandthreshs.swhist,'b')
    xwin = bz_RandomWindowInIntervals(SleepScoreMetrics(1,1).t_clus([1 end]),2000);

    subplot(4,1,3)
    plot(SleepScoreMetrics(9,10).t_clus,SleepScoreMetrics(9,10).broadbandSlowWave,'k')
    hold on
    plot(SleepScoreMetrics(1,10).t_clus,SleepScoreMetrics(1,10).broadbandSlowWave,'r')
    plot(SleepScoreMetrics(9,1).t_clus,SleepScoreMetrics(9,1).broadbandSlowWave,'b')
    %plot(SleepScoreMetrics(1,1).t_clus,SleepScoreMetrics(1,1).broadbandSlowWave,'g')
    xlim(xwin)
    subplot(4,1,4)
    plot(SleepScoreMetrics(9,10).t_clus,SleepScoreMetrics(9,10).thratio,'k')
    hold on
    plot(SleepScoreMetrics(1,10).t_clus,SleepScoreMetrics(1,10).thratio,'r')
    plot(SleepScoreMetrics(9,1).t_clus,SleepScoreMetrics(9,1).thratio,'b')
   % plot(SleepScoreMetrics(1,1).t_clus,SleepScoreMetrics(1,1).thratio,'g')
    xlim(xwin)
    NiceSave('SmoothWinOpt',figfolder,baseName,'includeDate',true)

    
