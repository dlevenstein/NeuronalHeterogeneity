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
load(fullfile(basePath,[baseName,'.SleepScoreLFP.LFP.mat']))
load(fullfile(basePath,[baseName,'.EMGFromLFP.LFP.mat']))


%%
wins = 1:12;
swins = 1:25;


%% Get the metrics (PSS, theta)
for ww = 1:length(wins)
    bz_Counter(ww,length(wins),'Window')
    parfor ss = 1:length(swins)
        
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
            continue
        end
    end
end


%%
examples = [2 10; 10 1; 10 10; 2 15];
 xwin = bz_RandomWindowInIntervals(SleepScoreMetrics(1,1).t_clus([1 end]),1000);

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
    title('Slow Wave')
subplot(3,3,2)
    imagesc(wins,swins,dipTH')
    axis xy
    colorbar
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    title('Theta Ratio')

    subplot(3,3,3)
    hold all
    for ee = 1:size(examples,1)
        plot(SleepScoreMetrics(examples(ee,1),examples(ee,2)).histsandthreshs.swhistbins,...
            SleepScoreMetrics(examples(ee,1),examples(ee,2)).histsandthreshs.swhist)
    end

   

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
    NiceSave('SmoothWinOpt',figfolder,baseName,'includeDate',true)

    
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