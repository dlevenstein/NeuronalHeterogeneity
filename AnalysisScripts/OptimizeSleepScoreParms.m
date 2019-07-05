function [dipmap] = OptimizeSleepScoreParms(basePath,figfolder)
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
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC/Achilles_10252013';
%basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX/YMV08_170922';
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%%
load(fullfile(basePath,[baseName,'.SleepScoreLFP.LFP.mat']))
load(fullfile(basePath,[baseName,'.EMGFromLFP.LFP.mat']))
SleepScoreLFP.params.SWweights = 'PSS';
SleepScoreLFP.params.SWfreqlist =[];

EMGFromLFP = EMGFromLFP;
%%
if ~exist('SleepScoreLFP','var')
    display(basePath)
    display(baseName)
end
%%
sessionInfo = bz_getSessionInfo(basePath);

%% Compare Bimodality IRASA vs no

%%
wins = 1:12;
swins = 1:25;


    %% Set up for parallel in cluster
    pc = parcluster('local');
    % store temporary files in the 'scratch' drive on the cluster, labeled by job ID
    pc.JobStorageLocation = strcat(getenv('SCRATCH'), '/', getenv('SLURM_JOB_ID'));
    % enable MATLAB to utilize the multiple cores allocated in the job script
    % SLURM_NTASKS_PER_NODE is a variable set in the job script by the flag --tasks-per-node
    % we use SLURM_NTASKS_PER_NODE - 1, because one of these tasks is the original MATLAB script itself
    parpool(pc, str2num(getenv('SLURM_NTASKS_PER_NODE'))-1);

%% Get the metrics (PSS, theta)
for ww = 1:length(wins)
    bz_Counter(ww,length(wins),'Window')
    
    
    %%
    parfor ss = 1:length(swins)
    %for ss = 1:length(swins)
        %SleepScoreLFP
[SleepScoreMetrics_IRASA(ww,ss),StatePlotMaterials_IRASA(ww,ss)] = ClusterStates_GetMetrics(...
                                           basePath,SleepScoreLFP,EMGFromLFP,true,...
                                           'window',wins(ww),'smoothfact',swins(ss),...
                                           'IRASA',true);
                                       
[SleepScoreMetrics(ww,ss),StatePlotMaterials(ww,ss)] = ClusterStates_GetMetrics(...
                                           basePath,SleepScoreLFP,EMGFromLFP,true,...
                                           'window',wins(ww),'smoothfact',swins(ss),...
                                           'IRASA',false);
    end
end
       
%% Test for bimodality
dipmap.SW = nan(length(wins),length(swins));
dipmap.TH = nan(length(wins),length(swins));
dipmap.SW_IRASA = nan(length(wins),length(swins));
dipmap.TH_IRASA = nan(length(wins),length(swins));

for ww = 1:length(wins)
    for ss = 1:length(swins)
        try
%         dipmap.SW(ww,ss) = bz_hartigansdiptest(sort(SleepScoreMetrics(ww,ss).broadbandSlowWave));
%         dipmap.TH(ww,ss) = bz_hartigansdiptest(sort(SleepScoreMetrics(ww,ss).thratio));
            dipmap.SW(ww,ss) = SleepScoreMetrics(ww,ss).SWdiptest;
            dipmap.TH(ww,ss) = SleepScoreMetrics(ww,ss).THdiptest;
        catch
        end
        
        try
       % dipmap.SW_IRASA(ww,ss) = bz_hartigansdiptest(sort(SleepScoreMetrics_IRASA(ww,ss).broadbandSlowWave));
        %dipmap.TH_IRASA(ww,ss) = bz_hartigansdiptest(sort(SleepScoreMetrics_IRASA(ww,ss).thratio));
            dipmap.SW_IRASA(ww,ss) = SleepScoreMetrics_IRASA(ww,ss).SWdiptest;
            dipmap.TH_IRASA(ww,ss) = SleepScoreMetrics_IRASA(ww,ss).THdiptest;
        catch  
        end
    end
end

dipmap.wins = wins;
dipmap.swins = swins;
%% Stuff for Saving

%%

%examples = [10 10; 2 15];
 xwin = bz_RandomWindowInIntervals(SleepScoreMetrics(1,1).t_clus([1 end]),2000);


figure

subplot(4,4,1)
    imagesc(wins,swins,dipmap.SW')
    hold all
    plot(wins(10),swins(10),'k+')
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.065])
    title('Bimodality: Linear Fit')


    
subplot(4,4,3)
    imagesc(wins,swins,dipmap.TH')
    hold all
        %plot(wins(2),swins(15),'r+')
        plot(wins(10),swins(10),'k+')
    
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.01])
    title('Bimodality: Theta IRASA')
    
subplot(4,4,4)
    imagesc(wins,swins,dipmap.TH_IRASA')
    hold all
        %plot(wins(2),swins(15),'r+')
        plot(wins(10),swins(10),'k+')
    
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.01])
    title('Bimodality: Theta')

subplot(4,4,2)
    imagesc(wins,swins,dipmap.SW_IRASA')
    hold all
        plot(wins(2),swins(15),'r+')
    
    xlabel('Win size (s)'); ylabel('Smooth window (s)')
    axis xy
    colorbar
    caxis([0 0.065])
    title('Bimodality: IRASA')
    
    
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

subplot(4,1,3)
    imagesc(SleepScoreMetrics_IRASA(2,15).t_clus,...
        log2(StatePlotMaterials_IRASA(2,15).swFFTfreqs),...
        log10(StatePlotMaterials_IRASA(2,15).swFFTspec))
    axis xy
        hold all
        plot(SleepScoreMetrics(10,10).t_clus,...
            bz_NormToRange(SleepScoreMetrics(10,10).broadbandSlowWave),'k')
        plot(SleepScoreMetrics_IRASA(2,15).t_clus,...
            bz_NormToRange(SleepScoreMetrics_IRASA(2,15).broadbandSlowWave),'r')

    xlim(xwin)
    SpecColorRange(log10(StatePlotMaterials(2,15).swFFTspec),[1.5 1.5])
    LogScale('y',2)
    bz_ScaleBar('s')
    ylabel('PSS (f, Hz)')



subplot(4,1,4)
    imagesc(SleepScoreMetrics_IRASA(2,15).t_clus,...
        log2(StatePlotMaterials_IRASA(2,15).thFFTfreqs),...
        (StatePlotMaterials_IRASA(2,15).thFFTspec))
    axis xy
        hold all
        plot(SleepScoreMetrics(10,10).t_clus,...
            bz_NormToRange(SleepScoreMetrics(10,10).thratio),'k')
        plot(SleepScoreMetrics_IRASA(2,15).t_clus,...
            bz_NormToRange(SleepScoreMetrics_IRASA(2,15).thratio),'r')

    xlim(xwin)
    SpecColorRange((StatePlotMaterials(2,15).thFFTspec),[1.5 1.5])
    LogScale('y',2)
    bz_ScaleBar('s')
    ylabel('Theta (f, Hz)')

 NiceSave('WinIRASAOptimization',figfolder,baseName)

