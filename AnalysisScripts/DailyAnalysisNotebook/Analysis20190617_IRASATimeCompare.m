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
basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX/YMV08_170922';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%%
load(fullfile(basePath,[baseName,'.SleepScoreLFP.LFP.mat']))
load(fullfile(basePath,[baseName,'.EMGFromLFP.LFP.mat']))

sessionInfo = bz_getSessionInfo(basePath);

%%
sessionInfo = bz_getSessionInfo(basePath);

  tic      
[SleepScoreMetrics(1),StatePlotMaterials(1)] = ClusterStates_GetMetrics(...
                                           basePath,SleepScoreLFP,EMGFromLFP,true,...
                                           'window',4,'smoothfact',15,...
                                           'IRASA',true);
     toc
     tic
[SleepScoreMetrics(2),StatePlotMaterials(2)] = ClusterStates_GetMetrics(...
                                           basePath,SleepScoreLFP,EMGFromLFP,true,...
                                           'window',2,'smoothfact',15,'IRASA',true);
  toc
  
  
  
  %% Compare Methods for channel selection
  scoretime = [SleepScoreLFP.t(1) SleepScoreLFP.t(end)];
  SWWeightsName = 'PSS';
  tic
  [SleepScoreLFP_test(1),PickChannelStats(1)] = PickSWTHChannel(basePath,scoretime,SWWeightsName,...
    false,false,false,false,0,0,...
    sessionInfo.badchannels,true,'noPrompts',true,'saveFiles',false,'SHOWFIG',true,...
                            'window',2,'smoothfact',15,'IRASA',true,'downsamplefactor',5);

     toc;tic                   
  [SleepScoreLFP_test(2),PickChannelStats(2)] = PickSWTHChannel(basePath,scoretime,SWWeightsName,...
    false,false,false,false,0,0,...
    sessionInfo.badchannels,true,'noPrompts',true,'saveFiles',false,'SHOWFIG',true,...
                            'window',2,'smoothfact',15,'IRASA',false,'downsamplefactor',5);
            toc;tic            
  [SleepScoreLFP_test(3),PickChannelStats(3)] = PickSWTHChannel(basePath,scoretime,SWWeightsName,...
    false,false,false,false,0,0,...
    sessionInfo.badchannels,true,'noPrompts',true,'saveFiles',false,'SHOWFIG',true,...
                            'window',10,'smoothfact',10,'IRASA',true,'downsamplefactor',5);
                   toc;tic     
  [SleepScoreLFP_test(4),PickChannelStats(4)] = PickSWTHChannel(basePath,scoretime,SWWeightsName,...
    false,false,false,false,0,0,...
    sessionInfo.badchannels,true,'noPrompts',true,'saveFiles',false,'SHOWFIG',true,...
                            'window',10,'smoothfact',10,'IRASA',false,'downsamplefactor',5);
toc

%% 

end