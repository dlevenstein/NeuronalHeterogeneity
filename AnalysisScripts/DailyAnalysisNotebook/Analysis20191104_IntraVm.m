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
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
basePath = '/mnt/NyuShare/valerm05/ihf/pvch6/pvch6_180607_164835/';
%basePath = '/mnt/NyuShare/valerm05/ihf/pvch5/pvch5_180518_115703/';

%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%%
intraname = fullfile(basePath,'dataCell.mat');
load(intraname)

%

%%
for cc = 1:length(spk2_cell)
    offset = (intraSpk{cc}(1)./fs_spk2)-intraSpk_timesIntan{cc}(1);
    
    spk2_timestamps{cc} = [1:length(spk2_cell{cc})]./fs_spk2-offset;
end
%%

figure
for whichcell = 1:3
    subplot(3,1,whichcell)
plot(spk2_timestamps{whichcell},spk2_cell{whichcell}(1,:))
hold on
plot(intraSpk_timesIntan{whichcell},ones(size(intraSpk_timesIntan{whichcell})),'r+')
end
%%
ExIn = cell(size(all_spikes.UID));
ExIn(1:end) = {'Ex'};
ExIn(all_spikes.UID==0)={'In'};
%%
ISIStats = bz_ISIStats(all_spikes,'showfig',true,'cellclass',ExIn);

%%
%For each point in time (Vm) - get the ISI it's in and time since last
%spike

%First: conditional distribution of Vm|ISI


