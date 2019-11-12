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
[spikes] = chasingSpikes(spk2_cell{2})

%%
ISIs = diff(spikes{1});
%%

thiscell = 2;
intra.timestamps = [1:length(spk2_cell{thiscell})]'./fs_spk2;
intra.Vm = spk2_cell{thiscell}(1,:)';
intra.spiketimes = intra.timestamps(spikes{1});
intra.ISIs = diff(intra.spiketimes);
intra.nextISI = interp1(intra.spiketimes(2:end),intra.ISIs,intra.timestamps,'next');
%%
ConditionalHist(log10(intra.nextISI),intra.Vm,'SHOWFIG',true,'Xbounds',[-2.5 1.25],...
    'Ybounds',[-75 -40])
xlabel('ISI (s)')
LogScale('x',10)
ylabel('V_m (mV)')
NiceSave('Vmtest',figfolder,baseName)
%%
figure
plot(intra.timestamps,intra.Vm)
hold on
plot(intra.spiketimes,zeros(size(intra.spiketimes)),'r+')
hold on
%plot(intra.timestamps,intra.nextISI)
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


