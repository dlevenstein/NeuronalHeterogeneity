function [ ] = AnalysisXXXXXXXX(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%%
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
regions = {'THAL','vCTX','fCTX','CA1'};
%regions = {'fCTX'};
%%
for rr = 1:length(regions)
    [ISIstats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);
    OccupancyStats.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'OccupancyStats','dataset',true,'catall',true,'baseNames',baseNames);

    numcells.(regions{rr}) = length(CellClass.(regions{rr}).UID);
    
end

%%
statenames = fieldnames(ISIstats.(regions{1}).summstats);
statecolors = {[0 0 0],[0 0 1],[1 0 0]};
numstates = length(statenames);

%%
for rr = 1:length(regions)
    for ss = 1:3
        OccupancyStats.(regions{rr}).(statenames{ss}).MTORate = (1./OccupancyStats.(regions{rr}).(statenames{ss}).median);
        OccupancyStats.(regions{rr}).(statenames{ss}).MTORatio = ...
            ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate./...
            (1./OccupancyStats.(regions{rr}).(statenames{ss}).median);
    end
end
%%
figure
for rr = 1:length(regions)
for ss =1:3
subplot(6,4,rr+(ss-1)*4)
scatter(log10(OccupancyStats.(regions{rr}).(statenames{ss}).MTORate),...
    log10(OccupancyStats.(regions{rr}).(statenames{ss}).MTORatio),0.5,...
     log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate))
 axis tight
 
 ylim([0 2])
 xlim([-2.75 1.75])
 LogScale('xy',10,'exp',true,'nohalf',true)
 
  caxis([-1.5 2])
  
 colorbar
 LogScale('c',10,'exp',true,'nohalf',true)
 if ss==1
     title(regions{rr})
 end
  if rr ==1
    ylabel('Activation Ratio') 
 end

subplot(6,4,rr+(ss-1)*4+12)
scatter(log10(OccupancyStats.(regions{rr}).(statenames{ss}).MTORate),...
    log10(OccupancyStats.(regions{rr}).(statenames{ss}).MTORatio),0.5,...
     (ISIstats.(regions{rr}).summstats.(statenames{ss}).meanCV2)) 
  axis tight
 ylim([0 2])
  xlim([-2.75 1.75])
 caxis([0.75 1.25])
 colorbar
 crameri('berlin','pivot',1)
 LogScale('xy',10,'exp',true,'nohalf',true)
 if ss ==3
     xlabel('Ground State Rate (MTO)')
 end
 if rr ==1
    ylabel('Activation Ratio') 
 end
end
end
NiceSave('RateCV2byGSAS',figfolder,[])

