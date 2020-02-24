function [ ] = AnalysisXXXXXXXX(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
datasetPath.BLA = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/GG_BLA';
datasetPath.PIR = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/GG_BLA';
regions = {'THAL','vCTX','fCTX','BLA','PIR','CA1'};
rnames =  {''    ,''    ,''    ,'bla','pir',''   };
%regions = {'fCTX'};
%%
for rr = 1:length(regions)
    [ISIstats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);
    numcells.(regions{rr}) = length(CellClass.(regions{rr}).UID);
    
    ISIstats.(regions{rr}) = rmfield( ISIstats.(regions{rr}),'allspikes');
    
    %Remove cells not in the proper region by removing their cell class!
    if ismember(rr,[4 5])
        inregion = cellfun(@(X) strcmp(X,rnames{rr}),ISIstats.(regions{rr}).cellinfo.regions);
        CellClass.(regions{rr}).label(~inregion)={[]};
        CellClass.(regions{rr}).pE(~inregion)=false;
        CellClass.(regions{rr}).pI(~inregion)=false;
    end
end


%%
statenames = {'NREMstate','WAKEstate','REMstate'};
%%
numAS = 3;
rr = 6;
ss = 1;
logtimebins = ISIstats.(regions{rr}).ISIhist.logbins(1,:);
logISIhist = ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(CellClass.(regions{rr}).pE,:)';
logISIhist = logISIhist./mode(diff(logtimebins));
[sharedfit,singlecell] = bz_FitISISharedGammaModes(logISIhist,logtimebins,'numAS',numAS,'showfig',true);

