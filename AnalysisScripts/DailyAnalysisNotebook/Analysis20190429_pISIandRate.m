function [ ] = Analysis20190429(basePath,figfolder)
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
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];

%%
datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC';
regions = {'fCTX','CA1'};


%%
for rr = 1:length(regions)
    ISIstats.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true);
    numcells.(regions{rr}) = length(CellClass.(regions{rr}).UID);
end

%%
statenames = fieldnames(ISIstats.(regions{1}).summstats);
statecolors = {'k','b','r'};
numstates = length(statenames);

%%
rr = 1;
ss = 1;

%%
for rr = 1:length(regions)
    for ss = 1:3
        ratecorr.(regions{rr}).(statenames{ss}) = corr(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(CellClass.(regions{rr}).pE,:),...
            log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).pE))','type','spearman');
    end
end

%%
figure
plot(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(CellClass.(regions{rr}).pE,25),...
    log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).pE)),'.')
%%
figure

    for rr = 1:length(regions)
        subplot(2,2,rr)
        hold on
        plot(ISIstats.(regions{rr}).ISIhist.logbins([1 end]),[0 0],'k')
        for ss = 1:3
            plot(ISIstats.(regions{rr}).ISIhist.logbins(1,:),ratecorr.(regions{rr}).(statenames{ss}),'color',statecolors{ss})
        end
        LogScale('x',10)
    end



