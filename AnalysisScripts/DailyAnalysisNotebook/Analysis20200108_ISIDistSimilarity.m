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
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
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
    numcells.(regions{rr}) = length(CellClass.(regions{rr}).UID);
end

%%
statenames = fieldnames(ISIstats.(regions{1}).summstats);
statecolors = {[0 0 0],[0 0 1],[1 0 0]};
numstates = length(statenames);
%% Comiple all histograms for tSNE
allISIhists.hists = [];
allISIhists.region = [];
allISIhists.celltype.pE = [];
allISIhists.celltype.pI = [];
allISIhists.state = [];

numPerregion = 400;
numISIthresh = 500;


for rr = 1:length(regions)
    for ss = 1:3
        whichcells = find((CellClass.(regions{rr}).pE | CellClass.(regions{rr}).pI) &...
            ~isnan(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(:,1))'); 
        whichcells = randsample(whichcells,numPerregion);
        allISIhists.region = [allISIhists.region rr.*ones(size(whichcells))];
        allISIhists.celltype.pE = [allISIhists.celltype.pE CellClass.(regions{rr}).pE(whichcells)];
        allISIhists.celltype.pI = [allISIhists.celltype.pI CellClass.(regions{rr}).pI(whichcells)];
        allISIhists.state = [allISIhists.state ss.*ones(size(whichcells))];
        allISIhists.hists = [allISIhists.hists ; ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(whichcells,:)];
    end
end
%% tSNE the histograms!
perplexity = 20;
no_dims = 2;
initial_dims = 30;
close all
tSNEmap = tsne(allISIhists.hists, allISIhists.celltype.pE, no_dims, initial_dims, perplexity);

%%

figure
subplot(3,3,[5 8])
hold on
plot(tSNEmap(allISIhists.state==2,1),tSNEmap(allISIhists.state==2,2),'b.')
plot(tSNEmap(allISIhists.state==1,1),tSNEmap(allISIhists.state==1,2),'k.')
plot(tSNEmap(allISIhists.state==3,1),tSNEmap(allISIhists.state==3,2),'r.')
axis tight
xticks(gca,[]);yticks(gca,[])
legend(statenames,'location','northoutside')

%legend('NREM','WAKE','REM')
subplot(3,3,[6 9])
hold on
for rr = 1:length(regions)
    plot(tSNEmap(allISIhists.region==rr,1),tSNEmap(allISIhists.region==rr,2),'.')
end
legend(regions,'location','northoutside')
axis tight
xticks(gca,[]);yticks(gca,[])
subplot(3,3,[4 7])
hold on
axis tight

xticks(gca,[]);yticks(gca,[])
plot(tSNEmap(allISIhists.celltype.pE==1,1),tSNEmap(allISIhists.celltype.pE==1,2),'.k')
plot(tSNEmap(allISIhists.celltype.pI==1,1),tSNEmap(allISIhists.celltype.pI==1,2),'.r')
legend(celltypes,'location','northoutside')
NiceSave('tSNE_Map',figfolder,'')
