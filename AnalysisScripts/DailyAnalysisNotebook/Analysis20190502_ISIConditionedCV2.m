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
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};



%% Restrict to state

state = states{2};
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);

ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end); nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);
%%
[ ISICV2 ] = cellfun(@(X,Y,Z,W) ConditionalHist( log10([X(W);Y(W)]),[Z(W);Z(W)],...
    'Xbounds',[-2.6 1],'numXbins',100,'Ybounds',[0 2],'numYbins',100),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.CV2,ISIStats.allspikes.instate,...
    'UniformOutput',false);

ISICV2 = cat(1,ISICV2{:});
ISICV2 = CollapseStruct( ISICV2,3);


for tt = 1:length(celltypes)
    CV2distbyPOP.(celltypes{tt}) = nanmean(ISICV2.XYprob(:,:,CellClass.(celltypes{tt})),3);
    meanCV2byPOP.(celltypes{tt}) = nanmean(ISICV2.meanYX(:,:,CellClass.(celltypes{tt})),3);
end

%%
figure
   
for tt = 1:length(celltypes)
subplot(3,2,tt*2-1)
    imagesc(ISICV2.Xbins(1,:,1),ISICV2.Ybins(1,:,1), (CV2distbyPOP.(celltypes{tt}))')
    hold on
    plot(ISICV2.Xbins(1,:,1),meanCV2byPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('x',10)
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Power')
    xlabel('ISI (s)');ylabel('CV2')
    title((celltypes{tt}))
    colorbar
    if tt ==1 
        %caxis([0 0.03])
        caxis([0 0.00025])
    elseif tt==2
         caxis([0 0.0004])
    end
    
end 
