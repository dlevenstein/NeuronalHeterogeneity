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

ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end); nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);
ISIStats.allspikes.CV2nm1 = cellfun(@(X) [nan; X(1:end-1)],...
    ISIStats.allspikes.CV2,'UniformOutput',false);

ISIStats.allspikes.meanCV2 = cellfun(@(prev,next) mean([prev next],2),...
    ISIStats.allspikes.CV2nm1,ISIStats.allspikes.CV2,'UniformOutput',false);

%% Restrict to state

for ss = 1:3
state = states{ss};
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);

%Both adjacent cpikes CV2
[ ISICV2 ] = cellfun(@(X,Y,Z,W) ConditionalHist( log10([X(W);X(W)]),[Y(W);Z(W)],...
    'Xbounds',[-2.6 1],'numXbins',75,'Ybounds',[0 2],'numYbins',75),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.CV2nm1,...
    ISIStats.allspikes.CV2,ISIStats.allspikes.instate,...
    'UniformOutput',false);

% %Both adjacent cpikes CV2 - mean normed ISI
% [ ISICV2 ] = cellfun(@(X,Y,Z,W) ConditionalHist( log10([X(W);X(W)]./mean([X(W);X(W)])),[Y(W);Z(W)],...
%     'Xbounds',[-3.5 1],'numXbins',75,'Ybounds',[0 2],'numYbins',75),...
%     ISIStats.allspikes.ISIs,ISIStats.allspikes.CV2nm1,...
%     ISIStats.allspikes.CV2,ISIStats.allspikes.instate,...
%     'UniformOutput',false);

%Mean of adj spikes cv2
% [ ISICV2 ] = cellfun(@(X,Y,W) ConditionalHist( log10(X(W)),Y(W),...
%     'Xbounds',[-2.6 1],'numXbins',75,'Ybounds',[0 2],'numYbins',75),...
%     ISIStats.allspikes.ISIs,ISIStats.allspikes.meanCV2,...
%     ISIStats.allspikes.instate,...
%     'UniformOutput',false);

ISICV2 = cat(1,ISICV2{:});
ISICV2 = CollapseStruct( ISICV2,3);


for tt = 1:length(celltypes)
    CV2distbyPOP.(state).(celltypes{tt}) = nanmean(ISICV2.pYX(:,:,CellClass.(celltypes{tt})),3);
    CV2distjointPOP.(state).(celltypes{tt}) = nanmean(ISICV2.XYprob(:,:,CellClass.(celltypes{tt})),3);
    meanCV2byPOP.(state).(celltypes{tt}) = nanmean(ISICV2.meanYX(:,:,CellClass.(celltypes{tt})),3);
end

end
%%
figure
for ss = 1:3   
for tt = 1:length(celltypes)
subplot(4,3,tt*3+ss-3)
    imagesc(ISICV2.Xbins(1,:,1),ISICV2.Ybins(1,:,1), (CV2distjointPOP.(states{ss}).(celltypes{tt}))')
    hold on
    plot(ISICV2.Xbins(1,:,1),meanCV2byPOP.(states{ss}).(celltypes{tt}),'w')
    axis xy
    LogScale('x',10)
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Power')
    xlabel('ISI (s)');
    %title((celltypes{tt}))
    colorbar
    if ss == 1
        ylabel({celltypes{tt},'CV2'})
    end
    if tt ==1 
        %caxis([0 0.03])
        caxis([0 0.0004])
        title(states{ss})
    elseif tt==2
         caxis([0 0.0007])
    end
    
    
subplot(4,3,tt*3+ss+6-3)
    imagesc(ISICV2.Xbins(1,:,1),ISICV2.Ybins(1,:,1), (CV2distbyPOP.(states{ss}).(celltypes{tt}))')
    hold on
    plot(ISICV2.Xbins(1,:,1),meanCV2byPOP.(states{ss}).(celltypes{tt}),'w')
    axis xy
    LogScale('x',10)
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Power')
    xlabel('ISI (s)');ylabel('CV2')
    title((celltypes{tt}))
    colorbar
    if tt ==1 
        caxis([0 0.025])
        %caxis([0 0.0004])
    elseif tt==2
         caxis([0 0.03])
    end
end 
end

NiceSave('CV2byISI',figfolder,baseName,'includeDate',true)

