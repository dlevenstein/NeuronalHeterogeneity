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
basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
%basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
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


position = bz_LoadBehavior( basePath,'position' );

%%
ISIStats.allspikes.position = cellfun(@(X) interp1(position.timestamps,position.position.lin,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);
ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end);nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);

%%
[ ISIbyPos ] = cellfun(@(X,Y,Z) ConditionalHist( [Z;Z],log10([X;Y]),...
    'Xbounds',[0 max(position.position.lin)],'numXbins',25,'Ybounds',[-3 2],'numYbins',125,'minX',15,...
    'conditionby',position.position.lin),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.position,...
    'UniformOutput',false);
ISIbyPos = cat(1,ISIbyPos{:});
ISIbyPos = CollapseStruct( ISIbyPos,3);
ISIbyPos.pYX = ISIbyPos.pYX./mode(diff(position.timestamps));
ISIbyPos.rate = sum(ISIbyPos.pYX,2);
%%
excell = [31 68 89 90 91 94];
excell = [31 89 90 94];
figure
for ee = 1:length(excell)

subplot(2,2,ee)
    imagesc(ISIbyPos.Xbins(:,:,excell(ee)),ISIbyPos.Ybins(:,:,excell(ee)),ISIbyPos.pYX(:,:,excell(ee))')
    hold on
    plot(ISIStats.allspikes.position{excell(ee)},log10(ISIStats.allspikes.ISIs{excell(ee)}),'k.','markersize',0.1)
    plot(ISIStats.allspikes.position{excell(ee)},log10(ISIStats.allspikes.ISInp1{excell(ee)}),'k.','markersize',0.1)
    plot(ISIStats.ISIhist.WAKEstate.log(excell(ee),:)*10+max(ISIbyPos.Xbins(:,:,excell(ee))),ISIStats.ISIhist.logbins,'k','linewidth',2)
    plot(ISIbyPos.Xbins(:,:,excell(ee)),log10(1./ISIbyPos.rate(:,:,excell(ee))),'k','linewidth',2)
    %
    axis tight
    ylim([-3 2])
    
    xlabel('Position (m)')
    title(excell(ee))
    box off
    ylabel('ISI (s)')
    LogScale('y',10,'exp',true)

% pause
% close all
end

    NiceSave('PlaceFields',figfolder,baseName,'includeDate',true)
