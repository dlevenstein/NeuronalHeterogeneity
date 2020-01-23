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
%basePath = pwd;
basePath = fullfile(reporoot,'Datasets/onProbox/AP_THAL/Mouse12-120807');
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis/HDCellex'];
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

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};

%%
%bz_GetLFP(1,'basePath',basePath);
%%
headdir.samplingRate = 39.06; %Hz
hdfilename = fullfile(basePath,[baseName,'.ang']);
headdir.data = importdata(hdfilename);
headdir.data(headdir.data==-1)=nan;
headdir.timestamps = [1:length(headdir.data)]'./headdir.samplingRate; 

%%
state = states{1};
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);
headdir.instate = InIntervals(headdir.timestamps,double(SleepState.ints.(state)));
%%
ISIStats.allspikes.position = cellfun(@(X) interp1(headdir.timestamps,headdir.data,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);
ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end);nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);
%%
[ ISIbyPos ] = cellfun(@(X,Y,Z,Q) ConditionalHist( [Z(Q);Z(Q)],log10([X(Q);Y(Q)]),...
    'Xbounds',[0 2*pi],'numXbins',20,'Ybounds',[-3 2],'numYbins',125,'minX',25),...%,'conditionby',headdir.data(headdir.instate)),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.position,ISIStats.allspikes.instate,...
    'UniformOutput',false);

ISIbyPos = cat(1,ISIbyPos{:});
ISIbyPos = CollapseStruct( ISIbyPos,3);
%ISIbyPos.pYX = ISIbyPos.pYX./mode(diff(T));
ISIbyPos.rate = sum(ISIbyPos.pYX,2);

%%
figure
for excell =1:spikes.numcells;
subplot(4,3,mod(excell,12)+1)
imagesc(ISIbyPos.Xbins(:,:,1),ISIbyPos.Ybins(:,:,1),ISIbyPos.pYX(:,:,excell)')
hold on
imagesc(ISIbyPos.Xbins(:,:,1)+2*pi,ISIbyPos.Ybins(:,:,1),ISIbyPos.pYX(:,:,excell)')
%axis xy
plot(ISIbyPos.Xbins(:,:,1),bz_NormToRange(-ISIbyPos.pX(:,:,excell),'ylim',[-0.3 0]),'r')
plot(ISIbyPos.Xbins(:,:,1)+2*pi,bz_NormToRange(-ISIbyPos.pX(:,:,excell),'ylim',[-0.3 0]),'r')
title(['UID: ',num2str(spikes.UID(excell))])
xlim([0 4*pi])
bz_piTickLabel('x')
xlabel('HD (rad)');ylabel('ISI (s)')
LogScale('y',10,'exp',true)
%drawnow
%pause
if (mod(excell,12)==0 && excell>0) || excell == spikes.numcells
NiceSave(['UID: ',num2str(spikes.UID(excell))],figfolder,baseName)
figure
end
end
