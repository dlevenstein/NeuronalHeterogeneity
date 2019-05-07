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


%% Load the LFP if needed

% lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID;
% downsamplefactor = 1;
% lfp = bz_GetLFP(lfpchan,...
%     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
% %Noralize the LFP
% lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);



%% Get SleepScoreParams at each spike time

ISIStats.allspikes.PSS =cellfun(@(X) ...
    interp1(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus,...
    SleepState.detectorinfo.detectionparms.SleepScoreMetrics.broadbandSlowWave,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);

ISIStats.allspikes.thetarat =cellfun(@(X) ...
    interp1(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus,...
    SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);

ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end); nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);



%% Restrict to state

state = states{1};
ISIStats.allspikes.instate_PSS = cellfun(@(X) InIntervals(X,double(SleepState.ints.ALL)),...
    ISIStats.allspikes.times,'UniformOutput',false);

ISIStats.allspikes.instate_notTheta = cellfun(@(X) InIntervals(X,double(SleepState.ints.NREMstate)),...
    ISIStats.allspikes.times,'UniformOutput',false);


%% COnditional Hists
[ ISIbyPSS ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(W);Z(W)],log10([X(W);Y(W)]),...
    'Xbounds',[0 1],'numXbins',60,'Ybounds',[-3 2],'numYbins',125,'minX',50),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.PSS,ISIStats.allspikes.instate_PSS,...
    'UniformOutput',false);
ISIbyPSS = cat(1,ISIbyPSS{:});
ISIbyPSS = CollapseStruct( ISIbyPSS,3);

[ ISIbytheta ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(~W);Z(~W)],log10([X(~W);Y(~W)]),...
    'Xbounds',[0 1],'numXbins',60,'Ybounds',[-3 2],'numYbins',125,'minX',50),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.thetarat,ISIStats.allspikes.instate_notTheta,...
    'UniformOutput',false);
ISIbytheta = cat(1,ISIbytheta{:});
ISIbytheta = CollapseStruct( ISIbytheta,3);


for tt = 1:length(celltypes)
    ISIbyPSS.pop.(celltypes{tt}) = nanmean(ISIbyPSS.pYX(:,:,CellClass.(celltypes{tt})),3);
    ISIbytheta.pop.(celltypes{tt}) = nanmean(ISIbytheta.pYX(:,:,CellClass.(celltypes{tt})),3);

end

%%
figure
   
for tt = 1:length(celltypes)
subplot(3,2,tt*2-1)
    imagesc(ISIbyPSS.Xbins(1,:,1),ISIbyPSS.Ybins(1,:,1), ISIbyPSS.pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('y',10)
    ylabel({(celltypes{tt}),'ISI (s)'});xlabel('PSS')
    %title((celltypes{tt}))
    colorbar
    if tt ==1 
        caxis([0 0.02])
    elseif tt==2
         caxis([0 0.03])
    end
    
end 


for tt = 1:length(celltypes)
subplot(3,2,tt*2)
    imagesc(ISIbytheta.Xbins(1,:,1),ISIbytheta.Ybins(1,:,1), ISIbytheta.pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('y',10)
    ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
    colorbar
    if tt ==1 
        caxis([0 0.02])
    elseif tt==2
         caxis([0 0.03])
    end
    
end 

NiceSave('ISIbyStateVars',figfolder,baseName,'includeDate',true)

