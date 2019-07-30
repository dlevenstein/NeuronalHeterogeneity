function [BShist,ISIbyPSS,ISIbytheta,statehists ] = SpikeStatsbyBrainState(basePath,figfolder)
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
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC/Achilles_10252013';
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
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
statecolors = {[0 0 0],[0 0 1],[1 0 0],[0.6 0.6 0.6]};

try
celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};


%% Load the LFP if needed

% lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID;
% downsamplefactor = 1;
% lfp = bz_GetLFP(lfpchan,...
%     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
% %Noralize the LFP
% lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

%% Normalize the brain state metrics
BSmetrics.timestamps = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus;
BSmetrics.EMG = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.EMG;
BSmetrics.PSS = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.broadbandSlowWave;
BSmetrics.thratio = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio;


for ss = 1:length(states)
   BSmetrics.instatetime.(states{ss}) = InIntervals(BSmetrics.timestamps,SleepState.ints.(states{ss}));
   
%    [N,X] = hist(BSmetrics.thratio(BSmetrics.instatetime.(states{ss})),20);
%    [~,BSmetrics.peaks.(states{ss}).thratio] = max(N);
%    BSmetrics.peaks.(states{ss}).thratio = X(BSmetrics.peaks.(states{ss}).thratio);
    BSmetrics.peaks.(states{ss}).thratio = median(BSmetrics.thratio(BSmetrics.instatetime.(states{ss})));
   
%    [N,X] = hist(BSmetrics.PSS(BSmetrics.instatetime.(states{ss})),20);
%    [~,BSmetrics.peaks.(states{ss}).PSS] = max(N);
%    BSmetrics.peaks.(states{ss}).PSS = X(BSmetrics.peaks.(states{ss}).PSS);
    BSmetrics.peaks.(states{ss}).PSS = median(BSmetrics.PSS(BSmetrics.instatetime.(states{ss})));

end

% BSmetrics.thratio = (BSmetrics.thratio-BSmetrics.peaks.WAKEstate.thratio)./...
%     (BSmetrics.peaks.REMstate.thratio-BSmetrics.peaks.WAKEstate.thratio);
% BSmetrics.PSS = (BSmetrics.PSS-BSmetrics.peaks.WAKEstate.PSS)./...
%     (BSmetrics.peaks.NREMstate.PSS-BSmetrics.peaks.WAKEstate.PSS);

%% BS metrics histogram by states
BShist.bins = linspace(0,1,50);
for ss = 1:length(states)
   BShist.(states{ss}).PSS = hist(BSmetrics.PSS(BSmetrics.instatetime.(states{ss})),BShist.bins);
   BShist.(states{ss}).PSS = BShist.(states{ss}).PSS./length(BSmetrics.timestamps);
   BShist.(states{ss}).thratio = hist(BSmetrics.thratio(BSmetrics.instatetime.(states{ss})),BShist.bins);
   BShist.(states{ss}).thratio = BShist.(states{ss}).thratio./length(BSmetrics.timestamps);
end


%% Get SleepScoreParams at each spike time

ISIStats.allspikes.PSS =cellfun(@(X) ...
    interp1(BSmetrics.timestamps,BSmetrics.PSS,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);

ISIStats.allspikes.thetarat =cellfun(@(X) ...
    interp1(BSmetrics.timestamps,BSmetrics.thratio,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);

ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end); nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);



%% Restrict to state

%state = states{1};
for ss = 1:length(states)
    ISIStats.allspikes.instate.(states{ss}) = cellfun(@(X) InIntervals(X,double(SleepState.ints.(states{ss}))),...
        ISIStats.allspikes.times,'UniformOutput',false);

ISIStats.allspikes.instate_PSS = cellfun(@(X) InIntervals(X,double(SleepState.ints.ALL)),...
    ISIStats.allspikes.times,'UniformOutput',false);

ISIStats.allspikes.instate_notTheta = cellfun(@(X) InIntervals(X,double(SleepState.ints.NREMstate)),...
    ISIStats.allspikes.times,'UniformOutput',false);
end


%% COnditional Hists
%[-0.6 1.5] (median normalize)
[ ISIbyPSS ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(W);Z(W)],log10([X(W);Y(W)]),...
    'Xbounds',[0 1],'numXbins',20,'Ybounds',[-3 2],'numYbins',125,'minX',20),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.PSS,ISIStats.allspikes.instate.ALL,...
    'UniformOutput',false);
ISIbyPSS = cat(1,ISIbyPSS{:});
ISIbyPSS = CollapseStruct( ISIbyPSS,3);

%[-1.5 1.75]
[ ISIbytheta ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(~W);Z(~W)],log10([X(~W);Y(~W)]),...
    'Xbounds',[0 1],'numXbins',20,'Ybounds',[-3 2],'numYbins',125,'minX',20),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.thetarat,ISIStats.allspikes.instate.NREMstate,...
    'UniformOutput',false);
ISIbytheta = cat(1,ISIbytheta{:});
ISIbytheta = CollapseStruct( ISIbytheta,3);

[ ISIbytheta.WAKE ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(W);Z(W)],log10([X(W);Y(W)]),...
    'Xbounds',[0 1],'numXbins',20,'Ybounds',[-3 2],'numYbins',125,'minX',20),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.thetarat,ISIStats.allspikes.instate.WAKEstate,...
    'UniformOutput',false);
ISIbytheta.WAKE = cat(1,ISIbytheta.WAKE{:});
ISIbytheta.WAKE = CollapseStruct( ISIbytheta.WAKE,3);

[ ISIbytheta.REM ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(W);Z(W)],log10([X(W);Y(W)]),...
    'Xbounds',[0 1],'numXbins',20,'Ybounds',[-3 2],'numYbins',125,'minX',20),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.thetarat,ISIStats.allspikes.instate.REMstate,...
    'UniformOutput',false);
ISIbytheta.REM = cat(1,ISIbytheta.REM{:});
ISIbytheta.REM = CollapseStruct( ISIbytheta.REM,3);


for tt = 1:length(celltypes)
    ISIbyPSS.pop.(celltypes{tt}) = nanmean(ISIbyPSS.pYX(:,:,CellClass.(celltypes{tt})),3);
    ISIbytheta.pop.(celltypes{tt}) = nanmean(ISIbytheta.pYX(:,:,CellClass.(celltypes{tt})),3);

    ISIbyPSS.celltypeidx.(celltypes{tt}) = CellClass.(celltypes{tt});
    ISIbytheta.celltypeidx.(celltypes{tt}) = CellClass.(celltypes{tt});
end
    if length(celltypes)==1
        ISIbyPSS.celltypeidx.pI = false(size(CellClass.pE));
        ISIbytheta.celltypeidx.pI = false(size(CellClass.pE));
    end
    
%% State variable histograms 
statehists.PSSbins = BShist.bins;
statehists.thetabins = BShist.bins;
statehists.EMGbins = BShist.bins;

for ss = 1:length(states)
    statehists.(states{ss}).PSS = hist3([BSmetrics.PSS(BSmetrics.instatetime.(states{ss})) BSmetrics.EMG(BSmetrics.instatetime.(states{ss}))],...
        {statehists.PSSbins,statehists.EMGbins});
    statehists.(states{ss}).theta = hist3([BSmetrics.thratio(BSmetrics.instatetime.(states{ss}))...
        BSmetrics.EMG(BSmetrics.instatetime.(states{ss}))],...
        {statehists.thetabins,statehists.EMGbins});
    statehists.(states{ss}).PSSvtheta = hist3([BSmetrics.PSS(BSmetrics.instatetime.(states{ss})) BSmetrics.thratio(BSmetrics.instatetime.(states{ss}))],...
        {statehists.PSSbins,statehists.thetabins});
end
for ss = 1:length(states)
    statehists.(states{ss}).PSS = statehists.(states{ss}).PSS./sum(statehists.ALL.PSS(:));
    statehists.(states{ss}).theta = statehists.(states{ss}).theta./sum(statehists.ALL.theta(:));
    statehists.(states{ss}).PSSvtheta = statehists.(states{ss}).PSSvtheta./sum(statehists.ALL.PSSvtheta(:));

end

%%
figure
   
for tt = 1:length(celltypes)
subplot(4,3,tt*3-2)
    imagesc(ISIbyPSS.Xbins(1,:,1),ISIbyPSS.Ybins(1,:,1), ISIbyPSS.pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    ylabel({(celltypes{tt}),'ISI (log(s))'});
    %title((celltypes{tt}))
   % colorbar
    if tt ==1 
        caxis([0 0.02])
    elseif tt==2
         caxis([0 0.03])
    end
    xlim(ISIbyPSS.Xbins(1,[1 end],1))
end 


for tt = 1:length(celltypes)
subplot(4,3,tt*3-1)
    imagesc(ISIbytheta.Xbins(1,:,1),ISIbytheta.Ybins(1,:,1), ISIbytheta.pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
   % colorbar
    if tt ==1 
        caxis([0 0.02])
    elseif tt==2
         caxis([0 0.03])
    end
    xlim(ISIbytheta.Xbins(1,[1 end],1))
end 

subplot(6,3,10)
    hold on
    for ss = 1:3
    plot(BShist.bins,BShist.(states{ss}).PSS,'color',statecolors{ss})
    end
    xlabel('PSS')
    xlim(ISIbyPSS.Xbins(1,[1 end],1))
    
    
subplot(6,3,11)
    hold on
    for ss = [1 3]
    plot(BShist.bins,BShist.(states{ss}).thratio,'color',statecolors{ss})
    end
    xlabel('Theta Ratio')
    xlim(ISIbytheta.Xbins(1,[1 end],1))

subplot(3,3,3)
    hold on
    %imagesc(statehists.PSSbins,statehists.thetabins,statehists.(states{ss}).PSSvtheta')
    for ss = 1:3
        plotcolor = cat(3,statecolors{ss}(1).*ones(size(statehists.(states{ss}).PSSvtheta))',...
            statecolors{ss}(2).*ones(size(statehists.(states{ss}).PSSvtheta))',...
            statecolors{ss}(3).*ones(size(statehists.(states{ss}).PSSvtheta))');
        h = image(statehists.PSSbins,statehists.thetabins,plotcolor);
        set(h,'AlphaData',200*statehists.(states{ss}).PSSvtheta')

        plot(BSmetrics.PSS(BSmetrics.instatetime.(states{ss})),...
            BSmetrics.thratio(BSmetrics.instatetime.(states{ss})),'.','color',statecolors{ss},...
            'markersize',0.1)
    end
    xlabel('PSS');ylabel('Theta Ratio')
    axis tight
    
    
subplot(3,3,7)
    hold on
    %imagesc(statehists.PSSbins,statehists.EMGbins,statehists.PSS')
    for ss = 1:3
        plotcolor = cat(3,statecolors{ss}(1).*ones(size(statehists.(states{ss}).PSS))',...
            statecolors{ss}(2).*ones(size(statehists.(states{ss}).PSS))',...
            statecolors{ss}(3).*ones(size(statehists.(states{ss}).PSS))');
        h = image(statehists.PSSbins,statehists.EMGbins,plotcolor);
        set(h,'AlphaData',200*statehists.(states{ss}).PSS')
        
        plot(BSmetrics.PSS(BSmetrics.instatetime.(states{ss})),...
            BSmetrics.EMG(BSmetrics.instatetime.(states{ss})),'.','color',statecolors{ss},...
            'markersize',1)
    end
    xlabel('PSS');ylabel('EMG')
    axis tight
    
    
subplot(3,3,8)
    hold on
    %imagesc(statehists.thetabins,statehists.EMGbins,statehists.theta')
    for ss = [1 3]
        plotcolor = cat(3,statecolors{ss}(1).*ones(size(statehists.(states{ss}).theta))',...
            statecolors{ss}(2).*ones(size(statehists.(states{ss}).theta))',...
            statecolors{ss}(3).*ones(size(statehists.(states{ss}).theta))');
        h = image(statehists.thetabins,statehists.EMGbins,plotcolor);
        set(h,'AlphaData',200*statehists.(states{ss}).theta')
        
        plot(BSmetrics.thratio(BSmetrics.instatetime.(states{ss})),...
            BSmetrics.EMG(BSmetrics.instatetime.(states{ss})),'.','color',statecolors{ss},...
            'markersize',1)
    end
    xlabel('Theta Ratio');ylabel('EMG')
    axis tight
    
    
NiceSave('ISIbyStateVars',figfolder,baseName)



