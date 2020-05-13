function [TH_ISIstats,ISIbythetaphase,ISIbytheta] = ThetaISIAnalysis(basePath,figfolder)
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
%reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
% basePath = '/Users/dl2820/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
% %basePath = pwd;
% %basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
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
BSmetrics.thratio = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio;

%%
for ss = 1:length(states)
   BSmetrics.instatetime.(states{ss}) = InIntervals(BSmetrics.timestamps,SleepState.ints.(states{ss}));
   

end


%% BS metrics histogram by states
BShist.bins = linspace(0,1,50);
for ss = 1:length(states)
   BShist.(states{ss}).thratio = hist(BSmetrics.thratio(BSmetrics.instatetime.(states{ss})),BShist.bins);
   BShist.(states{ss}).thratio = BShist.(states{ss}).thratio./length(BSmetrics.timestamps);
end


%% Get SleepScoreParams at each spike time

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

end

%%
ThetaIDX.threshs.hitheta = 0.65;
ThetaIDX.threshs.lotheta = 0.35;
ThetaIDX.timestamps = BSmetrics.timestamps;
ThetaIDX.states = zeros(size(ThetaIDX.timestamps));
ThetaIDX.states(BSmetrics.thratio>ThetaIDX.threshs.hitheta & BSmetrics.instatetime.WAKEstate) = 1;
ThetaIDX.states(BSmetrics.thratio<ThetaIDX.threshs.lotheta & BSmetrics.instatetime.WAKEstate) = 2;

ThetaIDX.statenames = {'hiTheta','loTheta'};
ThetaIDX.timestamps = BSmetrics.timestamps;

ThetaINT = bz_IDXtoINT(ThetaIDX);
%%
TH_ISIstats = bz_ISIStats(spikes,'ints',ThetaINT,'showfig',true,'cellclass',CellClass.label);
TH_ISIstats = rmfield(TH_ISIstats,'allspikes');

%% Get the Theta LFP

 lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID;
 downsamplefactor = 5;
 lfp = bz_GetLFP(lfpchan,...
     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%%
thetalfp = bz_Filter(lfp,'passband',[6 10]);
thetalfp.amp = NormToInt((thetalfp.amp),'mean',SleepState.ints.WAKEstate);

ISIStats.allspikes.thetaphase =cellfun(@(X) ...
    interp1(thetalfp.timestamps,thetalfp.phase,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);

%% COnditional Hists
%[-0.6 1.5] (median normalize)

%[-1.5 1.75]

[ ISIbythetaphase ] = cellfun(@(X,Y,Z,W,T) ConditionalHist( [Z(W&T>ThetaIDX.threshs.hitheta);Z(W&T>ThetaIDX.threshs.hitheta)],...
    log10([X(W&T>ThetaIDX.threshs.hitheta);Y(W&T>ThetaIDX.threshs.hitheta)]),...
    'Xbounds',[-pi pi],'numXbins',50,'Ybounds',[-3 2],'numYbins',125,'minX',100),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.thetaphase,ISIStats.allspikes.instate.WAKEstate,...
    ISIStats.allspikes.thetarat,...
    'UniformOutput',false);
ISIbythetaphase = cat(1,ISIbythetaphase{:});
ISIbythetaphase = CollapseStruct( ISIbythetaphase,3);


[ ISIbytheta ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(W);Z(W)],log10([X(W);Y(W)]),...
    'Xbounds',[0 1],'numXbins',30,'Ybounds',[-3 2],'numYbins',125,'minX',100),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.thetarat,ISIStats.allspikes.instate.WAKEstate,...
    'UniformOutput',false);
ISIbytheta = cat(1,ISIbytheta{:});
ISIbytheta = CollapseStruct( ISIbytheta,3);



for tt = 1:length(celltypes)
    ISIbytheta.pop.(celltypes{tt}) = nanmean(ISIbytheta.pYX(:,:,CellClass.(celltypes{tt})),3);
    ISIbythetaphase.pop.(celltypes{tt}) = nanmean(ISIbythetaphase.pYX(:,:,CellClass.(celltypes{tt})),3);

    ISIbythetaphase.celltypeidx.(celltypes{tt}) = CellClass.(celltypes{tt});
    ISIbytheta.celltypeidx.(celltypes{tt}) = CellClass.(celltypes{tt});
end
    if length(celltypes)==1
        ISIbytheta.celltypeidx.pI = false(size(CellClass.pE));
    end
    

%%
phasex = linspace(-pi,3*pi,100);
figure
for tt = 1:length(celltypes)
subplot(4,3,tt*3-2)
    imagesc(ISIbythetaphase.Xbins(1,:,1),ISIbythetaphase.Ybins(1,:,1), ISIbythetaphase.pop.(celltypes{tt})')
    hold on
    imagesc(ISIbythetaphase.Xbins(1,:,1)+2*pi,ISIbythetaphase.Ybins(1,:,1), ISIbythetaphase.pop.(celltypes{tt})')
    %axis xy
    plot(phasex,-cos(phasex),'k')
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    %axis xy
    %xlim([-pi 3*pi])
    LogScale('y',10,'exp',true)
    ylabel('ISI (s)');xlabel('Theta Phase')
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    xlim([-pi 3*pi])
end 
NiceSave('TH_Phase',figfolder,baseName)

%%
THlabels = {'hiThetastate','loThetastate'};
%%
figure
subplot(3,3,7)
    hold on
    for tt = 1:length(celltypes)
        plot(log10(TH_ISIstats.summstats.hiThetastate.meanrate(CellClass.(celltypes{tt}))),...
            log10(TH_ISIstats.summstats.loThetastate.meanrate(CellClass.(celltypes{tt}))),...
            '.','color',cellcolor{tt})
    end
    hold on
    UnityLine
    xlabel('hiTheta Rate');ylabel('loTheta Rate')
    

for tt = 1:length(celltypes) 
    
    subplot(4,4,tt+6)
        plot(TH_ISIstats.ISIhist.logbins,TH_ISIstats.meandists.hiThetastate.(celltypes{tt}).ISIdist,'k')
        hold on
        plot(TH_ISIstats.ISIhist.logbins,TH_ISIstats.meandists.loThetastate.(celltypes{tt}).ISIdist,'r')
        LogScale('x',10,'exp',true)
        title(celltypes{tt})
        box off 
        
    for ss = 1:length(THlabels)
        subplot(4,4,10+tt+(ss-1)*4)
            imagesc(TH_ISIstats.ISIhist.logbins,TH_ISIstats.ISIhist.logbins,...
            TH_ISIstats.meandists.(THlabels{ss}).(celltypes{tt}).Return)
            LogScale('xy',10,'exp',true)
            axis xy
    end
end


for tt = 1:length(celltypes)
subplot(4,3,tt*3-2)
    imagesc(ISIbytheta.Xbins(1,:,1),ISIbytheta.Ybins(1,:,1), ISIbytheta.pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('y',10,'exp',true)
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
    for ss = [1]
    plot(BShist.bins,BShist.(states{ss}).thratio,'color',statecolors{ss})
    end
    xlabel('Theta Ratio')
    xlim(ISIbytheta.Xbins(1,[1 end],1))

NiceSave('TH_ISIstats',figfolder,baseName)

end
