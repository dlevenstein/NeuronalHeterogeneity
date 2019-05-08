function [ISIbyRate,CV2byRate,ISIbyCspkRate,CV2byCspkRate ] = SpikeStatsbyRate(basePath,figfolder)
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
statecolors = {'k','b','r',[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};


%% Get Rate - constant time and ISI bins

%COnstant Time
spkmat = bz_SpktToSpkmat(spikes,'binsize',4,'dt',0.5);
spkmat.data = spkmat.data./spkmat.binsize;

%Constant Spike COunt
nspkintervals = 5;

cspkrate.meanISI = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.ISIs,'UniformOutput',false);
cspkrate.rate = cellfun(@(X) 1./X ,cspkrate.meanISI,'UniformOutput',false);
cspkrate.bindur = cellfun(@(X) movsum(X,nspkintervals),ISIStats.allspikes.ISIs,'UniformOutput',false);
cspkrate.stdISI = cellfun(@(X) movstd(X,nspkintervals),ISIStats.allspikes.ISIs,'UniformOutput',false);
cspkrate.CVISI = cellfun(@(X,Y) X./Y,cspkrate.stdISI,cspkrate.meanISI,'UniformOutput',false);
cspkrate.meanCV2 = cellfun(@(X) movmean(X,nspkintervals),ISIStats.allspikes.CV2,'UniformOutput',false);
cspkrate.times = cellfun(@(X) movmean(X,nspkintervals+1),ISIStats.allspikes.times,'UniformOutput',false);

%% Get Rate stuff at each spike
for cc = 1:spikes.numcells
    ISIStats.allspikes.rate{cc} = interp1(spkmat.timestamps,spkmat.data(:,cc),...
        ISIStats.allspikes.times{cc},'nearest');
    ISIStats.allspikes.cspkrate{cc} = interp1(cspkrate.times{cc},cspkrate.rate{cc},...
        ISIStats.allspikes.times{cc},'nearest');
end

ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end); nan],...
    ISIStats.allspikes.ISIs,'UniformOutput',false);
%% Restrict to state
for ss = 1:length(states)
state = states{ss};
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);
cspkrate.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    cspkrate.times,'UniformOutput',false);
spkmat.instats = InIntervals(spkmat.timestamps,SleepState.ints.(state));

%% Conditional distirbutions: ISI, CV2
maxrate = 40;
[ ISIbyRate.(state) ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z(W);Z(W)],log10([X(W);Y(W)]),...
    'Xbounds',[0 maxrate],'numXbins',(maxrate.*spkmat.binsize)+1,'Ybounds',[-3 2],'numYbins',125,'minX',50),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.rate,ISIStats.allspikes.instate,...
    'UniformOutput',false);
ISIbyRate.(state) = cat(1,ISIbyRate.(state){:});
ISIbyRate.(state) = CollapseStruct( ISIbyRate.(state),3);

[ CV2byRate.(state) ] = cellfun(@(X,Z,W) ConditionalHist((Z(W)),(X(W)),...
    'Xbounds',[0 maxrate],'numXbins',(maxrate.*spkmat.binsize)+1,'Ybounds',[0 2],'numYbins',75,'minX',50),...
    ISIStats.allspikes.CV2,...
    ISIStats.allspikes.rate,ISIStats.allspikes.instate,...
    'UniformOutput',false);
CV2byRate.(state) = cat(1,CV2byRate.(state){:});
CV2byRate.(state) = CollapseStruct( CV2byRate.(state),3);


[ ISIbyCspkRate.(state) ] = cellfun(@(X,Y,Z,W) ConditionalHist( log10([Z(W);Z(W)]),log10([X(W);Y(W)]),...
    'Xbounds',[-2 2.5],'numXbins',60,'Ybounds',[-3 2],'numYbins',125,'minX',40),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
    ISIStats.allspikes.cspkrate,ISIStats.allspikes.instate,...
    'UniformOutput',false);
ISIbyCspkRate.(state) = cat(1,ISIbyCspkRate.(state){:});
ISIbyCspkRate.(state) = CollapseStruct( ISIbyCspkRate.(state),3);

% [ normISIbyCspkRate.(state) ] = cellfun(@(X,Y,Z,W) ConditionalHist( log10([Z(W);Z(W)]),log10([X(W);Y(W)]./mean(X(W))),...
%     'Xbounds',[-2 2.5],'numXbins',60,'Ybounds',[-3 2],'numYbins',125,'minX',40),...
%     ISIStats.allspikes.ISIs,ISIStats.allspikes.ISInp1,...
%     ISIStats.allspikes.cspkrate,ISIStats.allspikes.instate,...
%     'UniformOutput',false);
% normISIbyCspkRate.(state) = cat(1,normISIbyCspkRate.(state){:});
% normISIbyCspkRate.(state) = CollapseStruct( normISIbyCspkRate.(state),3);


[ CV2byCspkRate.(state) ] = cellfun(@(X,Z,W) ConditionalHist(log10(Z(W)),(X(W)),...
    'Xbounds',[-2 2.5],'numXbins',60,'Ybounds',[0 2],'numYbins',75,'minX',40),...
    ISIStats.allspikes.CV2,...
    ISIStats.allspikes.cspkrate,ISIStats.allspikes.instate,...
    'UniformOutput',false);
CV2byCspkRate.(state) = cat(1,CV2byCspkRate.(state){:});
CV2byCspkRate.(state) = CollapseStruct( CV2byCspkRate.(state),3);

for tt = 1:length(celltypes)
    ISIbyRate.(state).pop.(celltypes{tt}) = nanmean(ISIbyRate.(state).XYprob(:,:,CellClass.(celltypes{tt})),3);
    CV2byRate.(state).pop.(celltypes{tt}) = nanmean(CV2byRate.(state).XYprob(:,:,CellClass.(celltypes{tt})),3);
    ISIbyCspkRate.(state).pop.(celltypes{tt}) = nanmean(ISIbyCspkRate.(state).XYprob(:,:,CellClass.(celltypes{tt})),3);
    CV2byCspkRate.(state).pop.(celltypes{tt}) = nanmean(CV2byCspkRate.(state).XYprob(:,:,CellClass.(celltypes{tt})),3);
end
end
%% Figure: rate

figure
for ss = 1:3
    state = states{ss};
for tt = 1:length(celltypes)
subplot(6,6,ss+(tt-1)*3)
    imagesc(ISIbyRate.(state).Xbins(1,:,1),ISIbyRate.(state).Ybins(1,:,1), ISIbyRate.(state).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('y',10)
    ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Rate')
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    
end

for tt = 1:length(celltypes)
subplot(6,6,ss+(tt-1)*3+6)
    imagesc(CV2byRate.(state).Xbins(1,:,1),CV2byRate.(state).Ybins(1,:,1), CV2byRate.(state).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    ylabel({(celltypes{tt}),'CV2'});xlabel('Rate')
    %title((celltypes{tt}))
    %colorbar
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    
end

for tt = 1:length(celltypes)
subplot(6,6,ss+(tt-1)*3+12)
    imagesc(CV2byRate.(state).Xbins(1,:,1),[0 1], squeeze(CV2byRate.(state).pX(:,:,ISIStats.sorts.(state).(['rate',celltypes{tt}])))')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    ylabel('Cell');xlabel('Rate')
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.03])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    
end 
end
NiceSave('ISIbyRate',figfolder,baseName,'includeDate',true)
%% Figure: CSPKRWATE
figure
for ss = 1:3
    state = states{ss};
for tt = 1:length(celltypes)
subplot(6,6,ss+(tt-1)*3)
    imagesc(ISIbyCspkRate.(state).Xbins(1,:,1),ISIbyCspkRate.(state).Ybins(1,:,1), ISIbyCspkRate.(state).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('y',10)
    if ss==1
    ylabel('ISI (s)');
    end
    %xlabel('Rate')
    %title((celltypes{tt}))
    %colorbar
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    
end 


for tt = 1:length(celltypes)
subplot(6,6,ss+(tt-1)*3+6)
    imagesc(CV2byCspkRate.(state).Xbins(1,:,1),CV2byCspkRate.(state).Ybins(1,:,1), CV2byCspkRate.(state).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    if ss==1
    ylabel('CV2');end
    %title((celltypes{tt}))
	%colorbar
    if tt ==1 
        caxis([0 0.6e-3])
    elseif tt==2
         caxis([0 0.003])
    end
    
end 

for tt = 1:length(celltypes)
subplot(6,6,ss+(tt-1)*3+12)
    imagesc(CV2byCspkRate.(state).Xbins(1,:,1),[1 sum(CellClass.(celltypes{tt}))], squeeze(CV2byCspkRate.(state).pX(:,:,ISIStats.sorts.(state).(['rate',celltypes{tt}])))')
    hold on
    plot(log10(ISIStats.summstats.(state).meanrate(ISIStats.sorts.(state).(['rate',celltypes{tt}]))),[1:sum(CellClass.(celltypes{tt}))],'w')
    axis xy
    %LogScale('y',10)
    if ss==1
    ylabel('Cell');
    end
    xlabel('Rate')
    %title((celltypes{tt}))
    %colorbar
%     if tt ==1 
%         caxis([0 0.03])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    
end 
end
NiceSave('ISIbyCspkRate',figfolder,baseName,'includeDate',true)

