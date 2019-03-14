function [ ] = Analysis20190304(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: Are slow rate fluctuations (Okun et al) changes in ground state
%rate, changes in activated state?
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
% 
% lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID;
% downsamplefactor = 1;
% lfp = bz_GetLFP(lfpchan,...
%     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
% %Noralize the LFP
% lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

%% Calculate long time scale rate matrix
dt = 6;
spkmat = bz_SpktToSpkmat(spikes,'binsize',60,'overlap',binsize./dt);
spkmat.data = spkmat.data./spkmat.binsize;


%%
state = states{2};
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);
spkmat.instate = InIntervals(spkmat.timestamps,SleepState.ints.(state));
%%
spkmat.data_meannorm = bsxfun(@(X,Y) X./Y,spkmat.data,mean(spkmat.data(spkmat.instate,:),1));

%%
for cc = 1:spikes.numcells
    ISIStats.allspikes.longrate{cc} = interp1(spkmat.timestamps,spkmat.data(:,cc),...
        ISIStats.allspikes.times{cc},'nearest');
    ISIStats.allspikes.normrate{cc} = interp1(spkmat.timestamps,spkmat.data_meannorm(:,cc),...
        ISIStats.allspikes.times{cc},'nearest');
end

%%
ISIStats.allspikes.normISI = cellfun(@(X,Y) X(Y)./mean(X(Y)),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.instate,'UniformOutput',false);


ISIStats.allspikes.logISI = cellfun(@(X) log10(X),...
    ISIStats.allspikes.normISI,'UniformOutput',false);
ISIStats.allspikes.logISInp1 = cellfun(@(X) [X(2:end);nan],...
    ISIStats.allspikes.logISI,'UniformOutput',false);

ISIStats.allspikes.bothISI = cellfun(@(X,Y) [X;Y],...
    ISIStats.allspikes.logISI,ISIStats.allspikes.logISInp1,'UniformOutput',false);
ISIStats.allspikes.bothlongrate = cellfun(@(X,Y) log10([X(Y);X(Y)]),...
    ISIStats.allspikes.longrate,ISIStats.allspikes.instate,'UniformOutput',false);
ISIStats.allspikes.bothlongrate_norm = cellfun(@(X,Y) log2([X(Y);X(Y)]),...
    ISIStats.allspikes.normrate,ISIStats.allspikes.instate,'UniformOutput',false);
%%
 [ISIbyRate] = ConditionalHist(ISIStats.allspikes.bothlongrate_norm,ISIStats.allspikes.bothISI,...
     'numYbins',50,'numXbins',15,'Xbounds',[-1.1 1.1],'Ybounds',[-3 1.5],'minX',50) ;
 
 for tt = 1:length(celltypes)
     ISIbyRate.classmean.(celltypes{tt}) = nanmean(ISIbyRate.pYX(:,:,CellClass.(celltypes{tt})),3);
     ISIbyRate.classmean.ratehist.(celltypes{tt}) = nanmean(ISIbyRate.Xhist(:,:,CellClass.(celltypes{tt})),3);

 end
 
 %%
 excell = randsample(spikes.numcells,1);
 figure
  for tt = 1:length(celltypes)
    subplot(3,3,tt)
        imagesc(ISIbyRate.Xbins,ISIbyRate.Ybins,ISIbyRate.classmean.(celltypes{tt})')
        hold on
        plot(ISIbyRate.Xbins([1 end]),[0 0],'w--')
        plot([0 0],ISIbyRate.Ybins([1 end]),'w--')
        title(celltypes{tt})
        %colorbar
       % plot(ISIbyRate.Xbins,log10(1./10.^ISIbyRate.Xbins),'r')
        axis xy
        LogScale('y',10)
        %LogScale('x',2)
        xlabel('Norm Rate (mean^-^1)');ylabel('Norm ISI (mean^-^1)')
        
    subplot(3,3,tt+3)
        plot(ISIbyRate.Xbins,ISIbyRate.classmean.ratehist.(celltypes{tt}))
        xlabel('Norm Rate')
        %axis 
        box off
        xlim(ISIbyRate.Xbins([1 end]))
  end
  
  
    subplot(3,3,7)
        imagesc(ISIbyRate.Xbins,ISIbyRate.Ybins,ISIbyRate.pYX(:,:,excell)')
        axis xy
        LogScale('y',10)
        
NiceSave(['ISIbyLongRate_',state],figfolder,baseName,'includeDate',true)
