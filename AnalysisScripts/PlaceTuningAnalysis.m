function [ISIbyPOS_norm,MutInfo ] = PlaceTuningAnalysis(basePath,figfolder)
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
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
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


position = bz_LoadBehavior( basePath,'position' );

%%
nantimes = isnan(position.position.lin);
position.data = position.position.lin(~(nantimes));
%Possible here: remove drops for more than a... second?
minjump = 2;
[ position.data ] = NanPadJumps( position.timestamps(~(nantimes)),position.data,minjump );
position.data = interp1(position.timestamps(~(nantimes)),position.data,position.timestamps);
[ISIbyPOS] = bz_ConditionalISI(spikes.times,position,...
    'ints',position.Epochs.MazeEpoch,...
    'showfig',false,'GammaFit',false,'numXbins',25,'numISIbins',100,...
    'normtype','none','Xwin',[0 3]);

%%
[~,ISIbyPOS.fieldpeak] = max(ISIbyPOS.Dist.SpikeRate,[],2);
ISIbyPOS.fieldpeak = ISIbyPOS.Dist.Xbins(ISIbyPOS.fieldpeak);

%% Compare ISI MI vs Rate MI
spkmat = bz_SpktToSpkmat(spikes.times,'binsize',1,'win',position.Epochs.MazeEpoch,'units','rate');
spkmat.pos = interp1(position.timestamps,position.data,spkmat.timestamps);

MutInfo.ISI = squeeze(ISIbyPOS.MutInf);
for cc = 1:spikes.numcells
    MutInfo.Rate(cc) = mutualinfo(spkmat.data(:,cc),spkmat.pos);
end

%%

for cc = 1:spikes.numcells
    position_norm = position;
    position_norm.data = position.data-ISIbyPOS.fieldpeak(:,:,cc);
[ISIbyPOS_norm(cc)] = bz_ConditionalISI(spikes.times{cc},position_norm,...
    'ints',position_norm.Epochs.MazeEpoch,...
    'showfig',false,'GammaFit',false,'numXbins',30,'numISIbins',100,...
    'normtype','none','Xwin',[-1 1],'minX',10);
end



MIthresh = 0.013;
ISIbyPOS_norm_mean = bz_CollapseStruct( ISIbyPOS_norm(squeeze(ISIbyPOS.MutInf)>MIthresh & MutInfo.Rate'>MIthresh),3,'mean',true);
%ISIbyPOS_norm_mean = bz_CollapseStruct( ISIbyPOS_norm(squeeze(ISIbyPOS.MutInf)>MIthresh),3,'mean',true);


ISIbyPOS_norm = bz_CollapseStruct( ISIbyPOS_norm,3,'justcat',true);


%%
MIforbest = ISIbyPOS.MutInf;
MIforbest(MutInfo.Rate<MIthresh) = nan;
[~,bestcell] =max(MIforbest);
[~,sortMI_ISI] = sort(ISIbyPOS.MutInf);
[~,sortMutInfo.Rate] = sort(MutInfo.Rate);
%%
figure
subplot(2,2,1)
    imagesc(ISIbyPOS_norm_mean.Dist.Xbins,[1 spikes.numcells],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortMI_ISI)))')
    ylabel('Sort by MIISI')
subplot(2,2,2)
    imagesc(ISIbyPOS_norm_mean.Dist.Xbins,[1 spikes.numcells],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortMutInfo.Rate)))')
    ylabel('Sort by MIRate')
    
subplot(2,2,3)
    imagesc(ISIbyPOS_norm_mean.Dist.Xbins,ISIbyPOS_norm_mean.Dist.Ybins,ISIbyPOS_norm_mean.Dist.pYX')
    hold on
    plot(ISIbyPOS_norm_mean.Dist.Xbins,-log10(ISIbyPOS_norm_mean.Dist.SpikeRate),'r')
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to PF Peak (m)')
    
subplot(2,2,4)
    plot(log10(MutInfo.Rate),squeeze(log10(ISIbyPOS.MutInf)),'.')
    xlabel('MI - rate');ylabel('MI - ISI')
    
NiceSave('PlaceCoding',figfolder,baseName)


%%
% [excell] = bz_ConditionalISI(spikes.times{bestcell},position,...
%     'ints',position_norm.Epochs.MazeEpoch,...
%     'showfig',true,'GammaFit',false,'numXbins',10,'numISIbins',100,'minX',0);
