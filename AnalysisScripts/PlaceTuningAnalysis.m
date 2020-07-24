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
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
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
nantimes = isnan(position.position.lin);
position.data = position.position.lin(~(nantimes));
%Possible here: remove drops for more than a... second?
minjump = 2;
numXbins = 80;
[ position.data ] = NanPadJumps( position.timestamps(~(nantimes)),position.data,minjump );
position.data = interp1(position.timestamps(~(nantimes)),position.data,position.timestamps);
[ISIbyPOS] = bz_ConditionalISI(spikes.times,position,...
    'ints',position.Epochs.MazeEpoch,...
    'showfig',false,'GammaFit',false,'numXbins',numXbins,'numISIbins',100,...
    'normtype','none','Xwin',[0 3],'Xbinoverlap',3);

%%
meanrate = nansum(ISIbyPOS.Dist.SpikeRate.*ISIbyPOS.Dist.pX,2);
meanrate_full = repmat(meanrate,1,numXbins,1);
MutInfo.SkaggsInf = squeeze(nansum(ISIbyPOS.Dist.SpikeRate.*log2(ISIbyPOS.Dist.SpikeRate./meanrate_full).*ISIbyPOS.Dist.pX,2));
MutInfo.SkaggsInf_sec = MutInfo.SkaggsInf .*squeeze(meanrate);


%%
[~,ISIbyPOS.fieldpeak] = max(ISIbyPOS.Dist.SpikeRate,[],2);
ISIbyPOS.fieldpeak = ISIbyPOS.Dist.Xbins(ISIbyPOS.fieldpeak);


%% Compare ISI MI vs Rate MI


binsizes = logspace(-2.5,1.5,25);
for bb = 1:length(binsizes)
    bz_Counter(bb,length(binsizes),'Bin')
    spkmat = bz_SpktToSpkmat(spikes.times,'dt',binsizes(bb),'binsize',binsizes(bb),...
        'win',position.Epochs.MazeEpoch,'units','rate');
    spkmat.pos = interp1(position.timestamps,position.data,spkmat.timestamps);

    for cc = 1:spikes.numcells
        MutInfo.Rate_BinCompare(cc,bb) = mutualinfo(spkmat.data(:,cc),spkmat.pos);
    end
end




%% Compare ISI MI vs Rate MI

MutInfo.ISI = squeeze(ISIbyPOS.MutInf);
spkmat = bz_SpktToSpkmat(spikes.times,'dt',0.3,'binsize',0.3,'win',position.Epochs.MazeEpoch,'units','rate');
spkmat.pos = interp1(position.timestamps,position.data,spkmat.timestamps);

for cc = 1:spikes.numcells
    MutInfo.Rate(cc) = mutualinfo(spkmat.data(:,cc),spkmat.pos);
    

end

%%
% figure
% subplot(2,2,1)
% plot(MutInfo.GSrate,log10(MutInfo.ISI),'.')
% 
% subplot(2,2,2)
% plot(MutInfo.GSweight,log10(MutInfo.ISI),'.')
%%
MIthresh = 0.05;
MIforbest = ISIbyPOS.MutInf;
MIforbest(MutInfo.Rate<MIthresh) = nan;
[~,bestcell] =max(MIforbest);
[~,sortMI_ISI] = sort(MutInfo.ISI);
[~,sortMutInfo.Rate] = sort(MutInfo.Rate);


%%

for cc = 1:spikes.numcells
    position_norm = position;
    position_norm.data = position.data-ISIbyPOS.fieldpeak(:,:,cc);
[ISIbyPOS_norm(cc)] = bz_ConditionalISI(spikes.times{cc},position_norm,...
    'ints',position_norm.Epochs.MazeEpoch,...
    'showfig',false,'GammaFit',false,'numXbins',80,'numISIbins',100,...
    'normtype','none','Xwin',[-1 1],'minX',10,'Xbinoverlap',3);
end



%MIthresh = 0.013;
ISIbyPOS_norm_mean = bz_CollapseStruct( ISIbyPOS_norm(squeeze(ISIbyPOS.MutInf)>MIthresh & MutInfo.Rate'>MIthresh),3,'mean',true);
%ISIbyPOS_norm_mean = bz_CollapseStruct( ISIbyPOS_norm(squeeze(ISIbyPOS.MutInf)>MIthresh),3,'mean',true);


ISIbyPOS_norm = bz_CollapseStruct( ISIbyPOS_norm,3,'justcat',true);

%% Get gamma mode
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');

MutInfo.GSrate = nan(1,spikes.numcells);
MutInfo.GSweight = nan(1,spikes.numcells);
for cc = 1:spikes.numcells
    cellUID(cc) = spikes.UID(cc);
    GFIDX = find(GammaFit.WAKEstate.cellstats.UID==cellUID(cc));
    if isempty(GFIDX)
        continue
    end
    cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
    MutInfo.GSrate(cc) = GammaFit.WAKEstate.sharedfit.GSlogrates(GFIDX);
    MutInfo.GSrate(cc) = GammaFit.WAKEstate.sharedfit.GSweights(GFIDX);
end

%%
MutInfo.cellclass = CellClass;
%%
figure
subplot(3,3,1)
    imagesc(ISIbyPOS_norm_mean.Dist.Xbins,[1 spikes.numcells],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortMI_ISI)))')
    ylabel('Sort by MIISI')
subplot(3,3,4)
    imagesc(ISIbyPOS_norm_mean.Dist.Xbins,[1 spikes.numcells],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortMutInfo.Rate)))')
    ylabel('Sort by MIRate')
    
subplot(3,3,7)
    imagesc(ISIbyPOS_norm_mean.Dist.Xbins,ISIbyPOS_norm_mean.Dist.Ybins,ISIbyPOS_norm_mean.Dist.pYX')
    hold on
    plot(ISIbyPOS_norm_mean.Dist.Xbins,-log10(ISIbyPOS_norm_mean.Dist.SpikeRate),'r')
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to PF Peak (m)')
    

subplot(3,3,9)
    scatter(log10(MutInfo.Rate),squeeze(log10(ISIbyPOS.MutInf)),5,(MutInfo.GSrate))
    xlabel('MI - rate');ylabel('MI - ISI')
    
subplot(3,3,6)
imagesc(log10(binsizes),[1 spikes.numcells],log10(MutInfo.Rate_BinCompare(sortMutInfo.Rate,:)))
LogScale('x',10)

subplot(3,3,3)
imagesc(log10(binsizes),[1 spikes.numcells],log10(MutInfo.Rate_BinCompare(sortMI_ISI,:)))
LogScale('x',10)

NiceSave('PlaceCoding',figfolder,baseName)


%%
% [excell] = bz_ConditionalISI(spikes.times{bestcell},position,...
%     'ints',position_norm.Epochs.MazeEpoch,...
%     'showfig',true,'GammaFit',false,'numXbins',10,'numISIbins',100,'minX',0);
