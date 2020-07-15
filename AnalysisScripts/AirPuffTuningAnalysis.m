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
basePath = '/Users/dl2820/Dropbox/research/Datasets/Rat08-20130713';
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


%position = bz_LoadBehavior( basePath,'position' );

%%
puff = LoadEvents([basePath,'/Rat08-20130713.puf.evt'])
load([basePath,'/runintervals.mat'])

%%
airpuff.dt = 0.01;
airpuff.timestamps = runintervals(2,1):airpuff.dt:runintervals(2,2);%Make timestamps inside the laps...
for tt = 1:length(airpuff.timestamps)
    [~,whichAP] = min(abs(airpuff.timestamps(tt)-puff.time));
    airpuff.data(tt) = (airpuff.timestamps(tt)-puff.time(whichAP));
end



%%
[ISIbyAP] = bz_ConditionalISI(spikes.times,airpuff,...
    'showfig',false,'GammaFit',false,'numXbins',20,'numISIbins',100,...
    'normtype','none','Xwin',[-10 10]);


%%
blacells = strcmp(cat(1,spikes.region{:}),'bla') & CellClass.pE';
popspikes = cat(1,spikes.times{blacells});

[popISIbyAP] = bz_ConditionalISI(popspikes,airpuff,...
    'showfig',true,'GammaFit',false,'numXbins',20,'numISIbins',100,...
    'normtype','none','Xwin',[-10 10],'figfolder',figfolder,'figname','popISI',...
    'basePath',basePath);

%%
ISIbyAP.blaMean = nanmean(ISIbyAP.Dist.pYX(:,:,blacells),3);
%%
figure
imagesc(ISIbyAP.Dist.Xbins(1,:,1),(ISIbyAP.Dist.Ybins(1,:,1)),ISIbyAP.blaMean')
LogScale('y',10)
bz_AddRightRateAxis
%%
MutInfo.ISI = squeeze(ISIbyAP.MutInf);
spkmat = bz_SpktToSpkmat(spikes.times,'dt',0.3,'binsize',0.3,'win',runintervals(2,:),'units','rate');
spkmat.AP = interp1(airpuff.timestamps,airpuff.data,spkmat.timestamps);

for cc = 1:spikes.numcells
    MutInfo.Rate(cc) = mutualinfo(spkmat.data(:,cc),spkmat.AP);
end

%%
[ISIbyAP_foract] = bz_ConditionalISI(spikes.times,airpuff,...
    'showfig',false,'GammaFit',false,'numXbins',4,'numISIbins',100,...
    'normtype','none','Xwin',[-10 10]);

APactivation_pre = (ISIbyAP_foract.Dist.SpikeRate(:,2,:)-ISIbyAP_foract.Dist.SpikeRate(:,1,:))./(ISIbyAP_foract.Dist.SpikeRate(:,2,:)+ISIbyAP_foract.Dist.SpikeRate(:,1,:));
APactivation_post = (ISIbyAP_foract.Dist.SpikeRate(:,3,:)-ISIbyAP_foract.Dist.SpikeRate(:,4,:))./(ISIbyAP_foract.Dist.SpikeRate(:,3,:)+ISIbyAP_foract.Dist.SpikeRate(:,4,:));

%%
%[~,sortMI]=sort(MutInfo.ISI);
[~,sortMI]=sort(MutInfo.Rate);
[~,sortPre]=sort(APactivation_pre);
[~,sortPost]=sort(APactivation_post);

%%
APthresh = APactivation_pre>0 & APactivation_post>0;

ISIbyAP.blaMeanAP = nanmean(ISIbyAP.Dist.pYX(:,:,squeeze(blacells)&squeeze(APthresh)),3);
%%
figure
subplot(2,2,1)
hist(squeeze(ISIbyAP.MutInf))

subplot(2,2,2)
imagesc(log10(squeeze(ISIbyAP.Dist.SpikeRate(:,:,sortPost)))');

subplot(2,2,3)
plot(squeeze(APactivation_pre),squeeze(APactivation_post),'k.')
hold on
plot(squeeze(APactivation_pre(blacells)),squeeze(APactivation_post(blacells)),'r.')

subplot(2,2,4)
imagesc(ISIbyAP.Dist.Xbins(1,:,1),(ISIbyAP.Dist.Ybins(1,:,1)),ISIbyAP.blaMeanAP')
LogScale('y',10)
bz_AddRightRateAxis
xlabel('time relative to airpuff (s)')

NiceSave('APCoding',figfolder,baseName)

%%
%load('airpuff.mat')
%load('clearunpos.mat')

%%

figure
plot(cleanrunpos(:,2),cleanrunpos(:,3),'.')
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
    'showfig',false,'GammaFit',false,'numXbins',30,'numISIbins',100,...
    'normtype','none','Xwin',[-1 1],'minX',10);
end



%MIthresh = 0.013;
ISIbyPOS_norm_mean = bz_CollapseStruct( ISIbyPOS_norm(squeeze(ISIbyPOS.MutInf)>MIthresh & MutInfo.Rate'>MIthresh),3,'mean',true);
%ISIbyPOS_norm_mean = bz_CollapseStruct( ISIbyPOS_norm(squeeze(ISIbyPOS.MutInf)>MIthresh),3,'mean',true);


ISIbyPOS_norm = bz_CollapseStruct( ISIbyPOS_norm,3,'justcat',true);

%% Get gamma mode
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');

MutInfo.GSrate = nan(1,spikes.numcells);
for cc = 1:spikes.numcells
    cellUID(cc) = spikes.UID(cc);
    GFIDX = find(GammaFit.WAKEstate.cellstats.UID==cellUID(cc));
    if isempty(GFIDX)
        continue
    end
    cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
    MutInfo.GSrate(cc) = GammaFit.WAKEstate.sharedfit.GSlogrates(GFIDX);
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
