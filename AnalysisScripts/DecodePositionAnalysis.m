function [ ] = DecodePositionAnalysis(basePath,figfolder)
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
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse12-120807';
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

%[celltypes,~,typeidx] = unique(CellClass.label);
%cellcolor = {'k','r'};

%%
try %(CA1)
    position = bz_LoadBehavior( basePath,'position' );

    nantimes = isnan(position.position.lin);
    position.data = position.position.lin(~(nantimes));
    %Possible here: remove drops for more than a... second?
    minjump = 2;
    [ position.data ] = NanPadJumps( position.timestamps(~(nantimes)),position.data,minjump );
    position.data = interp1(position.timestamps(~(nantimes)),position.data,position.timestamps,'nearest');

    codingEpochs = position.Epochs.MazeEpoch;
    
    region = 'CA1';
catch
    position.samplingRate = 39.06; %Hz
    hdfilename = fullfile(basePath,[baseName,'.ang']);
    position.data = importdata(hdfilename);
    position.data(position.data==-1)=nan;
    position.timestamps = [1:length(position.data)]'./position.samplingRate; 
    
    %Use only long wakes for decoding
    wakedur = diff(SleepState.ints.WAKEstate,[],2);
    codingEpochs =  SleepState.ints.WAKEstate(wakedur>100,:);
    
    
    %restrict to wake
    inwake = InIntervals(position.timestamps,codingEpochs);
    position.data = position.data(inwake);
    position.timestamps = position.timestamps(inwake);
    
    region = 'THAL';
end

%%
switch region
    case 'THAL'
        xbins = 25;
    case 'CA1'
        xbins = 50;
end


[ISIbyPOS] = bz_ConditionalISI(spikes.times,position,...
    'ints',codingEpochs,...
    'showfig',false,'GammaFit',false,'numXbins',xbins,'numISIbins',100,...
    'normtype','none','Xwin',[min(position.data) max(position.data)]);

[~,ISIbyPOS.fieldpeak] = max(ISIbyPOS.Dist.SpikeRate,[],2);
ISIbyPOS.fieldpeak = ISIbyPOS.Dist.Xbins(ISIbyPOS.fieldpeak);
[~,sortfieldpeak] = sort(ISIbyPOS.fieldpeak);



%% Compare ISI MI vs Rate MI
MutInfo.ISI = squeeze(ISIbyPOS.MutInf);
spkmat = bz_SpktToSpkmat(spikes.times,'dt',0.1,'binsize',1,'units','rate');
spkmat.pos = interp1(position.timestamps,position.data,spkmat.timestamps,'nearest');
spkmat.CodeEpochs = InIntervals(spkmat.timestamps,codingEpochs);

for cc = 1:spikes.numcells
    MutInfo.Rate(cc) = mutualinfo(spkmat.data(spkmat.CodeEpochs,cc),spkmat.pos(spkmat.CodeEpochs));
end


%%

%% Testing values for MI threshold and bin size
switch region
    case 'THAL'
        MIthresh = 0.1;
    case 'CA1'
        MIthresh = 0.03;
end

binsizes = logspace(-1,0.5,20);
for bb = 1:length(binsizes)
    bz_Counter(bb,length(binsizes),'bin')
binsize = binsizes(bb); %s (1)
dt = binsize./4;
spkmat = bz_SpktToSpkmat(spikes.times,'dt',dt,'binsize',binsize,'units','counts');
spkmat.CodeEpochs = InIntervals(spkmat.timestamps,codingEpochs);
spkmat.data = spkmat.data(spkmat.CodeEpochs,:);
spkmat.timestamps = spkmat.timestamps(spkmat.CodeEpochs);

rateMap = squeeze(ISIbyPOS.Dist.SpikeRate);
keep = find(sum(rateMap)>0 & CellClass.pE & MutInfo.Rate>MIthresh);
[Pr, prMax] = placeBayes(spkmat.data(:,keep), rateMap(:,keep)', binsize);
position.decoded = interp1(spkmat.timestamps,ISIbyPOS.Dist.Xbins(1,prMax,1),position.timestamps,'nearest');
position.decoderror = min([abs(position.decoded-position.data),abs(abs(position.decoded-position.data)-max(position.data))],[],2);
meanDecodError.binsize(bb) = nanmean(position.decoderror);
end


switch region
    case 'THAL'
        binsize = 0.32; %s (1)
    case 'CA1'
        binsize = 0.6; %s (1)
end
dt = binsize./4;
spkmat = bz_SpktToSpkmat(spikes.times,'dt',dt,'binsize',binsize,'units','counts');
spkmat.CodeEpochs = InIntervals(spkmat.timestamps,codingEpochs);
spkmat.data = spkmat.data(spkmat.CodeEpochs,:);
spkmat.timestamps = spkmat.timestamps(spkmat.CodeEpochs);

MIthreshs = logspace(-4,-0.5,20);
for bb = 1:length(binsizes)
    bz_Counter(bb,length(binsizes),'bin')

MIthresh = MIthreshs(bb);
rateMap = squeeze(ISIbyPOS.Dist.SpikeRate);
keep = find(sum(rateMap)>0 & CellClass.pE & MutInfo.Rate>MIthresh);
[Pr, prMax] = placeBayes(spkmat.data(:,keep), rateMap(:,keep)', binsize);
position.decoded = interp1(spkmat.timestamps,ISIbyPOS.Dist.Xbins(1,prMax,1),position.timestamps,'nearest');
position.decoderror = min([abs(position.decoded-position.data),abs(abs(position.decoded-position.data)-max(position.data))],[],2);
meanDecodError.MIthresh(bb) = nanmean(position.decoderror);
end

%%
figure
subplot(2,2,1)
plot(log10(binsizes),meanDecodError.binsize)
LogScale('x',10)
ylabel('Mean Decoder Error')
xlabel('Bin Size')
subplot(2,2,2)
plot(log10(MIthreshs),meanDecodError.MIthresh)
LogScale('x',10)
ylabel('Mean Decoder Error')
xlabel('MI Thresh')
NiceSave('Decoding_Optimize',figfolder,baseName)

subplot(2,2,3)
imagesc(squeeze(ISIbyPOS.Dist.SpikeRate(:,:,sortfieldpeak)))
% subplot(2,2,2)
% plot(log10(MutInfo.Rate),log10(MutInfo.ISI),'.')
subplot(2,2,4)
plot(squeeze(ISIbyPOS.fieldpeak),log10(MutInfo.Rate),'.')
NiceSave('Decoding_Optimize',figfolder,baseName)

%% Decode position using best bin size and MI threshhold

switch region
    case 'THAL'
        binsize = 0.32; %s (1)
        MIthresh = 0.1;
    case 'CA1'
        binsize = 0.6; %s (1)
        MIthresh = 0.03;
end

dt = binsize./4;

spkmat = bz_SpktToSpkmat(spikes.times,'dt',dt,'binsize',binsize,'units','counts');
%spkmat.CodeEpochs = InIntervals(spkmat.timestamps,codingEpochs);
%spkmat.data = spkmat.data(spkmat.CodeEpochs,:);
%spkmat.timestamps = spkmat.timestamps(spkmat.CodeEpochs);

rateMap = squeeze(ISIbyPOS.Dist.SpikeRate);
keep = find(sum(rateMap)>0 & CellClass.pE & MutInfo.Rate>MIthresh);
[Pr, prMax] = placeBayes(spkmat.data(:,keep), rateMap(:,keep)', binsize);
position.decoded = interp1(spkmat.timestamps,ISIbyPOS.Dist.Xbins(1,prMax,1),position.timestamps,'nearest');
position.decoderror = min([abs(position.decoded-position.data),abs(abs(position.decoded-position.data)-max(position.data))],[],2);
meanDecodError.offish = nanmean(position.decoderror);


%spkmat.position = interp1(position.timestamps,position.data,spkmat.timestamps);
%spkmat.position_decoded = interp1(position.timestamps,position.decoded,spkmat.timestamps);
%spkmat.position_error = interp1(position.timestamps,position.decoderror,spkmat.timestamps);
%
sortfieldpeak_keep = sortfieldpeak(ismember(sortfieldpeak,keep));


%%
% binsize_NREM = 0.1;
% spkmat_NREM = bz_SpktToSpkmat(spikes.times,'dt',0.1,'binsize',binsize,'units','counts');
% spkmat
% [Pr_NREM, prMax] = placeBayes(spkmat.data(:,keep), rateMap(:,keep)', binsize);

%%
viewwin = bz_RandomWindowInIntervals(codingEpochs,100,1);
figure
subplot(3,1,1)
imagesc(spkmat.timestamps,spikes.numcells,log10(spkmat.data(:,sortfieldpeak_keep))')
xlim(viewwin)
ylabel('Decoder cells')

subplot(3,1,2)
imagesc(spkmat.timestamps,ISIbyPOS.Dist.Xbins(1,:,1),Pr')
hold on
plot(position.timestamps,position.data,'r')
plot(position.timestamps,position.decoded,'o')
ylabel('Pos');xlabel('t')
xlim(viewwin)

subplot(3,1,3)
plot(position.timestamps,position.decoderror)
xlim(viewwin)

NiceSave('Decoding',figfolder,baseName)

