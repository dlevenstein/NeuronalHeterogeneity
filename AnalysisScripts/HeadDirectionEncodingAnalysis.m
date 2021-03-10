function [model_m,model_c,tuningcurve] = HeadDirectionEncodingAnalysis(basePath,figfolder)
% Date XX/XX/20XX %ISIbyHD_alignGam
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
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse12-120807';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse24-131213';
%basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true,'region','THAL');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
%ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

%[celltypes,~,typeidx] = unique(CellClass.label);
%cellcolor = {'k','r'};

%% Get gamma mode
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');
%%
headdir.samplingRate = 39.06; %Hz
hdfilename = fullfile(basePath,[baseName,'.ang']);
headdir.data = importdata(hdfilename);
headdir.data(headdir.data==-1)=nan;
headdir.timestamps = [1:length(headdir.data)]'./headdir.samplingRate; 

%%
% nantimes = isnan(headdir.position.lin);
% headdir.data = headdir.position.lin(~(nantimes));
%Possible here: remove drops for more than a... second?
%minjump = 2;
%[ headdir.data ] = NanPadJumps( headdir.timestamps(~(nantimes)),headdir.data,minjump );
%headdir.data = interp1(headdir.timestamps(~(nantimes)),headdir.data,headdir.timestamps);
numXbins = 100;
[ISIbyHD] = bz_ConditionalISI(spikes.times,headdir,...
    'ints',SleepState.ints.WAKEstate,...
    'showfig',false,'GammaFit',false,'numXbins',numXbins,'numISIbins',100,...
    'normtype','none','Xwin',[0 2.*pi],'Xbinoverlap',3);

%% Get Tuning Specificity (Skaggs)

meanrate = nansum(ISIbyHD.Dist.SpikeRate.*ISIbyHD.Dist.pX,2);
meanrate_full = repmat(meanrate,1,numXbins,1);
MutInfo.Skaggs = squeeze(nansum(ISIbyHD.Dist.SpikeRate.*log2(ISIbyHD.Dist.SpikeRate./meanrate_full).*ISIbyHD.Dist.pX,2));
MutInfo.Skaggs_sec = MutInfo.Skaggs .*squeeze(meanrate);


%%
[~,ISIbyHD.fieldpeak] = max(ISIbyHD.Dist.SpikeRate,[],2);
ISIbyHD.fieldpeak = ISIbyHD.Dist.Xbins(ISIbyHD.fieldpeak);

%%




%%
MutInfo.ISI = squeeze(ISIbyHD.MutInf)';
binsize = 0.05;
dt = 0.025;
spkmat = bz_SpktToSpkmat(spikes.times,'dt',dt,'binsize',binsize,'units','counts');
spkmat.pos = interp1(headdir.timestamps,headdir.data,spkmat.timestamps,'nearest');
spkmat.InWake = InIntervals(spkmat.timestamps,SleepState.ints.WAKEstate);

%Remove bins in which position changes too much


maxspikes = 15;
for cc = 1:spikes.numcells
    MutInfo.Rate(cc) = mutualinfo(spkmat.data(spkmat.InWake,cc),spkmat.pos(spkmat.InWake));
    CONDXY(cc) = ConditionalHist(spkmat.pos(spkmat.InWake),spkmat.data(spkmat.InWake,cc),...
         'numXbins',20,'Xbounds',[0 2*pi],'numYbins',maxspikes+1,'Ybounds',[0 maxspikes],...
         'Xbinoverlap',2);
end

%% Spike count histogram
MIthresh = 0.3;
HDcells = find(MutInfo.Skaggs>MIthresh);
[~,HDcells_GammaIDX] = intersect(GammaFit.WAKEstate.cellstats.UID,HDcells,'stable');

figure
for cc = 1:length(HDcells)
    subplot(4,3,cc)
        imagesc(CONDXY(1).Xbins,CONDXY(1).Ybins,CONDXY(HDcells(cc)).pYX')
        hold on
        %imagesc(CONDXY(1).Xbins+2.*pi,CONDXY(1).Ybins,CONDXY(HDcells(cc)).pYX')
        axis xy
        %xlim([0 4*pi])
        bz_piTickLabel('x')
        ColorbarWithAxis([0 0.4],'P[s|HD]')
        xlabel('Head Direction');ylabel('Spike Count')
        title(['UID: ',num2str(spikes.UID(HDcells(cc)))])
    
   if cc == 12
       break
   end
end
NiceSave('HDCells_SpikeCount',figfolder,baseName)
%% Example cell



%testcell = 29;
%testcell = 21;
testcell = randsample(HDcells,1);
FitEncodingModel_HD(spkmat.data(spkmat.InWake,testcell),spkmat.pos(spkmat.InWake),binsize);

subplot(2,3,4)
try
bz_PlotISIDistModes(GammaFit.WAKEstate,testcell)
catch
end
NiceSave(['EncodingModelFit_UID',num2str(testcell)],figfolder,baseName)

%% Run All cells
for cc = 1:length(HDcells)
cc
[model_m(cc),model_c(cc),tuningcurve(cc)] = FitEncodingModel_HD(spkmat.data(spkmat.InWake,HDcells(cc)),spkmat.pos(spkmat.InWake),binsize);
model_m(cc).GSRate = GammaFit.WAKEstate.sharedfit.GSlogrates(HDcells_GammaIDX(cc));
model_c(cc).GSRate = GammaFit.WAKEstate.sharedfit.GSlogrates(HDcells_GammaIDX(cc));
close all
end
%%
allcells_m = bz_CollapseStruct(model_m,'match','justcat',true);
allcells_c = bz_CollapseStruct(model_c,'match','justcat',true);
%%
figure
subplot(3,4,1)
BoxAndScatterPlot(allcells_m.BIC-allcells_c.BIC)
hold on
box off
UnityLine
ylabel('Modal-Continuous Model BIC')
subplot(3,3,2)
plot(allcells_m.parms.pAS_0,allcells_m.parms.pAS_pi,'.')
hold on
UnityLine
xlabel('p[AS|x=x_0]');ylabel('p[AS|x=x_0+pi]')

subplot(3,3,3)
plot(allcells_m.parms.rAS,allcells_m.parms.pAS_0,'.')
hold on
%UnityLine
xlabel('R_A_S');ylabel('p[AS|x=x_0]')

subplot(3,3,4)
plot(allcells_m.parms.rGS,allcells_m.parms.pAS_pi,'.')
hold on
%UnityLine
box off
xlabel('R_G_S');ylabel('p[AS|x=x_0+pi]')

subplot(3,3,5)
plot(log10(allcells_m.parms.rGS),GammaFit.WAKEstate.sharedfit.GSlogrates(HDcells_GammaIDX),'.')
LogScale('xy',10)
box off
xlabel('GS Rate: Encoding Model');ylabel('GS Rate: Gamma ISI Model')
NiceSave('EncodingModelStats',figfolder,baseName)
%%


