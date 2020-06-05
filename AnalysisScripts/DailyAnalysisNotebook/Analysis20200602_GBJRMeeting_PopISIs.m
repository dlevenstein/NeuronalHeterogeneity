function [PopConditional,PopConditional_phase] = ISIModulationAnalysis(basePath,figfolder)
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
% basePath = '/Users/dl2820/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
% %basePath = pwd;
% %basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
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
statecolors = {[0 0 0],[0 0 1],[1 0 0],[0.6 0.6 0.6]};

try
celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};


%%
for tt = 1:length(celltypes)
    popspikes.(celltypes{tt}) = cat(1,spikes.times{CellClass.(celltypes{tt})});
end

%% Synchrony: compare bins
numbins = 30;
binrange = [0.003 30];
mindt = 0.005; %minimum acceptable dt
% numbins = 10;
% binrange = [0.005 10];
% mindt = 0.005; %minimum acceptable dt
binsizes = logspace(log10(binrange(1)),log10(binrange(2)),numbins);
dt = binsizes./3;%max(binsizes/2,mindt);

%%

for bb = 1:numbins
    bz_Counter(bb,numbins,'Bin')
spikemat = bz_SpktToSpkmat(spikes,'dt',dt(bb),'binsize',binsizes(bb),'units','rate','bintype','boxcar','showfig',false);
spikemat.synch = sum(spikemat.data>0.5,2);
%%
ss = 1;
tt = 1;
synch.data = spikemat.synch + 0.1.*randn(size(spikemat.synch));
synch.timestamps = spikemat.timestamps;
    for ss = 1:3
    [PopConditionalISIDist.(states{ss}).(celltypes{tt})] = ConditionalISI(popspikes.(celltypes{tt}),synch,...
        'ints',SleepState.ints.(states{ss}),'ISIDist',true);
    popMI.(states{ss})(bb) = PopConditionalISIDist.(states{ss}).(celltypes{tt}).MutInf;
    end
end

%%
figure
subplot(2,2,1)
hold on
for ss = 1:3
plot(log10(binsizes),popMI.(states{ss}),'linewidth',2,'color',statecolors{ss})
end
axis tight
LogScale('x',10,'nohalf',true)
xlabel('Bin Timescale')
NiceSave('BinSizeMI',figfolder,baseName)
%% Specific bin - for examples etc
dt = 0.01;binsize = 0.03;
spikemat = bz_SpktToSpkmat(spikes,'dt',dt,'binsize',binsize,'units','rate','bintype','boxcar');
spikemat.synch = sum(spikemat.data>0.5,2);
%%
ss = 2;
tt = 1;
synch.data = spikemat.synch + 0.1.*randn(size(spikemat.synch));
synch.timestamps = spikemat.timestamps;
[PopConditionalISIDist.(states{ss}).(celltypes{tt})] = ConditionalISI(popspikes.(celltypes{tt}),synch,...
    'ints',SleepState.ints.(states{ss}),'ISIDist',true);



end