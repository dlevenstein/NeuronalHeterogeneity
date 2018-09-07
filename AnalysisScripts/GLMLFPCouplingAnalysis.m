function [  ] = GLMLFPCouplingAnalysis(basePath,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/GLMLFPCouplingAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
%SleepState.ints.ALL = [0 Inf];
statenames = fieldnames(SleepState.ints);
numstates = length(statenames);
statecolors = {'k','b','r'};
classnames = unique(CellClass.label);
numclasses = length(classnames);
classcolors = {'k','r'};
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
     'basepath',basePath,'downsample',2);
sessionInfo = bz_getSessionInfo(basePath);

%% Which cell to use?
usecell = 9; %Cell from jenn animal with lots of spikes (pI cell?)

%% get the wavelets (NREM only)

wavespec = bz_WaveSpec(lfp,'nfreqs',10,'frange',[1 100]);

%%
tic
[ GLMFP ] = GLMLFP(spikes.times(usecell),wavespec,...
    'intervals',SleepState.ints.NREMstate );
toc

%%

%% Simulate ISIs
for tt = 1:length(GLMFP.timestamps)
    GLMFP.spkmat(tt) = rand(1)<=GLMFP.predRate(tt);
end
GLMFP.spktimes = find(GLMFP.spkmat).*dt;
%% ISI
isis = diff(spktimes);
GLMFP.isis = diff(GLMFP.spktimes);

isidist.bins = linspace(-3,1,40);
isidist.hist = hist(log10(isis),isidist.bins);

isidist.histsim = hist(log10(GLMFP.isis),isidist.bins);
%%
figure
plot(GLMFP.powerbins,GLMFP.Rpower)
%%
xwin = SleepState.ints.NREMstate(1,1)+[0 5];
figure
subplot(2,2,3)
%plot(GLMFP.powerbins,GLMFP.Rpower,'k','linewidth',2)
imagesc(log10(wavespec.freqs),GLMFP.powerbins,GLMFP.Rpower)
colorbar
axis xy
axis tight
xlabel('Freq (Hz)');ylabel('Power')
LogScale('x',10)
subplot(4,1,2)
plot(GLMFP.timestamps,GLMFP.predRate,'k')
hold on
%plot(t,rate./1000,'r')
legend('Predicted','Actual')
ylabel('Rate (spk/ms)')
xlabel('t (s)')
xlim(xwin)

subplot(4,1,1)
plot(lfp.timestamps,lfp.data,'k')
hold on
plot(spikes.times{usecell},ones(size(spikes.times{usecell})),'k.')
ylabel('"LFP"')
xlabel('t (s)')
xlim(xwin)

% subplot(2,2,4)
% plot(isidist.bins,isidist.hist,'r','linewidth',2)
% hold on
% plot(isidist.bins,isidist.histsim,'k','linewidth',2)
% legend('Actual','Simulated','location','northwest')
% LogScale('x',10)
% xlabel('ISI (s)')

end