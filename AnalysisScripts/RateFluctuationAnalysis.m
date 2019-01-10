function [ ] = RateFluctuationAnalysis(basePath,figfolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% DEV
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/FRHetAndDynamics/AnalysisScripts/AnalysisFigs/RateFluctuationAnalysis';
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');


%%
timebins = logspace(-3,1.75,40);
for tt = 1:length(timebins)
    %%
    tt
dt = timebins(tt);
[spikemat,t] = SpktToSpkmat(spikes.times, [], dt);

statenames = {'NREMstate','REMstate','WAKEstate'};
numstates = length(statenames);

%%
%Which state to look at?
for ss = 1:numstates
    %%
    %ss=2;
    
%Find which timebins are during state of interest
[~,statebins] = RestrictInts(t,SleepState.ints.(statenames{ss}));
    
instatespikemat = spikemat(statebins,:);
covmat = cov(instatespikemat);
E = eig(covmat);
dimensionalty.(statenames{ss})(tt) = sum(E).^2./sum(E.^2);
%Coresponds to ~80% of variability (Gau et al 2017)
%Should really look at E only..... (E/I different, see Bittner et al 2017.)
%Non-negative?

%Do the dimensions at faster time scales fit in the dimensions at longer
%time scales? Heirarchy? ... (subspace etc). 
%Consistent between sleep and wake?

spikecountmean.(statenames{ss})(tt,:) = mean(instatespikemat);
spikecountCV.(statenames{ss})(tt,:) = std(instatespikemat)./mean(instatespikemat);

SPKhist.bins = 1:100;
SPKhist.(statenames{ss}).hist(:,:,tt) = hist(instatespikemat,SPKhist.bins);

[~,sorts.(statenames{ss}).rate(tt,:)] = sort(spikecountmean.(statenames{ss})(tt,:));
[~,sorts.(statenames{ss}).CV(tt,:)] = sort(spikecountCV.(statenames{ss})(tt,:));

 %[C,LAGS] = xcov(instatespikemat,100,'coeff'); 
 
 %%
%  cmap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
%  figure
%  imagesc(C)
%  colorbar
%  colormap(cmap)
%  caxis([-1 1])
%end

end
end


%% Dimensionality of the rate trajectories... 
%(compare states)
figure
    plot(log10(timebins),dimensionalty.NREMstate,'b','linewidth',2)
    hold on
    plot(log10(timebins),dimensionalty.REMstate,'r','linewidth',2)
    plot(log10(timebins),dimensionalty.WAKEstate,'k','linewidth',2)
    legend('NREM','REM','WAKE')
    LogScale('x',10)
    xlabel('Bin Size (s)');ylabel(['Dimensionality (',num2str(length(spikes.UID)),' neurons)'])
NiceSave('Dimensionality',figfolder,baseName)    
%%
figure
subplot(2,2,1)
    plot(log10(spikecountmean(:,CellClass.pE)),log10(spikecountCV(:,CellClass.pE)),'k.')
    hold on
    plot(log10(spikecountmean(:,CellClass.pI)),log10(spikecountCV(:,CellClass.pI)),'r.')
    LogScale('xy',10)
subplot(2,2,2)
    imagesc(SPKhist.bins,1:length(spikes.UID),log10(SPKhist.hist(:,sorts.rate))')
    axis xy
%Dimensionality
%SpikeCountDistribtion


%%
figure
subplot(2,2,1)
    plot(log10(timebins),log10(spikecountCV(:,CellClass.pE)),'k')
    hold on
    plot(log10(timebins),log10(spikecountCV(:,CellClass.pI)),'r')
    xlabel('Bin Size (s)')
    ylabel('Spike Count CV')
    LogScale('xy',10)
subplot(2,2,2)
    imagesc(log10(timebins),SPKhist.bins,log10(squeeze(SPKhist.hist(:,9,:))))
    LogScale('x',10)

%%
popvectcorr = corr(spikemat','type','spearman');

%%
figure
subplot(4,1,2)
plot(SleepState.idx.timestamps,SleepState.idx.states)
subplot(2,1,2)
imagesc(t,t,popvectcorr)
%colorbar