function [ output_args ] = UPAttractorAnalysis( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% DEV
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/FRHetAndDynamics/AnalysisScripts/AnalysisFigs';
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
SlowWaves = bz_LoadEvents(basePath,'SlowWaves');
%% Get rate vector

dt = 0.01; %because maximum dt with maximum dimensionality
[spikemat] = bz_SpktToSpkmat(spikes,'binsize',dt);

[~,spikemat.NREMtimes] = RestrictInts(spikemat.timestamps,SleepState.ints.NREMstate);
spikemat.NREMspikemat = spikemat.data(spikemat.NREMtimes,:);
%%

%Similarity in: cosine, moholonobis distance, pearsons corr
% simwindow = 10; %s
% simwindow_dt = simwindow./spikemat.dt;
% batchsize = 10000;
% numbatches = floor(length(spikemat.NREMspikemat)./(batchsize)); %missing the last start of batch...
% 
% nearbyidx = zeros(1,batchsize);
% nearbyidx(1:simwindow_dt) = 1;
% nearbyidx = toeplitz(nearbyidx);
% nearbyidx = logical(nearbyidx);
% nearbyidx(1:simwindow_dt./2,:) = false;
% nearbyidx(end-simwindow_dt./2:end,:) = false;
% %%
% covmat = nancov(spikemat.NREMspikemat);
% clear nearby
% for bb = 1:numbatches
%     bb
%     batchidx = [1+(bb-1).*batchsize:(bb.*batchsize)];
%     batchspikemat = spikemat.NREMspikemat(batchidx,:);
%     
%     corrdist = pdist(batchspikemat,'mahalanobis',covmat);
%     %mdist = mahal(batchspikemat,spikemat.NREMspikemat);
%     distmat = squareform(corrdist);
%     nearby(batchidx,:) = distmat(nearbyidx);
% end
% 
% %%
% mdist = mahal(spikemat.NREMspikemat,spikemat.NREMspikemat);
% %%
% corrdist = pdist(spikemat.data(spikemat.NREMtimes,:),'spearman');

%%
%%
covmat = nancov(spikemat.NREMspikemat);
simwin = 1;
simwin_dt = simwin./spikemat.dt;
%Note: might want to normalize each cell by it's mean rate...
corrtypes = {'spearman','mahalanobis','jaccard'};
for cc = 1:length(corrtypes)
    simmat.(corrtypes{cc}) = zeros(2.*simwin_dt+1,2.*simwin_dt+1,length(SlowWaves.timestamps));
end

for dd = 1:length(SlowWaves.timestamps)
    dd
    %Find window around slow wave
    SWidx = interp1(spikemat.timestamps,1:length(spikemat.timestamps),SlowWaves.timestamps(dd),'nearest');
    SWwinidx = (SWidx-simwin_dt):(simwin_dt+SWidx);
    
    %Calculate population vector correlations
    corrdist = pdist(spikemat.data(SWwinidx,:),'mahalanobis',covmat);
    simmat.mahalanobis(:,:,dd) = squareform(corrdist);
    
    corrdist = pdist(spikemat.data(SWwinidx,:),'spearman');
    simmat.spearman(:,:,dd) = squareform(corrdist);
    %corrdist = pdist(spikemat.data(SWwinidx,:),'cityblock');
    %corrdist = pdist(spikemat.data(SWwinidx,:),'correlation');
    %corrdist = pdist(spikemat.data(SWwinidx,:),'cosine');
    %corrdist = pdist(spikemat.data(SWwinidx,:)>0,'hamming');
    corrdist = pdist(spikemat.data(SWwinidx,:)>0,'jaccard');
    simmat.jaccard(:,:,dd) = squareform(corrdist);
    
end
%%
t_simwin = -simwin:spikemat.dt:simwin;
for cc = 1:length(corrtypes)
    simmat_mean.(corrtypes{cc}) = nanmean(simmat.(corrtypes{cc}),3);
end

figure
subplot(3,3,1)
    imagesc(t_simwin,t_simwin,1-simmat_mean.spearman)
    colorbar
    caxis([0 0.15])
    
subplot(3,3,2)
    imagesc(t_simwin,t_simwin,simmat_mean.mahalanobis)
    colorbar
   caxis([7 10])
    
subplot(3,3,3)
    imagesc(t_simwin,t_simwin,1-simmat_mean.jaccard)
    colorbar
    caxis([0 0.07])
    
end

