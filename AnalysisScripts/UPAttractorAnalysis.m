function [ output_args ] = UPAttractorAnalysis( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% DEV
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/FRHetAndDynamics/AnalysisScripts/AnalysisFigs/UPAttractorAnalysis';
figfolder = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/UPAttractorAnalysis';
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
SlowWaves = bz_LoadEvents(basePath,'SlowWaves');

%%
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
%% ISI stats related to SW time
SWwin = 1;
prewin = [-0.2 0];
postwin = [0 0.2];

for dd = 1:length(SlowWaves.timestamps)
    %Find Spikes near the slow waves
    dd
    %thisSWwin = SlowWaves.timestamps+SWwin.*[-1 1];
    reltime = cellfun(@(X) X-SlowWaves.timestamps(dd),ISIStats.allspikes.times,'UniformOutput',false);
    SWrelidx = cellfun(@(X) abs(X)<SWwin,reltime,'UniformOutput',false);
    SWrelCV2(dd,:) = cellfun(@(X,Y) X(Y),ISIStats.allspikes.CV2,SWrelidx,'UniformOutput',false);
    SWreltime(dd,:) = cellfun(@(X,Y) X(Y),reltime,SWrelidx,'UniformOutput',false);
    ISIn(dd,:) = cellfun(@(X,Y) X(Y),ISIStats.allspikes.ISIs,SWrelidx,'UniformOutput',false);
    ISInp1(dd,:) = cellfun(@(X,Y) X([false; Y(1:end-1)]),ISIStats.allspikes.ISIs,SWrelidx,'UniformOutput',false);

    
    
end

%%
for cc = 1:spikes.numcells
    SWrelCV2_all{cc} = cat(1,SWrelCV2{:,cc});
    SWreltime_all{cc} = cat(1,SWreltime{:,cc});
    ISIn_all{cc} = cat(1,ISIn{:,cc});
    ISInp1_all{cc} = cat(1,ISInp1{:,cc});
end

%%
binsize = 0.1;
binedges = -SWwin:binsize:SWwin;
[ bincenters,binmeans,binstd,binnum,binneddata ] = cellfun(@(X,Y) BinDataTimes(X,Y,binedges),SWrelCV2_all,SWreltime_all,'UniformOutput',false);
bincenters = bincenters{1};
binmeans = cat(2,binmeans{:});
binstd = cat(2,binstd{:});

%%
celltypes = unique(CellClass.label);
cellcolors = {'k','r'};
for tt = 1:length(celltypes)
SWCV2bypop.(celltypes{tt}).mean = mean(binmeans(:,CellClass.(celltypes{tt})),2);
SWCV2bypop.(celltypes{tt}).std = std(binmeans(:,CellClass.(celltypes{tt})),[],2);
end


%%
postSWreturnmaps = cellfun(@(X,Y,Z) hist3([log10(X(Z>0&Z<win)) log10(Y(Z>0&Z<win))],{ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins}),...
    ISIn_all,ISInp1_all,SWreltime_all,'UniformOutput',false);
postSWreturnmaps = cellfun(@(X) X./sum(X(:)),postSWreturnmaps,'UniformOutput',false);
postSWreturnmaps = cat(3,postSWreturnmaps{:});

preSWreturnmaps = cellfun(@(X,Y,Z) hist3([log10(X(Z<0&Z>-win)) log10(Y(Z<0&Z>-win))],{ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins}),...
    ISIn_all,ISInp1_all,SWreltime_all,'UniformOutput',false);
preSWreturnmaps = cellfun(@(X) X./sum(X(:)),preSWreturnmaps,'UniformOutput',false);
preSWreturnmaps = cat(3,preSWreturnmaps{:});

nonSWreturnmaps = cellfun(@(X,Y,Z) hist3([log10(X(Z<-win|Z>win)) log10(Y(Z<-win|Z>win))],{ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins}),...
    ISIn_all,ISInp1_all,SWreltime_all,'UniformOutput',false);
nonSWreturnmaps = cellfun(@(X) X./sum(X(:)),nonSWreturnmaps,'UniformOutput',false);
nonSWreturnmaps = cat(3,nonSWreturnmaps{:});


for tt = 1:length(celltypes)
    meanpostSWreturn.(celltypes{tt}) = mean(postSWreturnmaps(:,:,CellClass.(celltypes{tt})),3);
    meanpreSWreturn.(celltypes{tt}) = mean(preSWreturnmaps(:,:,CellClass.(celltypes{tt})),3);
    meannonSWreturn.(celltypes{tt}) = mean(nonSWreturnmaps(:,:,CellClass.(celltypes{tt})),3);
end

%%
figure
subplot(2,2,1)
    imagesc(bincenters,[1 spikes.numcells],binmeans(:,ISIStats.sorts.NREMstate.CV2byclass)')
    hold on
    plot([0 0],[1 spikes.numcells],'k')
    xlim([-1 1])
    colorbar('northoutside')
    caxis([0.9 1.4])
    ylabel('Cell - Sorted by <CV2>, type')
    
% subplot(2,2,2)
%     imagesc(bincenters,[1 spikes.numcells],binstd(:,ISIStats.sorts.NREMstate.CV2byclass)')
%     hold on
%     plot([0 0],[1 spikes.numcells],'k')
%     xlim([-1 1])
%     colorbar
%    % caxis([0.5 1.5])
   
subplot(4,2,5)
plot([0 0],[0.7 1.5],'k')
hold on
for tt = 1:length(celltypes)
    plot(bincenters,SWCV2bypop.(celltypes{tt}).mean,'o-','color',cellcolors{tt})
    hold on
    errorshade(bincenters,SWCV2bypop.(celltypes{tt}).mean,...
        SWCV2bypop.(celltypes{tt}).std,SWCV2bypop.(celltypes{tt}).std,cellcolors{tt},'scalar')
end
ylim([0.8 1.4])
xlabel('t - relative to SW');ylabel('<CV2>')


for tt = 1:length(celltypes)
    subplot(3,4,2+tt)
    colormap(gca,histcolors)
        imagesc(ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins,...
            meanpreSWreturn.(celltypes{tt})')
        LogScale('xy',10)
        %xlabel('ISI_n');ylabel('ISI_n+1')
        axis xy
        caxis([0 5e-3])
        title(celltypes{tt})
        %colorbar
        
    subplot(3,4,6+tt)
    colormap(gca,histcolors)
        imagesc(ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins,...
            meanpostSWreturn.(celltypes{tt})')
        LogScale('xy',10)
        %xlabel('Preceeding ISI');ylabel('Next ISI')
        axis xy
        %colorbar
        caxis([0 5e-3])
        
    subplot(3,4,10+tt)
    colormap(gca,histcolors)
        imagesc(ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins,...
            meannonSWreturn.(celltypes{tt})')
        LogScale('xy',10)
        %xlabel('Preceeding ISI');ylabel('Next ISI')
        axis xy
        caxis([0 5e-3])
        
end

NiceSave('CVandSW',figfolder,baseName)
%%
%cc = 10;

win = 0.15;
%[ bincenters,binmeans,binstd,binnum,binneddata ] = BinDataTimes(SWrelCV2_all{cc},SWreltime_all{cc},binedges );
excells = [randsample(find(CellClass.pE),1) randsample(find(CellClass.pI),1)];
histcolors = flipud(gray);

figure
for ee = 1:length(excells)
    thiscell = excells(ee);
    subplot(3,2,ee)
        plot(SWreltime_all{thiscell},SWrelCV2_all{thiscell},'.','color',cellcolors{ee},'markersize',0.5)
        hold on
        plot(SWreltime_all{thiscell}(SWreltime_all{thiscell}>0 & SWreltime_all{thiscell}<win),...
            SWrelCV2_all{thiscell}(SWreltime_all{thiscell}>0 & SWreltime_all{thiscell}<win),'g.','markersize',0.5)
        plot(bincenters,binmeans(:,thiscell),'o-','color',cellcolors{ee})
        errorshade(bincenters,binmeans(:,thiscell),binstd(:,thiscell),binstd(:,thiscell),cellcolors{ee},'scalar')
        xlim([-1 1])
    
    subplot(2,2,2+ee)
    colormap(histcolors)
        imagesc(ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins,ISIStats.ISIhist.NREMstate.return(:,:,thiscell))
        hold on
        plot(log10(ISIn_all{thiscell}(SWreltime_all{thiscell}<0 & SWreltime_all{thiscell}>-win)),...
            log10(ISInp1_all{thiscell}(SWreltime_all{thiscell}<0 & SWreltime_all{thiscell}>-win)),...
            '.','color',cellcolors{ee},'markersize',0.5)
       % hold on
        plot(log10(ISIn_all{thiscell}(SWreltime_all{thiscell}>0 & SWreltime_all{thiscell}<win)),...
            log10(ISInp1_all{thiscell}(SWreltime_all{thiscell}>0 & SWreltime_all{thiscell}<win)),...
            'g.','markersize',0.5)
        LogScale('xy',10)
        xlabel('Preceeding ISI');ylabel('Next ISI')
        axis xy
            
end



  
        
        
        
        
%%    Attractor change    
%% Get rate vector

binsize = 0.1; %because maximum dt with maximum dimensionality
[spikemat] = bz_SpktToSpkmat(spikes,'binsize',binsize,'overlap',6);

[~,spikemat.NREMtimes] = RestrictInts(spikemat.timestamps,SleepState.ints.NREMstate);
spikemat.NREMspikemat = spikemat.data(spikemat.NREMtimes,:);

spikemat.data = bsxfun(@(X,Y) X./Y,spikemat.data,mean(spikemat.NREMspikemat,1));


%%
%covmat = nancov(spikemat.NREMspikemat);
simwin = 1.5;
simwin_dt = simwin./spikemat.dt;
%Note: might want to normalize each cell by it's mean rate...
corrtypes = {'spearman','jaccard'};
for cc = 1:length(corrtypes)
    simmat.(corrtypes{cc}) = zeros(2.*simwin_dt+1,2.*simwin_dt+1,length(SlowWaves.timestamps));
end

SWidx = interp1(spikemat.timestamps,1:length(spikemat.timestamps),SlowWaves.timestamps,'nearest');
for dd = 1:length(SlowWaves.timestamps)
    dd
    %Find window around slow wave
    for tt = 1:length(celltypes)
    SWwinidx = (SWidx(dd)-simwin_dt):(simwin_dt+SWidx(dd));
    
    %Calculate population vector correlations
%     corrdist = pdist(spikemat.data(SWwinidx,:),'mahalanobis',covmat);
%     simmat.mahalanobis(:,:,dd) = squareform(corrdist);
    
    corrdist = pdist(spikemat.data(SWwinidx,CellClass.(celltypes{tt})),'spearman');
    simmat.(celltypes{tt}).spearman(:,:,dd) = squareform(corrdist);
    %corrdist = pdist(spikemat.data(SWwinidx,:),'cityblock');
    %corrdist = pdist(spikemat.data(SWwinidx,:),'correlation');
    %corrdist = pdist(spikemat.data(SWwinidx,:),'cosine');
    %corrdist = pdist(spikemat.data(SWwinidx,:)>0,'hamming');
    corrdist = pdist(spikemat.data(SWwinidx,CellClass.(celltypes{tt}))>0,'jaccard');
    simmat.(celltypes{tt}).jaccard(:,:,dd) = squareform(corrdist);
    end
    
    corrdist = pdist(spikemat.data(SWwinidx,:),'spearman');
    simmat.ALL.spearman(:,:,dd) = squareform(corrdist);
    %corrdist = pdist(spikemat.data(SWwinidx,:),'cityblock');
    %corrdist = pdist(spikemat.data(SWwinidx,:),'correlation');
    %corrdist = pdist(spikemat.data(SWwinidx,:),'cosine');
    %corrdist = pdist(spikemat.data(SWwinidx,:)>0,'hamming');
    corrdist = pdist(spikemat.data(SWwinidx,:)>0,'jaccard');
    simmat.ALL.jaccard(:,:,dd) = squareform(corrdist);
end
%%
t_simwin = -simwin:spikemat.dt:simwin;
types = {'pE','pI','ALL'};
for cc = 1:length(corrtypes)
    for tt=1:length(types)
    simmat_mean.(types{tt}).(corrtypes{cc}) = nanmean(simmat.(types{tt}).(corrtypes{cc}),3);
    end
end

tbounds = 0.7.*[-1 1];
figure
    subplot(3,3,1)
        imagesc(t_simwin,t_simwin,1-simmat_mean.pE.spearman)
        hold on
        plot([0 0],tbounds,'k');plot(tbounds,[0 0],'k');
        colorbar
        caxis([0 0.15])
        xlim(tbounds);ylim(tbounds)
        title('Pop. vector corr')
        xlabel('t (rel to SW, s)');ylabel({'pE','t (rel to SW, s)'})
 
    subplot(3,3,2)
        imagesc(t_simwin,t_simwin,1-simmat_mean.pE.jaccard)
        hold on
        plot([0 0],tbounds,'k');plot(tbounds,[0 0],'k');
        colorbar
        caxis([0 0.13])
        xlim(tbounds);ylim(tbounds)
        title('Cell Overlap')
        xlabel('t (rel to SW, s)');ylabel('t (rel to S, s)')
        
    subplot(3,3,4)
        imagesc(t_simwin,t_simwin,1-simmat_mean.pI.spearman)
        hold on
        plot([0 0],tbounds,'k');plot(tbounds,[0 0],'k');
        colorbar
        caxis([0 0.4])
        xlim(tbounds);ylim(tbounds)
        xlabel('t (rel to SW, s)');ylabel({'pI','t (rel to SW, s)'})
 
    subplot(3,3,5)
        imagesc(t_simwin,t_simwin,1-simmat_mean.pI.jaccard)
        hold on
        plot([0 0],tbounds,'k');plot(tbounds,[0 0],'k');
        colorbar
        caxis([0 0.5])
        xlim(tbounds);ylim(tbounds)
        xlabel('t (rel to SW, s)');ylabel('t (rel to SW, s)')

        
    subplot(3,3,7)
        imagesc(t_simwin,t_simwin,1-simmat_mean.ALL.spearman)
        hold on
        plot([0 0],tbounds,'k');plot(tbounds,[0 0],'k');
        colorbar
        caxis([0 0.3])
        xlim(tbounds);ylim(tbounds)
        xlabel('t (rel to SW, s)');ylabel({'ALL','t (rel to SW, s)'})
 
    subplot(3,3,8)
        imagesc(t_simwin,t_simwin,1-simmat_mean.ALL.jaccard)
        hold on
        plot([0 0],tbounds,'k');plot(tbounds,[0 0],'k');
        colorbar
        caxis([0 0.3])
        xlim(tbounds);ylim(tbounds)
        xlabel('t (rel to SW, s)');ylabel('t (rel to SW, s)')
NiceSave('UPStateSimilarity',figfolder,baseName)
end

