reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyRate'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
regions = {'THAL','vCTX','fCTX','CA1'};

%%
for rr = 1:length(regions)
    [ISIstats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);

    [SpikeStatsbyRateAll,baseNames] = GetMatResults(figfolder,'SpikeStatsbyRate','baseNames',baseNames);
    SpikeStatsbyRateAll = bz_CollapseStruct(SpikeStatsbyRateAll);
    %thisregion = 'fCTX';

    %%
    ISIbyCspkRate.(regions{rr}) = bz_CollapseStruct(SpikeStatsbyRateAll.ISIbyCspkRate,'match','justcat',true);
    CV2byCspkRate.(regions{rr}) = bz_CollapseStruct(SpikeStatsbyRateAll.CV2byCspkRate,'match','justcat',true);
    binCV2byCspkRate.(regions{rr}) = bz_CollapseStruct(SpikeStatsbyRateAll.binCV2byCspkRate,'match','justcat',true);
    ISIbyRate.(regions{rr}) = bz_CollapseStruct(SpikeStatsbyRateAll.ISIbyRate,'match','justcat',true);
    CV2byRate.(regions{rr}) = bz_CollapseStruct(SpikeStatsbyRateAll.CV2byRate,'match','justcat',true);
    Cellmeanrate.(regions{rr}) = bz_CollapseStruct(SpikeStatsbyRateAll.Cellmeanrate,'match','justcat',true);

    states = {'WAKEstate','NREMstate','REMstate'};
    statecolors = {'k','b','r'};

    celltypes = fieldnames(ISIbyCspkRate.(regions{rr}).celltypeidx);
    for ss = 1:length(states)
    for tt = 1:length(celltypes)
        ISIbyCspkRate.(regions{rr}).(states{ss}).pop.(celltypes{tt}) = nanmean(ISIbyCspkRate.(regions{rr}).(states{ss}).pYX(:,:,ISIbyCspkRate.(regions{rr}).celltypeidx.(celltypes{tt})),3);
        CV2byCspkRate.(regions{rr}).(states{ss}).pop.(celltypes{tt}) = nanmean(CV2byCspkRate.(regions{rr}).(states{ss}).pYX(:,:,CV2byCspkRate.(regions{rr}).celltypeidx.(celltypes{tt})),3);
        binCV2byCspkRate.(regions{rr}).(states{ss}).pop.(celltypes{tt}) = nanmean(binCV2byCspkRate.(regions{rr}).(states{ss}).pYX(:,:,binCV2byCspkRate.(regions{rr}).celltypeidx.(celltypes{tt})),3);

        ISIbyRate.(regions{rr}).(states{ss}).pop.(celltypes{tt}) = nanmean(ISIbyRate.(regions{rr}).(states{ss}).XYprob(:,:,ISIbyRate.(regions{rr}).celltypeidx.(celltypes{tt})),3);
        CV2byRate.(regions{rr}).(states{ss}).pop.(celltypes{tt}) = nanmean(CV2byRate.(regions{rr}).(states{ss}).XYprob(:,:,CV2byRate.(regions{rr}).celltypeidx.(celltypes{tt})),3);

        [~,Cellmeanrate.(regions{rr}).sorts.(states{ss}).(celltypes{tt})] = sort(Cellmeanrate.(regions{rr}).(states{ss}));
        Cellmeanrate.(regions{rr}).sorts.(states{ss}).(celltypes{tt})(ismember(Cellmeanrate.(regions{rr}).sorts.(states{ss}).(celltypes{tt}),find(~ISIbyCspkRate.(regions{rr}).celltypeidx.(celltypes{tt}))))=[];
    end
    end

end
%% Figure: CSPKRWATE
for tt=1:2
figure
for ss = 1:3
    state = states{ss};
%for tt = 1:length(celltypes)
subplot(4,3,ss)
    imagesc(ISIbyCspkRate.(regions{rr}).(state).Xbins(1,:,1),ISIbyCspkRate.(regions{rr}).(state).Ybins(1,:,1), ISIbyCspkRate.(regions{rr}).(state).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('xy',10)
    if ss==1
    ylabel('ISI (s)');
    end
    %xlabel('Rate')
    %title((celltypes{tt}))
    %colorbar
    if tt ==1 
        caxis([0 0.5e-1])
    elseif tt==2
         caxis([0 2.5e-3])
    end
    
%end 


%for tt = 1:length(celltypes)
subplot(4,3,ss+3)
    imagesc(CV2byCspkRate.(regions{rr}).(state).Xbins(1,:,1),CV2byCspkRate.(regions{rr}).(state).Ybins(1,:,1), CV2byCspkRate.(regions{rr}).(state).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('x',10)
    if ss==1
    ylabel('CV2');end
    %title((celltypes{tt}))
	%colorbar
    if tt ==1 
        caxis([0 0.04])
    elseif tt==2
         caxis([0 1e-3])
    end
    
%end 

%for tt = 1:length(celltypes)
subplot(4,3,ss+9)
    imagesc(CV2byCspkRate.(regions{rr}).(state).Xbins(1,:,1),[1 sum(ISIbyCspkRate.(regions{rr}).celltypeidx.(celltypes{tt}))],...
        squeeze(CV2byCspkRate.(regions{rr}).(state).pX(:,:,Cellmeanrate.(regions{rr}).sorts.(states{ss}).(celltypes{tt})))')
    hold on
    plot(log10(Cellmeanrate.(regions{rr}).(states{ss})(Cellmeanrate.(regions{rr}).sorts.(states{ss}).(celltypes{tt}))),[1:sum(ISIbyCspkRate.(regions{rr}).celltypeidx.(celltypes{tt}))],'w')
    axis xy
    LogScale('x',10)
    if ss==1
    ylabel('Cell');
    end
    xlabel('Rate')
    %title((celltypes{tt}))
    %colorbar
    if tt ==1 
        caxis([0 0.1])
    elseif tt==2
         caxis([0 0.15])
    end
    
%end 

%for tt = 1:length(celltypes)
subplot(4,3,ss+6)
    imagesc(binCV2byCspkRate.(regions{rr}).(state).Xbins(1,:,1),binCV2byCspkRate.(regions{rr}).(state).Ybins(1,:,1), binCV2byCspkRate.(regions{rr}).(state).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('x',10)
    if ss==1
    ylabel('<CV2>');end
    %title((celltypes{tt}))
	%colorbar
    if tt ==1 
        caxis([0 0.05])
    elseif tt==2
         caxis([0 0.05])
    end
%end
end 

 NiceSave(['ISIbycBinRate_',(celltypes{tt})],figfolder,thisregion)
end




%% Figure: SpikeRate
for tt=1:2
figure
for ss = 1:3
    state = states{ss};
%for tt = 1:length(celltypes)
subplot(4,3,ss)
    imagesc(ISIbyRate.(regions{rr}).(state).Xbins(1,:,1),ISIbyRate.(regions{rr}).(state).Ybins(1,:,1), ISIbyRate.(regions{rr}).(state).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    LogScale('y',10)
    if ss==1
    ylabel('ISI (s)');
    end
    %xlabel('Rate')
    %title((celltypes{tt}))
    %colorbar
    if tt ==1 
        xlim([0 10])
        caxis([0 3e-3])
    elseif tt==2
         caxis([0 00.5e-3])
    end
    
%end 


%for tt = 1:length(celltypes)
subplot(4,3,ss+3)
    imagesc(CV2byRate.(regions{rr}).(state).Xbins(1,:,1),CV2byRate.(regions{rr}).(state).Ybins(1,:,1), CV2byRate.(regions{rr}).(state).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    if ss==1
    ylabel('CV2');end
    %title((celltypes{tt}))
	%colorbar
    if tt ==1 
        xlim([0 10])
        caxis([0 3e-3])
    elseif tt==2
         caxis([0 5e-4])
    end
    
%end 

%for tt = 1:length(celltypes)
subplot(4,3,ss+6)
    imagesc(CV2byRate.(regions{rr}).(state).Xbins(1,:,1),[1 sum(ISIbyRate.(regions{rr}).celltypeidx.(celltypes{tt}))],...
        squeeze(CV2byRate.(regions{rr}).(state).pX(:,:,Cellmeanrate.(regions{rr}).sorts.(states{ss}).(celltypes{tt})))')
    hold on
    plot((Cellmeanrate.(regions{rr}).(states{ss})(Cellmeanrate.(regions{rr}).sorts.(states{ss}).(celltypes{tt}))),[1:sum(ISIbyRate.(regions{rr}).celltypeidx.(celltypes{tt}))],'w')
    axis xy
    %LogScale('y',10)
    if ss==1
    ylabel('Cell');
    end
    xlabel('Rate')
    %title((celltypes{tt}))
    %colorbar
     if tt ==1 
        xlim([0 10])
         caxis([0 0.3])
     elseif tt==2
          caxis([0 0.2])
     end
    
end 
 NiceSave(['ISIbyRate.(regions{rr})_',(celltypes{tt})],figfolder,thisregion)

end
