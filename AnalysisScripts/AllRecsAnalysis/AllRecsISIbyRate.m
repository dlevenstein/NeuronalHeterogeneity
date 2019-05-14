%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyRate'];

%datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/BW_CTX';
%datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC';
%regions = {'fCTX','CA1'};
%datasetPath = figfolder;

[SpikeStatsbyRateAll,baseNames] = GetMatResults(figfolder,'SpikeStatsbyRate','select',true);
SpikeStatsbyRateAll = bz_CollapseStruct(SpikeStatsbyRateAll);
thisregion = 'fCTX';

%%
ISIbyCspkRate = bz_CollapseStruct(SpikeStatsbyRateAll.ISIbyCspkRate,'match','justcat',true);
CV2byCspkRate = bz_CollapseStruct(SpikeStatsbyRateAll.CV2byCspkRate,'match','justcat',true);
binCV2byCspkRate = bz_CollapseStruct(SpikeStatsbyRateAll.binCV2byCspkRate,'match','justcat',true);
ISIbyRate = bz_CollapseStruct(SpikeStatsbyRateAll.ISIbyRate,'match','justcat',true);
CV2byRate = bz_CollapseStruct(SpikeStatsbyRateAll.CV2byRate,'match','justcat',true);
Cellmeanrate = bz_CollapseStruct(SpikeStatsbyRateAll.Cellmeanrate,'match','justcat',true);

states = {'WAKEstate','NREMstate','REMstate'};
statecolors = {'k','b','r'};

celltypes = fieldnames(ISIbyCspkRate.celltypeidx);
for ss = 1:length(states)
for tt = 1:length(celltypes)
    ISIbyCspkRate.(states{ss}).pop.(celltypes{tt}) = nanmean(ISIbyCspkRate.(states{ss}).pYX(:,:,ISIbyCspkRate.celltypeidx.(celltypes{tt})),3);
    CV2byCspkRate.(states{ss}).pop.(celltypes{tt}) = nanmean(CV2byCspkRate.(states{ss}).pYX(:,:,CV2byCspkRate.celltypeidx.(celltypes{tt})),3);
    binCV2byCspkRate.(states{ss}).pop.(celltypes{tt}) = nanmean(binCV2byCspkRate.(states{ss}).pYX(:,:,binCV2byCspkRate.celltypeidx.(celltypes{tt})),3);

    ISIbyRate.(states{ss}).pop.(celltypes{tt}) = nanmean(ISIbyRate.(states{ss}).XYprob(:,:,ISIbyRate.celltypeidx.(celltypes{tt})),3);
    CV2byRate.(states{ss}).pop.(celltypes{tt}) = nanmean(CV2byRate.(states{ss}).XYprob(:,:,CV2byRate.celltypeidx.(celltypes{tt})),3);

    [~,Cellmeanrate.sorts.(states{ss}).(celltypes{tt})] = sort(Cellmeanrate.(states{ss}));
    Cellmeanrate.sorts.(states{ss}).(celltypes{tt})(ismember(Cellmeanrate.sorts.(states{ss}).(celltypes{tt}),find(~ISIbyCspkRate.celltypeidx.(celltypes{tt}))))=[];
end
end


%% Figure: CSPKRWATE
for tt=1:2
figure
for ss = 1:3
    state = states{ss};
%for tt = 1:length(celltypes)
subplot(4,3,ss)
    imagesc(ISIbyCspkRate.(state).Xbins(1,:,1),ISIbyCspkRate.(state).Ybins(1,:,1), ISIbyCspkRate.(state).pop.(celltypes{tt})')
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
    imagesc(CV2byCspkRate.(state).Xbins(1,:,1),CV2byCspkRate.(state).Ybins(1,:,1), CV2byCspkRate.(state).pop.(celltypes{tt})')
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
    imagesc(CV2byCspkRate.(state).Xbins(1,:,1),[1 sum(ISIbyCspkRate.celltypeidx.(celltypes{tt}))],...
        squeeze(CV2byCspkRate.(state).pX(:,:,Cellmeanrate.sorts.(states{ss}).(celltypes{tt})))')
    hold on
    plot(log10(Cellmeanrate.(states{ss})(Cellmeanrate.sorts.(states{ss}).(celltypes{tt}))),[1:sum(ISIbyCspkRate.celltypeidx.(celltypes{tt}))],'w')
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
    imagesc(binCV2byCspkRate.(state).Xbins(1,:,1),binCV2byCspkRate.(state).Ybins(1,:,1), binCV2byCspkRate.(state).pop.(celltypes{tt})')
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
    imagesc(ISIbyRate.(state).Xbins(1,:,1),ISIbyRate.(state).Ybins(1,:,1), ISIbyRate.(state).pop.(celltypes{tt})')
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
    imagesc(CV2byRate.(state).Xbins(1,:,1),CV2byRate.(state).Ybins(1,:,1), CV2byRate.(state).pop.(celltypes{tt})')
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
    imagesc(CV2byRate.(state).Xbins(1,:,1),[1 sum(ISIbyRate.celltypeidx.(celltypes{tt}))],...
        squeeze(CV2byRate.(state).pX(:,:,Cellmeanrate.sorts.(states{ss}).(celltypes{tt})))')
    hold on
    plot((Cellmeanrate.(states{ss})(Cellmeanrate.sorts.(states{ss}).(celltypes{tt}))),[1:sum(ISIbyRate.celltypeidx.(celltypes{tt}))],'w')
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
 NiceSave(['ISIbyRate_',(celltypes{tt})],figfolder,thisregion)

end
