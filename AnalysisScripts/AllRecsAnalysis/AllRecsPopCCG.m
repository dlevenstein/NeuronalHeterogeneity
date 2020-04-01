%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/PopCCGAnalysis'];

% datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
% datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
% datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
% datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
regions = {'CA1'};

%%
for rr = 1:length(regions)
    
    [PopCCGAll,baseNames] = GetMatResults(figfolder,'PopCCGAnalysis');
    PopCCGAll = bz_CollapseStruct(PopCCGAll);
    %thisregion = 'fCTX';
    
    ISICCG.(regions{rr}) = bz_CollapseStruct(PopCCGAll.ISICCG,'match','justcat',true);
    popCCG.(regions{rr}) = bz_CollapseStruct(PopCCGAll.popCCG,'match','justcat',true);
    CellClass.(regions{rr}) = bz_CollapseStruct(PopCCGAll.CellClass,'match','justcat',true);
   % ISIStats.(regions{rr}) = bz_CollapseStruct(PopCCGAll.ISIStats,'match','justcat',true);
end

%%
states = {'WAKEstate','NREMstate','REMstate'};
celltypes = {'pE','pI'};
for ss = 1:3
for tt = 1:2
    for tt2 = 1:2
        ISICCG.(regions{rr}).(states{ss}).popmean.(celltypes{tt}).(celltypes{tt2}) = ...
            nanmean(ISICCG.(regions{rr}).(states{ss}).(celltypes{tt})(:,:,CellClass.(regions{rr}).(celltypes{tt2})),3);
    end
end
end

%%
figure
for tt = 1:2
    for tt2 = 1:2
        subplot(2,2,(tt-1)*2+tt2)
            imagesc(ISICCG.(regions{rr}).t_ccg,...
                ISIStats.ISIhist.logbins,ISICCG.(regions{rr}).(states{ss}).popmean.(celltypes{tt}).(celltypes{tt2})')
            hold on
            %plot(ISIStats.ISIhist.(state).log(cc,:),ISIStats.ISIhist.logbins,'k')
            LogScale('y',10)
            
            if tt == 1
                ColorbarWithAxis([0 2.5],'E Pop Rate')
            else 
                ColorbarWithAxis([0 30],'I Pop Rate')
            end
            xlabel(['t lag (s) - ',(celltypes{tt2})]);ylabel('ISI (s)')
            
    end
    
end





%%


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


%% Figure: CSPKRWATE
for rr = 1:length(regions)
for tt=1:2
figure
for ss = 1:3
    state = states{ss};
%for tt = 1:length(celltypes)
subplot(4,3,ss)
    s =imagesc(ISIbyCspkRate.(regions{rr}).(state).Xbins(1,:,1),ISIbyCspkRate.(regions{rr}).(state).Ybins(1,:,1),...
        ISIbyCspkRate.(regions{rr}).(state).pop.(celltypes{tt})');
    %hold on
       % alpha(s,single(ISIoccupancy.(regions{rr}).(state).mednormhist(:,sorts.(regions{rr}).(statenames{ss}).ratebyclass)'~=0))

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

 NiceSave(['ISIbycBinRate_',(celltypes{tt})],figfolder,(regions{rr}))
end
end



%% Figure: SpikeRate
for rr = 1:length(regions)
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
 NiceSave(['ISIbyRate_',(celltypes{tt})],figfolder,(regions{rr}))

end
end


%% Figure: Excitatory Cells Compare
tt = 1
figure
for rr = 1:length(regions)
    for ss = 1:length(states)
subplot(3,4,rr + (ss-1)*4)
    s = imagesc(CV2byCspkRate.(regions{rr}).(states{ss}).Xbins(1,:,1),[1 sum(ISIbyCspkRate.(regions{rr}).celltypeidx.(celltypes{tt}))],...
        squeeze(CV2byCspkRate.(regions{rr}).(states{ss}).pX(:,:,Cellmeanrate.(regions{rr}).sorts.(states{ss}).(celltypes{tt})))')
        alpha(s,single(squeeze(CV2byCspkRate.(regions{rr}).(states{ss}).pX(:,:,Cellmeanrate.(regions{rr}).sorts.(states{ss}).(celltypes{tt})))'~=0))

    hold on
    plot(log10(Cellmeanrate.(regions{rr}).(states{ss})(Cellmeanrate.(regions{rr}).sorts.(states{ss}).(celltypes{tt}))),[1:sum(ISIbyCspkRate.(regions{rr}).celltypeidx.(celltypes{tt}))],'w')
    axis xy
    LogScale('x',10)
    if rr==1
    ylabel({(states{ss}),'Cell'});
    end
    xlabel('Rate')
    %title((celltypes{tt}))
    %colorbar
    if ss == 1
        title((regions{rr}))
    end
    if tt ==1 
        caxis([0 0.1])
    elseif tt==2
         caxis([0 0.15])
    end
    
    set(gca,'yticklabel',[])
    
    
    end
    
end

 NiceSave(['RateDist_',(celltypes{tt})],figfolder,[])