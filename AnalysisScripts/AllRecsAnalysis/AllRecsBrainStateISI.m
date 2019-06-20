reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyBrainState'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
regions = {'vCTX','fCTX','CA1'};

for rr = 1:length(regions)
    [ISIstats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    BrainStateandSpikingAll = GetMatResults(figfolder,'SpikeStatsbyBrainState','baseNames',baseNames);
    BrainStateandSpikingAll = bz_CollapseStruct(BrainStateandSpikingAll);


%%
ISIbyPSS.(regions{rr}) = bz_CollapseStruct(BrainStateandSpikingAll.ISIbyPSS,'match','justcat',true);
ISIbytheta.(regions{rr}) = bz_CollapseStruct(BrainStateandSpikingAll.ISIbytheta,'match','justcat',true);
BShist.(regions{rr}) = bz_CollapseStruct(BrainStateandSpikingAll.BShist,1,'mean',true);
BShist.(regions{rr}).std = bz_CollapseStruct(BrainStateandSpikingAll.BShist,1,'std',true);
BShist_all.(regions{rr}) = bz_CollapseStruct(BrainStateandSpikingAll.BShist,1,'justcat',true);

celltypes = fieldnames(ISIbyPSS.(regions{rr}).celltypeidx);

for tt = 1:length(celltypes)
    ISIbyPSS.(regions{rr}).pop.(celltypes{tt}) = nanmean(ISIbyPSS.(regions{rr}).pYX(:,:,ISIbyPSS.(regions{rr}).celltypeidx.(celltypes{tt})),3);
    ISIbytheta.(regions{rr}).pop.(celltypes{tt}) = nanmean(ISIbytheta.(regions{rr}).pYX(:,:,ISIbytheta.(regions{rr}).celltypeidx.(celltypes{tt})),3);
end
end
%%
figure
for rr = 1:length(regions)
subplot(3,2,rr*2-1)
imagesc(BShist_all.(regions{rr}).bins(1,:),[0 1],BShist_all.(regions{rr}).ALL.PSS)
xlabel('PSS metric (Mean Norm)');ylabel(regions{rr})
subplot(3,2,rr*2)
imagesc(BShist_all.(regions{rr}).bins(1,:),[0 1],BShist_all.(regions{rr}).ALL.thratio)
xlabel('Theta metric (Mean Norm)');
end
%%
states = {'WAKEstate','NREMstate','REMstate'};
statecolors = {'k','b','r'};
%% ISI by PSS
figure
 for rr = 1:length(regions)
for tt = 1:length(celltypes)
subplot(4,3,tt*3-3+rr)
    imagesc(ISIbyPSS.(regions{rr}).Xbins(1,:,1),ISIbyPSS.(regions{rr}).Ybins(1,:,1), ISIbyPSS.(regions{rr}).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    if rr ==1
    ylabel({(celltypes{tt}),'ISI (s)'});
    LogScale('y',10,'exp',true)
    else
        set(gca,'yticklabel',[])
    end
    %title((celltypes{tt}))
   % colorbar
    if tt ==1 
        caxis([0 0.023])
        title(regions{rr})
    elseif tt==2
         caxis([0 0.032])
    end
    xlim(ISIbyPSS.(regions{rr}).Xbins(1,[1 end],1))
    if strcmp(regions{rr},'CA1')
        xlim([ISIbyPSS.(regions{rr}).Xbins(1,1,1) 1.9])
    end
end 



subplot(6,3,10+rr-1)
    hold on
    for ss = 1:3
        errorshade(BShist.(regions{rr}).bins,BShist.(regions{rr}).(states{ss}).PSS,...
            BShist.(regions{rr}).std.(states{ss}).PSS,BShist.(regions{rr}).std.(states{ss}).PSS,statecolors{ss},'scalar')
    plot(BShist.(regions{rr}).bins,BShist.(regions{rr}).(states{ss}).PSS,'color',statecolors{ss})
    end
    xlabel('PSS')
    
    box off
    axis tight
    xlim(ISIbyPSS.(regions{rr}).Xbins(1,[1 end],1))
        if strcmp(regions{rr},'CA1')
        xlim([ISIbyPSS.(regions{rr}).Xbins(1,1,1) 1.9])
    end
    


end
 NiceSave('ISIbyPSS',figfolder,regions{rr})
 
 
 %% ISI by Theta
figure
   
 for rr = 1:length(regions)
for tt = 1:length(celltypes)
subplot(4,3,tt*3-3+rr)
    imagesc(ISIbytheta.(regions{rr}).Xbins(1,:,1),ISIbytheta.(regions{rr}).Ybins(1,:,1), ISIbytheta.(regions{rr}).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
   % colorbar
    if rr ==1
    ylabel({(celltypes{tt}),'ISI (s)'});
    LogScale('y',10,'exp',true)
    else
        set(gca,'yticklabel',[])
    end
    if tt ==1 
        caxis([0 0.025])
        title(regions{rr})
    elseif tt==2
         caxis([0 0.035])
    end
    xlim(ISIbytheta.(regions{rr}).Xbins(1,[1 end],1))
end 

    
    
subplot(6,3,10+rr-1)
    hold on
    for ss = [1 3]
        errorshade(BShist.(regions{rr}).bins,BShist.(regions{rr}).(states{ss}).thratio,...
            BShist.(regions{rr}).std.(states{ss}).thratio,BShist.(regions{rr}).std.(states{ss}).thratio,statecolors{ss},'scalar')
    plot(BShist.(regions{rr}).bins,BShist.(regions{rr}).(states{ss}).thratio,'color',statecolors{ss})
    end
    xlabel('Theta Ratio')
        box off
    axis tight
    xlim(ISIbytheta.(regions{rr}).Xbins(1,[1 end],1))

 end
 NiceSave('ISIbyTheta',figfolder,regions{rr})
