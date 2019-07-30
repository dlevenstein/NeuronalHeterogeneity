reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyBrainState'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
regions = {'THAL','vCTX','fCTX','CA1'};

for rr = 1:length(regions)
    [ISIstats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    BrainStateandSpikingAll = GetMatResults(figfolder,'SpikeStatsbyBrainState','baseNames',baseNames);
    BrainStateandSpikingAll = bz_CollapseStruct(BrainStateandSpikingAll);



    ISIbyPSS.(regions{rr}) = bz_CollapseStruct(BrainStateandSpikingAll.ISIbyPSS,'match','justcat',true);
    ISIbytheta.(regions{rr}) = bz_CollapseStruct(BrainStateandSpikingAll.ISIbytheta,'match','justcat',true);
    BShist.(regions{rr}) = bz_CollapseStruct(BrainStateandSpikingAll.BShist,1,'mean',true);
    BShist.(regions{rr}).std = bz_CollapseStruct(BrainStateandSpikingAll.BShist,1,'std',true);
    BShist_all.(regions{rr}) = bz_CollapseStruct(BrainStateandSpikingAll.BShist,1,'justcat',true);

    statehists.(regions{rr}) = bz_CollapseStruct(BrainStateandSpikingAll.statehists,3,'mean',true);
    celltypes = fieldnames(ISIbyPSS.(regions{rr}).celltypeidx);
    
    %How many cells are contributing?
    nspkthresh = 50;
    ncellthresh = 200;
    %sum(ISIbytheta.(regions{rr}).Xhist>nspkthresh,3)
    %sum(ISIbyPSS.(regions{rr}).Xhist>nspkthresh,3)
    ISIbyPSS.(regions{rr}).pYX(sum(ISIbyPSS.(regions{rr}).Xhist>nspkthresh,3)<ncellthresh,:,:)=nan;
    ISIbytheta.(regions{rr}).pYX(sum(ISIbytheta.(regions{rr}).Xhist>nspkthresh,3)<ncellthresh,:,:)=nan;
    
    if rr ==4
       ISIbytheta.(regions{rr}).WAKE.pYX(sum(ISIbytheta.(regions{rr}).WAKE.Xhist>nspkthresh,3)<ncellthresh,:,:)=nan;
       ISIbytheta.(regions{rr}).REM.pYX(sum(ISIbytheta.(regions{rr}).REM.Xhist>nspkthresh,3)<ncellthresh,:,:)=nan;
    end
    
    for tt = 1:length(celltypes)
        ISIbyPSS.(regions{rr}).pop.(celltypes{tt}) = nanmean(ISIbyPSS.(regions{rr}).pYX(:,:,ISIbyPSS.(regions{rr}).celltypeidx.(celltypes{tt})),3);
        ISIbytheta.(regions{rr}).pop.(celltypes{tt}) = nanmean(ISIbytheta.(regions{rr}).pYX(:,:,ISIbytheta.(regions{rr}).celltypeidx.(celltypes{tt})),3);
        
        if rr ==4
        ISIbytheta.(regions{rr}).WAKE.pop.(celltypes{tt}) = nanmean(ISIbytheta.(regions{rr}).WAKE.pYX(:,:,ISIbytheta.(regions{rr}).celltypeidx.(celltypes{tt})),3);
        ISIbytheta.(regions{rr}).REM.pop.(celltypes{tt}) = nanmean(ISIbytheta.(regions{rr}).REM.pYX(:,:,ISIbytheta.(regions{rr}).celltypeidx.(celltypes{tt})),3);
        end
    end
end
%%
figure
for rr = 1:length(regions)
subplot(4,2,rr*2-1)
imagesc(BShist_all.(regions{rr}).bins(1,:),[0 1],BShist_all.(regions{rr}).ALL.PSS)
xlabel('PSS metric (Mean Norm)');ylabel(regions{rr})
subplot(4,2,rr*2)
imagesc(BShist_all.(regions{rr}).bins(1,:),[0 1],BShist_all.(regions{rr}).ALL.thratio)
xlabel('Theta metric (Mean Norm)');
end
%%
states = {'WAKEstate','NREMstate','REMstate'};
statecolors = {[0 0 0],[0 0 1],[1 0 0]};
%% ISI by PSS
figure
 for rr = 1:length(regions)
for tt = 1:length(celltypes)
subplot(4,4,tt*4-4+rr+8)
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
        %caxis([0 0.023])
        set(gca,'xticklabel',[])
    elseif tt==2
         %caxis([0 0.032])
         xlabel('PSS')
    end
        xlim([0 1]) 
    set(gca,'xtick',[0 1]);
    if strcmp(regions{rr},'CA1')
       % xlim([ISIbyPSS.(regions{rr}).Xbins(1,1,1) 1.9])
    end
end 



% subplot(6,3,10+rr-1)
%     hold on
%     for ss = 1:3
%         errorshade(BShist.(regions{rr}).bins,BShist.(regions{rr}).(states{ss}).PSS,...
%             BShist.(regions{rr}).std.(states{ss}).PSS,BShist.(regions{rr}).std.(states{ss}).PSS,statecolors{ss},'scalar')
%     plot(BShist.(regions{rr}).bins,BShist.(regions{rr}).(states{ss}).PSS,'color',statecolors{ss})
%     end
%     xlabel('PSS (peak norm)')
%     
%     box off
%     axis tight
%     xlim(ISIbyPSS.(regions{rr}).Xbins(1,[1 end],1))
%         if strcmp(regions{rr},'CA1') | strcmp(regions{rr},'vCTX')
%         %xlim([ISIbyPSS.(regions{rr}).Xbins(1,1,1) 1.9])
%         end
    
scale = 5;
imagescale = 200;
 subplot(2,4,rr)
    %crameri grayC
    hold on
    %imagesc(statehists.(regions{rr}).PSSbins,statehists.(regions{rr}).EMGbins,statehists.(regions{rr}).PSS')

    for ss = 1:3
        plotcolor = cat(3,statecolors{ss}(1).*ones(size(statehists.(regions{rr}).(states{ss}).PSS))',...
            statecolors{ss}(2).*ones(size(statehists.(regions{rr}).(states{ss}).PSS))',...
            statecolors{ss}(3).*ones(size(statehists.(regions{rr}).(states{ss}).PSS))');
        h = image(statehists.(regions{rr}).PSSbins,statehists.(regions{rr}).EMGbins,plotcolor);
        set(h,'AlphaData',imagescale*statehists.(regions{rr}).(states{ss}).PSS')
    end

    axis tight
    for ss = 1:3
        errorshade(BShist.(regions{rr}).bins,1+scale*BShist.(regions{rr}).(states{ss}).PSS,...
            scale*BShist.(regions{rr}).std.(states{ss}).PSS,scale*BShist.(regions{rr}).std.(states{ss}).PSS,statecolors{ss},'scalar')
    plot(BShist.(regions{rr}).bins,1+scale*BShist.(regions{rr}).(states{ss}).PSS,'color',statecolors{ss})
    end
    caxis([0 5e-3])
    set(gca,'xtick',[0 1])
    set(gca,'ytick',[0 1])
    title(regions{rr})
    if rr == 1
        ylabel('EMG')
    else
        set(gca,'yticklabel',[])
    end
   % set(gca,'xticklabel',[])
    axis tight
    %colorbar
    %xlabel('PSS')
    xlim([0 1]) 
    set(gca,'xtick',[0 1])
        

end
 NiceSave('ISIbyPSS',figfolder,[])
 
 
 %% ISI by Theta
figure
   
 for rr = 1:length(regions)
for tt = 1:length(celltypes)
subplot(4,4,tt*4-4+rr+8)
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
        %caxis([0 0.025])
        set(gca,'xticklabel',[])
        %title(regions{rr})
    elseif tt==2
         %caxis([0 0.035])
         xlabel('Theta')
    end
        xlim([0 1]) 
    set(gca,'xtick',[0 1])
end 

    
    

    
scale = 5;
 subplot(2,4,rr)
    %crameri grayC
    hold on
    for ss = [1 3]
        plotcolor = cat(3,statecolors{ss}(1).*ones(size(statehists.(regions{rr}).(states{ss}).theta))',...
            statecolors{ss}(2).*ones(size(statehists.(regions{rr}).(states{ss}).theta))',...
            statecolors{ss}(3).*ones(size(statehists.(regions{rr}).(states{ss}).theta))');
        h = image(statehists.(regions{rr}).thetabins,statehists.(regions{rr}).EMGbins,plotcolor);
        set(h,'AlphaData',200*statehists.(regions{rr}).(states{ss}).theta')
    end
    
    for ss = [1 3]
        errorshade(BShist.(regions{rr}).bins,1+scale*BShist.(regions{rr}).(states{ss}).thratio,...
            scale*BShist.(regions{rr}).std.(states{ss}).thratio,scale*BShist.(regions{rr}).std.(states{ss}).thratio,statecolors{ss},'scalar')
    plot(BShist.(regions{rr}).bins,1+scale*BShist.(regions{rr}).(states{ss}).thratio,'color',statecolors{ss})
    end
    %caxis([0 5e-3])
    title(regions{rr})
   % set(gca,'xticklabel',[])
    %xlabel('Theta');
    if rr == 1
        ylabel('EMG')
    else
        set(gca,'yticklabel',[])
    end
    axis tight
    %colorbar
        xlim([0 1]) 
    set(gca,'xtick',[0 1]);set(gca,'ytick',[0 1])
    
 end
 NiceSave('ISIbyTheta',figfolder,[])

 %% All State Variables
 figure
 %cmap = crameri('grayC');
 colormap(gcf,crameri('grayC'))
  for rr = 1:length(regions)
   subplot(3,4,rr)
    hold on
    for ss = 1:3
        plotcolor = cat(3,statecolors{ss}(1).*ones(size(statehists.(regions{rr}).(states{ss}).PSSvtheta))',...
            statecolors{ss}(2).*ones(size(statehists.(regions{rr}).(states{ss}).PSSvtheta))',...
            statecolors{ss}(3).*ones(size(statehists.(regions{rr}).(states{ss}).PSSvtheta))');
        h = image(statehists.(regions{rr}).PSSbins,statehists.(regions{rr}).thetabins,plotcolor);
        set(h,'AlphaData',imagescale*statehists.(regions{rr}).(states{ss}).PSSvtheta')
    end
    
    xlabel('PSS');ylabel('Theta')
    title(regions{rr})
        xlim([0 1]) ;ylim([0 1])
    set(gca,'xtick',[0 1]);set(gca,'ytick',[0 1])
    
    %xlim(ISIbyPSS.(regions{rr}).Xbins(1,[1 end],1)) 
    
  subplot(3,4,4+rr)
    hold on
    for ss = 1:3
        plotcolor = cat(3,statecolors{ss}(1).*ones(size(statehists.(regions{rr}).(states{ss}).PSS))',...
            statecolors{ss}(2).*ones(size(statehists.(regions{rr}).(states{ss}).PSS))',...
            statecolors{ss}(3).*ones(size(statehists.(regions{rr}).(states{ss}).PSS))');
        h = image(statehists.(regions{rr}).PSSbins,statehists.(regions{rr}).EMGbins,plotcolor);
        set(h,'AlphaData',imagescale*statehists.(regions{rr}).(states{ss}).PSS')
    end
    
    xlabel('PSS');ylabel('EMG')
        xlim([0 1]) ;ylim([0 1])
    set(gca,'xtick',[0 1]);set(gca,'ytick',[0 1])
    %xlim(ISIbyPSS.(regions{rr}).Xbins(1,[1 end],1)) 
        
 
  subplot(3,4,8+rr)
    hold on
    for ss = [1 3]
        plotcolor = cat(3,statecolors{ss}(1).*ones(size(statehists.(regions{rr}).(states{ss}).theta))',...
            statecolors{ss}(2).*ones(size(statehists.(regions{rr}).(states{ss}).theta))',...
            statecolors{ss}(3).*ones(size(statehists.(regions{rr}).(states{ss}).theta))');
        h = image(statehists.(regions{rr}).thetabins,statehists.(regions{rr}).EMGbins,plotcolor);
        set(h,'AlphaData',200*statehists.(regions{rr}).(states{ss}).theta')
    end
    
    xlabel('Theta Ratio');ylabel('EMG')
        xlim([0 1]) ;ylim([0 1])
    set(gca,'xtick',[0 1]);set(gca,'ytick',[0 1])
    %xlim(ISIbytheta.(regions{rr}).Xbins(1,[1 end],1))   
  end
  
   NiceSave('StateSpace',figfolder,[])

   
   %% HPC WAKE
   WAKEREM = {'WAKE','REM'};
   figure
   
 for rr = 4
     
     for ww = 1:2
for tt = 1:length(celltypes)
subplot(6,4,tt*4-4+rr+8+(ww-1)*8)
    imagesc(ISIbytheta.(regions{rr}).Xbins(1,:,1),ISIbytheta.(regions{rr}).Ybins(1,:,1), ISIbytheta.(regions{rr}).(WAKEREM{ww}).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
   %colorbar
   %caxis([0 0.022])
    if rr ==4
    ylabel({(celltypes{tt}),'ISI (s)'});
    LogScale('y',10,'exp',true)
    else
        set(gca,'yticklabel',[])
    end
    
    if tt ==1 
        %caxis([0 0.025])
        set(gca,'xticklabel',[])
        %title(regions{rr})
    elseif tt==2
         %caxis([0 0.035])
         xlabel('Theta')
    end
        xlim([0 1]) 
    set(gca,'xtick',[0 1])
end 
     end
    
    

    
scale = 5;
 subplot(4,4,rr)
    %crameri grayC
    hold on
    for ss = [1 3]
        plotcolor = cat(3,statecolors{ss}(1).*ones(size(statehists.(regions{rr}).(states{ss}).theta))',...
            statecolors{ss}(2).*ones(size(statehists.(regions{rr}).(states{ss}).theta))',...
            statecolors{ss}(3).*ones(size(statehists.(regions{rr}).(states{ss}).theta))');
        h = image(statehists.(regions{rr}).thetabins,statehists.(regions{rr}).EMGbins,plotcolor);
        set(h,'AlphaData',500*statehists.(regions{rr}).(states{ss}).theta')
    end
    
    for ss = [1 3]
        errorshade(BShist.(regions{rr}).bins,1+scale*BShist.(regions{rr}).(states{ss}).thratio,...
            scale*BShist.(regions{rr}).std.(states{ss}).thratio,scale*BShist.(regions{rr}).std.(states{ss}).thratio,statecolors{ss},'scalar')
    plot(BShist.(regions{rr}).bins,1+scale*BShist.(regions{rr}).(states{ss}).thratio,'color',statecolors{ss})
    end
    %caxis([0 5e-3])
    %colorbar
    title(regions{rr})
   % set(gca,'xticklabel',[])
    %xlabel('Theta');
    if rr == 4
        ylabel('EMG')
    else
        set(gca,'yticklabel',[])
    end
    axis tight
    %colorbar
        xlim([0 1]) 
    set(gca,'xtick',[0 1]);set(gca,'ytick',[0 1])
    
 end
 NiceSave('ISIbyTheta_WAKEREM',figfolder,[])