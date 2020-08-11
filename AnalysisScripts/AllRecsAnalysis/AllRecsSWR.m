reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ThetaISIAnalysis/Old_ThetafromSS'];
%Loading theta for cell class.. this is stupid. Add to SWR when adding
%modal fit

[ThetaALL,baseNames] = GetMatResults(figfolder,'ThetaISIAnalysis');
ThetaALL = bz_CollapseStruct(ThetaALL);
TH_ISIstats = bz_CollapseStruct(ThetaALL.TH_ISIstats,'match','justcat',true);
%%
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SharpWaveISIAnalysis'];


[SharpWaveALL,baseNames] = GetMatResults(figfolder,'SharpWaveISIAnalysis');
SharpWaveALL = bz_CollapseStruct(SharpWaveALL);

%%
PeriSWISIDist = bz_CollapseStruct(SharpWaveALL.PeriSWISIDist,3','justcat',true);
SW_ISIstats = bz_CollapseStruct(SharpWaveALL.SW_ISIstats,'match','justcat',true);

SW_ISIstats.cellinfo.celltype = TH_ISIstats.cellinfo.celltype;






%%
celltypes = {'pE','pI'};
swrlabels = {'SWR','iSWR'};
cellcolor = {'k','r'};
%THlabels = {'hiThetastate','loThetastate'};

for tt = 1:length(celltypes)
    PeriSWISIDist.pop.(celltypes{tt}).pYX = nanmean(PeriSWISIDist.cells.pYX(:,:,SW_ISIstats.cellinfo.celltype.(celltypes{tt})),3);
    PeriSWISIDist.pop.(celltypes{tt}).rate = nanmean((PeriSWISIDist.cells.rate(:,:,SW_ISIstats.cellinfo.celltype.(celltypes{tt}))),3);
    for ss = 1:length(swrlabels)
        SW_ISIstats.meandists.(swrlabels{ss}).(celltypes{tt}).ISIdist = nanmean(SW_ISIstats.ISIhist.(swrlabels{ss}).log(SW_ISIstats.cellinfo.celltype.(celltypes{tt}),:),1);
        SW_ISIstats.meandists.(swrlabels{ss}).(celltypes{tt}).return = nanmean(SW_ISIstats.ISIhist.(swrlabels{ss}).return(:,:,SW_ISIstats.cellinfo.celltype.(celltypes{tt})),3);
    end
end

%%
NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);

figure
for tt = 1:length(celltypes)
subplot(2,2,2.*(tt-1)+1)
    imagesc(PeriSWISIDist.pop.(celltypes{tt}).Xbins(1,:,1),...
        PeriSWISIDist.pop.(celltypes{tt}).Ybins(1,:,1),PeriSWISIDist.pop.(celltypes{tt}).pYX')
    hold on
    plot(PeriSWISIDist.pop.(celltypes{tt}).Xbins(1,:,1),...
        log10(1./PeriSWISIDist.pop.(celltypes{tt}).rate),cellcolor{tt},'linewidth',1)
    axis tight
    plot([0 0],ylim(gca),'w--')
    LogScale('y',10,'nohalf',true)
    xlabel('t (s) - relative to SW');ylabel('ISI (s)')
    title(celltypes{tt})
    xlim([-0.2 0.2])
    bz_AddRightRateAxis
    
    


    
    subplot(4,4,tt+6)
        plot(SW_ISIstats.ISIhist.logbins(1,:),SW_ISIstats.meandists.SWR.(celltypes{tt}).ISIdist,'k')
        hold on
        plot(SW_ISIstats.ISIhist.logbins(1,:),SW_ISIstats.meandists.iSWR.(celltypes{tt}).ISIdist,'r')
        LogScale('x',10,'exp',true)
        title(celltypes{tt})
        box off 
        
    for ss = 1:length(swrlabels)
        subplot(4,4,10+tt+(ss-1)*4)
            imagesc(SW_ISIstats.ISIhist.logbins(1,:),SW_ISIstats.ISIhist.logbins(1,:),...
            SW_ISIstats.meandists.(swrlabels{ss}).(celltypes{tt}).return)
            LogScale('xy',10,'exp',true)
            axis xy
            colormap(gca,NREMhistcolors)

    end
end

NiceSave('SW_ISIStats',figfolder,[])




%%
phasex = linspace(-pi,3*pi,100);

histcolors.loThetastate = flipud(gray);
histcolors.hiThetastate = makeColorMap([1 1 1],[0.8 0 0]);



figure

for tt = 1:2
subplot(4,3,tt*3-2+6)
    imagesc(ISIbythetaphase.Dist.Xbins(1,:,1),ISIbythetaphase.Dist.Ybins(1,:,1),...
        ISIbythetaphase.pop.(celltypes{tt})')
    hold on
    imagesc(ISIbythetaphase.Dist.Xbins(1,:,1)+2*pi,ISIbythetaphase.Dist.Ybins(1,:,1), ...
        ISIbythetaphase.pop.(celltypes{tt})')
    %axis xy
    plot(phasex,-cos(phasex),'k')
    
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    %axis xy
    %xlim([-pi 3*pi])
    LogScale('y',10,'exp',true)
    ylabel('ISI (s)');xlabel('Theta Phase')
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    xlim([-pi 3*pi])
    plot(xlim(gca),log10(1/8).*[1 1],'w--')
    bz_piTickLabel('x')
    bz_AddRightRateAxis
    %colormap(gca,histcolors.loThetastate)
end 

for tt = 1:2
subplot(4,3,tt*3-2)
    imagesc(ISIbySWR.Dist.Xbins(1,:,1),ISIbySWR.Dist.Ybins(1,:,1), ...
        ISIbySWR.pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    %axis xy
    LogScale('y',10,'exp',true)
    ylabel('ISI (s)');xlabel('Theta Power')
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    xlim([0 1])
    
    bz_AddRightRateAxis
    %colormap(gca,histcolors.loThetastate)
end 

for tt = 1:2
    
    subplot(4,4,tt+6)
        plot(TH_ISIstats.ISIhist.logbins(1,:),TH_ISIstats.meandists.hiThetastate.(celltypes{tt}).ISIdist,'r')
        hold on
        plot(TH_ISIstats.ISIhist.logbins(1,:),TH_ISIstats.meandists.loThetastate.(celltypes{tt}).ISIdist,'k')
        LogScale('x',10,'exp',true)
        title(celltypes{tt})
        box off 
        
    for ss = 1:length(THlabels)
        subplot(4,4,10+tt+(ss-1)*4)
        
            imagesc(TH_ISIstats.ISIhist.logbins(1,:),TH_ISIstats.ISIhist.logbins(1,:),...
            TH_ISIstats.meandists.(THlabels{ss}).(celltypes{tt}).Return)
            LogScale('xy',10,'exp',true)
            axis xy
            colormap(gca,histcolors.(THlabels{ss}))
            if ss==2
                xlabel('ISI (s)')
            end
            if tt ==1
                ylabel('ISI_n_+_1 (s)')
            end
    end
end

NiceSave('TH_ISIstats',figfolder,[])

%%
figure

for tt = 3:4
subplot(4,3,(tt-2)*3-2+6)
    imagesc(ISIbythetaphase.Dist.Xbins(1,:,1),ISIbythetaphase.Dist.Ybins(1,:,1),...
        ISIbythetaphase.pop.(celltypes{tt})')
    hold on
    imagesc(ISIbythetaphase.Dist.Xbins(1,:,1)+2*pi,ISIbythetaphase.Dist.Ybins(1,:,1), ...
        ISIbythetaphase.pop.(celltypes{tt})')
    %axis xy
    plot(phasex,-cos(phasex),'k')
    
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    %axis xy
    %xlim([-pi 3*pi])
    LogScale('y',10,'exp',true)
    ylabel('ISI (s)');xlabel('Theta Phase')
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    xlim([-pi 3*pi])
    plot(xlim(gca),log10(1/8).*[1 1],'w--')
    bz_piTickLabel('x')
    bz_AddRightRateAxis
    %colormap(gca,histcolors.loThetastate)
end 

for tt = 3:4
subplot(4,3,(tt-2)*3-2)
    imagesc(ISIbySWR.Dist.Xbins(1,:,1),ISIbySWR.Dist.Ybins(1,:,1), ...
        ISIbySWR.pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    %axis xy
    LogScale('y',10,'exp',true)
    ylabel('ISI (s)');xlabel('Theta Power')
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    xlim([0 1])
    
    bz_AddRightRateAxis
    %colormap(gca,histcolors.loThetastate)
end 

for tt = 3:4
    
    subplot(4,4,(tt-2)+6)
        plot(TH_ISIstats.ISIhist.logbins(1,:),TH_ISIstats.meandists.hiThetastate.(celltypes{tt}).ISIdist,'r')
        hold on
        plot(TH_ISIstats.ISIhist.logbins(1,:),TH_ISIstats.meandists.loThetastate.(celltypes{tt}).ISIdist,'k')
        LogScale('x',10,'exp',true)
        title(celltypes{tt})
        box off 
        
    for ss = 1:length(THlabels)
        subplot(4,4,10+(tt-2)+(ss-1)*4)
        
            imagesc(TH_ISIstats.ISIhist.logbins(1,:),TH_ISIstats.ISIhist.logbins(1,:),...
            TH_ISIstats.meandists.(THlabels{ss}).(celltypes{tt}).Return)
            LogScale('xy',10,'exp',true)
            axis xy
            colormap(gca,histcolors.(THlabels{ss}))
            if ss==2
                xlabel('ISI (s)')
            end
            if tt ==3
                ylabel('ISI_n_+_1 (s)')
            end
    end
end

NiceSave('TH_ISIstats_hilo',figfolder,[])

%%
GScolor = [0.6 0.4 0];

figure
%subplot(2,2,1)
   % hist(ThetaISImodes.GSModulation)
subplot(2,2,2)
   plot(mean(ThetaISImodes.GSlogRates,2),...
       ThetaISImodes.GSModulation,'.','color',GScolor)
hold on
for rr = 1:5
    scatter(ThetaISImodes.ASlogRates(:,rr),...
        ThetaISImodes.ASModulation(:,rr),20*ThetaISImodes.ASweight(:,rr)+0.00001,'k','filled')
end    

    
%    plot(ThetaISImodes.GSrate,...
%        ThetaISImodes.GSModulation,'.','color',GScolor)


    axis tight
    box off
        plot(xlim(gca),[0 0],'--','color',[0.5 0.5 0.5])
        LogScale('x',10)
        xlabel('Mode Rate (Hz)')
        ylabel('Weight-Power Corr')
    
    
%subplot(2,2,3)
  %  plot([1:10],log10(ThetaISImodes.GSlogRates),'.')
    
  
 NiceSave('ThetaMod',figfolder,[])