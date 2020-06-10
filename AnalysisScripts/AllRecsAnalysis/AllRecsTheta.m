reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ThetaISIAnalysis'];


[ThetaALL,baseNames] = GetMatResults(figfolder,'ThetaISIAnalysis');
ThetaALL = bz_CollapseStruct(ThetaALL);

%%
ISIbythetaphase = bz_CollapseStruct(ThetaALL.ISIbythetaphase,'match','justcat',true);
ISIbytheta = bz_CollapseStruct(ThetaALL.ISIbytheta,'match','justcat',true);
TH_ISIstats = bz_CollapseStruct(ThetaALL.TH_ISIstats,'match','justcat',true);
ThetaISImodes = bz_CollapseStruct(ThetaALL.ThetaISImodes,'match','justcat',true);


%%
celltypes = {'pE','pI'};
THlabels = {'hiThetastate','loThetastate'};

for tt = 1:length(celltypes)
    ISIbytheta.pop.(celltypes{tt}) = nanmean(ISIbytheta.Dist.pYX(:,:,TH_ISIstats.cellinfo.celltype.(celltypes{tt})),3);
    ISIbythetaphase.pop.(celltypes{tt}) = nanmean(ISIbythetaphase.Dist.pYX(:,:,TH_ISIstats.cellinfo.celltype.(celltypes{tt})),3);

    ISIbythetaphase.celltypeidx.(celltypes{tt}) = TH_ISIstats.cellinfo.celltype.(celltypes{tt});
    ISIbytheta.celltypeidx.(celltypes{tt}) = TH_ISIstats.cellinfo.celltype.(celltypes{tt});
    
    for hilo = 1:2
        TH_ISIstats.meandists.(THlabels{hilo}).(celltypes{tt}).ISIdist = nanmean(TH_ISIstats.ISIhist.(THlabels{hilo}).log(TH_ISIstats.cellinfo.celltype.(celltypes{tt}),:),1);
        TH_ISIstats.meandists.(THlabels{hilo}).(celltypes{tt}).Return = nanmean(TH_ISIstats.ISIhist.(THlabels{hilo}).return(:,:,TH_ISIstats.cellinfo.celltype.(celltypes{tt})),3);

    end
end




%%
phasex = linspace(-pi,3*pi,100);

figure

for tt = 1:length(celltypes)
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
end 

for tt = 1:length(celltypes)
subplot(4,3,tt*3-2)
    imagesc(ISIbytheta.Dist.Xbins(1,:,1),ISIbytheta.Dist.Ybins(1,:,1), ...
        ISIbytheta.pop.(celltypes{tt})')
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
    xlim(ISIbytheta.Dist.Xbins(1,[1 end],1))
end 

for tt = 1:length(celltypes) 
    
    subplot(4,4,tt+6)
        plot(TH_ISIstats.ISIhist.logbins(1,:),TH_ISIstats.meandists.hiThetastate.(celltypes{tt}).ISIdist,'k')
        hold on
        plot(TH_ISIstats.ISIhist.logbins(1,:),TH_ISIstats.meandists.loThetastate.(celltypes{tt}).ISIdist,'r')
        LogScale('x',10,'exp',true)
        title(celltypes{tt})
        box off 
        
    for ss = 1:length(THlabels)
        subplot(4,4,10+tt+(ss-1)*4)
            imagesc(TH_ISIstats.ISIhist.logbins(1,:),TH_ISIstats.ISIhist.logbins(1,:),...
            TH_ISIstats.meandists.(THlabels{ss}).(celltypes{tt}).Return)
            LogScale('xy',10,'exp',true)
            axis xy
    end
end

NiceSave('TH_ISIstats',figfolder,[])


%%
GScolor = [0.6 0.4 0];

figure
%subplot(2,2,1)
   % hist(ThetaISImodes.GSModulation)
subplot(2,2,2)
hold on
for rr = 1:5
    scatter(ThetaISImodes.ASlogRates(:,rr),...
        ThetaISImodes.ASModulation(:,rr),20*ThetaISImodes.ASweight(:,rr)+0.00001,'k','filled')
end    

    
  %  plot(ThetaISImodes.GSrate,...
  %      ThetaISImodes.GSModulation,'.')
   plot(mean(ThetaISImodes.GSlogRates,2),...
       ThetaISImodes.GSModulation,'.','color',GScolor)

    axis tight
    box off
        plot(xlim(gca),[0 0],'k--')
        LogScale('x',10)
        xlabel('Mode Rate (Hz)')
        ylabel('Weight-Power Corr')
    
    
%subplot(2,2,3)
  %  plot([1:10],log10(ThetaISImodes.GSlogRates),'.')
    
    %%
 NiceSave('ThetaMod',figfolder,[])