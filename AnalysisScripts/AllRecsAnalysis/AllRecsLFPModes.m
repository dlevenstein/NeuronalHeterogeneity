reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISIModeLFPAnalysis_HPC'];


[LFPModeALL,baseNames] = GetMatResults(figfolder,'ISIModeLFPAnalysis_HPC');
LFPModeALL = bz_CollapseStruct(LFPModeALL);


%%
ModalLFPModulation = bz_CollapseStruct(LFPModeALL.ModalLFPModulation,'match','justcat',true);


%%
celltypes = {'pE','pI'};
THlabels = {'hiThetastate','loThetastate'};

for tt = 1:length(celltypes)
    ISIbytheta.pop.(celltypes{tt}) = nanmean(ISIbytheta.Dist.pYX(:,:,TH_ISIstats.cellinfo.celltype.(celltypes{tt})),3);
    ISIbythetaphase.pop.(celltypes{tt}) = nanmean(ISIbythetaphase.shiftDist.pYX(:,:,TH_ISIstats.cellinfo.celltype.(celltypes{tt})),3);

    ISIbythetaphase.celltypeidx.(celltypes{tt}) = TH_ISIstats.cellinfo.celltype.(celltypes{tt});
    ISIbytheta.celltypeidx.(celltypes{tt}) = TH_ISIstats.cellinfo.celltype.(celltypes{tt});
    
    for hilo = 1:2
        TH_ISIstats.meandists.(THlabels{hilo}).(celltypes{tt}).ISIdist = nanmean(TH_ISIstats.ISIhist.(THlabels{hilo}).log(TH_ISIstats.cellinfo.celltype.(celltypes{tt}),:),1);
        TH_ISIstats.meandists.(THlabels{hilo}).(celltypes{tt}).Return = nanmean(TH_ISIstats.ISIhist.(THlabels{hilo}).return(:,:,TH_ISIstats.cellinfo.celltype.(celltypes{tt})),3);

    end
end





%%
ASfreq = repmat(ModalLFPModulation.freq(1,:)',[1,499,5]);
allASweight = repmat(ModalLFPModulation.ASweight,[150,1,1]);
%%
weightthresh = 0.01;
[ meanZ,N,Xbins,Ybins ] = ConditionalHist3( log10(ASfreq(allASweight>weightthresh)),...
    ModalLFPModulation.ASlogRates(allASweight>weightthresh),...
    ModalLFPModulation.ASModulation(allASweight>weightthresh),...
    'minXY',1,'numXbins',150,'numYbins',150);

%%
figure
imagesc(Xbins,Ybins,meanZ')
alpha((N'./max(N(:))).*2)
hold on
UnityLine

LogScale('xy',10)
xlabel('LFP Frequency (Hz)');ylabel('ISI Mode Rate (Hz)')

ColorbarWithAxis([-0.0 0.1],'Mode-Power Correlation','inclusive',{'<','>'})
crameri('vik','pivot',0)
axis xy
%%
numcells = length(ModalLFPModulation.UID);
figure
subplot(2,2,1)
% hold on
% for cc = 1:numcells
% for rr = 1:5
%     modeoccupancy = ModalLFPModulation.ASweight(1,cc,rr);
%     if (modeoccupancy<0.05)
%         continue
%     end
%     %showwhichmodes = (AllFConditionalISIModes(cc).AScorr_p(rr,:)<1); %& ...
%         %showwhichmodes = (modeoccupancy>0.05);
% 
%         
%     scatter(log10(ModalLFPModulation.freq(1,:)),...
%         ModalLFPModulation.ASlogRates(:,cc,rr),...
%         5*modeoccupancy,ModalLFPModulation.ASModulation(:,cc,rr),'filled')
% end
% end
% 
% axis tight
% UnityLine
% LogScale('xy',10)
% ColorbarWithAxis([-0.05 0.15],'Mode-Power Correlation','inclusive',{'<','>'})
% crameri('vik','pivot',0)
% xlabel('LFP frequency (Hz)');ylabel('ISI Mode Rate (Hz)')
imagesc(Xbins,Ybins,meanZ')
alpha((N'./max(N(:))).*2)
hold on
UnityLine

xlabel('LFP Frequency (Hz)');ylabel('ISI Mode Rate (Hz)')

ColorbarWithAxis([-0.0 0.1],'Mode-Power Correlation','inclusive',{'<','>'})
crameri('vik','pivot',0)
axis xy
xlim([0 2.5]);ylim([0 2.5])
LogScale('xy',10,'nohalf',true)

subplot(4,2,7)
hold on
plot(log10(ModalLFPModulation.freq(1,:)),nanmean(ModalLFPModulation.MutInf,2),'color','k','linewidth',1)
% errorshade(log10(ModalLFPModulation.freq(1,:)),nanmean(ModalLFPModulation.MutInf,2),...
%     nanstd(ModalLFPModulation.MutInf,[],2),nanstd(ModalLFPModulation.MutInf,[],2),'k','scalar')
LogScale('x',10)
xlabel('LFP Frequency (Hz)')
ylabel('MI[ISI;Power]')
colorbar
axis tight


subplot(4,2,8)
hold on

plot(log10(ModalLFPModulation.freq),-nanmean(ModalLFPModulation.GSModulation,2),'color','k','linewidth',2)
plot(xlim(gca),[0 0],'k--')
LogScale('x',10)
xlabel('LFP Frequency (Hz)')
ylabel('AR-Power Corr')
colorbar
axis tight


NiceSave('ModeLFPFreqMod',figfolder,[])






%%
phasex = linspace(-pi,3*pi,100);

histcolors.loThetastate = flipud(gray);
histcolors.hiThetastate = makeColorMap([1 1 1],[0.8 0 0]);



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
    bz_piTickLabel('x')
    bz_AddRightRateAxis
    %colormap(gca,histcolors.loThetastate)
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
    xlim([0 1])
    
    bz_AddRightRateAxis
    %colormap(gca,histcolors.loThetastate)
end 

for tt = 1:length(celltypes) 
    
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

    
%    plot(ThetaISImodes.GSrate,...
%        ThetaISImodes.GSModulation,'.','color',GScolor)
   plot(mean(ThetaISImodes.GSlogRates,2),...
       ThetaISImodes.GSModulation,'.','color',GScolor)

    axis tight
    box off
        plot(xlim(gca),[0 0],'--','color',[0.5 0.5 0.5])
        LogScale('x',10)
        xlabel('Mode Rate (Hz)')
        ylabel('Weight-Power Corr')
    
    
%subplot(2,2,3)
  %  plot([1:10],log10(ThetaISImodes.GSlogRates),'.')
    
  
 NiceSave('ThetaMod',figfolder,[])