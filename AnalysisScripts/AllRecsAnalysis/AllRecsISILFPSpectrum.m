reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISILFPSpectrumAnalysis']; 
%figfolder_BLA = [reporoot,'AnalysisScripts/AnalysisFigs/ISILFPSpectrumAnalysis_BLA']; 
%figfolder_PIR = [reporoot,'AnalysisScripts/AnalysisFigs/ISILFPSpectrumAnalysis_PIR']; 

[baseNames] = getDatasetBasenames();
regions = {'THAL','vCTX','fCTX','BLA','PIR','CA1'};
for rr = 1:length(regions)
    %Get baseNames
    if rr == 4 || rr == 5
        ISILFPSpectrumAnalysis_ALL = GetMatResults([figfolder,'_',regions{rr}],['ISILFPSpectrumAnalysis','_',regions{rr}]);
    else
        ISILFPSpectrumAnalysis_ALL = GetMatResults(figfolder,'ISILFPSpectrumAnalysis',...
            'baseNames',baseNames.(regions{rr}));
    end
    ISILFPSpectrumAnalysis_ALL = bz_CollapseStruct(ISILFPSpectrumAnalysis_ALL);
    
    MutInf.(regions{rr}) = bz_CollapseStruct(ISILFPSpectrumAnalysis_ALL.MutInf,'match','justcat',true);
    PSSConditionalISIDist.(regions{rr}) = bz_CollapseStruct(ISILFPSpectrumAnalysis_ALL.PSSConditionalISIDist,3,'justcat',true);
    PSSConditionalGamma.(regions{rr}) = bz_CollapseStruct(ISILFPSpectrumAnalysis_ALL.PSSConditionalGamma,'match','justcat');
    PSSConditionalGamma.(regions{rr}).modes = bz_CollapseStruct(PSSConditionalGamma.(regions{rr}).modes,3,'justcat',true);
    PSSConditionalGamma.(regions{rr}).dist = bz_CollapseStruct(PSSConditionalGamma.(regions{rr}).dist,3,'justcat',true);
    
    HiLowISIStats.(regions{rr}) = bz_CollapseStruct(ISILFPSpectrumAnalysis_ALL.HiLowISIStats,'match','justcat',true);
    GammaParms.(regions{rr}) = bz_CollapseStruct(ISILFPSpectrumAnalysis_ALL.GammaParms,'match','justcat',true);
    PSScorr.(regions{rr}) = bz_CollapseStruct(ISILFPSpectrumAnalysis_ALL.PSScorr,'match','justcat',true);
    PSShist.(regions{rr}).mean = bz_CollapseStruct(ISILFPSpectrumAnalysis_ALL.PSShist,3,'mean',true);
    PSShist.(regions{rr}).std = bz_CollapseStruct(ISILFPSpectrumAnalysis_ALL.PSShist,3,'std',true);
end

%%
states = {'WAKEstate','NREMstate','REMstate'};
celltypes = {'pE','pI'};
cellcolor = {'k','r'};
statecolors = {[0 0 0],[0 0 1],[1 0 0]};

for rr = 1:length(regions)
MutInf.(regions{rr}).freqs = logspace(log10(2),log10(128),150);
end
%%
for rr = 1:length(regions)
for ss = 1:3
    for tt = 1:length(celltypes)
        MeanMI.(regions{rr}).(states{ss}).(celltypes{tt}).Osci = ...
            nanmean((MutInf.(regions{rr}).(states{ss}).Osci(:,MutInf.(regions{rr}).CellClass.(celltypes{tt}))),2);
        MeanMI.(regions{rr}).(states{ss}).(celltypes{tt}).Power = ...
            nanmean((MutInf.(regions{rr}).(states{ss}).Power(:,MutInf.(regions{rr}).CellClass.(celltypes{tt}))),2);

        %Conditioned on spike AND power
        PSSConditionalISIDist.(regions{rr}).(states{ss}).(celltypes{tt}) = ...
            nanmean(PSSConditionalISIDist.(regions{rr}).(states{ss}).Dist.pYX(:,:,MutInf.(regions{rr}).CellClass.(celltypes{tt})),3);
        
        %Just conditioned on power
%         PSSConditionalISIDist.(states{ss}).(celltypes{tt}) = ...
%             nanmean(PSSConditionalISIDist.(states{ss}).Dist.XYprob(:,:,CellClass.(celltypes{tt})),3);
        
        PSSConditionalISIDist.(regions{rr}).(states{ss}).meanlograte.(celltypes{tt}) = ...
            nanmean(log10(PSSConditionalISIDist.(regions{rr}).(states{ss}).Dist.SpikeRate(:,:,MutInf.(regions{rr}).CellClass.(celltypes{tt}))),3);
    end
end
end

%%
hilow = {'HighPSSstate','LowPSSstate',};
for rr = 1:length(regions)
for ss = 1:3
    %ISIstatshist.(states{ss}) = HiLowISIStats.(states{ss}).ISIhist;
    for hl = 1:2
        meanISIhist.(regions{rr}).(states{ss}).(hilow{hl}).return =...
            mean(HiLowISIStats.(regions{rr}).(states{ss}).(hilow{hl}).return(:,:,MutInf.(regions{rr}).CellClass.pE),3);
        meanISIhist.(regions{rr}).(states{ss}).(hilow{hl}).logdist =...
            mean(HiLowISIStats.(regions{rr}).(states{ss}).(hilow{hl}).log(MutInf.(regions{rr}).CellClass.pE,:),1);
    end
end
end
meanISIhist.logbins = HiLowISIStats.(regions{rr}).(states{ss}).logbins(1,:);

%%
tt = 1;
for rr = 1:length(regions)
figure
for ss = 1:3

subplot(3,5,1+5.*(ss-1))
BoxAndScatterPlot({(MutInf.(regions{rr}).(states{ss}).PSS(MutInf.(regions{rr}).CellClass.(celltypes{tt})))},...
    'groupnumbers',-1)
if ss == 1
PSSlim = ylim(gca); PSSlim(1) = 0;
end
ylim(PSSlim)
box off
ylabel('MI - PSS')

subplot(3,5,[2:3]+5.*(ss-1))
% BoxAndScatterPlot({(MutInf.(regions{rr}).(states{ss}).PSS(MutInf.(regions{rr}).CellClass.(celltypes{tt})))},...
%     'groupnumbers',0,'labels','PSS')
plot(log2(MutInf.(regions{rr}).freqs),MeanMI.(regions{rr}).(states{ss}).(celltypes{tt}).Osci,'k')
hold on 
plot(log2(MutInf.(regions{rr}).freqs),MeanMI.(regions{rr}).(states{ss}).(celltypes{tt}).Power,'r:')

LogScale('x',2)

axis tight
ylim(PSSlim)
box off
ylabel('MI - Power')
xlabel('freq (Hz)')


subplot(3,5,[4:5]+5.*(ss-1))
imagesc(PSSConditionalISIDist.(regions{rr}).(states{ss}).Dist.Xbins(1,:,1),...
    PSSConditionalISIDist.(regions{rr}).(states{ss}).Dist.Ybins(1,:,1),...
    PSSConditionalISIDist.(regions{rr}).(states{ss}).(celltypes{tt})')
hold on
plot(PSSConditionalISIDist.(regions{rr}).(states{ss}).Dist.Xbins(1,:,1),...
    -PSSConditionalISIDist.(regions{rr}).(states{ss}).meanlograte.(celltypes{tt}),'r','linewidth',2)
LogScale('y',10,'exp',true,'nohalf',true)
xlabel('PSS (%ile)');ylabel('ISI (s)')
bz_AddRightRateAxis

end
NiceSave('ISIModSpectrum',figfolder,(regions{rr}))
end


%% Spectrum all recs

tt = 1;

figure
for rr = 1:length(regions)
for ss = 1:3

% subplot(3,5,1+5.*(ss-1))
% BoxAndScatterPlot({(MutInf.(regions{rr}).(states{ss}).PSS(MutInf.(regions{rr}).CellClass.(celltypes{tt})))},...
%     'groupnumbers',-1)

% ylim(PSSlim)
% box off
% ylabel('MI - PSS')

subplot(5,6,rr+(ss-1)*6)
BoxAndScatterPlot({(MutInf.(regions{rr}).(states{ss}).PSS(MutInf.(regions{rr}).CellClass.(celltypes{tt})))},...
    'groupnumbers',0.5,'labels','PSS','withpoints',false);
axis tight
plot(log2(MutInf.(regions{rr}).freqs),MeanMI.(regions{rr}).(states{ss}).(celltypes{tt}).Osci,'color',statecolors{ss},'linewidth',1)
plot(log2(MutInf.(regions{rr}).freqs),MeanMI.(regions{rr}).(states{ss}).(celltypes{tt}).Power,':','color',statecolors{ss})

LogScale('x',2)

%axis tight
if ss == 1
    title(regions{rr})
end
% PSSlim = ylim(gca); PSSlim(1) = 0;
% end
% ylim(PSSlim)
ylim([0 0.035]);xlim([0 log2(MutInf.(regions{rr}).freqs(end))])
box off
ylabel('MI - Power')
xlabel('freq (Hz)')

end
NiceSave('ISIModSpectrum',figfolder,'AllRegions')
end

%% Conditional spectrum all recs
tt = 1;

figure
for rr = 1:length(regions)
for ss = 1:2

% subplot(3,5,1+5.*(ss-1))
% BoxAndScatterPlot({(MutInf.(regions{rr}).(states{ss}).PSS(MutInf.(regions{rr}).CellClass.(celltypes{tt})))},...
%     'groupnumbers',-1)

% ylim(PSSlim)
% box off
% ylabel('MI - PSS')

subplot(4,6,rr*2+(-ss+1))
imagesc(PSSConditionalISIDist.(regions{rr}).(states{ss}).Dist.Xbins(1,:,1),...
    PSSConditionalISIDist.(regions{rr}).(states{ss}).Dist.Ybins(1,:,1),...
    PSSConditionalISIDist.(regions{rr}).(states{ss}).(celltypes{tt})')
hold on
LogScale('y',10,'exp',true,'nohalf',true)
xlabel('PSS (%ile)');
if ss == 2
    ylabel('ISI (s)')
    title(regions{rr})
else
    set(gca,'yticklabels',[])
    %bz_AddRightRateAxis
end

%


end
NiceSave('ConditionalISIPSS',figfolder,'AllRegions')
end



%%
for rr = 1:length(regions)
figure
for ss = 1:3
subplot(3,5,1+(ss-1)*5)
BoxAndScatterPlot({squeeze(-PSSConditionalGamma.(regions{rr}).modes.(states{ss}).GS_R)})
hold on
plot(xlim(gca),[0 0],'k--')
ylabel('AS Modulation')
box off

subplot(3,5,[2:3]+5.*(ss-1))
hold on
keepcells = PSSConditionalGamma.(regions{rr}).modes.(states{ss}).GScorr_p<=0.05;
%keepcells = true(size(keepcells))
for aa = 1:5
    %keepmodes = keepcells&mean(PSSConditionalGamma.modes.(states{ss}).ASweights(:,aa,:),1)>0.02;
    keepmodes = (PSSConditionalGamma.(regions{rr}).modes.(states{ss}).AScorr_p(:,aa,:))<=0.05;
    %keepmodes = true(size(keepmodes))
scatter(-PSSConditionalGamma.(regions{rr}).modes.(states{ss}).ASlogrates(1,aa,keepmodes),...
    log10(PSSConditionalGamma.(regions{rr}).modes.(states{ss}).ASCVs(1,aa,keepmodes)),...
    40*mean(PSSConditionalGamma.(regions{rr}).modes.(states{ss}).ASweights(:,aa,keepmodes),1)+eps,...
    squeeze(PSSConditionalGamma.(regions{rr}).modes.(states{ss}).AS_R(1,aa,keepmodes)),'filled')
end
% scatter(-PSSConditionalGamma.(regions{rr}).modes.(states{ss}).GSlogrates(1,1,keepcells),...
%     log10(PSSConditionalGamma.(regions{rr}).modes.(states{ss}).GSCVs(1,1,keepcells)),...
%     10,...
%     squeeze(PSSConditionalGamma.(regions{rr}).modes.(states{ss}).GS_R(1,1,keepcells)))
colorbar
axis tight
caxis([-0.1 0.1])
crameri('vik','pivot',0)
xlabel('Mean');ylabel('CV')


%tt = 1
subplot(3,5,[4:5]+5.*(ss-1))
imagesc(PSSConditionalISIDist.(regions{rr}).(states{ss}).Dist.Xbins(1,:,1),...
    PSSConditionalISIDist.(regions{rr}).(states{ss}).Dist.Ybins(1,:,1),...
    PSSConditionalISIDist.(regions{rr}).(states{ss}).pE')
end
NiceSave('GSASModPSS',figfolder,(regions{rr}))
end


%% AS/GS mod all regions

figure
for rr = 1:length(regions)
for ss = 1:3

subplot(4,6,rr+(ss-1)*6)
hold on
keepcells = PSSConditionalGamma.(regions{rr}).modes.(states{ss}).GScorr_p<=0.05;
%keepcells = true(size(keepcells))
for aa = 1:5
    %keepmodes = keepcells&mean(PSSConditionalGamma.modes.(states{ss}).ASweights(:,aa,:),1)>0.02;
    keepmodes = (PSSConditionalGamma.(regions{rr}).modes.(states{ss}).AScorr_p(:,aa,:))<=0.05;
    %keepmodes = true(size(keepmodes))
scatter(-PSSConditionalGamma.(regions{rr}).modes.(states{ss}).ASlogrates(1,aa,keepmodes),...
    log10(PSSConditionalGamma.(regions{rr}).modes.(states{ss}).ASCVs(1,aa,keepmodes)),...
    8*mean(PSSConditionalGamma.(regions{rr}).modes.(states{ss}).ASweights(:,aa,keepmodes),1)+eps,...
    squeeze(PSSConditionalGamma.(regions{rr}).modes.(states{ss}).AS_R(1,aa,keepmodes)),'filled')
end
scatter(-PSSConditionalGamma.(regions{rr}).modes.(states{ss}).GSlogrates(1,1,keepcells),...
    log10(PSSConditionalGamma.(regions{rr}).modes.(states{ss}).GSCVs(1,1,keepcells)),...
    2,...
    squeeze(PSSConditionalGamma.(regions{rr}).modes.(states{ss}).GS_R(1,1,keepcells)),'linewidth',0.1)
%colorbar

axis tight
caxis([-0.2 0.2])
crameri('vik','pivot',0)
%xlabel('Mean');

xlim([-3 1.7]);ylim([-2 0.5])
plot(xlim(gca),[0 0],'k--')
LogScale('x',10,'exp',true,'nohalf',true)
LogScale('y',10,'nohalf',true)
if rr == 1
    ylabel('CV')
else
    set(gca,'yticklabel',[])
end
if ss == 1
    title(regions{rr})
    set(gca,'xticklabel',[])
else
    xlabel('Mean (s)')
end


end
end
NiceSave('GSASModPSS',figfolder,'AllRegions')


%% In/Out Field: REturn Maps
histcolors = flipud(gray);
NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);
REMhistcolors = makeColorMap([1 1 1],[0.8 0 0]);
statecolormap = {histcolors,NREMhistcolors,REMhistcolors};

for rr =1:length(regions)
figure
for ss = 1:3

subplot(3,3,6+ss)
hold on
for hl = 1:2
    plot(meanISIhist.logbins,meanISIhist.(regions{rr}).(states{ss}).(hilow{hl}).logdist)
end
    axis tight
    box off
    LogScale('x',10,'exp',true,'nohalf',true)

%legend(cellISIStats.statenames{2:3},'location','southoutside')



for hl = 1:2
subplot(3,3,(hl-1)*3+ss)
    imagesc(meanISIhist.logbins,meanISIhist.logbins,meanISIhist.(regions{rr}).(states{ss}).(hilow{hl}).return)
    axis xy
    axis tight
    LogScale('xy',10,'exp',true,'nohalf',true)
    colormap(gca,statecolormap{ss})

end
%legend(cellISIStats.statenames{1:3})
end

NiceSave('HiLowPSSReturn',figfolder,(regions{rr}))
end


%% In/Out Field: REturn Maps ALL
histcolors = flipud(gray);
NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);
REMhistcolors = makeColorMap([1 1 1],[0.8 0 0]);
statecolormap = {histcolors,NREMhistcolors,REMhistcolors};


figure
for rr =1:length(regions)
for ss = 1:2

subplot(6,6,rr+(ss-1)*18)
hold on
for hl = 1:2
    plot(meanISIhist.logbins,meanISIhist.(regions{rr}).(states{ss}).(hilow{hl}).logdist,...
        'linewidth',1,'color',statecolormap{ss}(end./(1+(hl~=ss)),:))
end
    if ss == 1
        title(regions{rr})
    end
    set(gca,'yticklabel',[])
    set(gca,'xticklabel',[])
    axis tight
    box off
    %LogScale('x',10,'exp',true,'nohalf',true)

%legend(cellISIStats.statenames{2:3},'location','southoutside')



for hl = 1:2
subplot(6,6,(hl-1)*6+(ss-1)*18+6+rr)
    imagesc(meanISIhist.logbins,meanISIhist.logbins,meanISIhist.(regions{rr}).(states{ss}).(hilow{hl}).return)
    axis xy
    axis tight
    LogScale('xy',10,'exp',true,'nohalf',true)
    colormap(gca,statecolormap{ss})

    if hl == 1
        set(gca,'XTickLabels',[])
    end
end
%legend(cellISIStats.statenames{1:3})
end

NiceSave('HiLowPSSReturn',figfolder,'AllRegions')
end


%%
for rr = 1:length(regions)
figure
for ss = 1:3

    subplot(4,3,ss)
    hold on
    for tt = 1:2
        plot(log10(GammaParms.(regions{rr}).(states{ss}).meanrate(MutInf.(regions{rr}).CellClass.(celltypes{tt}))),...
            PSScorr.(regions{rr}).(states{ss})(MutInf.(regions{rr}).CellClass.(celltypes{tt})),'.','color',cellcolor{tt})
    end
    axis tight; box off
    plot(xlim(gca),[0 0],'k--')
    xlabel('Mean Rate (Hz)');ylabel('PSS-Rate Corr')
    
    subplot(4,3,ss+3)
        hold on
        ScatterWithLinFit(GammaParms.(regions{rr}).(states{ss}).GSrate,PSScorr.(regions{rr}).(states{ss}))
        axis tight; box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS Rate (Hz)');ylabel('PSS-Rate Corr')
        
    subplot(4,3,ss+6)
        hold on
        plot(-GammaParms.(regions{rr}).(states{ss}).GSmod(GammaParms.(regions{rr}).(states{ss}).sigmod),...
            PSScorr.(regions{rr}).(states{ss})(GammaParms.(regions{rr}).(states{ss}).sigmod),'k.')
        plot(-GammaParms.(regions{rr}).(states{ss}).GSmod(~GammaParms.(regions{rr}).(states{ss}).sigmod),...
            PSScorr.(regions{rr}).(states{ss})(~GammaParms.(regions{rr}).(states{ss}).sigmod),'.','color',[0.5 0.5 0.5])
        axis tight; box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('AS Mod(Hz)');ylabel('PSS-Rate Corr')
        
    subplot(4,3,ss+9)
        hold on
        plot(GammaParms.(regions{rr}).(states{ss}).GSrate(GammaParms.(regions{rr}).(states{ss}).sigmod),...
            -GammaParms.(regions{rr}).(states{ss}).GSmod(GammaParms.(regions{rr}).(states{ss}).sigmod),'k.')
        plot(GammaParms.(regions{rr}).(states{ss}).GSrate(~GammaParms.(regions{rr}).(states{ss}).sigmod),...
            -GammaParms.(regions{rr}).(states{ss}).GSmod(~GammaParms.(regions{rr}).(states{ss}).sigmod),'.','color',[0.5 0.5 0.5])
        axis tight; box off
        plot(xlim(gca),[0 0],'k--')
        ylabel('AS Mod(Hz)');ylabel('GS Rate')
        
end
NiceSave('PSSRateCorr',figfolder,(regions{rr}))
end

%%
BoxAndScatterPlot({PSScorr.(regions{:}).(states{1:2})(MutInf.(regions{:}).CellClass.(celltypes{2}))})

%%
clear allINcolors
for rr = 1:length(regions)
    for ss = 1:2
        allIN{rr+(ss-1)*length(regions)} = ...
            PSScorr.(regions{rr}).(states{ss})(MutInf.(regions{rr}).CellClass.(celltypes{2}));
        allINcolors(rr+(ss-1)*length(regions),:) = statecolors{ss};
    end
end
%%
figure
for rr = 1:length(regions)
for ss = 1:2

    subplot(4,6,(ss-1).*6+rr)
        hold on
        ScatterWithLinFit(GammaParms.(regions{rr}).(states{ss}).GSrate,PSScorr.(regions{rr}).(states{ss}),...
            'RemoveOutlierX',true,'color',statecolors{ss},'markersize',2)
        axis tight; box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS Rate (Hz)');ylabel('PSS-Rate Corr')
        if ss == 1
            title((regions{rr}))
        end
        LogScale('x',10,'exp',true,'nohalf',true)
        
end
subplot(4,2,7)

BoxAndScatterPlot(allIN,'colors',allINcolors,'labels',[regions,regions])
hold on
plot(xlim(gca),[0 0],'k--')
box off
ylabel('PSS-FR Mod (pI)')
end
NiceSave('PSSRateCorr',figfolder,'AllRegions')


%%
figure
for rr = 1:6
    subplot(4,6,rr)
hold on
for ss = 1:3
plot(PSShist.(regions{rr}).mean.bins,PSShist.(regions{rr}).mean.(states{ss}),'color',statecolors{ss},'linewidth',2)
errorshade(PSShist.(regions{rr}).mean.bins,PSShist.(regions{rr}).mean.(states{ss}),...
    PSShist.(regions{rr}).std.(states{ss}),PSShist.(regions{rr}).std.(states{ss}),statecolors{ss},'scalar')
end
axis tight
lim = ylim(gca);ylim([0 lim(2)])
title(regions{rr})
xlabel('PSS')
end
NiceSave('PSShist',figfolder,[])
