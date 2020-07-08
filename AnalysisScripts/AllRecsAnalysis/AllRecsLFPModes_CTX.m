reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop

regions = {'THAL','vCTX','fCTX','BLA','PIR','CA1'};

for rr = 1:3
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISIModeLFPAnalysis_',regions{rr}];
[LFPModeALL.(regions{rr}),baseNames] = GetMatResults(figfolder,['ISIModeLFPAnalysis_',regions{rr}]);
LFPModeALL.(regions{rr}) = bz_CollapseStruct(LFPModeALL.(regions{rr}));
ModalLFPModulation.(regions{rr}) = bz_CollapseStruct(LFPModeALL.(regions{rr}).ModalLFPModulation,'match','justcat',true);

end





%%
states = {'WAKEstate','NREMstate'};
%ss = 1;

for rr = 1:3
ModalLFPModulation.(regions{rr}).numcells = length(ModalLFPModulation.(regions{rr}).UID);
ModalLFPModulation.(regions{rr}).ASfreq = repmat(ModalLFPModulation.(regions{rr}).freq(1,:)',[1,ModalLFPModulation.(regions{rr}).numcells,5]);
%%
for ss = 1:2
    allASweight = repmat(ModalLFPModulation.(regions{rr}).(states{ss}).ASweight,[150,1,1]);
    
    weightthresh = 0.01;
    [ ModalLFPModulation.(regions{rr}).(states{ss}).meanASModulation,...
        ModalLFPModulation.(regions{rr}).(states{ss}).meanASModulation_N,...
        ModalLFPModulation.(regions{rr}).Xbins,ModalLFPModulation.(regions{rr}).Ybins ] = ...
        ConditionalHist3( log10(ModalLFPModulation.(regions{rr}).ASfreq(allASweight>weightthresh)),...
        ModalLFPModulation.(regions{rr}).(states{ss}).ASlogRates(allASweight>weightthresh),...
        ModalLFPModulation.(regions{rr}).(states{ss}).ASModulation(allASweight>weightthresh),...
        'minXY',5,'numXbins',75,'numYbins',75);
end

end
%%
for rr = 1:3
figure
for ss = 1:2
subplot(2,2,ss)
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
imagesc(ModalLFPModulation.(regions{rr}).Xbins,ModalLFPModulation.(regions{rr}).Ybins,...
    ModalLFPModulation.(regions{rr}).(states{ss}).meanASModulation')
alpha((ModalLFPModulation.(regions{rr}).(states{ss}).meanASModulation_N'./max(ModalLFPModulation.(regions{rr}).(states{ss}).meanASModulation_N(:))).*2)
hold on
UnityLine

xlabel('LFP Frequency (Hz)');ylabel('ISI Mode Rate (Hz)')
title(states{ss})
ColorbarWithAxis([-0.0 0.1],'Mode-Power Correlation','inclusive',{'<','>'})
crameri('vik','pivot',0)
axis xy
xlim([0 2.5]);ylim([0 2.5])
LogScale('xy',10,'nohalf',true)

subplot(4,2,7+(ss-1))
hold on
plot(log10(ModalLFPModulation.(regions{rr}).freq(1,:)),nanmean(ModalLFPModulation.(regions{rr}).(states{ss}).MutInf,2),'color','k','linewidth',1)
% errorshade(log10(ModalLFPModulation.freq(1,:)),nanmean(ModalLFPModulation.(states{ss}).MutInf,2),...
%     nanstd(ModalLFPModulation.MutInf,[],2),nanstd(ModalLFPModulation.(states{ss}).MutInf,[],2),'k','scalar')
LogScale('x',10)
xlabel('LFP Frequency (Hz)')
ylabel('MI[ISI;Power]')
colorbar
axis tight


subplot(4,2,5+(ss-1))
hold on

plot(log10(ModalLFPModulation.(regions{rr}).freq),-nanmean(ModalLFPModulation.(regions{rr}).(states{ss}).GSModulation,2),'color','k','linewidth',2)
plot(xlim(gca),[0 0],'k--')
LogScale('x',10)
%xlabel('LFP Frequency (Hz)')
ylabel('AR-Power Corr')
colorbar
axis tight
end

NiceSave('ModeLFPFreqMod',figfolder,(regions{rr}))
end


