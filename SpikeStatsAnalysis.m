function [ISIbins,ISIhist,summstats,...
    ISIreturn,sorts] = SpikeStatsAnalysis(basePath,figfolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%

%% DEV
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/FRHetAndDynamics/AnalysisScripts/AnalysisFigs';
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');

%LFP for plot
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath);

%%
statenames = {'NREMstate','REMstate','WAKEstate'};
statecolors = {'b','r','k'};
numstates = length(statenames);


[ ISIstats ] = bz_ISIStats( spikes,'states',SleepState.ints,...
    'cellclass',CellClass.label,'figfolder',figfolder);


%%
numtimewins = 30;
numjits = 10;
for ss = 1:numstates
    CV2_jitt.(statenames{ss}) = zeros(length(spikes.times),numtimewins,numjits);
end
timebins = logspace(-3,1.75,numtimewins);
for tt = 1:numtimewins
    tt
    for jj = 1:numjits
    [spiketimes_jitt] = JitterSpiketimes(spikes.times,timebins(tt));
    jitspikes = spikes;
    jitspikes.times = spiketimes_jitt;
    [ ISIstats_jitt(tt) ] = bz_ISIStats( jitspikes,'states',SleepState.ints,'showfig',false );
    for ss = 1:numstates
        CV2_jitt.(statenames{ss})(:,tt,jj) = ISIstats_jitt(tt).summstats.(statenames{ss}).meanCV2;
    end
    end
end
%%
    for ss = 1:numstates
        meanCV2_jitt.(statenames{ss}) = mean(CV2_jitt.(statenames{ss}),3);
        stdCV2_jitt.(statenames{ss}) = std(CV2_jitt.(statenames{ss}),[],3);
        difffrom1 = abs(meanCV2_jitt.(statenames{ss})-1);
        sigCV2_jitt.(statenames{ss}) = difffrom1>2.*stdCV2_jitt.(statenames{ss});
        
        meanCV2_jitt.pE.(statenames{ss}) = mean(meanCV2_jitt.(statenames{ss})(CellClass.pE,:),1);
        meanCV2_jitt.pI.(statenames{ss}) = mean(meanCV2_jitt.(statenames{ss})(CellClass.pI,:),1);
        
        stdCV2_jitt.pE.(statenames{ss}) = std(meanCV2_jitt.(statenames{ss})(CellClass.pE,:),[],1);
        stdCV2_jitt.pI.(statenames{ss}) = std(meanCV2_jitt.(statenames{ss})(CellClass.pI,:),[],1);
    end

    


%%
figure
subplot(2,2,1)
hold on
for ss=1:numstates
    %errorshade(log10(timebins),meanCV2_jitt.pE.(statenames{ss}),...
    %    stdCV2_jitt.pE.(statenames{ss}),...
    %    stdCV2_jitt.pE.(statenames{ss}),statecolors{ss},'scalar')
    plot(log10(timebins),meanCV2_jitt.pE.(statenames{ss}),statecolors{ss},'linewidth',2)
    xlabel('Jitter Window (s)');ylabel('<CV2>')
    title('pE Cells')
end
LogScale('x',10)

subplot(2,2,2)
hold on
for ss=1:numstates
    %errorshade(log10(timebins),meanCV2_jitt.pE.(statenames{ss}),...
    %    stdCV2_jitt.pE.(statenames{ss}),...
    %    stdCV2_jitt.pE.(statenames{ss}),statecolors{ss},'scalar')
    plot(log10(timebins),meanCV2_jitt.pI.(statenames{ss}),statecolors{ss},'linewidth',2)
    xlabel('Jitter Window (s)');ylabel('<CV2>')
    title('pI Cells')
end
LogScale('x',10)


bwcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);

for ss=1:numstates
    subplot(2,3,ss+3)
        imagesc(log10(timebins),[0 1],meanCV2_jitt.(statenames{ss})(ISIstats.sorts.(statenames{ss}).CV2byclass,:))
        colorbar
        colormap(bwcolormap)
        caxis([0.5 1.5])
        LogScale('x',10)
        title((statenames{ss}))
        xlabel('Jitter Window (s)');ylabel('Cell, Sorted By Rate/Type')
end

NiceSave(['CV2jitter'],figfolder,baseName);
%%
exdt = [0.1,1,10];
exjit = [13,19,26]; %find this instead of hard code
for dd = 1:length(exdt)
    dt = exdt(dd);
    [spkmat(dd).spkmat, spkmat(dd).timestamps] = SpktToSpkmat(spikes.times,...
        [],dt);
    spkmat(dd).ratemat = spkmat(dd).spkmat./dt;
end

%%
cellnum = 10;
[ twin ] = bz_RandomWindowInIntervals( SleepState.ints.NREMstate,20 );

figure
subplot(2,2,1)
    hold on
    for ss=1:numstates
        errorshade(log10(timebins),meanCV2_jitt.(statenames{ss})(cellnum,:),...
            stdCV2_jitt.(statenames{ss})(cellnum,:),...
            stdCV2_jitt.(statenames{ss})(cellnum,:),statecolors{ss},'scalar')
        plot(log10(timebins),meanCV2_jitt.(statenames{ss})(cellnum,:),statecolors{ss})
        plot([log10(timebins(1)) log10(1./ISIstats.summstats.(statenames{ss}).meanrate(cellnum))],...
            ISIstats.summstats.(statenames{ss}).meanCV2(cellnum).*[1 1],'--','color',statecolors{ss})
        plot(log10(timebins(~sigCV2_jitt.(statenames{ss})(cellnum,:))),~sigCV2_jitt.(statenames{ss})(cellnum,~sigCV2_jitt.(statenames{ss})(cellnum,:)),'ko')
        xlabel('Time Window (s)');ylabel('<CV2>')
    end
    title(['Cell ',num2str(cellnum)])
    LogScale('x',10)
    
subplot(2,1,2)
    %bar(spkmat(3).timestamps,spkmat(3).ratemat(:,cellnum))
    
    bar(spkmat(2).timestamps,spkmat(2).spkmat(:,cellnum),'facecolor','w','edgecolor','k')
    hold on
   % bar(spkmat(1).timestamps,spkmat(1).spkmat(:,cellnum))
    plot(twin,[1 1]-2,'k--')
    plot(twin,[2 2]-2,'k-')
    plot(spikes.times{cellnum},10.*ones(size(spikes.times{cellnum})),'k.')
    plot(ISIstats_jitt(exjit(1)).allspikes.times{cellnum},...
        9.*ones(size(ISIstats_jitt(exjit(1)).allspikes.times{cellnum})),'r.')
    plot(ISIstats_jitt(exjit(2)).allspikes.times{cellnum},...
        8.*ones(size(ISIstats_jitt(exjit(2)).allspikes.times{cellnum})),'g.')
    plot(ISIstats_jitt(exjit(3)).allspikes.times{cellnum},...
        7*ones(size(ISIstats_jitt(exjit(3)).allspikes.times{cellnum})),'b.')
    %plot(spikes.times{cellnum},10.*ones(size(spikes.times{cellnum})),'k.')
    %plot(spikes.times{cellnum},10.*ones(size(spikes.times{cellnum})),'k.')
    plot(ISIstats.allspikes.times{cellnum},ISIstats.allspikes.CV2{cellnum}-2,'k.-')
    %plot(ISIstats_jitt(exjit(1)).allspikes.times{cellnum},ISIstats_jitt(exjit(1)).allspikes.CV2{cellnum}-2,'r.-')
    %plot(ISIstats_jitt(exjit(2)).allspikes.times{cellnum},ISIstats_jitt(exjit(2)).allspikes.CV2{cellnum}-2,'g.-')
    %plot(ISIstats_jitt(exjit(3)).allspikes.times{cellnum},ISIstats_jitt(exjit(3)).allspikes.CV2{cellnum}-2,'b.-')
    xlim(twin);ylim([-2 11])
    
subplot(4,6,4)
    imagesc(ISIstats.ISIhist.logbins,ISIstats.ISIhist.logbins,...
        ISIstats.ISIhist.NREMstate.return(:,:,cellnum)')
    LogScale('xy',10)
    axis xy
    xlabel('ISI_n');ylabel('ISI_n_+_1')
    
for dd = 1:length(exdt)
    subplot(4,6,9+dd)
        imagesc(ISIstats_jitt(exjit(dd)).ISIhist.logbins,ISIstats_jitt(exjit(dd)).ISIhist.logbins,...
            ISIstats_jitt(exjit(dd)).ISIhist.NREMstate.return(:,:,cellnum)')
        LogScale('xy',10)
        axis xy
        set(gca,'ytick',[]);set(gca,'xtick',[])
end

NiceSave(['CV2jitter_Cell',num2str(cellnum)],figfolder,baseName);


%%
twin
figure
plot(spikes.times{1},ones(size(spikes.times{1})),'.')
hold on
plot(spiketimes_jitt{1},ones(size(spiketimes_jitt{1})),'.')

%%
figure
for ss=1:numstates
    subplot(2,2,ss)
        plot(log2(summstats.(statenames{ss}).ISICV(CellClass.pE)),...
            log2(summstats.(statenames{mod(ss,numstates)+1}).ISICV(CellClass.pE)),'k.')
        hold on
        plot(log2(summstats.(statenames{ss}).ISICV(CellClass.pI)),...
            log2(summstats.(statenames{mod(ss,numstates)+1}).ISICV(CellClass.pI)),'r.')
        plot(log2([1 4]),log2([1 4]),'k')
        xlabel([statenames{ss},' CV']);ylabel([statenames{mod(ss,numstates)+1},' CV'])
        LogScale('xy',2)
end

figure
for ss=1:numstates
    subplot(2,2,ss)
        plot((summstats.(statenames{ss}).meanCV2(CellClass.pE)),...
            (summstats.(statenames{mod(ss,numstates)+1}).meanCV2(CellClass.pE)),'k.')
        hold on
        plot((summstats.(statenames{ss}).meanCV2(CellClass.pI)),...
            (summstats.(statenames{mod(ss,numstates)+1}).meanCV2(CellClass.pI)),'r.')
        plot(([0.5 2]),([0.5 2]),'k')
        xlabel([statenames{ss},' CV2']);ylabel([statenames{mod(ss,numstates)+1},' CV2'])
       % LogScale('xy',2)
end

%%
% cc=21;
% figure
% plot(log10(ISIs{cc}),CV2{cc},'.')
% LogScale('x',10)
% hold on
% plot(get(gca,'xlim'),summstats.(statenames{ss}).meanCV2(cc).*[1 1],'r')
%%
% figure
% colormap(histcolors)
% for ss = 1:numstates
%     
% subplot(2,3,ss)
%     imagesc((ISIbins.loginv),[1 numcells],ISIhist.(statenames{ss}).loginv(sorts.(statenames{ss}).rate,:))
%     hold on
%     plot(log10(summstats.(statenames{ss}).meanrate(sorts.(statenames{ss}).rate)),[1:numcells],'k','LineWidth',2)
%     LogScale('x',10)
%     xlabel('1./ISI (Hz)')
%     xlim(ISIbins.loginv([1 end]))
%    % colorbar
%     %legend('Mean Firing Rate (Hz)','location','southwest')
%     ylabel('Cell (Sorted by Mean FR)')
%     caxis([0 0.1])
%     title((statenames{ss}))
%     
% subplot(2,3,ss+3)
%     imagesc((ISIbins.norm),[1 numcells],ISIhist.(statenames{ss}).norm(sorts.(statenames{ss}).rate,:))
%     hold on
%     plot([0 0],[1 numcells],'k','LineWidth',2)
%     LogScale('x',10)
%     xlabel('ISI^-^1/meanrate')
%     xlim(ISIbins.norm([1 end]))
%    % colorbar
%     %legend('Mean Firing Rate (Hz)','location','southwest')
%     ylabel('Cell (Sorted by Mean FR)')
%     caxis([0 0.1])
%     title([(statenames{ss}),': Norm'])
%     
% end
%     
% NiceSave('FRDist_all',figfolder,baseName);

%%

figure
colormap(histcolors)
ff=0;
for cc = 1:numcells
    cellnum = sortrate.NREMstate(cc);
    subplot(6,7,mod(cc-1,42)+1)
    imagesc((ISIbins.log),(ISIbins.log),squeeze(ISIreturn.NREMstate.log(cellnum,:,:)))
    hold on
    plot(log10(1./summstats.NREMstate.meanrate(cellnum)),log10(1./summstats.NREMstate.meanrate(cellnum)),'k+')
    axis xy
    LogScale('xy',10)
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    caxis([0 0.003])
    xlim(ISIbins.log([1 end]));ylim(ISIbins.log([1 end]))
    xlabel(['FR: ',num2str(round(summstats.NREMstate.meanrate(cellnum),2)),'Hz'])
    
    if mod(cc,42) == 0 || cc ==numcells
        ff= ff+1;
        NiceSave(['invISImap',num2str(ff)],figfolder,baseName);
        figure
        colormap(histcolors)
    end
end

%%
ISICVdiff = (summstats.NREMstate.ISICV)-(summstats.WAKEstate.ISICV);
figure
plot(log10(summstats.NREMstate.meanrate(CellClass.pE)),ISICVdiff(CellClass.pE),'.')
hold on
plot(log10(summstats.NREMstate.meanrate(CellClass.pI)),ISICVdiff(CellClass.pI),'r.')
plot(get(gca,'xlim'),[0 0],'k')
%%
% figure
% colormap(histcolors)
% ff=0;
% for cc = 1:numcells
%     cellnum = sortrate.NREMpacket(cc);
%     subplot(6,6,mod(cc-1,36)+1)
%     imagesc((ISIbins.log),(ISIbins.log),squeeze(ISIreturn.NREMpacket.log(cellnum,:,:)))
%     hold on
%     plot(log10(1./meanrate.NREMpacket(cellnum)),log10(1./meanrate.NREMpacket(cellnum)),'r+')
%     axis xy
%     LogScale('xy',10)
%     set(gca,'ytick',[]);set(gca,'xtick',[]);
%     caxis([0 0.003])
%     
%     if mod(cc,36) == 0 | cc ==numcells
%         ff= ff+1;
%         NiceSave(['ISImap',num2str(ff)],figfolder,baseName);
%         figure
%         colormap(histcolors)
%     end
% end

%%
end
    %%
%     
% subplot(2,2,2)
%     imagesc((normbins),[1 numcells],normcv2hist(sortrate,:))
%     hold on
%     plot([0 0],[0 numcells+1],'k','LineWidth',2)
%     %LogScale('x',10)
%     xlabel('1./CV2 (Hz)')
%     %xlim(isibins([1 end]))
%     colorbar
%     %legend('Mean Firing Rate (Hz)','location','southwest')
%     ylabel('Cell (Sorted by Mean FR)')
%     caxis([0 0.1])
%     
% NiceSave('FRDist_',figfolder,baseName);
%%
% dt = 1;
% overlap = 10;
% winsize = dt*overlap;
% [ratemat,t] = SpktToSpkmat(Se,[],dt,overlap);
% ratemat = ratemat./winsize;
% 
% 
% 
% for ss = 1:numstates
%     stateIDX.(statenames{ss}) = INTtoIDX(StateIntervals.(statenames{ss}),length(ratemat(:,1)),1/dt);
%     meanFR.(statenames{ss}) = mean(ratemat(stateIDX.(statenames{ss}),:));
%     [~,sort.(statenames{ss}).rate] = sort(meanFR.(statenames{ss}));
% end
% 
% 
% 
% %%
% %bins = unique(NREMratemat(:));
% bins = 0:0.1:10;
% for ss = 1:numstates
%     FRdist.(statenames{ss}) = hist(ratemat(stateIDX.(statenames{ss}),:),bins);
%     FRdist.(statenames{ss}) = bsxfun(@(X,Y) X./Y,FRdist.(statenames{ss}),sum(FRdist.(statenames{ss})));
% end    
% % [FRdist.REM] = hist(ratemat(REMidx,:),bins);
% % [FRdist.WAKE] = hist(ratemat(WAKEidx,:),bins);
% % FRdist = structfun(@(A) bsxfun(@(X,Y) X./Y,A,sum(A)),FRdist,'UniformOutput',false);
% 
% %%
% ff = figure;
% for distidx = 1:numstates
%     for sortidx = 1:numstates
%         subplot(numstates,numstates,distidx+(sortidx-1)*numstates)
%             imagesc(bins,1:length(Se),log10(FRdist.(statenames{distidx})(:,sortrate.(statenames{sortidx})))')
%             ColorbarWithAxis([-3 -0.3],'log p(FR)')
%            % colorbar
%             xlabel('Firing Rate (Hz)');ylabel(['Cell - Sorted by Mean FR, ',(statenames{sortidx})])
%             title(['Single Cell Firing Rate Distributions - ',(statenames{distidx})])
%             axis xy
%     end
% end
%     saveas(ff,[figfolder,baseName,'_FRHet'],'jpeg')
%     