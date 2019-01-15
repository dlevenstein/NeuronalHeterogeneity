function [ jitterCV2,ISIstats ] = SpikeStatsAnalysis(basePath,figfolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%

%% DEV
repoRoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity'; %desktop

basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/mnt/NyuShare/Buzsakilabspace/Datasets/GrosmarkAD/Gatsby/Gatsby_08022013';

%figfolder = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs';
figfolder = [repoRoot,'/AnalysisScripts/AnalysisFigs/SpikeStatsAnalysis'];

%figfolder = '/mnt/data1/Dropbox/research/Current Projects/FRHET_temp/SpikeStatsAnalysis';

%%

baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath,'noPrompts',true);

%%
statenames = {'NREMstate','WAKEstate','REMstate'};
statecolors = {'b','k','r'};
numstates = length(statenames);


[ ISIstats ] = bz_ISIStats( spikes,'ints',SleepState.ints,...
    'cellclass',CellClass.label,'figfolder',figfolder,'shuffleCV2',true,...
    'savecellinfo',true,'basePath',basePath,'forceRedetect',true);



%% CV/CV2 as a f'n of rate
figure

for ss=1:numstates
    subplot(3,3,ss)
        plot(log10(ISIstats.summstats.(statenames{ss}).meanrate(CellClass.pE)),...
            log10(ISIstats.summstats.(statenames{mod(ss,numstates)+1}).meanrate(CellClass.pE)),'k.')
        hold on
        plot(log10(ISIstats.summstats.(statenames{ss}).meanrate(CellClass.pI)),...
            log10(ISIstats.summstats.(statenames{mod(ss,numstates)+1}).meanrate(CellClass.pI)),'r.')
        plot(log10([0.03 30]),log10([0.03 30]),'k')
        xlabel([statenames{ss},' Rate']);ylabel([statenames{mod(ss,numstates)+1},' Rate'])
        LogScale('xy',10)
end

for ss=1:numstates
    subplot(3,3,ss+3)
        plot(log2(ISIstats.summstats.(statenames{ss}).ISICV(CellClass.pE)),...
            log2(ISIstats.summstats.(statenames{mod(ss,numstates)+1}).ISICV(CellClass.pE)),'k.')
        hold on
        plot(log2(ISIstats.summstats.(statenames{ss}).ISICV(CellClass.pI)),...
            log2(ISIstats.summstats.(statenames{mod(ss,numstates)+1}).ISICV(CellClass.pI)),'r.')
        plot(log2([1 6]),log2([1 6]),'k')
        xlabel([statenames{ss},' CV']);ylabel([statenames{mod(ss,numstates)+1},' CV'])
        LogScale('xy',2)
end

for ss=1:numstates
    subplot(3,3,ss+6)
        plot((ISIstats.summstats.(statenames{ss}).meanCV2(CellClass.pE)),...
            (ISIstats.summstats.(statenames{mod(ss,numstates)+1}).meanCV2(CellClass.pE)),'k.')
        hold on
        plot((ISIstats.summstats.(statenames{ss}).meanCV2(CellClass.pI)),...
            (ISIstats.summstats.(statenames{mod(ss,numstates)+1}).meanCV2(CellClass.pI)),'r.')
        plot(([0.5 2]),([0.5 2]),'k')
        xlabel([statenames{ss},' CV2']);ylabel([statenames{mod(ss,numstates)+1},' CV2'])
       % LogScale('xy',2)
end

NiceSave(['StateComparison'],figfolder,baseName);


%%
%histcolors = [makeColorMap([1 1 1],[0.8 0.2 0.2],[0.8 0 0]);makeColorMap([0.8 0 0],[0.3 0 0])];
histcolors = flipud(gray);
cellnos = ISIstats.sorts.NREMstate.CV2([1,round(end/2),end]);
sortmetric = nan(size(ISIstats.summstats.NREMstate.meanCV2));
sortmetric(cellnos) = [1 2 3];

figure
colormap(histcolors)
%Return Map Samples
for ss = 1:numstates
subplot(5,2,2.*(ss-1)+2)

    [ twin ] = bz_RandomWindowInIntervals( SleepState.ints.(statenames{ss}),6 );
    bz_MultiLFPPlot(lfp,'timewin',twin,'spikes',spikes,'plotcells',cellnos,...
        'sortmetric',sortmetric)
    ylabel(statenames{ss})
end

%Time series samples
for cc = 1:length(cellnos)
    for ss=1:numstates
        subplot(5,6,(ss-1).*6+cc)
            imagesc(ISIstats.ISIhist.logbins,ISIstats.ISIhist.logbins,...
                ISIstats.ISIhist.(statenames{ss}).return(:,:,cellnos(cc))')
            hold on
            plot(log10(1./ISIstats.summstats.(statenames{ss}).meanrate(cellnos(cc))),...
                log10(1./ISIstats.summstats.(statenames{ss}).meanrate(cellnos(cc))),'r+')
            LogScale('xy',10)
            set(gca,'Ytick',[]);set(gca,'Xtick',[])
            axis xy
            %xlabel('ISI_n');ylabel('ISI_n_+_1')
            title(['CV2: ',num2str(round(ISIstats.summstats.(statenames{ss}).meanCV2(cellnos(cc)),2))])
    end
end
NiceSave(['CV2Examples'],figfolder,baseName);

%% All the return Maps

figure
colormap(histcolors)
ff=0;
for cc = 1:spikes.numcells
    cellnum = ISIstats.sorts.NREMstate.CV2byclass(cc);   %%sortrate.NREMstate(cc);
    subplot(6,7,mod(cc-1,42)+1)
    imagesc((ISIstats.ISIhist.logbins),(ISIstats.ISIhist.logbins),(ISIstats.ISIhist.NREMstate.return(:,:,cellnum)))
    hold on
    plot(log10(1./ISIstats.summstats.NREMstate.meanrate(cellnum)),log10(1./ISIstats.summstats.NREMstate.meanrate(cellnum)),'k+')
    axis xy
    LogScale('xy',10)
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    caxis([0 0.003])
    xlim(ISIstats.ISIhist.logbins([1 end]));ylim(ISIstats.ISIhist.logbins([1 end]))
    %xlabel(['FR: ',num2str(round(ISIstats.summstats.NREMstate.meanrate(cellnum),2)),'Hz'])
    xlabel(['CV2: ',num2str(round(ISIstats.summstats.NREMstate.meanCV2(cellnum),2))])
    if mod(cc,42) == 0 || cc ==spikes.numcells
        ff= ff+1;
        NiceSave(['ISIreturnmap',num2str(ff)],figfolder,baseName);
        figure
        colormap(histcolors)
    end
end
close
%% Measuring Similarity between ISI maps (same cell different state)
%Need to run PCA denoise on all maps
%     linearizedreturn1 = reshape(ISIstats.ISIhist.NREMstate.return,[],spikes.numcells);
%     linearizedreturn2 =reshape(ISIstats.ISIhist.WAKEstate.return,[],spikes.numcells);
% X = [reshape(ISIstats.ISIhist.NREMstate.return,[],spikes.numcells) ,...
%     reshape(ISIstats.ISIhist.WAKEstate.return,[],spikes.numcells),...
%     reshape(ISIstats.ISIhist.REMstate.return,[],spikes.numcells)]';
% 
% [Xdn, sigma, npars, u, vals, v] = bz_PCAdenoise(X);
% 
% cleanmaps = reshape(Xdn',...
%     size(ISIstats.ISIhist.NREMstate.return,1),size(ISIstats.ISIhist.NREMstate.return,2),...
%     size(X,1));
% 
% %%
% cellnum = 3;
% figure
% subplot(2,2,1)
% imagesc(ISIstats.ISIhist.NREMstate.return(:,:,cellnum))
% subplot(2,2,2)
% imagesc(cleanmaps(:,:,cellnum))
%Should do this with the PCA denoised maps!
%for ss=1:numstates
%     linearizedreturn1 = reshape(ISIstats.ISIhist.NREMstate.return,[],spikes.numcells);
%     linearizedreturn2 =reshape(ISIstats.ISIhist.WAKEstate.return,[],spikes.numcells);
%     returnmapsimilarity = corr(linearizedreturn1,linearizedreturn2,'type','spearman');
% 
%     
%     NWsimilarity = diag(returnmapsimilarity);
%     
% %     figure
% %     imagesc(returnmapsimilarity(ISIstats.sorts.NREMstate.ratebyclass,ISIstats.sorts.NREMstate.ratebyclass))
% %     colorbar
% hist(NWsimilarity)
% figure
% plot(log10(ISIstats.summstats.NREMstate.meanrate),NWsimilarity,'.')

%%
% figure
%         plot((ISIstats.summstats.(statenames{ss}).meanCV2(CellClass.pE)),...
%             (ISIstats.summstats.(statenames{ss}).ISICV(CellClass.pE)),'k.')
%         hold on
%         plot((ISIstats.summstats.(statenames{ss}).meanCV2(CellClass.pI)),...
%             (ISIstats.summstats.(statenames{ss}).ISICV(CellClass.pI)),'r.')

%% Jitter control
numtimewins = 30;
numjits = 20;
for ss = 1:numstates
    CV2_jitt.(statenames{ss}) = zeros(length(spikes.times),numtimewins,numjits);
end
jitterCV2.timebins = logspace(-3,1.75,numtimewins);
for tt = 1:numtimewins
    tt
    for jj = 1:numjits
    [spiketimes_jitt] = JitterSpiketimes(spikes.times,jitterCV2.timebins(tt));
    jitspikes = spikes;
    jitspikes.times = spiketimes_jitt;
    [ ISIstats_jitt(tt) ] = bz_ISIStats( jitspikes,'ints',SleepState.ints,'showfig',false );
    for ss = 1:numstates
        CV2_jitt.(statenames{ss})(:,tt,jj) = ISIstats_jitt(tt).summstats.(statenames{ss}).meanCV2;
    end
    end
end
%%
    for ss = 1:numstates
        jitterCV2.mean.(statenames{ss}) = mean(CV2_jitt.(statenames{ss}),3);
        jitterCV2.std.(statenames{ss}) = std(CV2_jitt.(statenames{ss}),[],3);
        difffrom1 = abs(jitterCV2.mean.(statenames{ss})-1);
        sigCV2_jitt.(statenames{ss}) = difffrom1>2.*jitterCV2.std.(statenames{ss});
        
        jitterCV2.mean.pE.(statenames{ss}) = mean(jitterCV2.mean.(statenames{ss})(CellClass.pE,:),1);
        jitterCV2.mean.pI.(statenames{ss}) = mean(jitterCV2.mean.(statenames{ss})(CellClass.pI,:),1);
        
        jitterCV2.std.pE.(statenames{ss}) = std(jitterCV2.mean.(statenames{ss})(CellClass.pE,:),[],1);
        jitterCV2.std.pI.(statenames{ss}) = std(jitterCV2.mean.(statenames{ss})(CellClass.pI,:),[],1);
    end

    


%%
figure
subplot(2,2,1)
hold on
for ss=1:numstates
    %errorshade(log10(timebins),meanCV2_jitt.pE.(statenames{ss}),...
    %    stdCV2_jitt.pE.(statenames{ss}),...
    %    stdCV2_jitt.pE.(statenames{ss}),statecolors{ss},'scalar')
    plot(log10(jitterCV2.timebins),jitterCV2.mean.pE.(statenames{ss}),statecolors{ss},'linewidth',2)
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
    plot(log10(jitterCV2.timebins),jitterCV2.mean.pI.(statenames{ss}),statecolors{ss},'linewidth',2)
    xlabel('Jitter Window (s)');ylabel('<CV2>')
    title('pI Cells')
end
LogScale('x',10)


bwcolormap = [makeColorMap([0 0 0.2],[0 0 0.9],[1 1 1]);makeColorMap([1 1 1],[0.9 0 0],[0.2 0 0])];

for ss=1:numstates
    subplot(2,3,ss+3)
        imagesc(log10(jitterCV2.timebins),[0 1],jitterCV2.mean.(statenames{ss})(ISIstats.sorts.(statenames{ss}).CV2byclass,:))
        colorbar
        colormap(bwcolormap)
        caxis([0.5 1.5])
        LogScale('x',10)
        title((statenames{ss}))
        xlabel('Jitter Window (s)');ylabel('Cell, Sorted By <CV2>/Type')
end

NiceSave(['CV2jitter'],figfolder,baseName);
%%
% exdt = [0.1,1,10];
% exjit = [13,19,26]; %find this instead of hard code
% for dd = 1:length(exdt)
%     dt = exdt(dd);
%     [spkmat(dd).spkmat, spkmat(dd).timestamps] = SpktToSpkmat(spikes.times,...
%         [],dt);
%     spkmat(dd).ratemat = spkmat(dd).spkmat./dt;
% end
% 
% 
% %%
% cellnum = 48;
% [ twin ] = bz_RandomWindowInIntervals( SleepState.ints.NREMstate,20 );
% 
% figure
% subplot(2,2,1)
%     hold on
%     for ss=1:numstates
%         errorshade(log10(timebins),meanCV2_jitt.(statenames{ss})(cellnum,:),...
%             stdCV2_jitt.(statenames{ss})(cellnum,:),...
%             stdCV2_jitt.(statenames{ss})(cellnum,:),statecolors{ss},'scalar')
%         plot(log10(timebins),meanCV2_jitt.(statenames{ss})(cellnum,:),statecolors{ss})
%         plot([log10(timebins(1)) log10(1./ISIstats.summstats.(statenames{ss}).meanrate(cellnum))],...
%             ISIstats.summstats.(statenames{ss}).meanCV2(cellnum).*[1 1],'--','color',statecolors{ss})
%         plot(log10(timebins(~sigCV2_jitt.(statenames{ss})(cellnum,:))),~sigCV2_jitt.(statenames{ss})(cellnum,~sigCV2_jitt.(statenames{ss})(cellnum,:)),'ko')
%         xlabel('Time Window (s)');ylabel('<CV2>')
%     end
%     title(['Cell ',num2str(cellnum)])
%     LogScale('x',10)
%     
% subplot(3,1,3)
%     %bar(spkmat(3).timestamps,spkmat(3).ratemat(:,cellnum))
%     
%     bar(spkmat(2).timestamps,spkmat(2).spkmat(:,cellnum),'facecolor','w','edgecolor','k')
%     hold on
%    % bar(spkmat(1).timestamps,spkmat(1).spkmat(:,cellnum))
%     plot(twin,[1 1]-2,'k--')
%     plot(twin,[2 2]-2,'k-')
%     plot(spikes.times{cellnum},10.*ones(size(spikes.times{cellnum})),'k.')
%     plot(ISIstats_jitt(exjit(1)).allspikes.times{cellnum},...
%         9.*ones(size(ISIstats_jitt(exjit(1)).allspikes.times{cellnum})),'r.')
%     plot(ISIstats_jitt(exjit(2)).allspikes.times{cellnum},...
%         8.*ones(size(ISIstats_jitt(exjit(2)).allspikes.times{cellnum})),'g.')
%     plot(ISIstats_jitt(exjit(3)).allspikes.times{cellnum},...
%         7*ones(size(ISIstats_jitt(exjit(3)).allspikes.times{cellnum})),'b.')
%     %plot(spikes.times{cellnum},10.*ones(size(spikes.times{cellnum})),'k.')
%     %plot(spikes.times{cellnum},10.*ones(size(spikes.times{cellnum})),'k.')
%     plot(ISIstats.allspikes.times{cellnum},ISIstats.allspikes.CV2{cellnum}-2,'k.-')
%     %plot(ISIstats_jitt(exjit(1)).allspikes.times{cellnum},ISIstats_jitt(exjit(1)).allspikes.CV2{cellnum}-2,'r.-')
%     %plot(ISIstats_jitt(exjit(2)).allspikes.times{cellnum},ISIstats_jitt(exjit(2)).allspikes.CV2{cellnum}-2,'g.-')
%     %plot(ISIstats_jitt(exjit(3)).allspikes.times{cellnum},ISIstats_jitt(exjit(3)).allspikes.CV2{cellnum}-2,'b.-')
%     xlim(twin);ylim([-2 11])
%     
%     for ss=1:numstates
%         subplot(4,6,3+ss)
%             imagesc(ISIstats.ISIhist.logbins,ISIstats.ISIhist.logbins,...
%                 ISIstats.ISIhist.(statenames{ss}).return(:,:,cellnum)')
%             hold on
%             plot(log10(1./ISIstats.summstats.(statenames{ss}).meanrate(cellnum)),...
%                 log10(1./ISIstats.summstats.(statenames{ss}).meanrate(cellnum)),'r+')
%             LogScale('xy',10)
%             axis xy
%             xlabel('ISI_n');ylabel('ISI_n_+_1')
%             title({(statenames{ss}),['<CV2> = ',num2str(ISIstats.summstats.(statenames{ss}).meanCV2(cellnum))]})
%     end
%     
% for dd = 1:length(exdt)
%     subplot(4,6,9+dd)
%         imagesc(ISIstats_jitt(exjit(dd)).ISIhist.logbins,ISIstats_jitt(exjit(dd)).ISIhist.logbins,...
%             ISIstats_jitt(exjit(dd)).ISIhist.NREMstate.return(:,:,cellnum)')
%         LogScale('xy',10)
%         axis xy
%         set(gca,'ytick',[]);set(gca,'xtick',[])
% end
% 
% NiceSave(['CV2jitter_Cell',num2str(cellnum)],figfolder,baseName);






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
% ISICVdiff = (summstats.NREMstate.ISICV)-(summstats.WAKEstate.ISICV);
% figure
% plot(log10(summstats.NREMstate.meanrate(CellClass.pE)),ISICVdiff(CellClass.pE),'.')
% hold on
% plot(log10(summstats.NREMstate.meanrate(CellClass.pI)),ISICVdiff(CellClass.pI),'r.')
% plot(get(gca,'xlim'),[0 0],'k')
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