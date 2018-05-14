function [ISIbins,ISIhist,summstats,...
    ISIreturn,sorts] = SpikeStatsA(basePath,figfolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%   INPUTS
%       spikes
%
%       (options)
%       'states'
%       'savecellinfo'
%       'basePath'
%       'figfolder'
%
%   OUTPUTS
%       
%
%DLevenstein
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
numstates = length(statenames);
%% ISI and CV2 statistics
numcells = length(spikes.UID);

%Calculate ISI and CV2 for allspikes
spikes.ISIs = cellfun(@diff,spikes.times,'UniformOutput',false);
spikes.CV2 = cellfun(@(X) 2.*abs(X(2:end)-X(1:end-1))./(X(2:end)+X(1:end-1)),spikes.ISIs ,'UniformOutput',false);
%spikes.estrate = cellfun(@(X) 2./(X(2:end)+X(1:end-1)),spikes.ISIs ,'UniformOutput',false);
%%
%Which state to look at?
for ss = 1:numstates
%ss=1;
%ss = 1;

%Find which spikes are during state of interest
[statespiketimes,statespikes] = cellfun(@(X) RestrictInts(X,SleepState.ints.(statenames{ss})),...
    spikes.times,'UniformOutput',false);
CV2 = cellfun(@(X,Y) X(Y(2:end-1)),spikes.CV2,statespikes,'Uniformoutput',false);
ISIs = cellfun(@(X,Y) X(Y(2:end-1)),spikes.ISIs,statespikes,'Uniformoutput',false);

%Summary Statistics
summstats.(statenames{ss}).meanISI = cellfun(@(X) mean(X),ISIs);
summstats.(statenames{ss}).meanrate = 1./summstats.(statenames{ss}).meanISI;
summstats.(statenames{ss}).ISICV = cellfun(@(X) std(X)./mean(X),ISIs);
summstats.(statenames{ss}).meanCV2 = cellfun(@(X) mean(X),CV2);


%% CV2 for shuffle (shows that CV2 is not much meaningful?)
% numshuffle = 100;
% for sh = 1:numshuffle
%     ISIs_shuffle = cellfun(@(X) shuffle(X),ISIs,'UniformOutput',false);
%     CV2_shuffle = cellfun(@(X) 2.*abs(X(2:end)-X(1:end-1))./(X(2:end)+X(1:end-1)),...
%         ISIs_shuffle ,'UniformOutput',false);
%     meanshuffle(sh,:) = cellfun(@(X) mean(X),CV2_shuffle);
% end
% shufflemean = mean(meanshuffle);
% shufflestd = std(meanshuffle);
% 
% %CV2_reltoshuff = (summstats.(statenames{ss}).meanCV2-shufflemean)./shufflestd;
% 
% figure
% subplot(2,2,1)
% plot([log10(summstats.(statenames{ss}).meanrate);log10(summstats.(statenames{ss}).meanrate)],...
%     [shufflemean;summstats.(statenames{ss}).meanCV2],'color',0.7.*[1 1 1],'linewidth',0.5)
% hold on
% plot([log10(summstats.(statenames{ss}).meanrate);log10(summstats.(statenames{ss}).meanrate)],...
%     [shufflemean-shufflestd;shufflemean+shufflestd],'color',0.7.*[1 1 1],'linewidth',3)
% plot(log10(summstats.(statenames{ss}).meanrate),summstats.(statenames{ss}).meanCV2,'.r','markersize',10)
% xlabel('FR (Hz)');ylabel('<CV2>')
% NiceSave(['CV2_',(statenames{ss})],figfolder,baseName);

%%
%Set up all the bins and matrices
numbins = 60;
ISIbins.lin = linspace(0,10,numbins);
ISIbins.log = linspace(log10(0.001),log10(200),numbins);
ISIhist.(statenames{ss}).lin = zeros(numcells,numbins);
ISIhist.(statenames{ss}).log = zeros(numcells,numbins);
normcv2hist = zeros(numcells,numbins);

ISIreturn.(statenames{ss}).log = zeros(numcells,numbins,numbins);

%Calculate all the histograms: ISI, log(ISI), 1/ISI, log(1/ISI)
for cc = 1:numcells
    numspks(cc) = length(ISIs{cc});
    
    %Calculate ISI histograms
    ISIhist.(statenames{ss}).lin(cc,:) = hist(ISIs{cc},ISIbins.lin);
    ISIhist.(statenames{ss}).log(cc,:) = hist(log10(ISIs{cc}),ISIbins.log);
    
    %Normalize histograms to number of spikes
    ISIhist.(statenames{ss}).lin(cc,:) = ISIhist.(statenames{ss}).lin(cc,:)./numspks(cc);
    ISIhist.(statenames{ss}).log(cc,:) = ISIhist.(statenames{ss}).log(cc,:)./numspks(cc);
    
    %Calculate Return maps
    if numspks(cc)>1
    ISIreturn.(statenames{ss}).log(cc,:,:) = hist3(log10([ISIs{cc}(1:end-1) ISIs{cc}(2:end)]),{ISIbins.log,ISIbins.log});
    end
    ISIreturn.(statenames{ss}).log(cc,:,:) = ISIreturn.(statenames{ss}).log(cc,:,:)./numspks(cc);
  
end

%Sortings
[~,sorts.(statenames{ss}).rate]=sort(summstats.(statenames{ss}).meanrate);
[~,sorts.(statenames{ss}).ISICV]=sort(summstats.(statenames{ss}).ISICV);
[~,sorts.(statenames{ss}).CV2]=sort(summstats.(statenames{ss}).meanCV2);

sorts.(statenames{ss}).rateE = intersect(sorts.(statenames{ss}).rate,...
                                    find(CellClass.pE),'stable');
sorts.(statenames{ss}).rateI = intersect(sorts.(statenames{ss}).rate,...
                                    find(CellClass.pI),'stable');
sorts.(statenames{ss}).rateEI = [sorts.(statenames{ss}).rateE sorts.(statenames{ss}).rateI];

sorts.(statenames{ss}).ISICVE = intersect(sorts.(statenames{ss}).ISICV,...
                                    find(CellClass.pE),'stable');
sorts.(statenames{ss}).ISICVI = intersect(sorts.(statenames{ss}).ISICV,...
                                    find(CellClass.pI),'stable');
sorts.(statenames{ss}).ISICVEI = [sorts.(statenames{ss}).ISICVE sorts.(statenames{ss}).ISICVI];

sorts.(statenames{ss}).CV2E = intersect(sorts.(statenames{ss}).CV2,...
                                    find(CellClass.pE),'stable');
sorts.(statenames{ss}).CV2I = intersect(sorts.(statenames{ss}).CV2,...
                                    find(CellClass.pI),'stable');
sorts.(statenames{ss}).CV2EI = [sorts.(statenames{ss}).CV2E sorts.(statenames{ss}).CV2I];


%% ACG
%[ccg,t] = CCG(statespiketimes,[],<options>)

%%
figure
    subplot(2,2,1)
        plot(log10(summstats.(statenames{ss}).meanrate(CellClass.pE)),...
            log2(summstats.(statenames{ss}).ISICV(CellClass.pE)),'k.')
        hold on
        plot(log10(summstats.(statenames{ss}).meanrate(CellClass.pI)),...
            log2(summstats.(statenames{ss}).ISICV(CellClass.pI)),'r.')
        LogScale('x',10);LogScale('y',2);
        xlabel('Mean Rate (Hz)');ylabel('ISI CV')
        title(statenames{ss})
        box off
        
    subplot(2,2,3)
        plot(log10(summstats.(statenames{ss}).meanrate(CellClass.pE)),...
            (summstats.(statenames{ss}).meanCV2(CellClass.pE)),'k.')
        hold on
        plot(log10(summstats.(statenames{ss}).meanrate(CellClass.pI)),...
            (summstats.(statenames{ss}).meanCV2(CellClass.pI)),'r.')
        plot(get(gca,'xlim'),[1 1])
        LogScale('x',10);
        xlabel('Mean Rate (Hz)');ylabel('ISI <CV2>')
        title(statenames{ss})
        box off
        
    subplot(2,2,2)
        imagesc((ISIbins.log),[1 numcells],...
            ISIhist.(statenames{ss}).log(sorts.(statenames{ss}).rateEI,:))
        hold on
        plot(log10(1./(summstats.(statenames{ss}).meanrate(sorts.(statenames{ss}).rateEI))),[1:numcells],'k.','LineWidth',2)
        plot(ISIbins.log([1 end]),sum(CellClass.pE).*[1 1]+0.5,'r')
        LogScale('x',10)
        xlabel('ISI (s)')
        xlim(ISIbins.log([1 end]))
        colorbar
      %  legend('1/Mean Firing Rate (s)','location','southeast')
        ylabel('Cell (Sorted by FR, Type)')
        %legend('1/Mean Firing Rate (s)','location','southeast')
        caxis([0 0.1])
        title('ISI Distribution (Log Scale)')
        
    subplot(2,2,4)
        imagesc((ISIbins.log),[1 numcells],...
            ISIhist.(statenames{ss}).log(sorts.(statenames{ss}).ISICVEI,:))
        hold on
       % plot(log10(1./(summstats.(statenames{ss}).meanrate(sorts.(statenames{ss}).rateEI))),[1:numcells],'k.','LineWidth',2)
        plot(ISIbins.log([1 end]),sum(CellClass.pE).*[1 1]+0.5,'r')
        LogScale('x',10)
        xlabel('ISI (s)')
        xlim(ISIbins.log([1 end]))
        colorbar
      %  legend('1/Mean Firing Rate (s)','location','southeast')
        ylabel('Cell (Sorted by CV2, Type)')
        %legend('1/Mean Firing Rate (s)','location','southeast')
        caxis([0 0.1])
        title('ISI Distribution (Log Scale)')
        
%     subplot(2,2,3)
%         plot(log2(summstats.(statenames{ss}).meanrate(CellClass.pE)),...
%             log2(summstats.(statenames{ss}).meanCV2(CellClass.pE)),'k.')
%         hold on
%         plot(log2(summstats.(statenames{ss}).meanrate(CellClass.pI)),...
%             log2(summstats.(statenames{ss}).meanCV2(CellClass.pI)),'r.')
%         LogScale('xy',2)
%         xlabel('Mean Rate (Hz)');ylabel('Mean CV2')
%         title(statenames{ss})
%         box off


NiceSave(['ISIstats_',(statenames{ss})],figfolder,baseName);


%exneurons %top/bottom 25 %ile rate/ISICV
%%
% exwindur = 4; %s
% STATEtimepoints = Restrict(lfp.timestamps,double(SleepState.ints.(statenames{ss})));
% samplewin = STATEtimepoints(randi(length(STATEtimepoints))) + [0 exwindur];
% %%
% figure
% bz_MultiLFPPlot( lfp,'spikes',spikes,'timewin',samplewin,...
%     'sortmetric',summstats.(statenames{ss}).meanCV2,...
%     'cellgroups',{CellClass.pI,CellClass.pE})


end
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