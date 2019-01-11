

datasetPath = '/mnt/proraidDL/Database/BWCRCNS';
ISIstats = bz_LoadCellinfo(datasetPath,'ISIStats','dataset',true,'catall',true);
%%
CellClass = bz_LoadCellinfo(datasetPath,'CellClass','dataset',true,'catall',true);
%
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsAnalysis'];
%%
statenames = fieldnames(ISIstats.summstats);
statecolors = {'k','b','r'};
numstates = length(statenames);
numcells = length(CellClass.UID);

%%
sorttypes = {'rate','ISICV','CV2'};
%Make the cell-type specific sortings
for ss = 1:length(statenames)
    [~,sorts.(statenames{ss}).rate]=sort(ISIstats.summstats.(statenames{ss}).meanrate);
    [~,sorts.(statenames{ss}).ISICV]=sort(ISIstats.summstats.(statenames{ss}).ISICV);
    [~,sorts.(statenames{ss}).CV2]=sort(ISIstats.summstats.(statenames{ss}).meanCV2);
    
    classnames = unique(CellClass.label);
    numclasses = length(classnames);
    for cl = 1:numclasses
        inclasscells{cl} = strcmp(classnames{cl},CellClass.label);
        
        for tt = 1:length(sorttypes)
        sorts.(statenames{ss}).([sorttypes{tt},classnames{cl}]) = ...
            intersect(sorts.(statenames{ss}).(sorttypes{tt}),find(inclasscells{cl}),'stable');
        
        if cl==1
            sorts.(statenames{ss}).([sorttypes{tt},'byclass'])=[];
        end
        sorts.(statenames{ss}).([sorttypes{tt},'byclass']) = ...
            [sorts.(statenames{ss}).([sorttypes{tt},'byclass']) sorts.(statenames{ss}).([sorttypes{tt},classnames{cl}])];
        end
            
    end  
end

%Calculate mean ISI dists by state and cell type
meanISIhist.logbins = ISIstats.ISIhist.logbins(1,:);
for ss = 1:length(statenames)
   for cc = 1:length(classnames)
       meanISIhist.(statenames{ss}).(classnames{cc}) = ...
           nanmean(ISIstats.ISIhist.(statenames{ss}).log(CellClass.(classnames{cc}),:),1);
       meanISIhist.std.(statenames{ss}).(classnames{cc}) = ...
           nanstd(ISIstats.ISIhist.(statenames{ss}).log(CellClass.(classnames{cc}),:),[],1);

   end
end
%% CV2-rate correlation
for cl = 1:numclasses 
    for ss = 1:length(statenames)
    [rateCV2corr.(statenames{ss}).(classnames{cl}).rho,rateCV2corr.(statenames{ss}).(classnames{cl}).p]=...
        corr(log10(ISIstats.summstats.(statenames{ss}).meanrate(CellClass.(classnames{cl})))',...
            ISIstats.summstats.(statenames{ss}).meanCV2(CellClass.(classnames{cl}))',...
            'type','spearman','rows','complete')
    end
end
%%
%excells = 981, 869 559 613 513 585
excells = [552 281 356 932];
histcolors = flipud(gray);
figure
for ss = 1:length(statenames)
subplot(3,3,ss)
    plot(log10(ISIstats.summstats.(statenames{ss}).meanrate(CellClass.pE)),...
        ISIstats.summstats.(statenames{ss}).meanCV2(CellClass.pE),'k.','markersize',4)
    hold on
    plot(log10(ISIstats.summstats.(statenames{ss}).meanrate(CellClass.pI)),...
        ISIstats.summstats.(statenames{ss}).meanCV2(CellClass.pI),'r.','markersize',4)
    plot(log10(ISIstats.summstats.(statenames{ss}).meanrate(excells)),...
        ISIstats.summstats.(statenames{ss}).meanCV2(excells),...
        'o','color',[0.1 0.7 0],'markersize',5,'LineWidth',2)
    LogScale('x',10)
    xlim([-2.2 1.7]); ylim([0.4 1.6])
    plot(get(gca,'xlim'),[1 1],'k')
    title(statenames{ss})
    xlabel('FR (Hz)');ylabel('<CV2>')
    
    
    
    %cc=1;
    for cc = 1:length(excells)
        subplot(6,6,((cc+1).*6)+2.*ss-0.5)
        colormap(gca,histcolors)
            imagesc(ISIstats.ISIhist.(statenames{ss}).return(:,:,excells(cc)))
            set(gca,'ytick',[]);set(gca,'xtick',[])
            axis xy
    end
end

NiceSave('RateandCV2',figfolder,[])

%%
figure
for ss = 1:length(statenames)
subplot(2,3,ss)
colormap(histcolors)
   % subplot(2,3,4)
        imagesc((ISIstats.ISIhist.logbins(1,:)),[1 numcells],...
            ISIstats.ISIhist.(statenames{ss}).log(sorts.(statenames{ss}).ratebyclass,:))
        hold on
        plot(log10(1./(ISIstats.summstats.(statenames{ss}).meanrate(sorts.(statenames{ss}).ratebyclass))),...
            [1:numcells],'k.','markersize',1)
        plot(ISIstats.ISIhist.logbins([1 end]),sum(inclasscells{1}).*[1 1]+0.5,'r')
        LogScale('x',10)
        xlabel('ISI (s)')
        xlim(ISIstats.ISIhist.logbins([1 end]))
        %colorbar
      %  legend('1/Mean Firing Rate (s)','location','southeast')
        ylabel('Cell (Sorted by FR, Type)')
        %legend('1/Mean Firing Rate (s)','location','southeast')
        caxis([0 0.1])
        %title('ISI Distribution (Log Scale)')
        title(statenames{ss})
        
        
    subplot(2,3,ss+3)
        imagesc((ISIstats.ISIhist.logbins(1,:)),[1 numcells],...
            ISIstats.ISIhist.(statenames{ss}).log(sorts.(statenames{ss}).CV2byclass,:))
        hold on
        plot(log10(1./(ISIstats.summstats.(statenames{ss}).meanrate(sorts.(statenames{ss}).CV2byclass))),...
            [1:numcells],'k.','markersize',1)
        plot(ISIstats.ISIhist.logbins([1 end]),sum(inclasscells{1}).*[1 1]+0.5,'r')
        LogScale('x',10)
        xlabel('ISI (s)')
        xlim(ISIstats.ISIhist.logbins([1 end]))
        title(statenames{ss})
        %colorbar
      %  legend('1/Mean Firing Rate (s)','location','southeast')
        ylabel('Cell (Sorted by CV2, Type)')
        %legend('1/Mean Firing Rate (s)','location','southeast')
        caxis([0 0.1])
end

NiceSave('ISIdistssorted',figfolder,[])

%%
%Get 3 random E cells and 1 I cell
%Sort the E cells by rate
exE = randsample(find(CellClass.pE),3);
rates = ISIstats.summstats.NREMstate.meanrate(exE);
[~,sortedrateidx] = sort(rates);
excell = [exE(sortedrateidx) randsample(find(CellClass.pI),1)]
figure
for ee = 1:4
    subplot(8,4,4.*ee)
    hold on
    for ss = 1:2 
        plot(ISIstats.ISIhist.logbins(1,:),ISIstats.ISIhist.(statenames{ss}).log(excell(ee),:),...
            'color',statecolors{ss})
        plot(log10(1./ISIstats.summstats.(statenames{ss}).meanrate(excell(ee))).*[1 1],...
            [0 0.12],'color',statecolors{ss});
    end
    box off
    xlim([-3 2.25]);ylim([0 0.12])
    LogScale('x',10)
end


for ss = 1:2
    subplot(3,3,ss*3-2)
    colormap(histcolors)
       % subplot(2,3,4)
            imagesc((ISIstats.ISIhist.logbins(1,:)),[1 numcells],...
                ISIstats.ISIhist.(statenames{ss}).log(sorts.(statenames{ss}).ratebyclass,:))
            hold on
            plot(log10(1./(ISIstats.summstats.(statenames{ss}).meanrate(sorts.(statenames{ss}).ratebyclass))),...
                [1:numcells],'k.','markersize',1)
            plot(ISIstats.ISIhist.logbins([1 end]),sum(inclasscells{1}).*[1 1]+0.5,'r')
            LogScale('x',10)
            xlabel('ISI (s)')
            xlim(ISIstats.ISIhist.logbins([1 end]))
            %colorbar
          %  legend('1/Mean Firing Rate (s)','location','southeast')
            ylabel('Cell (Sorted by FR, Type)')
            %legend('1/Mean Firing Rate (s)','location','southeast')
            caxis([0 0.1])
            %title('ISI Distribution (Log Scale)')
            title(statenames{ss})
end

for cc = 1:length(classnames)
    subplot(8,4,20+4.*cc)
    hold on
    for ss = 1:2
        errorshade(meanISIhist.logbins,meanISIhist.(statenames{ss}).(classnames{cc}),...
            meanISIhist.std.(statenames{ss}).(classnames{cc}),meanISIhist.std.(statenames{ss}).(classnames{cc}),...
            statecolors{ss},'scalar')
    end
    for ss = 1:2
        
        plot(meanISIhist.logbins,meanISIhist.(statenames{ss}).(classnames{cc}),...
            'color',statecolors{ss},'linewidth',1)
    end
    
        axis tight
        xlim([-3 2.25])
        %ylim([0 0.087])
        LogScale('x',10)
end
NiceSave('ISIfig',figfolder,[])

%%
plotstates = {'WAKEstate','REMstate','WAKEstate'};
plotstates2 = {'NREMstate','NREMstate','REMstate'};
%% CV/CV2 by state
figure

for ss=1:numstates
    subplot(3,3,ss)
        plot(log10(ISIstats.summstats.(plotstates{ss}).meanrate(CellClass.pE)),...
            log10(ISIstats.summstats.(plotstates2{ss}).meanrate(CellClass.pE)),...
            'k.','markersize',4)
        hold on
        plot(log10(ISIstats.summstats.(plotstates{ss}).meanrate(CellClass.pI)),...
            log10(ISIstats.summstats.(plotstates2{ss}).meanrate(CellClass.pI)),...
            'r.','markersize',4)
        plot(log10([0.03 30]),log10([0.03 30]),'k')
        xlabel([plotstates{ss},' Rate']);ylabel([plotstates2{ss},' Rate'])
        LogScale('xy',10)
end

for ss=1:numstates
    subplot(3,3,ss+3)
        plot(log2(ISIstats.summstats.(plotstates{ss}).ISICV(CellClass.pE)),...
            log2(ISIstats.summstats.(plotstates2{ss}).ISICV(CellClass.pE)),'k.','markersize',4)
        hold on
        plot(log2(ISIstats.summstats.(plotstates{ss}).ISICV(CellClass.pI)),...
            log2(ISIstats.summstats.(plotstates2{ss}).ISICV(CellClass.pI)),'r.','markersize',4)
        plot(log2([1 6]),log2([1 6]),'k')
        xlabel([plotstates{ss},' CV']);ylabel([plotstates2{ss},' CV'])
        LogScale('xy',2)
end

for ss=1:numstates
    subplot(3,3,ss+6)
        plot((ISIstats.summstats.(plotstates{ss}).meanCV2(CellClass.pE)),...
            (ISIstats.summstats.(plotstates2{ss}).meanCV2(CellClass.pE)),'k.','markersize',4)
        hold on
        plot((ISIstats.summstats.(plotstates{ss}).meanCV2(CellClass.pI)),...
            (ISIstats.summstats.(plotstates2{ss}).meanCV2(CellClass.pI)),'r.','markersize',4)
        plot(([0 2]),([0 2]),'k')
        
        xlabel([plotstates{ss},' CV2']);ylabel([plotstates2{ss},' CV2'])
        xlim([0.4 1.6]);ylim([0.4 1.6])
       % LogScale('xy',2)
end

NiceSave('ISIstatsbystate',figfolder,[])
%%
figure
colormap(histcolors)
ff=0;
for cc = 1:numcells
    cellnum = sorts.NREMstate.CV2byclass(cc);   %%sortrate.NREMstate(cc);
    subplot(6,7,mod(cc-1,42)+1)
    imagesc((ISIstats.ISIhist.NREMstate.return(:,:,cellnum)))
    hold on
    plot(log10(1./ISIstats.summstats.NREMstate.meanrate(cellnum)),log10(1./ISIstats.summstats.NREMstate.meanrate(cellnum)),'k+')
    axis xy
    LogScale('xy',10)
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    %caxis([0 0.003])
    %xlim(ISIstats.ISIhist.logbins([1 end]));ylim(ISIstats.ISIhist.logbins([1 end]))
    %xlabel(['FR: ',num2str(round(ISIstats.summstats.NREMstate.meanrate(cellnum),2)),'Hz'])
    title([num2str(round(ISIstats.summstats.NREMstate.meanCV2(cellnum),2))])
    if mod(cc,42) == 0 || cc ==numcells
        ff= ff+1;
        NiceSave(['ISIreturnmap',num2str(ff)],figfolder,[]);
        figure
        colormap(histcolors)
    end
end
close

%% Manually Classify ISI Types


%% %% Measuring Similarity between ISI maps (same cell different state)
% %Need to run PCA denoise on all maps
%     linearizedreturn1 = reshape(ISIstats.ISIhist.NREMstate.return,[],numcells);
%     linearizedreturn2 =reshape(ISIstats.ISIhist.WAKEstate.return,[],numcells);
% X = [reshape(ISIstats.ISIhist.NREMstate.return,[],numcells) ,...
%     reshape(ISIstats.ISIhist.WAKEstate.return,[],numcells),...
%     reshape(ISIstats.ISIhist.REMstate.return,[],numcells)]';
% 
% X(isnan(X))=0;
% [Xdn, sigma, npars, u, vals, v] = bz_PCAdenoise(X);
% 
% cleanmaps = reshape(Xdn',...
%     size(ISIstats.ISIhist.NREMstate.return,1),size(ISIstats.ISIhist.NREMstate.return,2),...
%     size(X,1));
% 
% %%
% cellnum = 20;
% figure
% subplot(2,2,1)
% imagesc(ISIstats.ISIhist.NREMstate.return(:,:,cellnum))
% subplot(2,2,2)
% imagesc(cleanmaps(:,:,cellnum))
% %%
% %Should do this with the PCA denoised maps!
% %for ss=1:numstates
%     linearizedreturn1 = reshape(ISIstats.ISIhist.NREMstate.return,[],numcells);
%     linearizedreturn2 =reshape(ISIstats.ISIhist.WAKEstate.return,[],numcells);
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
