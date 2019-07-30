reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
regions = {'THAL','vCTX','fCTX','CA1'};
%regions = {'fCTX'};
%%
for rr = 1:length(regions)
    [ISIstats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);
    numcells.(regions{rr}) = length(CellClass.(regions{rr}).UID);
end

%%
statenames = fieldnames(ISIstats.(regions{1}).summstats);
statecolors = {[0 0 0],[0 0 1],[1 0 0]};
numstates = length(statenames);


%%
sorttypes = {'rate','ISICV','CV2'};
numperciles = 6;
%Make the cell-type specific sortings
for rr = 1:length(regions)
    for ss = 1:length(statenames)
        [~,sorts.(regions{rr}).(statenames{ss}).rate]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate);
        [~,sorts.(regions{rr}).(statenames{ss}).ISICV]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).ISICV);
        [~,sorts.(regions{rr}).(statenames{ss}).CV2]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanCV2);

        
        %Check for empty cell class entries
        noclass = cellfun(@isempty,CellClass.(regions{rr}).label);
        sorts.(regions{rr}).numclassycells = sum(~noclass);
        %cellclass(noclass)={'none'};
        classnames = unique(CellClass.(regions{rr}).label(~noclass));
        numclasses = length(classnames);

        for cl = 1:numclasses
            inclasscells.(regions{rr}){cl} = ...
                strcmp(classnames{cl},CellClass.(regions{rr}).label);

            for tt = 1:length(sorttypes)
            sorts.(regions{rr}).(statenames{ss}).([sorttypes{tt},classnames{cl}]) = ...
                intersect(sorts.(regions{rr}).(statenames{ss}).(sorttypes{tt}),...
                find(inclasscells.(regions{rr}){cl}),'stable');

                if cl==1
                    sorts.(regions{rr}).(statenames{ss}).([sorttypes{tt},'byclass'])=[];
                   
                end
            sorts.(regions{rr}).(statenames{ss}).([sorttypes{tt},'byclass']) = ...
                [sorts.(regions{rr}).(statenames{ss}).([sorttypes{tt},'byclass']),...
                sorts.(regions{rr}).(statenames{ss}).([sorttypes{tt},classnames{cl}])];
            end

        end
        
    end
end

%Add pE sextiles to classes


%Calculate mean ISI dists by state and cell type
meanISIhist.logbins = ISIstats.(regions{1}).ISIhist.logbins(1,:);
meanCV2hist.bins = ISIstats.(regions{1}).CV2hist.bins(1,:);

for rr = 1:length(regions)
    for ss = 1:length(statenames)
        
        noclass = cellfun(@isempty,CellClass.(regions{rr}).label);
        classnames = unique(CellClass.(regions{rr}).label(~noclass));
        percilenames = {};
        [percidx,edg] = discretize(1:length(sorts.(regions{rr}).(statenames{ss}).ratepE),...
            linspace(1,length(sorts.(regions{rr}).(statenames{ss}).ratepE),numperciles+1));
        for pp = 1:numperciles
            classnames = [classnames,['P',num2str(pp)]];
            percilenames = [percilenames ['P',num2str(pp)]];
            CellClass.(regions{rr}).(['P',num2str(pp)]) = sorts.(regions{rr}).(statenames{ss}).ratepE(percidx==pp);
        end
        
        
       for cc = 1:length(classnames)
           meanISIhist.(regions{rr}).(statenames{ss}).(classnames{cc}) = ...
               nanmean(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(CellClass.(regions{rr}).(classnames{cc}),:),1);
           meanISIhist.(regions{rr}).std.(statenames{ss}).(classnames{cc}) = ...
               nanstd(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(CellClass.(regions{rr}).(classnames{cc}),:),[],1);

           meannormISIhist.(regions{rr}).(statenames{ss}).(classnames{cc}) = ...
               nanmean(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).meannorm(CellClass.(regions{rr}).(classnames{cc}),:),1);
           
           meanreturnhist.(regions{rr}).(statenames{ss}).(classnames{cc}) = ...
               nanmean(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).return(:,:,CellClass.(regions{rr}).(classnames{cc})),3);
           meanreturnhist.(regions{rr}).std.(statenames{ss}).(classnames{cc}) = ...
               nanstd(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).return(:,:,CellClass.(regions{rr}).(classnames{cc})),[],3);
            
           meanCV2hist.(regions{rr}).(statenames{ss}).(classnames{cc}) = ...
               nanmean(ISIstats.(regions{rr}).CV2hist.(statenames{ss})(CellClass.(regions{rr}).(classnames{cc}),:),1);
           meanCV2hist.(regions{rr}).std.(statenames{ss}).(classnames{cc}) = ...
               nanstd(ISIstats.(regions{rr}).CV2hist.(statenames{ss})(CellClass.(regions{rr}).(classnames{cc}),:),[],1);
           
           meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).log = ...
               squeeze(nanmean(ISIstats.(regions{rr}).Jointhist.(statenames{ss}).log(CellClass.(regions{rr}).(classnames{cc}),:,:),1));
           meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).norm = ...
               squeeze(nanmean(ISIstats.(regions{rr}).Jointhist.(statenames{ss}).norm(CellClass.(regions{rr}).(classnames{cc}),:,:),1));
       
       end
    end
end
%% CV2-rate correlation
for rr = 1:length(regions)
    for cl = 1:numclasses 
        if ~any(CellClass.(regions{rr}).(classnames{cl}))
            rateCV2corr.(regions{rr}).(statenames{ss}).(classnames{cl}).rho = [];
            rateCV2corr.(regions{rr}).(statenames{ss}).(classnames{cl}).p = [];
            continue
        end
        for ss = 1:length(statenames)
        [rateCV2corr.(regions{rr}).(statenames{ss}).(classnames{cl}).rho,rateCV2corr.(regions{rr}).(statenames{ss}).(classnames{cl}).p]=...
            corr(log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).(classnames{cl})))',...
                ISIstats.(regions{rr}).summstats.(statenames{ss}).meanCV2(CellClass.(regions{rr}).(classnames{cl}))',...
                'type','spearman','rows','complete');
        end
    end
end


%%
%excells = 981, 869 559 613 513 585
%excells = [281 552 356 932];
exE = randsample(find(CellClass.(regions{rr}).pE),3);
rates = ISIstats.(regions{rr}).summstats.NREMstate.meanrate(exE);
[~,sortedrateidx] = sort(rates);
excells = [exE(sortedrateidx) randsample(find(CellClass.(regions{rr}).pI),1)];

histcolors = flipud(gray);
figure
for rr = 1:length(regions)
    for ss = 1:3
    subplot(4,4,ss+(rr-1).*4)
        plot(log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).pE)),...
            ISIstats.(regions{rr}).summstats.(statenames{ss}).meanCV2(CellClass.(regions{rr}).pE),'k.','markersize',4)
        hold on
        plot(log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).pI)),...
            ISIstats.(regions{rr}).summstats.(statenames{ss}).meanCV2(CellClass.(regions{rr}).pI),'r.','markersize',4)
%         plot(log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(excells)),...
%             ISIstats.(regions{rr}).summstats.(statenames{ss}).meanCV2(excells),...
%             'o','color',[0.1 0.7 0],'markersize',5,'LineWidth',2)
        xlim([-2.2 1.8]); ylim([0.4 1.6])
        LogScale('x',10)
        plot(get(gca,'xlim'),[1 1],'k')
        if rr == 1
        title(statenames{ss})
        end
        if rr==3
            xlabel('FR (Hz)');
        else
            set(gca,'xticklabel',[])
        end
        if ss==1
        ylabel({regions{rr},'<CV2>'})
        else
            set(gca,'yticklabel',[])
        end




    end
end
NiceSave('RateandCV2',figfolder,[])


%%

figure
for rr = 1:length(regions)
    for ss = 1:3

    subplot(4,4,ss+(rr-1).*4)
        plot(log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).pE)),...
            log2(ISIstats.(regions{rr}).summstats.(statenames{ss}).ISICV(CellClass.(regions{rr}).pE)),'k.','markersize',4)
        hold on
        plot(log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).pI)),...
            log2(ISIstats.(regions{rr}).summstats.(statenames{ss}).ISICV(CellClass.(regions{rr}).pI)),'r.','markersize',4)
%         plot(log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(excells)),...
%             log2(ISIstats.(regions{rr}).summstats.(statenames{ss}).ISICV(excells)),...
%             'o','color',[0.1 0.7 0],'markersize',5,'LineWidth',2)

        xlim([-2.2 1.8]); ylim([-1 5])
        LogScale('x',10);LogScale('y',2);
        plot(get(gca,'xlim'),[0 0],'k')
        if rr == 1
            title(statenames{ss})
        end
        if rr==3
            xlabel('FR (Hz)');
        else
            set(gca,'xticklabel',[])
        end
        if ss==1
            ylabel({regions{rr},'CV'})
        else
            set(gca,'yticklabel',[])
        end
    end
end
NiceSave('RateandCV',figfolder,[])


%%
%Get 3 random E cells and 1 I cell
%Sort the E cells by rate

%excell = excells;
NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);
REMhistcolors = makeColorMap([1 1 1],[0.8 0 0]);
statecolormap = {histcolors,NREMhistcolors,REMhistcolors};




figure
for rr = 1:length(regions)
for ss = 1:3
    subplot(3,4,ss*4-3+(rr-1))
    colormap(gca,statecolormap{ss})

       % subplot(2,3,4)
            imagesc((ISIstats.(regions{rr}).ISIhist.logbins(1,:)),[1 sorts.(regions{rr}).numclassycells],...
                ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(sorts.(regions{rr}).(statenames{ss}).ratebyclass,:))
            hold on
            plot(log10(1./(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(sorts.(regions{rr}).(statenames{ss}).ratebyclass))),...
                [1:sorts.(regions{rr}).numclassycells],'k.','markersize',1)
            plot(ISIstats.(regions{rr}).ISIhist.logbins([1 end]),sum(inclasscells.(regions{rr}){1}).*[1 1]+0.5,'r')
            
            plot(meanISIhist.logbins,-meanISIhist.(regions{rr}).(statenames{ss}).pE*5000+...
                sum(inclasscells.(regions{rr}){1})+0.5,...
                'color',statecolors{ss},'linewidth',2)
            
            plot(meanISIhist.logbins,-meanISIhist.(regions{rr}).(statenames{ss}).pI*2000+...
                sum(sorts.(regions{rr}).numclassycells),...
                'color',statecolors{ss},'linewidth',2)
            
            xlim(ISIstats.(regions{rr}).ISIhist.logbins([1 end]))
            xlim([-3 1.9])
            LogScale('x',10,'exp',true)
            if ss==3
                xlabel('ISI (s)')
            else
                set(gca,'xticklabels',[])
            end
            %colorbar
          %  legend('1/Mean Firing Rate (s)','location','southeast')
          if rr ==1
            ylabel({statenames{ss},'Cell'})
          end
            set(gca,'yticklabel',[])
            %legend('1/Mean Firing Rate (s)','location','southeast')
            caxis([0 0.1])
            %title('ISI Distribution (Log Scale)')
            if ss==1
                title(regions{rr})
            end

                
end



% subplot(8,4,4.*cc-3+(rr-1)) %Mean ISIHist
%     hold on
%     for ss = 1:2
%         errorshade(meanISIhist.logbins,meanISIhist.(regions{rr}).(statenames{ss}).(classnames{cc}),...
%             meanISIhist.(regions{rr}).std.(statenames{ss}).(classnames{cc}),meanISIhist.(regions{rr}).std.(statenames{ss}).(classnames{cc}),...
%             statecolors{ss},'scalar')
%     end
%     for ss = 1:2
%         
%         plot(meanISIhist.logbins,meanISIhist.(regions{rr}).(statenames{ss}).(classnames{cc}),...
%             'color',statecolors{ss},'linewidth',2)
%     end
%     
%         axis tight
%         yrange = get(gca,'ylim');ylim([0 yrange(2)])
%         xlim([-3 2.25])
%         %ylim([0 0.087])
%         LogScale('x',10)
%         set(gca,'ytick',[])
%         if rr == 1
%             ylabel(classnames{cc})
%         end
%         if cc ==2
%             xlabel('ISI (s)')
%         else
%             set(gca,'xticklabel',[])
%             title(regions{rr})
%         end


end

NiceSave('ISIDists',figfolder,[])

%%
figure
for rr = 1:length(regions)


for cc = 1:length(classnames)
    
        
	for ss = 1:3
        subplot(6,6,rr+(ss-1)*6+(cc-1)*18)    
        colormap(gca,statecolormap{ss})

            imagesc(ISIstats.(regions{rr}).ISIhist.logbins(1,:),...
                ISIstats.(regions{rr}).ISIhist.logbins(1,:),...
                meanreturnhist.(regions{rr}).(statenames{ss}).(classnames{cc}))
            axis xy
            set(gca,'ytick',[]);set(gca,'xtick',[]);
            if ss==1 & cc==1
                title(regions{rr})
            elseif ss==3 
                if cc ==2
                xlabel('ISI_n (s)')
                end
                set(gca,'xtick',[-2:1]);
                LogScale('x',10,'exp',true)
            end
            if rr==1 
                ylabel('ISI_n_+_1 (s)')
                set(gca,'ytick',[-2:1]);
                LogScale('y',10,'exp',true)
            end
            %LogScale('xy',10)
	end

end

end

NiceSave('ISIReturnMap',figfolder,[])
%% CV2 figure

figure
for rr = 1:length(regions)
for ss = 1:3
    subplot(3,4,ss*4-3+(rr-1))
    colormap(gca,statecolormap{ss})

       % subplot(2,3,4)
            imagesc((ISIstats.(regions{rr}).CV2hist.bins(1,:)),[1 sorts.(regions{rr}).numclassycells],...
                ISIstats.(regions{rr}).CV2hist.(statenames{ss})(sorts.(regions{rr}).(statenames{ss}).ratebyclass,:))
            hold on
            plot((ISIstats.(regions{rr}).summstats.(statenames{ss}).meanCV2(sorts.(regions{rr}).(statenames{ss}).ratebyclass)),...
                [1:sorts.(regions{rr}).numclassycells],'k.','markersize',1)
            plot(ISIstats.(regions{rr}).CV2hist.bins([1 end]),sum(inclasscells.(regions{rr}){1}).*[1 1]+0.5,'r')
            
            plot(meanCV2hist.bins,-meanCV2hist.(regions{rr}).(statenames{ss}).pE*7000+...
                sum(inclasscells.(regions{rr}){1})+0.5,...
                'color',statecolors{ss},'linewidth',2)
            
            plot(meanCV2hist.bins,-meanCV2hist.(regions{rr}).(statenames{ss}).pI*4000+...
                sum(numcells.(regions{rr})),...
                'color',statecolors{ss},'linewidth',2)
            
            xlim(ISIstats.(regions{rr}).CV2hist.bins([1 end]))
            %LogScale('x',10)
            if ss==3
                xlabel('CV2 (s)')
            else
                set(gca,'xticklabels',[])
            end
            %colorbar
          %  legend('1/Mean Firing Rate (s)','location','southeast')
          if rr ==1
            ylabel({statenames{ss},'Cell'})
          end
            set(gca,'yticklabel',[])
            %legend('1/Mean Firing Rate (s)','location','southeast')
            caxis([0 0.1])
            %title('ISI Distribution (Log Scale)')
            if ss==1
                title(regions{rr})
            end

                
end

end

NiceSave('CV2fig',figfolder,[])

%% CV2 figure

figure
for rr = 1:length(regions)

for cc = 1:length(classnames)
	for ss = 1:3
        subplot(6,4,rr+(ss-1)*4+(cc-1)*12)    
        %colormap(gca,statecolormap{ss})

            imagesc(meanISIhist.logbins,ISIstats.(regions{rr}).CV2hist.bins(1,:),...
                meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).log')
            hold on
            plot(meanISIhist.logbins,meanISIhist.(regions{rr}).(statenames{ss}).(classnames{cc})*25,...
                'color',statecolors{ss},'linewidth',1)
            
            axis xy
            set(gca,'ytick',[]);set(gca,'xtick',[]);
            if ss==1 & cc==1
                title(regions{rr})
            elseif ss==3 
                if cc ==2
                xlabel('ISI (s)')
                end
                set(gca,'xtick',[-2:1]);
                LogScale('x',10)
            end
            if rr==1 
                ylabel('CV2')
                set(gca,'ytick',[0 1 2]);
            end
            ylim([0 2]);
            xlim([-2.5 1.7])
            xlim(ISIstats.(regions{rr}).ISIhist.logbins([1 end]))
            
            caxis([0 max([meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).log(:,1);0])])
            switch cc
                case 1
                   % caxis([0.5e-4 1.05e-3])
                case 2
                   % caxis([0.5e-4 1.8e-3])
            end
            
    end
end
end

NiceSave('JointCV2ISI',figfolder,[])


%% FR %Ile FIgure
figure
for rr = 1:length(regions)
for ss = 1:3
    pcolor = makeColorMap([0.7 0.7 0.7],statecolors{ss},numperciles);

    subplot(5,4,(rr-1)+(ss-1)*4+1)
        hold on
        for cc = 1:length(percilenames)
            plot(meanISIhist.logbins,meanISIhist.(regions{rr}).(statenames{ss}).(percilenames{cc}),...
                'linewidth',1,'color',pcolor(cc,:))
        end
        axis tight
        if ss==1
            title(regions{rr})
        end
        if rr == 1
            ylabel('p(ISI)');
        end
        xlim([-3 1.9])
        LogScale('x',10,'exp',true)
        set(gca,'ytick',[])
            if ss==3
                xlabel('ISI (s)')
            else
                set(gca,'xticklabels',[])
            end
end
end


for rr=1:length(regions)
    subplot(4,4,rr+12)
        plot(log10(ISIstats.(regions{rr}).summstats.(plotstates{ss}).meanrate(CellClass.(regions{rr}).pE)),...
            log10(ISIstats.(regions{rr}).summstats.(plotstates2{ss}).meanrate(CellClass.(regions{rr}).pE)),...
            'k.','markersize',2)
        hold on
        plot(log10(ISIstats.(regions{rr}).summstats.(plotstates{ss}).meanrate(CellClass.(regions{rr}).pI)),...
            log10(ISIstats.(regions{rr}).summstats.(plotstates2{ss}).meanrate(CellClass.(regions{rr}).pI)),...
            'r.','markersize',2)
        %plot(log10([0.03 30]),log10([0.03 30]),'k')
        xlabel([plotstates{ss},' Rate']);ylabel([plotstates2{ss},' Rate'])
        %axis tight
        title(regions{rr})
        xlim([-2 2]);ylim([-2 2])
        LogScale('xy',10)
        UnityLine('linetype','-')
end


NiceSave('Percentiles',figfolder,[])
%%

for ss = 1:3
   figure
for rr = 1:length(regions)
   for cc = 1:length(percilenames)
        subplot(length(percilenames),4,(cc-1)*4+rr)    
        %colormap(gca,statecolormap{ss})

            imagesc(meanISIhist.logbins,ISIstats.(regions{rr}).CV2hist.bins(1,:),...
                meanJointhist.(regions{rr}).(statenames{ss}).(percilenames{cc}).log')
            hold on
            plot(meanISIhist.logbins,meanISIhist.(regions{rr}).(statenames{ss}).(percilenames{cc})*25,...
                'color',statecolors{ss},'linewidth',1)
            
            axis xy
            xlim([-3 1.9])
            set(gca,'ytick',[]);set(gca,'xticklabel',[]);
            if cc==1
                title(regions{rr})
            end
            if cc==length(percilenames) 
                xlabel('ISI (s)')
        
                LogScale('x',10,'exp',true)
            end
            if rr==1 
                ylabel('CV2')
                set(gca,'ytick',[0 1 2]);
            end
            ylim([0 2]);
            xlim([-2.5 1.7])
            xlim(ISIstats.(regions{rr}).ISIhist.logbins([1 end]))
            
            caxis([0 max([meanJointhist.(regions{rr}).(statenames{ss}).(percilenames{cc}).log(:,1);0])])

    end
end

NiceSave(['Percentiles_',(statenames{ss})],figfolder,[])
end
%% Interneuron Figure

figure
for rr = 2:length(regions)
for ss = 1:3
    subplot(6,4,ss*4-3+(rr-1))
    colormap(gca,statecolormap{ss})

       % subplot(2,3,4)
            imagesc((ISIstats.(regions{rr}).ISIhist.logbins(1,:)),[1 length(sorts.(regions{rr}).(statenames{ss}).ratepI)],...
                ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(sorts.(regions{rr}).(statenames{ss}).ratepI,:))
            hold on
            plot(log10(1./(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(sorts.(regions{rr}).(statenames{ss}).ratepI))),...
                [1:length(sorts.(regions{rr}).(statenames{ss}).ratepI)],'k.','markersize',1)
            
            
            plot(meanISIhist.logbins,-meanISIhist.(regions{rr}).(statenames{ss}).pI*2000+...
                length(sorts.(regions{rr}).(statenames{ss}).ratepI),...
                'color',statecolors{ss},'linewidth',2)
            
            xlim(ISIstats.(regions{rr}).ISIhist.logbins([1 end]))
            xlim([-3 1.9])
            LogScale('x',10,'exp',true)
            if ss==3
                xlabel('ISI (s)')
            else
                set(gca,'xticklabels',[])
            end
            %colorbar
          %  legend('1/Mean Firing Rate (s)','location','southeast')
          if rr ==1
            ylabel({statenames{ss},'Cell'})
          end
            set(gca,'yticklabel',[])
            %legend('1/Mean Firing Rate (s)','location','southeast')
            caxis([0 0.1])
            %title('ISI Distribution (Log Scale)')
            if ss==1
                title(regions{rr})
            end

                
end
end
NiceSave('InterNeuron',figfolder,[])
%% CTX WAKE Figure

actrange = [0.01 0.05 ; 0.05 0.12; 0.08 0.2];
actrange = [0.01 0.2 ; 0.01 0.2; 0.01 0.2];
for rr = 1:3
    actStateBins = 10.^(ISIstats.(regions{rr}).ISIhist.logbins(1,:)) > actrange(rr,1) & ...
        10.^(ISIstats.(regions{rr}).ISIhist.logbins(1,:)) < actrange(rr,2);
    actStatePower = max(ISIstats.(regions{rr}).ISIhist.WAKEstate.log(:,actStateBins),[],2);
        
        [~,sortAct]=...
            sort(actStatePower);
        
        sorts.(regions{rr}).actWAKE = ...
            intersect(sortAct,...
            find(inclasscells.(regions{rr}){1}),'stable');

end

figure
for rr = 1:3

	ss = 1; cc = 1;
        subplot(6,4,(rr-1)*4+1)    
        %colormap(gca,statecolormap{ss})

            imagesc(meanISIhist.logbins,ISIstats.(regions{rr}).CV2hist.bins(1,:),...
                meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).log')
            hold on
            plot(meanISIhist.logbins,meanISIhist.(regions{rr}).(statenames{ss}).(classnames{cc})*25,...
                'color',statecolors{ss},'linewidth',1)
            
            axis xy
            %set(gca,'ytick',[]);%set(gca,'xtick',[]);
             ylabel({(regions{rr}),'CV2'})
             set(gca,'xticklabel',[]);
             set(gca,'ytick',[0 1 2]);
            xlim([-2.75 1.75])
            if rr ==3
                xlabel('ISI (s)')
                %set(gca,'xtick',[-2:1]);
                LogScale('x',10,'exp',true)
            end
            if rr==1 
                title('WAKE')
            end
            ylim([0 2]);
            
            %xlim(ISIstats.(regions{rr}).ISIhist.logbins([1 end]))
            
            caxis([0 max([meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).log(:,1);0])])


    subplot(5,4,(rr-1)*4+3) 
    colormap(gca,statecolormap{ss})

       % subplot(2,3,4)
            imagesc((ISIstats.(regions{rr}).ISIhist.logbins(1,:)),[1 sum(inclasscells.(regions{rr}){1})],...
                ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(sorts.(regions{rr}).actWAKE,:))
            hold on
%             plot(log10(1./(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(sorts.(regions{rr}).(statenames{ss}).ratebyclass))),...
%                 [1:sorts.(regions{rr}).numclassycells],'k.','markersize',1)

            
            plot(meanISIhist.logbins,-meanISIhist.(regions{rr}).(statenames{ss}).pE*5000+...
                sum(inclasscells.(regions{rr}){1})+0.5,...
                'color',statecolors{ss},'linewidth',2)
            

            
            xlim(ISIstats.(regions{rr}).ISIhist.logbins([1 end]))
            xlim([-3 1.9])
            xlim([-2.75 1.75])
            LogScale('x',10,'exp',true)
            if rr==3
                xlabel('ISI (s)')
            else
                set(gca,'xticklabels',[])
            end
            %colorbar
          %  legend('1/Mean Firing Rate (s)','location','southeast')
         
            ylabel({(regions{rr}),'(Sort by', '10-200ms Peak)'})
          
            set(gca,'yticklabel',[])
            %legend('1/Mean Firing Rate (s)','location','southeast')
            caxis([0 0.1])
            %title('ISI Distribution (Log Scale)')
%             if ss==1
%                 title(regions{rr})
%             end

                
end
NiceSave('ThCTX_WAKE',figfolder,[])

%% Normalized ISI figure
figure
for rr = 1:length(regions)

for ss = 1:3
    subplot(3,4,ss*4-3+(rr-1))
    colormap(gca,statecolormap{ss})

       % subplot(2,3,4)
            imagesc((ISIstats.(regions{rr}).ISIhist.logbins(1,:)),[1 sorts.(regions{rr}).numclassycells],...
                ISIstats.(regions{rr}).ISIhist.(statenames{ss}).meannorm(sorts.(regions{rr}).(statenames{ss}).ratebyclass,:))
            hold on
            plot(zeros(1,sorts.(regions{rr}).numclassycells),...
                [1:sorts.(regions{rr}).numclassycells],'k.','markersize',1)
            plot(ISIstats.(regions{rr}).ISIhist.logbins([1 end]),sum(inclasscells.(regions{rr}){1}).*[1 1]+0.5,'r')
            
            plot(meanISIhist.logbins,-meannormISIhist.(regions{rr}).(statenames{ss}).pE*5000+...
                sum(inclasscells.(regions{rr}){1})+0.5,...
                'color',statecolors{ss},'linewidth',2)
            
            plot(meanISIhist.logbins,-meannormISIhist.(regions{rr}).(statenames{ss}).pI*2000+...
                sum(sorts.(regions{rr}).numclassycells),...
                'color',statecolors{ss},'linewidth',2)
            
            xlim([-3 1.5])
            LogScale('x',10,'exp',true)
            if ss==3
                xlabel('norm ISI log(mean^-^1)')
            else
                set(gca,'xticklabels',[])
            end
            %colorbar
          %  legend('1/Mean Firing Rate (s)','location','southeast')
          if rr ==1
            ylabel({statenames{ss},'Cell'})
          end
            set(gca,'yticklabel',[])
            %legend('1/Mean Firing Rate (s)','location','southeast')
            caxis([0 0.1])
            %title('ISI Distribution (Log Scale)')
            if ss==1
                title(regions{rr})
            end

                
end

end

NiceSave('ISIdist_meannorm',figfolder,[])


%% Normalized ISI figure: CV2
figure
for rr = 1:length(regions)


for cc = 1:length(classnames)
	for ss = 1:3
        subplot(6,4,rr+(ss-1)*4+(cc-1)*12)    
        %colormap(gca,statecolormap{ss})

            imagesc(meanISIhist.logbins,ISIstats.(regions{rr}).CV2hist.bins(1,:),...
                meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).norm')
            hold on
            plot(meanISIhist.logbins,meannormISIhist.(regions{rr}).(statenames{ss}).(classnames{cc})*25,...
                'color',statecolors{ss},'linewidth',1)
            
            axis xy
            ylim([0 2]);
            set(gca,'ytick',[]);%set(gca,'xtick',[]);
            if ss==1 &rr==1
                title(classnames{cc})
            elseif ss==3 
                if rr ==2
                xlabel('Norm ISI (log(mean^-^1))')
                end
                %set(gca,'xtick',[-2:1]);
                %LogScale('x',10)
            end
            if cc==1 
                ylabel('CV2')
                set(gca,'ytick',[0 1 2]);
            end
            
            xlim([-3 1])
            
            switch cc
                case 1
                    switch rr
                        case 1
                            caxis([0.5e-4 1.2e-3])
                        case 2
                            if ss == 2; caxis([0.5e-4 0.8e-3])
                            else; caxis([0.5e-4 0.6e-3])
                            end
                    end
                case 2
                    caxis([0.5e-4 1.8e-3])
            end
             caxis([0 max([meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).norm(:,1);0])])
            
    end
end
end

NiceSave('CV2fig_meannorm',figfolder,[])
%% Example Figure
figure
for ee = 1:4
    subplot(6,4,4.*ee-2)
    hold on
    for ss = 1:2 
        plot(ISIstats.(regions{rr}).ISIhist.logbins(1,:),ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(excells(ee),:),...
            'color',statecolors{ss},'linewidth',2)
        plot(log10(1./ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(excells(ee))).*[1 1],...
            [0 0.12],'color',statecolors{ss});
    end
    box off
    xlim([-3 2.25]);ylim([0 0.125])
    LogScale('x',10)
        if ee==4
            xlabel('ISI (s)')
        end
        
	for ss = 1:2
        subplot(6,6,6.*ee+ss-2)    
        colormap(histcolors)
        if ss==2
            colormap(gca,NREMhistcolors)
        end

            imagesc(ISIstats.(regions{rr}).ISIhist.logbins(1,:),ISIstats.(regions{rr}).ISIhist.logbins(1,:),...
                ISIstats.(regions{rr}).ISIhist.(statenames{ss}).return(:,:,excells(ee)))
            hold on
            plot(log10(1./ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(excells(ee))),...
                log10(1./ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(excells(ee))),'r+')

            axis xy
            set(gca,'ytick',[]);set(gca,'xtick',[]);
    end

end
NiceSave('ISexamples',figfolder,[])


%%
plotstates = {'WAKEstate','REMstate','WAKEstate'};
plotstates2 = {'NREMstate','NREMstate','REMstate'};
%% CV/CV2 by state

for rr = 1:length(regions)
figure
suptitle(regions{rr})
%Rate
for ss=1:3
    subplot(4,4,ss)
        plot(log10(ISIstats.(regions{rr}).summstats.(plotstates{ss}).meanrate(CellClass.(regions{rr}).pE)),...
            log10(ISIstats.(regions{rr}).summstats.(plotstates2{ss}).meanrate(CellClass.(regions{rr}).pE)),...
            'k.','markersize',2)
        hold on
        plot(log10(ISIstats.(regions{rr}).summstats.(plotstates{ss}).meanrate(CellClass.(regions{rr}).pI)),...
            log10(ISIstats.(regions{rr}).summstats.(plotstates2{ss}).meanrate(CellClass.(regions{rr}).pI)),...
            'r.','markersize',2)
        %plot(log10([0.03 30]),log10([0.03 30]),'k')
        xlabel([plotstates{ss},' Rate']);ylabel([plotstates2{ss},' Rate'])
        %axis tight
        xlim([-2 2]);ylim([-2 2])
        LogScale('xy',10)
        UnityLine('linetype','-')
end

%CV
for ss=1:3
    subplot(4,4,ss+4)
        plot(log2(ISIstats.(regions{rr}).summstats.(plotstates{ss}).ISICV(CellClass.(regions{rr}).pE)),...
            log2(ISIstats.(regions{rr}).summstats.(plotstates2{ss}).ISICV(CellClass.(regions{rr}).pE)),...
            'k.','markersize',2)
        hold on
        plot(log2(ISIstats.(regions{rr}).summstats.(plotstates{ss}).ISICV(CellClass.(regions{rr}).pI)),...
            log2(ISIstats.(regions{rr}).summstats.(plotstates2{ss}).ISICV(CellClass.(regions{rr}).pI)),...
            'r.','markersize',2)
        %plot(log2([1 6]),log2([1 6]),'k')
        xlabel([plotstates{ss},' CV']);ylabel([plotstates2{ss},' CV'])
        xlim([-0.3 5]);        ylim([-0.3 5])
        LogScale('xy',2)
        UnityLine('linetype','-')
end


%CV2
for ss=1:3
    subplot(4,4,ss+8)
        plot((ISIstats.(regions{rr}).summstats.(plotstates{ss}).meanCV2(CellClass.(regions{rr}).pE)),...
            (ISIstats.(regions{rr}).summstats.(plotstates2{ss}).meanCV2(CellClass.(regions{rr}).pE)),...
            'k.','markersize',2)
        hold on
        plot((ISIstats.(regions{rr}).summstats.(plotstates{ss}).meanCV2(CellClass.(regions{rr}).pI)),...
            (ISIstats.(regions{rr}).summstats.(plotstates2{ss}).meanCV2(CellClass.(regions{rr}).pI)),...
            'r.','markersize',2)
        %plot(([0 2]),([0 2]),'k')
        
        xlabel([plotstates{ss},' CV2']);ylabel([plotstates2{ss},' CV2'])
        
       % LogScale('xy',2)
       xlim([0.4 1.8]);ylim([0.4 1.8])
        UnityLine('linetype','-')
end

NiceSave(['ISIstatsbystate',regions{rr}],figfolder,[])
end
%%


%% All return maps
% figure
% colormap(histcolors)
% ff=0;
% for cc = 1:numcells.(regions{rr})
%     cellnum = sorts.(regions{rr}).NREMstate.CV2byclass(cc);   %%sortrate.NREMstate(cc);
%     subplot(6,7,mod(cc-1,42)+1)
%     imagesc((ISIstats.ISIhist.NREMstate.return(:,:,cellnum)))
%     hold on
%     plot(log10(1./ISIstats.summstats.NREMstate.meanrate(cellnum)),log10(1./ISIstats.summstats.NREMstate.meanrate(cellnum)),'k+')
%     axis xy
%     LogScale('xy',10)
%     set(gca,'ytick',[]);set(gca,'xtick',[]);
%     %caxis([0 0.003])
%     %xlim(ISIstats.ISIhist.logbins([1 end]));ylim(ISIstats.ISIhist.logbins([1 end]))
%     %xlabel(['FR: ',num2str(round(ISIstats.summstats.NREMstate.meanrate(cellnum),2)),'Hz'])
%     title([num2str(round(ISIstats.summstats.NREMstate.meanCV2(cellnum),2))])
%     if mod(cc,42) == 0 || cc ==numcells.(regions{rr})
%         ff= ff+1;
%         NiceSave(['ISIreturnmap',num2str(ff)],figfolder,[]);
%         figure
%         colormap(histcolors)
%     end
% end
% close

%% Manually Classify ISI Types


%% %% Measuring Similarity between ISI maps (same cell different state)
% %Need to run PCA denoise on all maps
%     linearizedreturn1 = reshape(ISIstats.ISIhist.NREMstate.return,[],numcells.(regions{rr}));
%     linearizedreturn2 =reshape(ISIstats.ISIhist.WAKEstate.return,[],numcells.(regions{rr}));
% X = [reshape(ISIstats.ISIhist.NREMstate.return,[],numcells.(regions{rr})) ,...
%     reshape(ISIstats.ISIhist.WAKEstate.return,[],numcells.(regions{rr})),...
%     reshape(ISIstats.ISIhist.REMstate.return,[],numcells.(regions{rr}))]';
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
%     linearizedreturn1 = reshape(ISIstats.ISIhist.NREMstate.return,[],numcells.(regions{rr}));
%     linearizedreturn2 =reshape(ISIstats.ISIhist.WAKEstate.return,[],numcells.(regions{rr}));
%     returnmapsimilarity = corr(linearizedreturn1,linearizedreturn2,'type','spearman');
% 
%     
%     NWsimilarity = diag(returnmapsimilarity);
%     
% %     figure
% %     imagesc(returnmapsimilarity(ISIstats.sorts.(regions{rr}).NREMstate.ratebyclass,ISIstats.sorts.(regions{rr}).NREMstate.ratebyclass))
% %     colorbar
% hist(NWsimilarity)
% figure
% plot(log10(ISIstats.summstats.NREMstate.meanrate),NWsimilarity,'.')
