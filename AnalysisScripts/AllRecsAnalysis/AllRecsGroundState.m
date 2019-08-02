reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/GroundStateAnalysis'];

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
    
    GroundStateAll = GetMatResults(figfolder,'GroundStateAnalysis','baseNames',baseNames);
    GroundStateAll = bz_CollapseStruct(GroundStateAll);
   
    ISIoccupancy.(regions{rr}) = bz_CollapseStruct(GroundStateAll.ISIoccupancy,'match',...
        'justcat',true );
    OccupancyStats.(regions{rr}) = bz_CollapseStruct(GroundStateAll.OccupancyStats,'match',...
        'justcat',true );
    normISIhist.(regions{rr}) = bz_CollapseStruct(GroundStateAll.normISIhist,'match',...
        'justcat',true );
end

%%
statenames = fieldnames(ISIstats.(regions{1}).summstats);
statecolors = {[0 0 0],[0 0 1],[1 0 0]};
numstates = length(statenames);


%% Sorts for plot
sorttypes = {'rate','medISI','MTOrat'};
%tt =1
%Make the cell-type specific sortings
%sorttypes = {'rate','ISICV','CV2'};
%Make the cell-type specific sortings
for rr = 1:length(regions)
    for ss = 1:3
        
        OccupancyStats.(regions{rr}).(statenames{ss}).MTORatio = ...
            ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate./...
            (1./OccupancyStats.(regions{rr}).(statenames{ss}).median);
        
        [~,sorts.(regions{rr}).(statenames{ss}).rate]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate);
        [~,sorts.(regions{rr}).(statenames{ss}).ISICV]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).ISICV);
        [~,sorts.(regions{rr}).(statenames{ss}).CV2]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanCV2);
        [~,sorts.(regions{rr}).(statenames{ss}).medISI]=...
            sort(1./OccupancyStats.(regions{rr}).(statenames{ss}).median);
        [~,sorts.(regions{rr}).(statenames{ss}).MTOrat]=...
            sort(OccupancyStats.(regions{rr}).(statenames{ss}).MTORatio);
        
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


%% Calculate mean ISI/occupancy dists by state and cell type
meanISIhist.logbins = ISIstats.(regions{1}).ISIhist.logbins(1,:);
meannormISIhist.bins = normISIhist.(regions{1}).bins(1,:);
meannormOcc.bins = ISIoccupancy.(regions{rr}).logbins(1,:);
numperciles = 6;

for rr = 1:length(regions)
    for ss = 1:3
        
        noclass = cellfun(@isempty,CellClass.(regions{rr}).label);
        %classnames = unique(CellClass.(regions{rr}).label(~noclass));
        classnames = {'pE','pI'};
        percilenames = {};
        [percidx,edg] = discretize(1:length(sorts.(regions{rr}).(statenames{ss}).medISIpE),...
            linspace(1,length(sorts.(regions{rr}).(statenames{ss}).medISIpE),numperciles+1));
        for pp = 1:numperciles
            classnames = [classnames,['P',num2str(pp)]];
            percilenames = [percilenames ['P',num2str(pp)]];
            CellClass.(regions{rr}).(['P',num2str(pp)]) = sorts.(regions{rr}).(statenames{ss}).medISIpE(percidx==pp);
            
            meanpercmedISI.(regions{rr}).(statenames{ss})(pp) = mean(1./OccupancyStats.(regions{rr}).(statenames{ss}).median(CellClass.(regions{rr}).(['P',num2str(pp)])));
        end
        
        
       for cc = 1:length(classnames)
           meanISIhist.(regions{rr}).(statenames{ss}).(classnames{cc}) = ...
               nanmean(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(CellClass.(regions{rr}).(classnames{cc}),:),1);
           meanISIhist.(regions{rr}).std.(statenames{ss}).(classnames{cc}) = ...
               nanstd(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(CellClass.(regions{rr}).(classnames{cc}),:),[],1);

           meannormISIhist.(regions{rr}).(statenames{ss}).(classnames{cc}) = ...
               nanmean(normISIhist.(regions{rr}).(statenames{ss}).mednorm(CellClass.(regions{rr}).(classnames{cc}),:),1);
           meannormISIhist.(regions{rr}).std.(statenames{ss}).(classnames{cc}) = ...
               nanstd(normISIhist.(regions{rr}).(statenames{ss}).mednorm(CellClass.(regions{rr}).(classnames{cc}),:),[],1);
           
           meannormOcc.(regions{rr}).(statenames{ss}).(classnames{cc}) = ...
               nanmean(ISIoccupancy.(regions{rr}).(statenames{ss}).mednormhist(:,CellClass.(regions{rr}).(classnames{cc})),2);
           meannormOcc.(regions{rr}).std.(statenames{ss}).(classnames{cc}) = ...
               nanstd(ISIoccupancy.(regions{rr}).(statenames{ss}).mednormhist(:,CellClass.(regions{rr}).(classnames{cc})),[],2);
           
           meanreturnhist.(regions{rr}).(statenames{ss}).(classnames{cc}) = ...
               nanmean(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).return(:,:,CellClass.(regions{rr}).(classnames{cc})),3);
           meanreturnhist.(regions{rr}).std.(statenames{ss}).(classnames{cc}) = ...
               nanstd(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).return(:,:,CellClass.(regions{rr}).(classnames{cc})),[],3);
           
           meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).log = ...
               squeeze(nanmean(ISIstats.(regions{rr}).Jointhist.(statenames{ss}).log(CellClass.(regions{rr}).(classnames{cc}),:,:),1));
           meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).norm = ...
               squeeze(nanmean(ISIstats.(regions{rr}).Jointhist.(statenames{ss}).norm(CellClass.(regions{rr}).(classnames{cc}),:,:),1));
       
       end
    end
end


%% Figure: Rate Sort Occupancy
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,4,(ss-1)*4+rr)
    s = imagesc(ISIoccupancy.(regions{rr}).logbins(1,:),[1 length(sorts.(regions{rr}).(statenames{ss}).ratebyclass)],...
        (ISIoccupancy.(regions{rr}).(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).ratebyclass))');
    alpha(s,single(ISIoccupancy.(regions{rr}).(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).ratebyclass)'~=0))

    hold on
    plot(log10(1./ISIstats.(regions{rr}).summstats.(state).meanrate(sorts.(regions{rr}).(statenames{ss}).ratebyclass)),...
        [1:length(sorts.(regions{rr}).(statenames{ss}).ratebyclass)],'.')
    LogScale('x',10)
    caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('ISI')
          if rr ==1
            ylabel({statenames{ss},'Cell'})
          end
            set(gca,'yticklabel',[])
    if ss==1
        title(regions{rr})
    end
    
    LogScale('x',10,'exp',true)
% subplot(2,1,2)
%     imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
%         ISIoccupancy.(state).hist(:,ISIStats.sorts.(state).ratebyclass)')
end
end

NiceSave('ISIOccupancy_ratesort',figfolder,[])

%% Median Occupancy Sort
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,4,(ss-1)*4+rr)
    s = imagesc(ISIoccupancy.(regions{rr}).logbins(1,:),[1 length(sorts.(regions{rr}).(statenames{ss}).medISIbyclass)],...
        (ISIoccupancy.(regions{rr}).(statenames{ss}).loghist(:,sorts.(regions{rr}).(statenames{ss}).medISIbyclass))');
    alpha(s,single(ISIoccupancy.(regions{rr}).(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).medISIbyclass)'~=0))

    hold on
    plot(log10(1./ISIstats.(regions{rr}).summstats.(state).meanrate(sorts.(regions{rr}).(statenames{ss}).medISIbyclass)),...
        [1:length(sorts.(regions{rr}).(statenames{ss}).medISIbyclass)],'.')
    LogScale('x',10)
    caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('ISI')
          if rr ==1
            ylabel({statenames{ss},'Cell'})
          end
            set(gca,'yticklabel',[])
    if ss==1
        title(regions{rr})
    end
    
    LogScale('x',10,'exp',true)
% subplot(2,1,2)
%     imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
%         ISIoccupancy.(state).hist(:,ISIStats.sorts.(state).ratebyclass)')
end
end

NiceSave('ISIOccupancy_MedOccupancysort',figfolder,[])
%% Mean-Normalized
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,4,(ss-1)*4+rr)
    s = imagesc(ISIoccupancy.(regions{rr}).logbins(1,:),[1 length(sorts.(regions{rr}).(statenames{ss}).ratebyclass)],...
        (ISIoccupancy.(regions{rr}).(state).normhist(:,sorts.(regions{rr}).(statenames{ss}).ratebyclass))');
    alpha(s,single(ISIoccupancy.(regions{rr}).(state).normhist(:,sorts.(regions{rr}).(statenames{ss}).ratebyclass)'~=0))

    hold on
    plot(0*log10(1./ISIstats.(regions{rr}).summstats.(state).meanrate(sorts.(regions{rr}).(statenames{ss}).ratebyclass)),...
        [1:length(sorts.(regions{rr}).(statenames{ss}).medISIbyclass)],'.')
    LogScale('x',10)
    caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('norm ISI (mean^-^1)')
          if rr ==1
            ylabel({statenames{ss},'Cell'})
          end
            set(gca,'yticklabel',[])
    if ss==1
        title(regions{rr})
    end
    
    LogScale('x',10,'exp',true)
% subplot(2,1,2)
%     imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
%         ISIoccupancy.(state).hist(:,ISIStats.sorts.(state).ratebyclass)')
end
end
NiceSave('ISIOccupancy_MeanNorm',figfolder,[])


%% MedianOccupancy-Normalized
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,4,(ss-1)*4+rr)
    s = imagesc(ISIoccupancy.(regions{rr}).logbins(1,:),[1 length(sorts.(regions{rr}).(statenames{ss}).ratebyclass)],...
        (ISIoccupancy.(regions{rr}).(state).mednormhist(:,sorts.(regions{rr}).(statenames{ss}).ratebyclass))');
    alpha(s,single(ISIoccupancy.(regions{rr}).(state).mednormhist(:,sorts.(regions{rr}).(statenames{ss}).ratebyclass)'~=0))

    hold on
    plot(0*log10(1./ISIstats.(regions{rr}).summstats.(state).meanrate(sorts.(regions{rr}).(statenames{ss}).ratebyclass)),...
        [1:length(sorts.(regions{rr}).(statenames{ss}).medISIbyclass)],'.')
    LogScale('x',10)
    caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('norm ISI (medOcc)')
          if rr ==1
            ylabel({statenames{ss},'Cell'})
          end
            set(gca,'yticklabel',[])
    if ss==1
        title(regions{rr})
    end
    
    LogScale('x',10,'exp',true)
% subplot(2,1,2)
%     imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
%         ISIoccupancy.(state).hist(:,ISIStats.sorts.(state).ratebyclass)')
end
end
NiceSave('ISIOccupancy_MedNorm',figfolder,[])
%% Median Occupancy Sort
histcolors = flipud(gray);
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
                ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(sorts.(regions{rr}).(statenames{ss}).medISIbyclass,:))
            hold on
            plot(log10((OccupancyStats.(regions{rr}).(statenames{ss}).median(sorts.(regions{rr}).(statenames{ss}).medISIbyclass))),...
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
            ylabel({statenames{ss},'Cell, sorted by MTO Rate'})
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
NiceSave('ISIDist_MedOccupancysort',figfolder,[])

%% MedOcc-normalized ISI
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,4,(ss-1)*4+rr)
colormap(gca,statecolormap{ss})
    s = imagesc(normISIhist.(regions{rr}).bins(1,:),[1 length(sorts.(regions{rr}).(statenames{ss}).medISIbyclass)],...
        (normISIhist.(regions{rr}).(statenames{ss}).mednorm(sorts.(regions{rr}).(statenames{ss}).medISIbyclass,:)));
    %alpha(s,single(ISIoccupancy.(regions{rr}).(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).medISIbyclass)'~=0))

    hold on
    plot(0*log10(OccupancyStats.(regions{rr}).(statenames{ss}).median(sorts.(regions{rr}).(statenames{ss}).medISIbyclass)),...
        [1:length(sorts.(regions{rr}).(statenames{ss}).medISIbyclass)],'k.','markersize',4)
    LogScale('x',10)
    %caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('norm ISI (medOcc)')
          if rr ==1
            ylabel({statenames{ss},'Cell, sorded by MTORate'})
          end
            set(gca,'yticklabel',[])
    if ss==1
        title(regions{rr})
    end
    caxis([0 0.1])
    LogScale('x',10,'exp',true)
% subplot(2,1,2)
%     imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
%         ISIoccupancy.(state).hist(:,ISIStats.sorts.(state).ratebyclass)')
end
end

NiceSave('normISIDist',figfolder,[])


%% Mean Distributions

%ISI (MedOcc Sextiles)
%NormISI
%Occupancy: NormISI

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

classcolors = {'k','r'};
subplot(3,4,8+rr)
hold on
for cc=1:2
plot(log10(ISIstats.(regions{rr}).summstats.WAKEstate.meanrate(CellClass.(regions{rr}).(classnames{cc}))),...
    log10(1./OccupancyStats.(regions{rr}).WAKEstate.median(CellClass.(regions{rr}).(classnames{cc}))),'.','color',classcolors{cc})

end
xlim([-4 2]);ylim([-4 2])
UnityLine
LogScale('xy',10,'exp',true)
xlabel('Mean Rate');ylabel('MTO Rate')
end

NiceSave('ISIdistMedOccPercile',figfolder,[])

%% GSRate by state
plotstates = {'WAKEstate','REMstate','WAKEstate'};
plotstates2 = {'NREMstate','NREMstate','REMstate'};



figure
%suptitle(regions{rr})
%Rate
for rr = 1:length(regions)
for ss=1:3
    subplot(4,4,(ss-1)*4+rr)
    hold on
    for cc=1:2
        plot(log10(1./OccupancyStats.(regions{rr}).(plotstates{ss}).median(CellClass.(regions{rr}).(classnames{cc}))),...
            log10(1./OccupancyStats.(regions{rr}).(plotstates2{ss}).median(CellClass.(regions{rr}).(classnames{cc}))),...
            'k.','markersize',2,'color',classcolors{cc})
        
    end

        %plot(log10([0.03 30]),log10([0.03 30]),'k')
        xlabel([plotstates{ss},' MTO Rate']);ylabel([plotstates2{ss},' MTO Rate'])
        %axis tight
        xlim([-2.5 2]);ylim([-2.5 2])
        LogScale('xy',10,'exp',true)
        UnityLine('linetype','-')
end
end
%%
%%
figure
for rr = 1:length(regions)
for ss = 1:3
    pcolor = makeColorMap([0.7 0.7 0.7],statecolors{ss},numperciles);

    subplot(5,4,(rr-1)+(ss-1)*4+1)
        hold on
        for cc = 1:length(percilenames)
            plot(meannormISIhist.bins,meannormISIhist.(regions{rr}).(statenames{ss}).(percilenames{cc}),...
                'linewidth',1,'color',pcolor(cc,:))
        end
        axis tight
        if ss==1
            title(regions{rr})
        end
        if rr == 1
            ylabel('p(ISI)');
        end
        %xlim([-3 1.9])
        LogScale('x',10,'exp',true)
        set(gca,'ytick',[])
            if ss==3
                xlabel('norm ISI (medOcc)')
            else
                set(gca,'xticklabels',[])
            end
            
subplot(5,4,(rr-1)+17)
 hold on
        for ss = 1:2
           
            plot(meannormISIhist.bins,meannormISIhist.(regions{rr}).(statenames{ss}).pE,...
                'linewidth',1,'color',statecolors{ss})
        end
        axis tight
        if ss==1
            title(regions{rr})
        end
        if rr == 1
            ylabel('P(ISI)');
        end
        xlim([-4 2.5])
        LogScale('x',10,'exp',true)
        set(gca,'ytick',[])
            %if ss==3
                xlabel('norm ISI (medOcc)')
            %else
             %   set(gca,'xticklabels',[])
            %end    
end
end
NiceSave('normISIdistMedOccPercile',figfolder,[])

%%
figure
for rr = 1:length(regions)
for ss = 1:3
    pcolor = makeColorMap([0.7 0.7 0.7],statecolors{ss},numperciles);

            
    subplot(5,4,(rr-1)+(ss-1)*4+1)
        hold on
        for cc = 1:length(percilenames)
            plot(meannormOcc.bins,meannormOcc.(regions{rr}).(statenames{ss}).(percilenames{cc}),...
                'linewidth',1,'color',pcolor(cc,:))
        end
        axis tight
        if ss==1
            title(regions{rr})
        end
        if rr == 1
            ylabel('p_t(ISI)');
        end
        xlim([-3 1.9])
        LogScale('x',10,'exp',true)
        set(gca,'ytick',[])
            if ss==3
                xlabel('norm ISI (medOcc)')
            else
                set(gca,'xticklabels',[])
            end
end

subplot(5,4,(rr-1)+17)
 hold on
        for ss = 1:2
           
            plot(meannormOcc.bins,meannormOcc.(regions{rr}).(statenames{ss}).pE,...
                'linewidth',1,'color',statecolors{ss})
        end
        axis tight
        if ss==1
            title(regions{rr})
        end
        if rr == 1
            ylabel('P_t(ISI)');
        end
        xlim([-4 2.5])
        LogScale('x',10,'exp',true)
        set(gca,'ytick',[])
            %if ss==3
                xlabel('norm ISI (medOcc)')
           % else
              %  set(gca,'xticklabels',[])
           % end

end

NiceSave('occISIdistMedOccPercile',figfolder,[])


%% Excells
figure
for rr =1:4
[~,excell(1)] = find(log10(OccupancyStats.(regions{rr}).WAKEstate.MTORatio)<0.3 & CellClass.(regions{rr}).pE,...
    1,'last')
[~,excell(2)] = find(log10(OccupancyStats.(regions{rr}).WAKEstate.MTORatio)>0.8 & CellClass.(regions{rr}).pE,...
    1,'last')
%excell = [108 112];


for ss = 1:3
    for ee = 1:2
        log10(OccupancyStats.(regions{rr}).(statenames{ss}).MTORatio(excell(ee)))
    subplot(6,4,(ss-1)*4+rr+(ee-1)*12)
    hold on
plot((ISIstats.(regions{rr}).ISIhist.logbins(1,:)),...
                ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(excell(ee),:),'linewidth',2,'color',statecolors{ss})
plot((ISIoccupancy.(regions{rr}).logbins(1,:)),...
                ISIoccupancy.(regions{rr}).(statenames{ss}).loghist(:,excell(ee)),':','linewidth',1,'color',statecolors{ss})
            
            box off
            plot(log10(OccupancyStats.(regions{rr}).(statenames{ss}).median(excell(ee))),0,'+')
            plot(log10(1./ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(excell(ee))),0,'+')
            xlim([-2.75 2.25])
            set(gca,'yticklabel',[])
            %ylabel('P(ISI), P_t(ISI)')
                        LogScale('x',10,'exp',true)
            
                xlabel('ISI (s)')
    end
end
end

NiceSave('ExampleISIs',figfolder,[])
