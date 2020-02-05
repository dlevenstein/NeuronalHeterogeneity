reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/GroundStateAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
datasetPath.BLA = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/GG_BLA';
datasetPath.PIR = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/GG_BLA';
regions = {'THAL','vCTX','fCTX','BLA','PIR','CA1'};
rnames =  {''    ,''    ,''    ,'bla','pir',''   };
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
    
    %Remove cells not in the proper region by removing their cell class!
    if ismember(rr,[4 5])
        inregion = cellfun(@(X) strcmp(X,rnames{rr}),ISIstats.(regions{rr}).cellinfo.regions);
        CellClass.(regions{rr}).label(~inregion)={[]};
        CellClass.(regions{rr}).pE(~inregion)=false;
        CellClass.(regions{rr}).pI(~inregion)=false;
    end
    clear GroundStateAll
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
numperciles = 5;

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
               squeeze(nanmean(normISIhist.(regions{rr}).(statenames{ss}).jointCV2(CellClass.(regions{rr}).(classnames{cc}),:,:),1));
       
       end
    end
end


%% Figure: Rate Sort Occupancy
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,6,(ss-1)*6+rr)
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

%% Figure: Rate Sort Occupancy Eonly
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,6,(ss-1)*6+rr)
    s = imagesc(ISIoccupancy.(regions{rr}).logbins(1,:),[1 length(sorts.(regions{rr}).(statenames{ss}).ratepE)],...
        (ISIoccupancy.(regions{rr}).(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).ratepE))');
    alpha(s,single(ISIoccupancy.(regions{rr}).(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).ratepE)'~=0))

    hold on
    plot(log10(1./ISIstats.(regions{rr}).summstats.(state).meanrate(sorts.(regions{rr}).(statenames{ss}).ratepE)),...
        [1:length(sorts.(regions{rr}).(statenames{ss}).ratepE)],'.')
    LogScale('x',10)
    caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('ISI (s)')
          if rr ==1
            ylabel({statenames{ss},'Cell, sorted by FR'})
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

NiceSave('ISIOccupancy_ratesort_E',figfolder,[])

%% Median Occupancy Sort
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,6,(ss-1)*6+rr)
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


%% Median Occupancy Sort - Eonly
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,6,(ss-1)*6+rr)
    s = imagesc(ISIoccupancy.(regions{rr}).logbins(1,:),[1 length(sorts.(regions{rr}).(statenames{ss}).medISIpE)],...
        (ISIoccupancy.(regions{rr}).(statenames{ss}).loghist(:,sorts.(regions{rr}).(statenames{ss}).medISIpE))');
    alpha(s,single(ISIoccupancy.(regions{rr}).(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).medISIpE)'~=0))

    hold on
    plot(log10(1./ISIstats.(regions{rr}).summstats.(state).meanrate(sorts.(regions{rr}).(statenames{ss}).medISIpE)),...
        [1:length(sorts.(regions{rr}).(statenames{ss}).medISIpE)],'.')
    LogScale('x',10)
    caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('ISI')
          if rr ==1
            ylabel({statenames{ss},'Cell, sorted by MTORate'})
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

NiceSave('ISIOccupancy_MedOccupancysort_E',figfolder,[])
%% Mean-Normalized
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,6,(ss-1)*6+rr)
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

subplot(3,6,(ss-1)*6+rr)
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
    subplot(3,6,ss*6-5+(rr-1))
    colormap(gca,statecolormap{ss})

       % subplot(2,3,4)
            imagesc((ISIstats.(regions{rr}).ISIhist.logbins(1,:)),[1 sorts.(regions{rr}).numclassycells],...
                ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(sorts.(regions{rr}).(statenames{ss}).medISIbyclass,:))
            hold on
            plot(log10((OccupancyStats.(regions{rr}).(statenames{ss}).median(sorts.(regions{rr}).(statenames{ss}).medISIbyclass))),...
                [1:sorts.(regions{rr}).numclassycells],'.','markersize',1,'color',[0.6 0.4 0])
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


%% Median Occupancy Sort - E only
histcolors = flipud(gray);
NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);
REMhistcolors = makeColorMap([1 1 1],[0.8 0 0]);
statecolormap = {histcolors,NREMhistcolors,REMhistcolors};



figure
for rr = 1:length(regions)
for ss = 1:3
    subplot(3,6,ss*6-5+(rr-1))
    colormap(gca,statecolormap{ss})

       % subplot(2,3,4)
            imagesc((ISIstats.(regions{rr}).ISIhist.logbins(1,:)),[1 length(sorts.(regions{rr}).(statenames{ss}).medISIpE)],...
                ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(sorts.(regions{rr}).(statenames{ss}).medISIpE,:))
            hold on
            plot(log10((OccupancyStats.(regions{rr}).(statenames{ss}).median(sorts.(regions{rr}).(statenames{ss}).medISIpE))),...
                [1:length(sorts.(regions{rr}).(statenames{ss}).medISIpE)],'.','markersize',1,'color',[0.6 0.4 0])

            
            plot(meanISIhist.logbins,-bz_NormToRange(meanISIhist.(regions{rr}).(statenames{ss}).pE,0.3)+length(sorts.(regions{rr}).(statenames{ss}).medISIpE),...
                'color',statecolors{ss},'linewidth',2)
            
            
            xlim(ISIstats.(regions{rr}).ISIhist.logbins([1 end]))
            xlim([-3 1.9])
            LogScale('x',10,'exp',true)
            if ss==3
                xlabel('ISI (s)')
            else
                set(gca,'xticklabels',[])
            end
            set(gca,'ytick',[])
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
NiceSave('ISIDist_MedOccupancysort_E',figfolder,[])
%% MedOcc-normalized ISI
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,6,(ss-1)*6+rr)
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

%% MedOcc-normalized ISI - Eonly
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,6,(ss-1)*6+rr)
colormap(gca,statecolormap{ss})
    s = imagesc(normISIhist.(regions{rr}).bins(1,:),[1 length(sorts.(regions{rr}).(statenames{ss}).medISIpE)],...
        (normISIhist.(regions{rr}).(statenames{ss}).mednorm(sorts.(regions{rr}).(statenames{ss}).medISIpE,:)));
    %alpha(s,single(ISIoccupancy.(regions{rr}).(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).medISIbyclass)'~=0))

    hold on
    plot(0*log10(OccupancyStats.(regions{rr}).(statenames{ss}).median(sorts.(regions{rr}).(statenames{ss}).medISIpE)),...
        [1:length(sorts.(regions{rr}).(statenames{ss}).medISIpE)],'k.','markersize',4)
    LogScale('x',10)
    %caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('norm ISI (MTO)')
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

NiceSave('normISIDist_E',figfolder,[])
%% MedOcc-normalized ISI
figure
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,6,(ss-1)*6+rr)
colormap(gca,statecolormap{ss})
    s = imagesc(normISIhist.(regions{rr}).bins(1,:),[1 length(sorts.(regions{rr}).(statenames{ss}).MTOratpE)],...
        (normISIhist.(regions{rr}).(statenames{ss}).mednorm(sorts.(regions{rr}).(statenames{ss}).MTOratpE,:)));
    %alpha(s,single(ISIoccupancy.(regions{rr}).(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).medISIbyclass)'~=0))

    hold on
    plot(0*log10(OccupancyStats.(regions{rr}).(statenames{ss}).median(sorts.(regions{rr}).(statenames{ss}).MTOratpE)),...
        [1:length(sorts.(regions{rr}).(statenames{ss}).MTOratpE)],'k.','markersize',4)
    LogScale('x',10)
    %caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('norm ISI (medOcc)')
          if rr ==1
            ylabel({statenames{ss},'Cell, sorded by MTORatio'})
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

NiceSave('normISIDist_MTOratio',figfolder,[])
%% Mean Distributions

%ISI (MedOcc Sextiles)
%NormISI
%Occupancy: NormISI

%% FR %Ile FIgure
figure
for rr = 1:length(regions)
for ss = 1:3
    pcolor = makeColorMap([0.7 0.7 0.7],statecolors{ss},numperciles);

    subplot(5,6,(rr-1)+(ss-1)*6+1)
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
subplot(3,6,12+rr)
hold on
for cc=1:2
plot(log10(1./OccupancyStats.(regions{rr}).WAKEstate.median(CellClass.(regions{rr}).(classnames{cc}))),...
    log10(ISIstats.(regions{rr}).summstats.WAKEstate.meanrate(CellClass.(regions{rr}).(classnames{cc}))),...
    '.','color',classcolors{cc})

end
xlim([-4 2]);ylim([-4 2])
UnityLine
LogScale('xy',10,'exp',true)
xlabel('MTO Rate');ylabel('Mean Rate')
end

NiceSave('ISIdistMedOccPercile',figfolder,[])

%%
figure
for rr = 1:length(regions)
for ss = 1:3
classcolors = {'k','r'};
subplot(3,6,(ss-1)*6+rr)
hold on
for cc=1
plot(log10(1./OccupancyStats.(regions{rr}).(statenames{ss}).median(CellClass.(regions{rr}).(classnames{cc}))),...
    log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).(classnames{cc}))),...
    '.','color',classcolors{cc})

end
xlim([-3 2]);ylim([-3 2])
UnityLine
LogScale('xy',10,'exp',true)
if ss == 3
xlabel('MTO Rate (Hz)');
end
if rr == 1
    ylabel({statenames{ss},'Mean Rate (Hz)'})
end
end
end
NiceSave('MTOandMeanRate',figfolder,[])

%%
for ss = 1:3
   figure
for rr = 1:length(regions)
   for cc = 1:length(percilenames)
       
       [~,idx] = min(abs(10.^meanISIhist.logbins - 1./meanpercmedISI.(regions{rr}).(statenames{ss})(cc)));
       
        subplot(length(percilenames),6,(cc-1)*6+rr)    
        colormap(gca,statecolormap{ss})

            imagesc(meanISIhist.logbins,meanISIhist.logbins,...
                meanreturnhist.(regions{rr}).(statenames{ss}).(percilenames{cc})')
            hold on
            plot(log10(1./meanpercmedISI.(regions{rr}).(statenames{ss})(cc)),log10(1./meanpercmedISI.(regions{rr}).(statenames{ss})(cc)),'+','color',[0.6 0.4 0])
            plot(meanISIhist.logbins,bz_NormToRange(meanISIhist.(regions{rr}).(statenames{ss}).(percilenames{cc}),0.3),...
                'color',statecolors{ss},'linewidth',0.5)
            UnityLine
            axis xy
            %xlim([-3 1.9])
            set(gca,'ytick',[]);set(gca,'xticklabel',[]);
            if cc==1
                title(regions{rr})
            end
            if cc==length(percilenames) 
                xlabel('ISI (s)')
        
                LogScale('x',10,'exp',true)
            end
            if rr==1 
                ylabel('ISI_n_+_1 (s)')
        
               % LogScale('y',10,'exp',true)
            end
            %ylim([0 2]);
            %xlim([-2.5 1.7])
            xlim(ISIstats.(regions{rr}).ISIhist.logbins([1 end]))
            ylim(ISIstats.(regions{rr}).ISIhist.logbins([1 end]))
            
           %caxis([0 1.75*meanreturnhist.(regions{rr}).(statenames{ss}).(percilenames{cc})(idx,idx)])

    end
end
NiceSave(['PercentilesReturn_',(statenames{ss})],figfolder,[])

end
%% Box plot comparing regions

regnames = repmat(regions,2,1);

figure
%for ss = 1:2
subplot(3,1,1)
BoxAndScatterPlot({log10(1./OccupancyStats.(regions{1}).(statenames{1}).median(CellClass.(regions{1}).(classnames{1}))),...
    log10(1./OccupancyStats.(regions{1}).(statenames{2}).median(CellClass.(regions{1}).(classnames{1}))),...
    log10(1./OccupancyStats.(regions{2}).(statenames{1}).median(CellClass.(regions{2}).(classnames{1}))),...
    log10(1./OccupancyStats.(regions{2}).(statenames{2}).median(CellClass.(regions{2}).(classnames{1}))),...
    log10(1./OccupancyStats.(regions{3}).(statenames{1}).median(CellClass.(regions{3}).(classnames{1}))),...
    log10(1./OccupancyStats.(regions{3}).(statenames{2}).median(CellClass.(regions{3}).(classnames{1}))),...
    log10(1./OccupancyStats.(regions{4}).(statenames{1}).median(CellClass.(regions{4}).(classnames{1}))),...
    log10(1./OccupancyStats.(regions{4}).(statenames{2}).median(CellClass.(regions{4}).(classnames{1}))),...
    log10(1./OccupancyStats.(regions{5}).(statenames{1}).median(CellClass.(regions{5}).(classnames{1}))),...
    log10(1./OccupancyStats.(regions{5}).(statenames{2}).median(CellClass.(regions{5}).(classnames{1}))),...
    log10(1./OccupancyStats.(regions{6}).(statenames{1}).median(CellClass.(regions{6}).(classnames{1}))),...
    log10(1./OccupancyStats.(regions{6}).(statenames{2}).median(CellClass.(regions{6}).(classnames{1})))},...
    'colors',repmat([0 0 0;0 0 1],6,1),...
    'labels',regnames(:))
ylim([-3 1.9])
box off
ylabel('GS Rate')
LogScale('y',10)


subplot(3,1,2)
BoxAndScatterPlot({log10(OccupancyStats.(regions{1}).(statenames{1}).MTORatio(CellClass.(regions{1}).(classnames{1}))),...
    log10(OccupancyStats.(regions{1}).(statenames{2}).MTORatio(CellClass.(regions{1}).(classnames{1}))),...
    log10(OccupancyStats.(regions{2}).(statenames{1}).MTORatio(CellClass.(regions{2}).(classnames{1}))),...
    log10(OccupancyStats.(regions{2}).(statenames{2}).MTORatio(CellClass.(regions{2}).(classnames{1}))),...
    log10(OccupancyStats.(regions{3}).(statenames{1}).MTORatio(CellClass.(regions{3}).(classnames{1}))),...
    log10(OccupancyStats.(regions{3}).(statenames{2}).MTORatio(CellClass.(regions{3}).(classnames{1}))),...
    log10(OccupancyStats.(regions{4}).(statenames{1}).MTORatio(CellClass.(regions{4}).(classnames{1}))),...
    log10(OccupancyStats.(regions{4}).(statenames{2}).MTORatio(CellClass.(regions{4}).(classnames{1}))),...
    log10(OccupancyStats.(regions{5}).(statenames{1}).MTORatio(CellClass.(regions{5}).(classnames{1}))),...
    log10(OccupancyStats.(regions{5}).(statenames{2}).MTORatio(CellClass.(regions{5}).(classnames{1}))),...
    log10(OccupancyStats.(regions{6}).(statenames{1}).MTORatio(CellClass.(regions{6}).(classnames{1}))),...
    log10(OccupancyStats.(regions{6}).(statenames{2}).MTORatio(CellClass.(regions{6}).(classnames{1})))},...
    'colors',repmat([0 0 0;0 0 1],6,1),...
    'labels',regnames(:))
box off
ylabel('Activation Ratio')
ylim([0 2])
%end
NiceSave('MTORegins',figfolder,[])
%%
figure
for rr = 1:length(regions)
for ss =1:3
subplot(6,6,rr+(ss-1)*6)
scatter(log10(1./OccupancyStats.(regions{rr}).(statenames{ss}).median(CellClass.(regions{rr}).pE)),...
    log10(OccupancyStats.(regions{rr}).(statenames{ss}).MTORatio(CellClass.(regions{rr}).pE)),0.5,...
     log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).pE)))
 axis tight
 
 ylim([0 2])
 xlim([-2.75 1.75])
 LogScale('xy',10,'exp',true,'nohalf',true)
 
  caxis([-1.5 2])
  
 colorbar
 LogScale('c',10,'exp',true,'nohalf',true)
 if ss==1
     title(regions{rr})
 end
  if rr ==1
    ylabel('Activation Ratio') 
 end

subplot(6,6,rr+(ss-1)*6+18)
scatter(log10(1./OccupancyStats.(regions{rr}).(statenames{ss}).median(CellClass.(regions{rr}).pE)),...
    log10(OccupancyStats.(regions{rr}).(statenames{ss}).MTORatio(CellClass.(regions{rr}).pE)),0.5,...
     (ISIstats.(regions{rr}).summstats.(statenames{ss}).meanCV2(CellClass.(regions{rr}).pE))) 
  axis tight
 ylim([0 2])
  xlim([-2.75 1.75])
 caxis([0.75 1.25])
 colorbar
 crameri('berlin','pivot',1)
 LogScale('xy',10,'exp',true,'nohalf',true)
 if ss ==3
     xlabel('Ground State Rate (MTO)')
 end
 if rr ==1
    ylabel('Activation Ratio') 
 end
end
end
NiceSave('RateCV2byGSAS',figfolder,[])
%% GSRate by state
plotstates = {'WAKEstate','REMstate','WAKEstate'};
plotstates2 = {'NREMstate','NREMstate','REMstate'};



figure
%suptitle(regions{rr})
%Rate
for rr = 1:length(regions)
for ss=1:3
    subplot(4,6,(ss-1)*6+rr)
    hold on
    for cc=1:2
        plot(log10(1./OccupancyStats.(regions{rr}).(plotstates{ss}).median(CellClass.(regions{rr}).(classnames{cc}))),...
            log10(1./OccupancyStats.(regions{rr}).(plotstates2{ss}).median(CellClass.(regions{rr}).(classnames{cc}))),...
            'k.','markersize',2,'color',classcolors{cc})
        
    end

        %plot(log10([0.03 30]),log10([0.03 30]),'k')
        xlabel([plotstates{ss},' MTO Rate']);ylabel([plotstates2{ss},' MTO Rate'])
        %axis tight
        xlim([-3 2]);ylim([-3 2])
        LogScale('xy',10,'exp',true)
        UnityLine('linetype','-')
                if ss==1
            title(regions{rr})
        end
end
end

NiceSave('MTOAcrossStates',figfolder,[])

%% GSRate by state
plotstates = {'WAKEstate','REMstate','WAKEstate'};
plotstates2 = {'NREMstate','NREMstate','REMstate'};



figure
%suptitle(regions{rr})
%Rate
for rr = 1:length(regions)
for ss=1:3
    subplot(4,6,(ss-1)*6+rr)
    hold on
    for cc=1
        plot(log10(1./OccupancyStats.(regions{rr}).(plotstates{ss}).median(CellClass.(regions{rr}).(classnames{cc}))),...
            log10(1./OccupancyStats.(regions{rr}).(plotstates2{ss}).median(CellClass.(regions{rr}).(classnames{cc}))),...
            'k.','markersize',2,'color',classcolors{cc})
        
    end

        %plot(log10([0.03 30]),log10([0.03 30]),'k')
        xlabel([plotstates{ss},' MTO Rate']);ylabel([plotstates2{ss},' MTO Rate'])
        %axis tight
        xlim([-3 2]);ylim([-3 2])
        LogScale('xy',10,'exp',true)
        UnityLine('linetype','-')
                if ss==1
            title(regions{rr})
        end
end
end

NiceSave('MTOAcrossStates_Eonly',figfolder,[])
%% Simulate: Poisson MTO-norm ISI/occupancy dists
dt = 0.005;
[ s ] = PoissonRateSpikeBins(1,dt,10000000);
timestamps = [1:length(s)]*dt;
stimes = timestamps(s);
ISIs = diff(stimes);
stimes(1) = [];

ISIrate = interp1(stimes,ISIs,timestamps,'next');

logbins = linspace(-3,1,101);
poissISIoccupancy = hist(log10(ISIrate),logbins);

MTO = nanmedian(ISIrate);

MTOnormoccupancy = hist(log10(ISIrate./MTO),logbins);
MTOnormoccupancy = MTOnormoccupancy./sum(MTOnormoccupancy);
MTOnormISI = hist(log10(ISIs./MTO),meannormISIhist.bins);
[~,normbin] = min(abs(meannormISIhist.bins-0));
MTOnormISI = MTOnormISI./sum(MTOnormISI);
%%
figure
subplot(2,2,1)
hist(log10(ISIs))
subplot(2,2,2)
plot(logbins,MTOnormoccupancy,'k')
subplot(2,2,3)
plot(logbins,MTOnormISI,'k')
%%
figure
for rr = 1:length(regions)
for ss = 1:3
    pcolor = makeColorMap([0.7 0.7 0.7],statecolors{ss},numperciles);

    subplot(5,6,(rr-1)+(ss-1)*6+1)
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
                xlabel('norm ISI (MTO^-^1)')
            else
                set(gca,'xticklabels',[])
            end
            
subplot(5,6,(rr-1)+25)
 hold on
        for ss = 1:2
 
            plot(meannormISIhist.bins,meannormISIhist.(regions{rr}).(statenames{ss}).pE,...
                'linewidth',1,'color',statecolors{ss})
        end
        plot(meannormISIhist.bins,MTOnormISI,'k--')
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
                xlabel('norm ISI (MTO^-^1)')
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

            
    subplot(5,6,(rr-1)+(ss-1)*6+1)
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
                xlabel('norm ISI (MTO^-^1)')
            else
                set(gca,'xticklabels',[])
            end
end

subplot(5,6,(rr-1)+25)
 hold on
        for ss = 1:2
           
            plot(meannormOcc.bins,meannormOcc.(regions{rr}).(statenames{ss}).pE,...
                'linewidth',1,'color',statecolors{ss})
        end
        plot(logbins,MTOnormoccupancy,'k')
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
                xlabel('norm ISI (MTO^-^1)')
           % else
              %  set(gca,'xticklabels',[])
           % end

end

NiceSave('occISIdistMedOccPercile',figfolder,[])


%% Excells
figure
for rr =1:length(regions)
[~,excell(1)] = find(log10(OccupancyStats.(regions{rr}).WAKEstate.MTORatio)<0.3 & CellClass.(regions{rr}).pE,...
    1,'last')
[~,excell(2)] = find(log10(OccupancyStats.(regions{rr}).WAKEstate.MTORatio)>0.8 & CellClass.(regions{rr}).pE,...
    1,'last')
%excell = [108 112];


for ss = 1:3
    for ee = 1:2
        log10(OccupancyStats.(regions{rr}).(statenames{ss}).MTORatio(excell(ee)))
    subplot(6,6,(ss-1)*6+rr+(ee-1)*18)
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

%%
for ss = 1:3
   figure
for rr = 1:length(regions)
   for cc = 1:length(percilenames)
       
       [~,idx] = min(abs(10.^meanISIhist.logbins - 1./meanpercmedISI.(regions{rr}).(statenames{ss})(cc)));
       
        subplot(length(percilenames),6,(cc-1)*6+rr)    
        %colormap(gca,statecolormap{ss})

            imagesc(meanISIhist.logbins,ISIstats.(regions{rr}).CV2hist.bins(1,:),...
                meanJointhist.(regions{rr}).(statenames{ss}).(percilenames{cc}).log')
            hold on
            plot(log10(1./meanpercmedISI.(regions{rr}).(statenames{ss})(cc)),ISIstats.(regions{rr}).CV2hist.bins(1,1),'r+')
            plot(meanISIhist.logbins,bz_NormToRange(meanISIhist.(regions{rr}).(statenames{ss}).(percilenames{cc}),0.3),...
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
            
            if ss==2
                caxis([0 2.7*meanJointhist.(regions{rr}).(statenames{ss}).(percilenames{cc}).log(idx,1)])
            else
                caxis([0 3.1*meanJointhist.(regions{rr}).(statenames{ss}).(percilenames{cc}).log(idx,1)])
                if cc == length(percilenames) & rr~=4
                caxis([0 max([meanJointhist.(regions{rr}).(statenames{ss}).(percilenames{cc}).log(10:end,1);0])])
                end
            end
            %caxis([0 max([meanJointhist.(regions{rr}).(statenames{ss}).(percilenames{cc}).log(10:end,1);0])])
            %crameri tokyo

    end
end

NiceSave(['MTOPercentiles_',(statenames{ss})],figfolder,[])
end

%%
figure
for rr = 1:length(regions)


for cc = 1:2
	for ss = 1:3
        
        [~,idx] = min(abs(meannormISIhist.bins - 0));
        
        
        subplot(6,6,rr+(ss-1)*6+(cc-1)*18)    
        %colormap(gca,statecolormap{ss})

            imagesc(meannormISIhist.bins,ISIstats.(regions{rr}).CV2hist.bins(1,:),...
                meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).norm')
            hold on
            plot(meannormISIhist.bins,meannormISIhist.(regions{rr}).(statenames{ss}).(classnames{cc})*25,...
                'color',statecolors{ss},'linewidth',1)
            
            axis xy
            ylim([0 2]);
            set(gca,'ytick',[]);%set(gca,'xtick',[]);

            if ss==3 
            
                xlabel('Norm ISI (MTO^-^1)')
            end
                %set(gca,'xtick',[-2:1]);
                %LogScale('x',10)
            
            if cc==1 
                ylabel('CV2')
                set(gca,'ytick',[0 1 2]);
            end
            
            
        if ss==1
            title(regions{rr})
        end
        if rr == 1
            ylabel('CV2');
        end            %xlim([-3 1])
            LogScale('x',10,'exp',true)
%             switch cc
%                 case 1
%                     switch rr
%                         case 1
%                             caxis([0.5e-4 1.2e-3])
%                         case 2
%                             if ss == 2; caxis([0.5e-4 0.8e-3])
%                             else; caxis([0.5e-4 0.6e-3])
%                             end
%                     end
%                 case 2
%                     caxis([0.5e-4 1.8e-3])
%             end
        %colorbar
        
             caxis([0 max([meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).norm(:,2);0])])
            %caxis([0 max([2*meanJointhist.(regions{rr}).(statenames{ss}).(classnames{cc}).norm(idx,2);0])])
    end
end
end
NiceSave('JointDist_MTONorm',figfolder,[])

%%
