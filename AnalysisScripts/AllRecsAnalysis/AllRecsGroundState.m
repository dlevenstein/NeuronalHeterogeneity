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
sorttypes = {'rate','medISI'};
%tt =1
%Make the cell-type specific sortings
%sorttypes = {'rate','ISICV','CV2'};
%Make the cell-type specific sortings
for rr = 1:length(regions)
    for ss = 1:3
        [~,sorts.(regions{rr}).(statenames{ss}).rate]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate);
        [~,sorts.(regions{rr}).(statenames{ss}).ISICV]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).ISICV);
        [~,sorts.(regions{rr}).(statenames{ss}).CV2]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanCV2);
        [~,sorts.(regions{rr}).(statenames{ss}).medISI]=...
            sort(1./OccupancyStats.(regions{rr}).(statenames{ss}).median);
        
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

%% Rate Sort
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
%colormap(cmap)
for rr = 1:length(regions)
for ss = 1:3
        state = statenames{ss};

subplot(3,4,(ss-1)*4+rr)
colormap(gca,statecolormap{ss})
    s = imagesc(ISIstats.(regions{rr}).ISIhist.logbins(1,:),[1 length(sorts.(regions{rr}).(statenames{ss}).medISIbyclass)],...
        (ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(sorts.(regions{rr}).(statenames{ss}).medISIbyclass,:)));
    %alpha(s,single(ISIoccupancy.(regions{rr}).(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).medISIbyclass)'~=0))

    hold on
    plot(log10(OccupancyStats.(regions{rr}).(statenames{ss}).median(sorts.(regions{rr}).(statenames{ss}).medISIbyclass)),...
        [1:length(sorts.(regions{rr}).(statenames{ss}).medISIbyclass)],'k.','markersize',4)
    LogScale('x',10)
    %caxis([0 0.05])
    %ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('ISI')
          if rr ==1
            ylabel({statenames{ss},'Cell'})
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

NiceSave('ISIDist_MedOccupancysort',figfolder,[])

%%
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
            ylabel({statenames{ss},'Cell, sorded by medOcc'})
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

NiceSave('ISIDist_MedOccupancysort',figfolder,[])