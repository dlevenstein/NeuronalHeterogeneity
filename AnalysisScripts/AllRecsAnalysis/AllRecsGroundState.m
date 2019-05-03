reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/GroundStateAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC';

[GroundStateAll,baseNames] = GetMatResults(figfolder,'GroundStateAnalysis','select',true);
GroundStateAll = bz_CollapseStruct(GroundStateAll);
thisregion = 'CA1';
%%
for rr = 1:length(regions)
    ISIstats.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true);
    numcells.(regions{rr}) = length(CellClass.(regions{rr}).UID);
end
%%

ISIoccupancy = bz_CollapseStruct(GroundStateAll.ISIoccupancy,'match',...
    'justcat',true );


%% Sorts for plot
sorttypes = {'rate'};
tt =1
%Make the cell-type specific sortings
sorttypes = {'rate','ISICV','CV2'};
%Make the cell-type specific sortings
for rr = 2
    for ss = 1:length(statenames)
        [~,sorts.(regions{rr}).(statenames{ss}).rate]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate);
        [~,sorts.(regions{rr}).(statenames{ss}).ISICV]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).ISICV);
        [~,sorts.(regions{rr}).(statenames{ss}).CV2]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanCV2);

        classnames = unique(CellClass.(regions{rr}).label);
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

%%
figure
%colormap(cmap)
for ss = 1:3
        state = statenames{ss};

subplot(3,2,ss*2-1)
    s = imagesc(ISIoccupancy.logbins(1,:),[1 numcells.(regions{rr})],...
        (ISIoccupancy.(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).ratebyclass))');
    alpha(s,single(ISIoccupancy.(state).loghist(:,sorts.(regions{rr}).(statenames{ss}).ratebyclass)'~=0))

    hold on
    plot(log10(1./ISIstats.(regions{rr}).summstats.(state).meanrate(sorts.(regions{rr}).(statenames{ss}).ratebyclass)),...
        [1:numcells.(regions{rr})],'.')
    LogScale('x',10)
    ColorbarWithAxis([0 0.05],'P_t(log(ISI))')
    xlabel('ISI')
    ylabel(state)
% subplot(2,1,2)
%     imagesc(ISIoccupancy.bins,[1 spikes.numcells],...
%         ISIoccupancy.(state).hist(:,ISIStats.sorts.(state).ratebyclass)')
end

NiceSave('ISIOccupancy',figfolder,[])

%%
