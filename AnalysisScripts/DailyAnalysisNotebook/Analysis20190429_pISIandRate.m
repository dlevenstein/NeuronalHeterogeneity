function [ ] = Analysis20190429(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];

%%
datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC';
regions = {'fCTX','CA1'};


%%
for rr = 1:length(regions)
    ISIstats.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true);
    numcells.(regions{rr}) = length(CellClass.(regions{rr}).UID);
end

%%
statenames = fieldnames(ISIstats.(regions{1}).summstats);
statecolors = {'k','b','r'};
numstates = length(statenames);

%%
rr = 1;
ss = 1;

%%
for rr = 1:length(regions)
    for ss = 1:3
        [ratecorr.(regions{rr}).(statenames{ss}) corrsig.(regions{rr}).(statenames{ss})] = corr(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(CellClass.(regions{rr}).pE,:),...
            (ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).pE))','type','spearman','rows','complete');
    end
end


%% Sorts for plot
sorttypes = {'rate'};
%Make the cell-type specific sortings
for rr = 1:length(regions)
    for ss = 1:length(statenames)
        [~,sorts.(regions{rr}).(statenames{ss}).rate]=...
            sort(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate);

        %classnames = unique(CellClass.(regions{rr}).label);
        classnames = {'pE'};
        numclasses = length(classnames);
        for cl = 1:numclasses
            inclasscells.(regions{rr}){cl} = ...
                strcmp(classnames{cl},CellClass.(regions{rr}).label);

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

%Calculate mean ISI dists by state and cell type
meanISIhist.logbins = ISIstats.(regions{1}).ISIhist.logbins(1,:);
for rr = 1:length(regions)
    for ss = 1:length(statenames)
       for cc = 1:length(classnames)
           meanISIhist.(regions{rr}).(statenames{ss}).(classnames{cc}) = ...
               nanmean(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(CellClass.(regions{rr}).(classnames{cc}),:),1);
           meanISIhist.(regions{rr}).std.(statenames{ss}).(classnames{cc}) = ...
               nanstd(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(CellClass.(regions{rr}).(classnames{cc}),:),[],1);

       end
    end
end


%%
figure
plot(ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(CellClass.(regions{rr}).pE,25),...
    log10(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).pE)),'.')
%%
figure

    for rr = 1:length(regions)
        subplot(2,2,rr)
        hold on
        plot(ISIstats.(regions{rr}).ISIhist.logbins([1 end]),[0 0],'k')
        for ss = 1:3
            plot(ISIstats.(regions{rr}).ISIhist.logbins(1,:),ratecorr.(regions{rr}).(statenames{ss}),'color',statecolors{ss})
        end
        LogScale('x',10)
    end

    
    %%
    
    NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);
REMhistcolors = makeColorMap([1 1 1],[0.8 0 0]);
statecolormap = {histcolors,NREMhistcolors,REMhistcolors};




figure
for rr = 1:length(regions)
for ss = 1:3
    subplot(4,3,ss*3-2+(rr-1))
    colormap(gca,statecolormap{ss})

       % subplot(2,3,4)
            imagesc((ISIstats.(regions{rr}).ISIhist.logbins(1,:)),[1 sum(CellClass.(regions{rr}).pE)],...
                ISIstats.(regions{rr}).ISIhist.(statenames{ss}).log(sorts.(regions{rr}).(statenames{ss}).ratebyclass,:))
            hold on
            plot(log10(1./(ISIstats.(regions{rr}).summstats.(statenames{ss}).meanrate(sorts.(regions{rr}).(statenames{ss}).ratebyclass))),...
                [1:sum(CellClass.(regions{rr}).pE)],'k.','markersize',1)
            
            plot(meanISIhist.logbins,-meanISIhist.(regions{rr}).(statenames{ss}).pE*5000+...
                sum(inclasscells.(regions{rr}){1})+0.5,...
                'color',statecolors{ss},'linewidth',2)
            

            
            xlim(ISIstats.(regions{rr}).ISIhist.logbins([1 end]))
            LogScale('x',10)
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

for rr = 1:length(regions)
    subplot(4,3,rr+9)
    hold on
    plot(ISIstats.(regions{rr}).ISIhist.logbins([1 end]),[0 0],'k')
    for ss = 1:3
        plot(ISIstats.(regions{rr}).ISIhist.logbins(1,:),ratecorr.(regions{rr}).(statenames{ss}).*(corrsig.(regions{rr}).(statenames{ss})<0.05),...
            'color',statecolors{ss},'linewidth',1)
        %plot(ISIstats.(regions{rr}).ISIhist.logbins(1,:),ratecorr.(regions{rr}).(statenames{ss}),'color',statecolors{ss})
        
    end
    LogScale('x',10)
    xlim((ISIstats.(regions{rr}).ISIhist.logbins(1,[1 end])))
end

    NiceSave('pISIandCellRate',figfolder,baseName,'includeDate',true)

