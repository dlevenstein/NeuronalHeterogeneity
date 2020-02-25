reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SharedGammaModeFitAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
datasetPath.BLA = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/GG_BLA';
datasetPath.PIR = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/GG_BLA';
regions = {'THAL','vCTX','fCTX','BLA','PIR','CA1'};
rnames =  {''    ,''    ,''    ,'bla','pir',''   };
regioncolors = crameri('batlow',length(regions));
celltypes = {'pE','pI'};
cellcolor = {'k','r'};
statenames = {'NREMstate','WAKEstate','REMstate'};

for rr = 1:length(regions)
    disp(['Loading ',regions{rr}])
    %[ISIStats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    %CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);

    %[GammaFitAll,baseNames] = bz_LoadAnalysisResults(datasetPath.(regions{rr}),'SharedGammaModeFitAnalysis','dataset',true);
   % CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);
   % ISIStats.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true,'baseNames',baseNames);
   % ISIStats.(regions{rr}) = rmfield( ISIStats.(regions{rr}),'allspikes');
    
    %PopActivityAll = GetMatResults(figfolder,'SpikeStatsbyPopActivityAnalysis','baseNames',baseNames);
    %GammaFitAll = bz_CollapseStruct(GammaFitAll);
    %GammaFit.(regions{rr}) = bz_CollapseStruct(GammaFitAll.GammaFit,'match','justcat',true);
    [GammaFit.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'GammaFit','dataset',true,'catall',true);
    %Remove cells not in the proper region by removing their cell class!
    for ss = 1:3
    if ismember(rr,[4 5])
        
        GammaFit.(regions{rr}).(statenames{ss}).inregion = cellfun(@(X) strcmp(X,rnames{rr}),GammaFit.(regions{rr}).(statenames{ss}).cellstats.region);
%         CellClass.(regions{rr}).label(~inregion)={[]};
%         CellClass.(regions{rr}).pE(~inregion)=false;
%         CellClass.(regions{rr}).pI(~inregion)=false;
%         %clear ISIStats
        
    elseif rr==1
        GammaFit.(regions{rr}).(statenames{ss}).inregion = true(size(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSweights));
    else
        GammaFit.(regions{rr}).(statenames{ss}).inregion = true(size(GammaFit.(regions{rr}).(statenames{ss}).cellstats.meanrate));
    end
    end
    clear GammaFitAll
end

%%

GammaFit.BLA.NREMstate.singlecell.GSCVs

%%
figure
for ss = 1:3
for rr = 1:length(regions)
    subplot(length(regions),3,(rr-1)*3+ss)
    hist(log10(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASweights(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASweights>0)))
end
end
%%
weightthresh = 0.05; %perc of spikes
GScolor = [0.6 0.4 0];
close all
for cc = 1
figure
for ss = 1:3
    for rr = 1:length(regions)
        GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSweights(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSweights<weightthresh) = nan;
        GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSlogrates(isnan(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSweights))=nan;
        GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASweights(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASweights<weightthresh) = nan;
        GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASlogrates(isnan(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASweights))=nan;

%subplot(length(regions),3,(rr-1)*3+ss)
subplot(5,length(regions),(ss-1)*length(regions)+rr)
hold on
        
%     plot(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSlogrates(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
%         log10(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSCVs(GammaFit.(regions{rr}).(statenames{ss}).inregion)),...
%         '.','color',GScolor,'markersize',0.5)
%     plot(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASlogrates(GammaFit.(regions{rr}).(statenames{ss}).inregion,:),...
%         log10(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASCVs(GammaFit.(regions{rr}).(statenames{ss}).inregion,:)),...
%         '.','color',cellcolor{cc},'markersize',0.5)
    scatter(-GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSlogrates(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        log10(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSCVs(GammaFit.(regions{rr}).(statenames{ss}).inregion)),...
        GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSweights(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        GScolor,'filled')
   for aa = 1:5
    scatter(-GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASlogrates(GammaFit.(regions{rr}).(statenames{ss}).inregion,aa),...
        log10(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASCVs(GammaFit.(regions{rr}).(statenames{ss}).inregion,aa)),...
        GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASweights(GammaFit.(regions{rr}).(statenames{ss}).inregion,aa),...
        cellcolor{cc},'filled')
   end
    plot(xlim(gca),[0 0],'k--')
    ylim([-2 1])
    %ylim([0 4])
    xlim([-3 1.9])
     LogScale('y',10)
     LogScale('x',10,'exp',true)

if rr == 1
    ylabel({(statenames{ss}),'CV'})
else
    set(gca,'ytick',[])
end
if ss == 1
    title((regions{rr}))
    
elseif ss ==3
    xlabel('Mean ISI (s)');
end


    end 

end
NiceSave(['AllISImodes',(celltypes{cc})],figfolder,[])
end

%%

figure
for ss = 1:3
    for rr = 1:length(regions)
        subplot(length(regions),3,(rr-1)*3+ss)
            hist(ISIfits.(regions{rr}).(statenames{ss}).Nmodes(CellClass.(regions{rr}).(celltypes{cc})))
    end
end

%%
close all
%Nmodes = max(ISIfits.(regions{rr}).(statenames{ss}).Nmodes);
%Nmodes = 5;
figure
for ss = 1:3
    for rr = 1:length(regions)
subplot(length(regions),3,(rr-1)*3+ss)
hold on
for cc = 1
    scatter(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSlogrates(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        log10(GammaFit.(regions{rr}).(statenames{ss}).cellstats.meanrate(GammaFit.(regions{rr}).(statenames{ss}).inregion)),...
        GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSweights(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        GScolor,'filled')
   for aa = 1:5
    scatter(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASlogrates(GammaFit.(regions{rr}).(statenames{ss}).inregion,aa),...
        log10(GammaFit.(regions{rr}).(statenames{ss}).cellstats.meanrate(GammaFit.(regions{rr}).(statenames{ss}).inregion)),...
        GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASweights(GammaFit.(regions{rr}).(statenames{ss}).inregion,aa),...
        cellcolor{cc},'filled')
   end
end
    axis tight
    yrange = ylim(gca);
    UnityLine
    %ylim(log10([min(ISIStats.(regions{rr}).summstats.(statenames{ss}).meanrate) max(ISIStats.(regions{rr}).summstats.(statenames{ss}).meanrate)]))
    
    if rr == 1
    title((statenames{ss}))
    elseif rr ==length(regions)
        xlabel('Mode Rate (Hz)');
    end
    if ss == 1
        ylabel({(regions{rr}),' Cell Rate (Hz)'})
    end

    %axis tight
    xlim([-2 2.5])
    ylim(yrange)
    LogScale('xy',10)
    end
end
NiceSave(['ISImodeandRate'],figfolder,[])

%%
close all
figure
for ss = 1:3
    for rr = 1:length(regions)
subplot(length(regions),3,(rr-1)*3+ss)
hold on
for cc = 1
    scatter(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSlogrates(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        log10(GammaFit.(regions{rr}).(statenames{ss}).cellstats.meanrate(GammaFit.(regions{rr}).(statenames{ss}).inregion)),...
        0.5,GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSweights(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        'filled')
   for aa = 1:5
    scatter(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASlogrates(GammaFit.(regions{rr}).(statenames{ss}).inregion,aa),...
        log10(GammaFit.(regions{rr}).(statenames{ss}).cellstats.meanrate(GammaFit.(regions{rr}).(statenames{ss}).inregion)),...
        0.5,GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASweights(GammaFit.(regions{rr}).(statenames{ss}).inregion,aa),...
        'filled')
   end
end
    %colorbar
    axis tight
    yrange = ylim(gca);
    UnityLine
    %ylim(log10([min(ISIStats.(regions{rr}).summstats.(statenames{ss}).meanrate) max(ISIStats.(regions{rr}).summstats.(statenames{ss}).meanrate)]))
    if rr == 1
    title((statenames{ss}))
    elseif rr ==length(regions)
        xlabel('Mode Rate (Hz)');
    end
    if ss == 1
        ylabel({(regions{rr}),' Cell Rate (Hz)'})
    end
    %axis tight
    %xlim([-2 1])
    ylim(yrange)
    LogScale('xy',10)
    end
end
NiceSave(['GSModeandRate_pE'],figfolder,[])





%%
close all
figure
for ss = 1:3
    for rr = 1:length(regions)
subplot(length(regions),3,(rr-1)*3+ss)
    scatter(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSlogrates(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        1-GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSweights(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        2,log10(GammaFit.(regions{rr}).(statenames{ss}).cellstats.meanrate(GammaFit.(regions{rr}).(statenames{ss}).inregion)),...
        'filled')
    xlabel('GS rate');ylabel('Total AS weight')
    end
end
