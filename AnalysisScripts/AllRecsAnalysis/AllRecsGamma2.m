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
statenames = {'WAKEstate','NREMstate','REMstate'};

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
        
    %elseif rr==1
    %    GammaFit.(regions{rr}).(statenames{ss}).inregion = true(size(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSweights));
    else
        GammaFit.(regions{rr}).(statenames{ss}).inregion = true(size(GammaFit.(regions{rr}).(statenames{ss}).cellstats.meanrate));
    end
    end
    clear GammaFitAll
end


%%
weightthresh = 0.02; %perc of spikes
figure
for ss = 1:3
for rr = 1:length(regions)
    subplot(length(regions),3,(rr-1)*3+ss)
    histogram(log10(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASweights(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASweights>0)),...
        linspace(-3,0,50))
    hold on
    plot(log10(weightthresh).*[1 1],ylim(gca),'r--')
    
end
end
figure
for ss = 1:3
for rr = 1:length(regions)
    GammaFit.(regions{rr}).(statenames{ss}).singlecell.numAS = sum(GammaFit.(regions{rr}).(statenames{ss}).singlecell.ASweights(GammaFit.(regions{rr}).(statenames{ss}).inregion,:)>weightthresh,2);
    subplot(length(regions),3,(rr-1)*3+ss)
    histogram(GammaFit.(regions{rr}).(statenames{ss}).singlecell.numAS,0:5,'facecolor',regioncolors(rr,:))
    axis tight
    box off
    set(gca,'yticklabel',[])
    if ss == 1
        ylabel({(regions{rr}),'% Cells'})
    end
    if rr == 1
        title((statenames{ss}))
        set(gca,'xticklabel',[])
    elseif rr == length(regions)
        xlabel('# Activated States')
    else
        set(gca,'xticklabel',[])
    end
end
end
NiceSave(['NumModes'],figfolder,[])
%%
figure
for ss = 1:3
    for rr = 1:length(regions)
        try
subplot(6,length(regions),(ss-1)*length(regions)+rr)
    plot([GammaFit.(regions{rr}).WAKEstate.singlecell.GSlogrates(GammaFit.(regions{rr}).WAKEstate.cellstats.NW & GammaFit.(regions{rr}).WAKEstate.inregion)],...
        [GammaFit.(regions{rr}).NREMstate.singlecell.GSlogrates(GammaFit.(regions{rr}).NREMstate.cellstats.NW& GammaFit.(regions{rr}).NREMstate.inregion)],'.');
    hold on
    UnityLine
    xlabel('WAKE ');ylabel('NREM')
    title('GS Rate')

subplot(6,length(regions),(ss-1)*length(regions)+rr+3*length(regions))
    plot(1-[GammaFit.(regions{rr}).WAKEstate.singlecell.GSweights(GammaFit.(regions{rr}).WAKEstate.cellstats.NW& GammaFit.(regions{rr}).WAKEstate.inregion)],...
        1-[GammaFit.(regions{rr}).NREMstate.singlecell.GSweights(GammaFit.(regions{rr}).NREMstate.cellstats.NW & GammaFit.(regions{rr}).NREMstate.inregion )],'.');
    hold on
    UnityLine
    xlabel('WAKE');ylabel('NREM')
    title('AS Ratio')
        catch
            continue
        end
    end
end
%%

figure
subplot(3,3,1)
for ss = 1:3
    for rr = 1:length(regions)
        hold on
        try
    scatter([GammaFit.(regions{rr}).WAKEstate.singlecell.GSlogrates(GammaFit.(regions{rr}).WAKEstate.cellstats.NW & GammaFit.(regions{rr}).WAKEstate.inregion)],...
        [GammaFit.(regions{rr}).NREMstate.singlecell.GSlogrates(GammaFit.(regions{rr}).NREMstate.cellstats.NW& GammaFit.(regions{rr}).NREMstate.inregion)],...
        1,regioncolors(rr,:),'filled');
        catch
            continue
        end
    end
end
    axis tight
    UnityLine
    xlabel('WAKE ');ylabel('NREM')
    LogScale('xy',10)
    title('GS Rate')

        
subplot(3,3,2)
for ss = 1:3
    hold on
    for rr = 1:length(regions)
        try
    scatter(1-[GammaFit.(regions{rr}).WAKEstate.singlecell.GSweights(GammaFit.(regions{rr}).WAKEstate.cellstats.NW& GammaFit.(regions{rr}).WAKEstate.inregion)],...
        1-[GammaFit.(regions{rr}).NREMstate.singlecell.GSweights(GammaFit.(regions{rr}).NREMstate.cellstats.NW & GammaFit.(regions{rr}).NREMstate.inregion )],...
        1,regioncolors(rr,:),'filled');
            catch
            continue
        end
        end
end
    axis tight
    UnityLine
    xlabel('WAKE');ylabel('NREM')
    title('AS Ratio')
    
    
subplot(3,3,3)
for ss = 1:3
    for rr = 1:length(regions)
        hold on

    scatter(log10([GammaFit.(regions{rr}).WAKEstate.singlecell.GSCVs(GammaFit.(regions{rr}).WAKEstate.cellstats.NW & GammaFit.(regions{rr}).WAKEstate.inregion)]),...
        log10([GammaFit.(regions{rr}).NREMstate.singlecell.GSCVs(GammaFit.(regions{rr}).NREMstate.cellstats.NW& GammaFit.(regions{rr}).NREMstate.inregion)]),...
        1,regioncolors(rr,:),'filled');

    end
end
    %axis tight
    ylim([-0.5 0.7]);xlim([-0.5 0.7])
    plot([0 0],ylim(gca),'k--')
    plot(ylim(gca),[0 0],'k--')
    UnityLine
    
    xlabel('WAKE ');ylabel('NREM')
    LogScale('xy',10)
    title('GS CV')

    
subplot(3,3,7)
for ss = 1:3
    hold on
    for rr = 1:length(regions)
        try
    scatter([GammaFit.(regions{rr}).NREMstate.singlecell.GSlogrates(GammaFit.(regions{rr}).NREMstate.cellstats.NW& GammaFit.(regions{rr}).NREMstate.inregion)],...
        (1-[GammaFit.(regions{rr}).NREMstate.singlecell.GSweights(GammaFit.(regions{rr}).NREMstate.cellstats.NW & GammaFit.(regions{rr}).NREMstate.inregion)] )-...
        (1-[GammaFit.(regions{rr}).WAKEstate.singlecell.GSweights(GammaFit.(regions{rr}).WAKEstate.cellstats.NW& GammaFit.(regions{rr}).WAKEstate.inregion)]),...
        1,regioncolors(rr,:),'filled');
            catch
            continue
        end
        end
end
    axis tight
    plot(xlim(gca),[0 0],'k--')
    xlabel('GS Rate');ylabel('Change in AR')
    title('AS Ratio')
    
    
subplot(3,3,8)
for ss = 1:3
    hold on
    for rr = 1:length(regions)
    scatter([GammaFit.(regions{rr}).NREMstate.singlecell.GSlogrates(GammaFit.(regions{rr}).NREMstate.cellstats.NW& GammaFit.(regions{rr}).NREMstate.inregion)],...
        ([GammaFit.(regions{rr}).NREMstate.singlecell.GSlogrates(GammaFit.(regions{rr}).NREMstate.cellstats.NW & GammaFit.(regions{rr}).NREMstate.inregion)] )-...
        ([GammaFit.(regions{rr}).WAKEstate.singlecell.GSlogrates(GammaFit.(regions{rr}).WAKEstate.cellstats.NW& GammaFit.(regions{rr}).WAKEstate.inregion)]),...
        1,regioncolors(rr,:),'filled');

        end
        end
end
    axis tight
    plot(xlim(gca),[0 0],'k--')
    xlabel('GS Rate');ylabel('Change in GS')
    title('AS Ratio')
NiceSave(['GAARcArossStates'],figfolder,[])

%%

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
close all
%Nmodes = max(ISIfits.(regions{rr}).(statenames{ss}).Nmodes);
%Nmodes = 5;
figure
for ss = 1:3
    for rr = 1:length(regions)
%subplot(length(regions),3,(rr-1)*3+ss)
subplot(5,length(regions),(ss-1)*length(regions)+rr)
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

    %axis tight
    xlim([-1.5 2.5])
    ylim(yrange)
    LogScale('xy',10,'exp',true,'nohalf',true)
    if ss == 1
    title(regions{rr})
    end
    
    if ss ==3
        xlabel('Mode Rate (Hz)');
    else
        set(gca,'xticklabels',[])
    end
    if rr == 1
        ylabel({(statenames{ss}),' Cell Rate (Hz)'})
    else
        set(gca,'yticklabels',[])
    end

    end
end
NiceSave(['ISImodeandRate'],figfolder,[])

%%
close all
figure
for ss = 1:3
    for rr = 1:length(regions)
%subplot(length(regions),3,(rr-1)*3+ss)
subplot(5,length(regions),(ss-1)*length(regions)+rr)
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
    LogScale('xy',10,'exp',true,'nohalf',true)
    if ss == 1
    title(regions{rr})
    end
    
    if ss ==3
        xlabel('Mode Rate (Hz)');
    else
        set(gca,'xticklabels',[])
    end
    if rr == 1
        ylabel({(statenames{ss}),' Cell Rate (Hz)'})
    else
        set(gca,'yticklabels',[])
    end
    %axis tight
    %xlim([-2 1])
    ylim(yrange)
    end
end
NiceSave(['GSModeandRate_pE'],figfolder,[])





%%
close all
figure
for ss = 1:3
    for rr = 1:length(regions)
%subplot(length(regions),3,(rr-1)*3+ss)
subplot(5,length(regions),(ss-1)*length(regions)+rr)
    scatter(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSlogrates(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        1-GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSweights(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        2,log10(GammaFit.(regions{rr}).(statenames{ss}).cellstats.meanrate(GammaFit.(regions{rr}).(statenames{ss}).inregion)),...
        'filled')
    xlabel('GS rate');ylabel('Total AS weight')
    colorbar
    LogScale('c',10,'exp',true)
    end
end

subplot(3,3,7)
for ss = 1:2
    for rr = 1:length(regions)
        hold on
    scatter(GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSlogrates(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        1-GammaFit.(regions{rr}).(statenames{ss}).singlecell.GSweights(GammaFit.(regions{rr}).(statenames{ss}).inregion),...
        0.5,log10(GammaFit.(regions{rr}).(statenames{ss}).cellstats.meanrate(GammaFit.(regions{rr}).(statenames{ss}).inregion)),...
        'filled')
    end
end


ColorbarWithAxis([-1.5 1.5],'Mean FR')
xlabel('GS Rate')
LogScale('c',10,'nohalf',true)
ylabel('Total AS Weight')
LogScale('x',10,'exp',true)
NiceSave(['GSASandRate_pE'],figfolder,[])

%%
regnames = repmat(regions,2,1);
scolors = repmat([0 0 0;0 0 1],length(regions),1);

figure
%for ss = 1:2
subplot(3,2,1)
BoxAndScatterPlot({GammaFit.(regions{1}).(statenames{1}).singlecell.GSlogrates(GammaFit.(regions{1}).(statenames{1}).inregion),...
    GammaFit.(regions{1}).(statenames{2}).singlecell.GSlogrates(GammaFit.(regions{1}).(statenames{2}).inregion),...
    GammaFit.(regions{2}).(statenames{1}).singlecell.GSlogrates(GammaFit.(regions{2}).(statenames{1}).inregion),...
    GammaFit.(regions{2}).(statenames{2}).singlecell.GSlogrates(GammaFit.(regions{2}).(statenames{2}).inregion),...
    GammaFit.(regions{3}).(statenames{1}).singlecell.GSlogrates(GammaFit.(regions{3}).(statenames{1}).inregion),...
    GammaFit.(regions{3}).(statenames{2}).singlecell.GSlogrates(GammaFit.(regions{3}).(statenames{2}).inregion),...
    GammaFit.(regions{4}).(statenames{1}).singlecell.GSlogrates(GammaFit.(regions{4}).(statenames{1}).inregion),...
    GammaFit.(regions{4}).(statenames{2}).singlecell.GSlogrates(GammaFit.(regions{4}).(statenames{2}).inregion),...
    GammaFit.(regions{5}).(statenames{1}).singlecell.GSlogrates(GammaFit.(regions{5}).(statenames{1}).inregion),...
    GammaFit.(regions{5}).(statenames{2}).singlecell.GSlogrates(GammaFit.(regions{5}).(statenames{2}).inregion),...
    GammaFit.(regions{6}).(statenames{1}).singlecell.GSlogrates(GammaFit.(regions{6}).(statenames{1}).inregion),...
    GammaFit.(regions{6}).(statenames{2}).singlecell.GSlogrates(GammaFit.(regions{6}).(statenames{2}).inregion)},...
    'colors',scolors,...
    'labels',regnames(:))
    %plot(xlim(gca),[0.5 0.5],'k--')
%ylim([0 1])
box off
ylabel('GS rate')
LogScale('y',10)
%end

subplot(3,2,2)
BoxAndScatterPlot({log10(GammaFit.(regions{1}).(statenames{1}).sharedfit.GSCVs(GammaFit.(regions{1}).(statenames{1}).inregion)),...
    log10(GammaFit.(regions{1}).(statenames{2}).sharedfit.GSCVs(GammaFit.(regions{1}).(statenames{2}).inregion)),...
    log10(GammaFit.(regions{2}).(statenames{1}).sharedfit.GSCVs(GammaFit.(regions{2}).(statenames{1}).inregion)),...
    log10(GammaFit.(regions{2}).(statenames{2}).sharedfit.GSCVs(GammaFit.(regions{2}).(statenames{2}).inregion)),...
    log10(GammaFit.(regions{3}).(statenames{1}).sharedfit.GSCVs(GammaFit.(regions{3}).(statenames{1}).inregion)),...
    log10(GammaFit.(regions{3}).(statenames{2}).sharedfit.GSCVs(GammaFit.(regions{3}).(statenames{2}).inregion)),...
    log10(GammaFit.(regions{4}).(statenames{1}).sharedfit.GSCVs(GammaFit.(regions{4}).(statenames{1}).inregion)),...
    log10(GammaFit.(regions{4}).(statenames{2}).sharedfit.GSCVs(GammaFit.(regions{4}).(statenames{2}).inregion)),...
    log10(GammaFit.(regions{5}).(statenames{1}).sharedfit.GSCVs(GammaFit.(regions{5}).(statenames{1}).inregion)),...
    log10(GammaFit.(regions{5}).(statenames{2}).sharedfit.GSCVs(GammaFit.(regions{5}).(statenames{2}).inregion)),...
    log10(GammaFit.(regions{6}).(statenames{1}).sharedfit.GSCVs(GammaFit.(regions{6}).(statenames{1}).inregion)),...
    log10(GammaFit.(regions{6}).(statenames{2}).sharedfit.GSCVs(GammaFit.(regions{6}).(statenames{2}).inregion))},...
    'colors',scolors,...
    'labels',regnames(:))
    plot(xlim(gca),[0 0],'k--')
ylim([-0.5 0.7])
LogScale('y',10)
box off
ylabel('GS CV')

subplot(3,2,3)
BoxAndScatterPlot({GammaFit.(regions{1}).(statenames{1}).singlecell.GSweights(GammaFit.(regions{1}).(statenames{1}).inregion),...
    GammaFit.(regions{1}).(statenames{2}).singlecell.GSweights(GammaFit.(regions{1}).(statenames{2}).inregion),...
    GammaFit.(regions{2}).(statenames{1}).singlecell.GSweights(GammaFit.(regions{2}).(statenames{1}).inregion),...
    GammaFit.(regions{2}).(statenames{2}).singlecell.GSweights(GammaFit.(regions{2}).(statenames{2}).inregion),...
    GammaFit.(regions{3}).(statenames{1}).singlecell.GSweights(GammaFit.(regions{3}).(statenames{1}).inregion),...
    GammaFit.(regions{3}).(statenames{2}).singlecell.GSweights(GammaFit.(regions{3}).(statenames{2}).inregion),...
    GammaFit.(regions{4}).(statenames{1}).singlecell.GSweights(GammaFit.(regions{4}).(statenames{1}).inregion),...
    GammaFit.(regions{4}).(statenames{2}).singlecell.GSweights(GammaFit.(regions{4}).(statenames{2}).inregion),...
    GammaFit.(regions{5}).(statenames{1}).singlecell.GSweights(GammaFit.(regions{5}).(statenames{1}).inregion),...
    GammaFit.(regions{5}).(statenames{2}).singlecell.GSweights(GammaFit.(regions{5}).(statenames{2}).inregion),...
    GammaFit.(regions{6}).(statenames{1}).singlecell.GSweights(GammaFit.(regions{6}).(statenames{1}).inregion),...
    GammaFit.(regions{6}).(statenames{2}).singlecell.GSweights(GammaFit.(regions{6}).(statenames{2}).inregion)},...
    'colors',scolors,...
    'labels',regnames(:))
    plot(xlim(gca),[0.5 0.5],'k--')
ylim([0 1])
%box off
ylabel('p_G_S')
NiceSave('PGS',figfolder,[])

