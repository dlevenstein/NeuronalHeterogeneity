reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISIModesbyPopActivityAnalysis']; 

[baseNames] = getDatasetBasenames();
regions = {'THAL','vCTX','fCTX','BLA','PIR','CA1'};
for rr = 1:6
    %Get baseNames
    if rr == 4 || rr == 5
        ISIPop_ALL = GetMatResults([figfolder,'_',regions{rr}],['ISIModesbyPopActivityAnalysis','_',regions{rr}]);
    else
        ISIPop_ALL = GetMatResults(figfolder,'ISIModesbyPopActivityAnalysis',...
            'baseNames',baseNames.(regions{rr}));
    end
    ISIPop_ALL = bz_CollapseStruct(ISIPop_ALL);
    
    PopCorr.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.PopCorr,'match','justcat',true);
    MUAConditionalISIDist.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.MUAConditionalISIDist,3,'justcat',true);
    %MUAConditionalISIDist_all.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.MUAConditionalISIDist_all,'match','justcat');
    MUAConditionalISIDist_gamma.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.MUAConditionalISIDist_gamma,'match','justcat');
    MUAConditionalISIDist_gamma.(regions{rr}).modes = bz_CollapseStruct(MUAConditionalISIDist_gamma.(regions{rr}).modes,3,'justcat',true);
    MUAConditionalISIDist_gamma.(regions{rr}).dist = bz_CollapseStruct(MUAConditionalISIDist_gamma.(regions{rr}).dist,3,'justcat',true);

    HiLowISIStats.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.ISIdists,'match','justcat',true);
end
%%
synchrate = {'rate','synch'};
statenames = {'WAKEstate','NREMstate','REMstate'};
celltypes = {'pE','pI'};
cellcolor = {'k','r'};
statecolors = {[0 0 0],[0 0 1],[1 0 0]};

cellthresh.pE = 20;
cellthresh.pI = 5;

hilow = {'HighPopRatestate','LowPopRatestate',};
%%
for sr = 1:2
for rr = 1:6
for ss = 1:3
    for tt = 1:length(celltypes)
        
        if ~isfield(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(synchrate{sr}),(celltypes{tt}))
            continue
        end
        for tt2 = 1:length(celltypes)
            usecells = PopCorr.(regions{rr}).CellClass.(celltypes{tt2}) & ...
                PopCorr.(regions{rr}).cellcount.(celltypes{tt})>cellthresh.(celltypes{tt});
            
            
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).Dist.pYX(:,:,...
                (usecells)),3);
            
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.SpikeRate = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).Dist.SpikeRate(:,:,...
                (usecells)),3);
            
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).Dist.Xbins(:,:,...
                (usecells)),3);
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Ybins = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).Dist.Ybins(:,:,...
                (usecells)),3);
            
            for hl = 1:2
            meanISIhist.(regions{rr}).(statenames{ss}).(celltypes{tt}).(hilow{hl}).(celltypes{tt2}).return =...
                mean(HiLowISIStats.(regions{rr}).(statenames{ss}).(celltypes{tt}).(hilow{hl}).return(:,:,usecells),3);
            meanISIhist.(regions{rr}).(statenames{ss}).(celltypes{tt}).(hilow{hl}).(celltypes{tt2}).logdist =...
                mean(HiLowISIStats.(regions{rr}).(statenames{ss}).(celltypes{tt}).(hilow{hl}).log(:,:,usecells),3);
            end
        end
        
    end
end
end

end

meanISIhist.logbins = HiLowISIStats.(regions{rr}).(statenames{ss}).(celltypes{tt}).logbins(1,:,1);



%%
for sr = 1:2
figure
for rr = 1:6
for tt = 1:2
    for ss = 1:3
        if ~isfield(PopCorr.(regions{rr}).(statenames{ss}).(synchrate{sr}),(celltypes{tt}))
            continue
        end
        
        usecells = PopCorr.(regions{rr}).cellcount.(celltypes{tt})>cellthresh.(celltypes{tt});
        
        subplot(6,6,(tt-1)*18+(ss-1)*6+rr)
        plot(PopCorr.(regions{rr}).(statenames{ss}).GSrate(usecells),...
            PopCorr.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt})(usecells),'k.','markersize',1)
        hold on
        plot(log10(PopCorr.(regions{rr}).(statenames{ss}).meanRate(PopCorr.(regions{rr}).CellClass.pI&usecells)),...
            PopCorr.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt})(PopCorr.(regions{rr}).CellClass.pI&usecells),...
            'r.','markersize',1)
        axis tight
        box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS/Mean Rate (Hz)');%ylabel([(celltypes{tt}),' Corr'])
        LogScale('x',10,'exp',true,'nohalf',true)
        if ss == 1
            title(regions{rr})
        end
        if rr ==1
            ylabel([(celltypes{tt}),' ',(synchrate{sr}),' Corr'])
        end
%         subplot(4,4,(tt-1)*4+ss+8)
%         plot(PopCorr.(regions{rr}).(statenames{ss}).GSCV,PopCorr.(regions{rr}).(statenames{ss}).(celltypes{tt}),'.')
%         hold on
%         axis tight
%         box off
%         plot(xlim(gca),[0 0],'k--')
%         xlabel('GS CV');ylabel([(celltypes{tt}),' Corr'])
%         %LogScale('x',10,'exp',false,'nohalf',true)
        
    end 
    
end
end
NiceSave(['MUACorrandGSRate_',(synchrate{sr})],figfolder,[])
end
%%
%% 
for sr = 1:2
for tt2 = 1:length(celltypes)
    for tt = 1:length(celltypes)
figure
for rr = 1:6
for ss = 1:3
     %Pop
        if ~isfield( MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}),(celltypes{tt})) | ...
            ~isfield(MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}),(celltypes{tt2})) | ...
            all(isnan(MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX))
            continue
        end
     %Ref

    subplot(4,6,(ss-1)*6+rr)
        imagesc(MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Ybins,...
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX')
        hold on
        plot(MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            -log10(MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.SpikeRate),...
            'r','LineWidth',2)
           
        
        LogScale('y',10,'exp',true,'nohalf',true)
        if rr == 1
            ylabel({(statenames{ss}),[(celltypes{tt2}),' ISI (s)']})
        elseif rr == 6
            set (gca,'yticklabels',[])
            bz_AddRightRateAxis
        else
            set (gca,'yticklabels',[])
        end
        
        if ss == 3
            xlabel([(celltypes{tt}),' ',(synchrate{sr})]);
        else
        set (gca,'xticklabels',[])
        end
        if ss == 1
            title(regions{rr})  
        end
        
    end
    end 

NiceSave(['ISIby',(synchrate{sr}),'_',(celltypes{tt}),'pop_',(celltypes{tt2}),'cells'],figfolder,[])
end
end
end


%%
figure
for rr = 1:6
for tt = 1:length(celltypes)
for ss = 1:2
subplot(4,6,rr+(tt-1)*12+6.*(ss-1))
hold on
keepcells = MUAConditionalISIDist_gamma.(regions{rr}).modes.(statenames{ss}).(celltypes{tt}).GScorr_p<=0.05;
%keepcells = true(size(keepcells))
for aa = 1:5
    %keepmodes = keepcells&mean(PSSConditionalGamma.modes.(states{ss}).ASweights(:,aa,:),1)>0.02;
    keepmodes = (MUAConditionalISIDist_gamma.(regions{rr}).modes.(statenames{ss}).(celltypes{tt}).AScorr_p(:,aa,:))<=0.05&...
        mean(MUAConditionalISIDist_gamma.(regions{rr}).modes.(statenames{ss}).(celltypes{tt}).ASweights(:,aa,:),1)>0.02;
    %keepmodes = true(size(keepmodes))
scatter(-MUAConditionalISIDist_gamma.(regions{rr}).modes.(statenames{ss}).(celltypes{tt}).ASlogrates(1,aa,keepmodes),...
    log10(MUAConditionalISIDist_gamma.(regions{rr}).modes.(statenames{ss}).(celltypes{tt}).ASCVs(1,aa,keepmodes)),...
    10*mean(MUAConditionalISIDist_gamma.(regions{rr}).modes.(statenames{ss}).(celltypes{tt}).ASweights(:,aa,keepmodes),1)+eps,...
    squeeze(MUAConditionalISIDist_gamma.(regions{rr}).modes.(statenames{ss}).(celltypes{tt}).AS_R(1,aa,keepmodes)),'filled')
end
scatter(-MUAConditionalISIDist_gamma.(regions{rr}).modes.(statenames{ss}).(celltypes{tt}).GSlogrates(1,1,keepcells),...
    log10(MUAConditionalISIDist_gamma.(regions{rr}).modes.(statenames{ss}).(celltypes{tt}).GSCVs(1,1,keepcells)),...
    5,...
    squeeze(MUAConditionalISIDist_gamma.(regions{rr}).modes.(statenames{ss}).(celltypes{tt}).GS_R(1,1,keepcells)),'linewidth',0.1)
%colorbar
axis tight
caxis([-0.2 0.2])
crameri('vik','pivot',0)
xlabel('Mean ISI');
if rr == 1
    ylabel({[(celltypes{tt}),' Pop'],(statenames{ss}),'CV'})
end

LogScale('x',10,'exp',true,'nohalf',true)
LogScale('y',10,'nohalf',true)
if ss == 1
    title(regions{rr})
end
end
end
end
NiceSave('ASModPopRate',figfolder,[])


%%
%% In/Out Field: REturn Maps ALL
histcolors = flipud(gray);
NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);
REMhistcolors = makeColorMap([1 1 1],[0.8 0 0]);
statecolormap = {histcolors,NREMhistcolors,REMhistcolors};

for tt = 1:2
    for tt2 = 1:2
figure
for rr =1:length(regions)
for ss = 1:2

subplot(6,6,rr+(ss-1)*18)
hold on
for hl = 1:2
    plot(meanISIhist.logbins,meanISIhist.(regions{rr}).(statenames{ss}).(celltypes{tt}).(hilow{hl}).(celltypes{tt2}).logdist,...
        'linewidth',1,'color',statecolormap{ss}(end./hl,:))
end
    if ss == 1
        title(regions{rr})
    end
    set(gca,'yticklabel',[])
    set(gca,'xticklabel',[])
    axis tight
    box off
    %LogScale('x',10,'exp',true,'nohalf',true)

%legend(cellISIStats.statenames{2:3},'location','southoutside')



for hl = 1:2
subplot(6,6,(hl-1)*6+(ss-1)*18+6+rr)
    imagesc(meanISIhist.logbins,meanISIhist.logbins,meanISIhist.(regions{rr}).(statenames{ss}).(celltypes{tt}).(hilow{hl}).(celltypes{tt2}).return)
    axis xy
    axis tight
    LogScale('xy',10,'exp',true,'nohalf',true)
    colormap(gca,statecolormap{ss})

    if hl == 1
        set(gca,'XTickLabels',[])
    end
end
%legend(cellISIStats.statenames{1:3})
end

%
end

NiceSave(['HiLowPopReturn_',(celltypes{tt}),'pop_',(celltypes{tt2}),'cells'],figfolder,[])
    end
end