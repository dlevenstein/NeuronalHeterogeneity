reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/HeadDirectionTuningAnalysis'];


[HDALL,baseNames] = GetMatResults(figfolder,'HeadDirectionTuningAnalysis');
HDALL = bz_CollapseStruct(HDALL);

%%
ISIbyHD_align = bz_CollapseStruct(HDALL.ISIbyHD_align,3,'justcat',true);
ISIbyHD_alignGam = bz_CollapseStruct(HDALL.ISIbyHD_alignGam,3,'justcat',true);
MutInfo = bz_CollapseStruct(HDALL.MutInfo,'match','justcat',true);

%%
%MutInfo.Skaggs = MutInfo.SkaggsInf;
MIkinds = {'Skaggs','Rate','ISI'};

MIthresh.Rate = 0.02;
MIthresh.ISI = 0.02;
MIthresh.Skaggs = 1;

for kk = 1:3
    tunedcells.(MIkinds{kk}) = MutInfo.(MIkinds{kk})>MIthresh.(MIkinds{kk});
    MeanPlaceField.(MIkinds{kk}).pISI = nanmean(ISIbyHD_align.Dist.pYX(:,:,tunedcells.(MIkinds{kk})),3);
    MeanPlaceField.(MIkinds{kk}).Rate = nanmean(ISIbyHD_align.Dist.SpikeRate(:,:,tunedcells.(MIkinds{kk})),3);
    MeanPlaceField.(MIkinds{kk}).pGS = nanmean(ISIbyHD_alignGam.GammaModes.GSweights(:,:,tunedcells.(MIkinds{kk})),3);

    [~,sortMutInfo.(MIkinds{kk})] = sort(MutInfo.(MIkinds{kk}));
    [~,sortAR.(MIkinds{kk})] = sort(MutInfo.GSweight);
    sortAR.(MIkinds{kk}) = sortAR.(MIkinds{kk})(ismember(sortAR.(MIkinds{kk}),find(tunedcells.(MIkinds{kk}))));
end

numcells = length(MutInfo.GSrate);

%% Figure Information Metrics
figure
for kk = 1:3
subplot(3,2,1+(kk-1)*2)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),[0 numcells],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMutInfo.(MIkinds{kk}))))')
    hold on
    plot(xlim(gca),numcells-sum(tunedcells.(MIkinds{kk})).*[1 1],'r--','linewidth',1)
    ylabel(['Sort by I ',(MIkinds{kk})])
    bz_piTickLabel('x')
end

subplot(3,2,2)
scatter(log10(MutInfo.Skaggs),log10(MutInfo.Rate),3,MutInfo.GSweight,'filled')
axis tight
hold on
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'r--')
plot(xlim(gca),log10(MIthresh.Rate).*[1 1],'r--')
xlabel('I Skaggs');ylabel('MI Rate')
LogScale('xy',10)
colorbar

subplot(3,2,4)
scatter(log10(MutInfo.Skaggs),log10(MutInfo.ISI),3,MutInfo.GSweight,'filled')
axis tight
hold on
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'r--')
plot(xlim(gca),log10(MIthresh.ISI).*[1 1],'r--')
xlabel('I Skaggs');ylabel('MI ISI')
LogScale('xy',10)
colorbar

subplot(3,2,6)
scatter(log10(MutInfo.Rate),log10(MutInfo.ISI),3,MutInfo.GSweight,'filled')
axis tight
hold on
plot(xlim(gca),log10(MIthresh.ISI).*[1 1],'r--')
plot(log10(MIthresh.Rate).*[1 1],ylim(gca),'r--')
ylabel('MI ISI');xlabel('MI Rate') 
LogScale('xy',10)
colorbar

NiceSave('InfoNetrics',figfolder,[])
%% Figure: Information Metrics and GS/AS
figure
for kk = 1:3
subplot(3,2,1+(kk-1)*2)
scatter(log10(MutInfo.(MIkinds{kk})),MutInfo.GSrate,5,MutInfo.GSweight,'filled')
axis tight
hold on
plot(log10(MIthresh.(MIkinds{kk})).*[1 1],ylim(gca),'r--')
box off
xlabel(['I ',(MIkinds{kk})]);ylabel('GS Rate (HZ)')
LogScale('xy',10)
colorbar

subplot(3,2,2+(kk-1)*2)
scatter(log10(MutInfo.(MIkinds{kk})),MutInfo.GSweight,5,MutInfo.GSrate,'filled')
hold on
plot(log10(MIthresh.(MIkinds{kk})).*[1 1],ylim(gca),'r--')
box off
axis tight
if kk==3
    plot(xlim(gca).*[0 1]+log10(MIthresh.(MIkinds{kk})).*[1 0],0.5.*[1 1],'k--')
end
caxis([-0.5 1.25])
colorbar
LogScale('x',10)
LogScale('c',10)
xlabel(['I ',(MIkinds{kk})]);ylabel('GS Weight')
end

NiceSave('HDGSAS',figfolder,[])

%%

figure
for kk = 1:3
subplot(3,3,kk+3)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(MIkinds{kk}).pISI')
    hold on
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(MIkinds{kk}).pISI')
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(MIkinds{kk}).pISI')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1),-log10(MeanPlaceField.(MIkinds{kk}).Rate),'r')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,-log10(MeanPlaceField.(MIkinds{kk}).Rate),'r')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,-log10(MeanPlaceField.(MIkinds{kk}).Rate),'r')
    LogScale('y',10,'nohalf',true,'exp',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to HD Peak (m)')
    xlim([-1.5*pi 1.5.*pi])
    bz_piTickLabel('x')
    
    
subplot(3,3,6+(kk))
    plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1),MeanPlaceField.(MIkinds{kk}).pGS,'k')
    hold on
    plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1)+2*pi,MeanPlaceField.(MIkinds{kk}).pGS,'k')
    plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1)-2*pi,MeanPlaceField.(MIkinds{kk}).pGS,'k')
    box off
    xlim([-1.5*pi 1.5.*pi])
    ylim([0.1 0.5])
    bz_piTickLabel('x')
    ylabel('pGS')
    

subplot(3,3,kk)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),[0 length(sortAR.(MIkinds{kk}))],...
        squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortAR.(MIkinds{kk}))))')
    hold on
    %ylabel(['Sort by I ',(MIkinds{kk})])
    bz_piTickLabel('x')
    title([(MIkinds{kk}),'-tuned cells (',num2str(length(sortAR.(MIkinds{kk}))),')'])
    ylabel('Sorted by AR')

    
end  
    
NiceSave('HDCoding',figfolder,[])


%% Groups

hilowAR = 0.5;
groups = {'NonHiGS','TunedHiAR','NonLoGS','TunedLoAR'};
tunedcells.NonHiGS = MutInfo.GSweight>0.95 & MutInfo.GSrate>-0.5;
tunedcells.NonLoGS = MutInfo.GSweight>0.95 & MutInfo.GSrate<-0.5;
tunedcells.TunedHiAR = tunedcells.ISI & MutInfo.GSweight<hilowAR;
tunedcells.TunedLoAR = tunedcells.ISI & MutInfo.GSweight>hilowAR;
%tunedcells.TunedLoAR = ~tunedcells.ISI & MutInfo.GSweight<0.2;


for kk = 1:4
    %tunedcells.(groups{kk}) = MutInfo.(groups{kk})>MIthresh.(groups{kk});
    MeanPlaceField.(groups{kk}).pISI = nanmean(ISIbyHD_align.Dist.pYX(:,:,tunedcells.(groups{kk})),3);
    MeanPlaceField.(groups{kk}).Rate = nanmean(ISIbyHD_align.Dist.SpikeRate(:,:,tunedcells.(groups{kk})),3);
    MeanPlaceField.(groups{kk}).pGS = nanmean(ISIbyHD_alignGam.GammaModes.GSweights(:,:,tunedcells.(groups{kk})),3);

    %[~,sortMutInfo.(groups{kk})] = sort(MutInfo.(groups{kk}));
    %[~,sortAR.(groups{kk})] = sort(MutInfo.GSweight);
    %sortAR.(groups{kk}) = sortAR.(groups{kk})(ismember(sortAR.(groups{kk}),find(tunedcells.(groups{kk}))));
end
%Non-activating (pGS>0.95)
%Divide into high and low GS rate

%ISI-tuned
%Divide into high AR, low AR, non-activating?
%%
figure
ScatterWithLinFit(MutInfo.GSweight(tunedcells.ISI),MutInfo.peakwidth(tunedcells.ISI))
box off
xlabel('GS Weight');ylabel('Peak Width')
%%
figure
for kk = 1:4
    subplot(3,2,kk)
        imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(groups{kk}).pISI')
        hold on
        imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(groups{kk}).pISI')
        imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(groups{kk}).pISI')
        plot(ISIbyHD_align.Dist.Xbins(1,:,1),-log10(MeanPlaceField.(groups{kk}).Rate),'r')
        plot(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,-log10(MeanPlaceField.(groups{kk}).Rate),'r')
        plot(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,-log10(MeanPlaceField.(groups{kk}).Rate),'r')
        LogScale('y',10,'nohalf',true)
        ylabel('ISI (s)')
        bz_AddRightRateAxis
        xlabel('Position relative to HD Peak (m)')
        xlim([-1.5*pi 1.5.*pi])
        bz_piTickLabel('x')


    subplot(3,2,mod(kk+1,2)+5)
        plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1),MeanPlaceField.(groups{kk}).pGS,'k')
        hold on
        plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1)+2*pi,MeanPlaceField.(groups{kk}).pGS,'k')
        plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1)-2*pi,MeanPlaceField.(groups{kk}).pGS,'k')
        box off
        xlim([-1.5*pi 1.5.*pi])
        bz_piTickLabel('x')
    
end  
NiceSave('ARGroups',figfolder,[])