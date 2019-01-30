%%
repoRoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity'; %desktop
basePath = '/mnt/NyuShare/Buzsakilabspace/Datasets/GrosmarkAD/Gatsby/Gatsby_08022013';
figfolder = [repoRoot,'/AnalysisScripts/AnalysisFigs/ReturnMapExample'];
baseName = bz_BasenameFromBasepath(basePath);
%%
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
SleepState = bz_LoadStates(basePath,'SleepState');
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID,...
    'basepath',basePath,'noPrompts',true);

%%

ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,SleepState.ints.NREMstate),...
    ISIStats.allspikes.times,'UniformOutput',false);
ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end);nan],ISIStats.allspikes.ISIs,...
    'UniformOutput',false);


%%
% thresh_n = cellfun(@(X,Y) bz_BimodalThresh(log10(X(Y)),'diptest',false,'startbins',25),...
%     ISIStats.allspikes.ISIs,ISIStats.allspikes.instate,'UniformOutput',false);
thresh_n = -1.5;
ISIStats.allspikes.noverthresh = cellfun(@(X) log10(X)>thresh_n, ISIStats.allspikes.ISIs,...
    'UniformOutput',false);
ISIStats.allspikes.np1overthresh = cellfun(@(X) log10(X)>thresh_n, ISIStats.allspikes.ISInp1,...
    'UniformOutput',false);

ISIStats.allspikes.zone1 = cellfun(@(X,Y) X&Y, ISIStats.allspikes.noverthresh,...
    ISIStats.allspikes.np1overthresh,'UniformOutput',false);
ISIStats.allspikes.zone2 = cellfun(@(X,Y) ~X&~Y, ISIStats.allspikes.noverthresh,...
    ISIStats.allspikes.np1overthresh,'UniformOutput',false);
ISIStats.allspikes.zone3 = cellfun(@(X,Y) X&~Y, ISIStats.allspikes.noverthresh,...
    ISIStats.allspikes.np1overthresh,'UniformOutput',false);
ISIStats.allspikes.zone4 = cellfun(@(X,Y) ~X&Y, ISIStats.allspikes.noverthresh,...
    ISIStats.allspikes.np1overthresh,'UniformOutput',false);
%
%figure
%[tes1,~,~,~,tes2] = bz_BimodalThresh(log10(ISIStats.allspikes.ISIs{excell}(ISIStats.allspikes.instate{excell})))

ISIStats.CV2hist.bins = linspace(0,2,20);
zones = {'zone1','zone2','zone3','zone4'};
ISIStats.CV2hist.ALL = cellfun(@(X,Y,Z) hist(X(Y),ISIStats.CV2hist.bins),...
    ISIStats.allspikes.CV2,ISIStats.allspikes.instate,...
    'UniformOutput',false);
for zz = 1:length(zones)
	zonespks = cellfun(@(Y,Z) (Y&Z),ISIStats.allspikes.instate,ISIStats.allspikes.(zones{zz}),...
        'UniformOutput',false);
    ISIStats.CV2hist.(zones{zz}) = cellfun(@(X,Y) hist(X(Y),ISIStats.CV2hist.bins),...
        ISIStats.allspikes.CV2,zonespks,'UniformOutput',false);
    exspks.(zones{zz}) = cellfun(@(X) randsample(find(X),1),zonespks);
end


%%
excell = 15;
%excell = 27;

histcolors = flipud(gray);


figure
subplot(3,4,1)
    colormap(gca,histcolors)
    imagesc(ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins,...
        ISIStats.ISIhist.NREMstate.return(:,:,excell))
    hold on
    %plot((thresh_n).*[1 1],ylim(gca),'r--')
    %plot(xlim(gca),thresh_n.*[1 1],'r--')
    xlabel('ISI_n');ylabel('ISI_n_+_1')
    xlim([-3 1]);    ylim([-3 1])

    axis xy
    LogScale('xy',10)
    
subplot(6,4,9)
    bar(ISIStats.ISIhist.logbins,ISIStats.ISIhist.NREMstate.log(excell,:),'facecolor','k')
    axis tight
    box off
    LogScale('x',10)
    xlim([-3 1]); 
    
subplot(6,4,17)
    hist(ISIStats.allspikes.ISIs{excell}(ISIStats.allspikes.instate{excell} & ...
        ISIStats.allspikes.ISIs{excell}<4),50)
    axis tight
    box off
    
subplot(6,4,21)
    hist(ISIStats.allspikes.ISIs{excell}(ISIStats.allspikes.instate{excell} & ...
        ISIStats.allspikes.ISIs{excell}<0.04),100)
    axis tight
    box off
    xlabel('ISI (s)')

NiceSave('ISIhist',figfolder,baseName)
%%
figure
    
subplot(3,3,1)
    scatter(log10(ISIStats.allspikes.ISIs{excell}(ISIStats.allspikes.instate{excell})),...
        log10(ISIStats.allspikes.ISInp1{excell}(ISIStats.allspikes.instate{excell})),0.1,...
        ISIStats.allspikes.CV2{excell}(ISIStats.allspikes.instate{excell}))
    %colorbar
    caxis([0 2])
    xlim(ISIStats.ISIhist.logbins([1 end]));
    ylim(ISIStats.ISIhist.logbins([1 end]));
    LogScale('xy',10)
    
    
subplot(6,3,3)
    bar(ISIStats.CV2hist.bins,ISIStats.CV2hist.ALL{excell},'k')
    hold on
%     for zz = 1:length(zones)
%         plot(ISIStats.CV2hist.bins,ISIStats.CV2hist.(zones{zz}){excell},'linewidth',2)
% 
%     end
    xlabel('CV2')
    ylabel('# Spikes')
    axis tight
    title(['<CV2>: ',num2str(ISIStats.summstats.NREMstate.meanCV2(excell))])
    box off
    
    
    zonecolor = {[0.5 0.5 0.5],'b','r'};
for zz = 1:3
    subplot(6,3,3*zz+7)
    bar(ISIStats.CV2hist.bins,ISIStats.CV2hist.(zones{zz}){excell},'facecolor',zonecolor{zz})
    if zz ==3
    xlabel('CV2')
    end
    ylabel('# Spikes')
    axis tight
    box off
end


%for zz = 1:4
    
zz = 2;
zonespks = cellfun(@(Y,Z) (Y&Z),ISIStats.allspikes.instate,ISIStats.allspikes.(zones{zz}),...
    'UniformOutput',false);
exspks.(zones{zz}) = cellfun(@(X) randsample(find(X),1),zonespks);

xwin = ISIStats.allspikes.times{excell}(exspks.(zones{zz})(excell))+[-0.1 0.5];
inwinspks = ISIStats.allspikes.times{excell}>(xwin(1)-10) & ISIStats.allspikes.times{excell}<(xwin(2)+10);

subplot(4,2,6)
bz_MultiLFPPlot(lfp,'timewin',xwin)
bz_ScaleBar('s')
xlabel('')

subplot(4,2,8)
plot(ISIStats.allspikes.times{excell}(inwinspks),ISIStats.allspikes.CV2{excell}(inwinspks),'.-')
hold on
plot([ISIStats.allspikes.times{excell}(inwinspks)*[1 1]]',...
    [ones(size(ISIStats.allspikes.times{excell}(inwinspks)))*[2.1 2.4]]','k')
plot(xwin,[2 2],'k')
plot(xwin,[2 2],'k--')

xlim(xwin)
ylim([0 2.5])
%bz_ScaleBar('s')
box off
ylabel('CV2')
set(gca,'xtick',[])
%end
NiceSave('CV2Hist',figfolder,baseName)
