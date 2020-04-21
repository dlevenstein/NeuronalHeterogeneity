
CVs = [0.05, 0.5, 1, 2];

%%
GSASparms.GSlogrates = 0;
GSASparms.GSCVs = 1;
GSASparms.GSweights =1;
GSASparms.ASlogrates = 0;
GSASparms.ASCVs = 1;
GSASparms.ASweights =0;

figure
[allISIdist] = GSASmodel(GSASparms,'sample',1,1);
spiketimes = cumsum(exp(allISIdist(1:20)));
plot(spiketimes,ones(size(spiketimes)),'k.')
%%
logtbins = linspace(-8,3,500);
figure
clear cc
for cc = 1:length(CVs)
    
subplot(6,4,cc)
    plot(logtbins,LogGamma(0,CVs(cc),1,logtbins'),'k','linewidth',2)
    box off
    axis tight
    set(gca,'ytick',[])

    xlabel('Log t')
    title(['CV = ',num2str(CVs(cc))])
    
    
    GSASparms.GSCVs = CVs(cc);
[allISIdist] = GSASmodel(GSASparms,'sample',1,1);
spiketimes = cumsum(exp(allISIdist(1:50)));
linrange = 15;

subplot(6,4,cc+8)
    isihist = hist(exp(allISIdist),linspace(0,linrange,100));
    bar(linspace(0,linrange,100),isihist,'facecolor','k','EdgeColor','none','barwidth',1)
    box off
    axis tight
    set(gca,'ytick',[])
    set(gca,'xtick',1)

subplot(6,4,cc+12)
    plot(spiketimes,ones(size(spiketimes)),'k.','markersize',0.01)
    xlim([0 linrange])
    box off
    set(gca,'ytick',[])
    set(gca,'xtick',[])

end

NiceSave('LogGamma',pwd,[])
