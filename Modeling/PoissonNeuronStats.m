%Poisson Neuron Statistics
savefolder = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/Modeling/Figures';

meanrate = [0.1 0.3 1 3];
%meanrate = sort(FRhet.meanrate.NREMpacket);
numcells = length(meanrate);
numbins = 60;
ISIbins = linspace(log10(0.001),log10(200),numbins);
ISIhist = zeros(numcells,numbins);

for cc = 1:numcells
    neuronrate = meanrate(cc); %Hz
    simdur = 100000; %s

    %Generate random spikes at constant rate
    numspks = round(simdur*neuronrate);
    spiketimes = simdur*rand(numspks,1);
    spiketimes = sort(spiketimes);

    ISIs = diff(spiketimes);

    ISIhist(cc,:) = hist(log10(ISIs),ISIbins);
    ISIhist(cc,:) = ISIhist(cc,:)./numspks;
    
    ISIreturn(cc,:,:) = hist3(log10([ISIs(1:end-1) ISIs(2:end)]),{ISIbins,ISIbins});
    ISIreturn(cc,:,:) = ISIreturn(cc,:,:)./numspks;

end
%%
figure
subplot(2,2,2)
    plot(log10(ISIs(1:end-1)),log10(ISIs(2:end)),'.')
    hold on
    plot(log10(1./neuronrate),log10(1./neuronrate),'k+')
    LogScale('xy',10)
subplot(2,2,1)
    hist(log10(ISIs),40)
    
    %%
figure
    histcolors = makeColorMap([1 1 1],[0.8 0 0]);
    colormap(histcolors)

    subplot(2,3,1)
    imagesc((ISIbins),[1 numcells],ISIhist)
    hold on
    plot(log10(1./(meanrate)),[1:numcells],'k','LineWidth',2)
    LogScale('x',10)
    xlabel('ISI (s)')
    xlim(ISIbins([1 end]))
    %colorbar
  %  legend('1/Mean Firing Rate (s)','location','southeast')
    ylabel('Cell (Sorted by Mean FR)')
    caxis([0 0.1])
    title('ISI Distribution (Log Scale)')
    
    subplot(6,3,10)
    bar(ISIbins,mean(ISIhist,1))
    xlim(ISIbins([1 end]))
    LogScale('x',10)
    axis tight
NiceSave('ISIdistribution',savefolder,[]);
%%
%% Detailed Example Return Maps
numex = 4;
randomcells = randi(numcells,1,numex);
randomcells = sort(randomcells);

figure
colormap(histcolors)
for cc = 1:numex
    cellnum = randomcells(cc);

    subplot(4,2,2*cc-1)
        bar(ISIbins,(ISIhist(cellnum,:)),'facecolor',histcolors(end,:))
        hold on
        axis tight
        plot(log10(1./meanrate(cellnum)).*[1 1],get(gca,'ylim'),'k','LineWidth',2)
        LogScale('x',10)
        xlim(ISIbins([1 end]))
        xlabel('ISI (s)')
       ylabel(['FR: ',num2str(round(meanrate(cellnum),2)),'Hz'])
        set(gca,'yticklabel',[]);%set(gca,'xticklabel',[]);
        
    subplot(4,4,4*cc-1)    
        imagesc((ISIbins),(ISIbins),squeeze(ISIreturn(cellnum,:,:)))
        hold on
        plot(log10(1./meanrate(cellnum)),log10(1./meanrate(cellnum)),'k+')
        axis xy
        LogScale('xy',10)
        xlabel('ISI_n (s)');ylabel('ISI_n_+_1 (s)')
       % set(gca,'ytick',[]);set(gca,'xtick',[]);
        caxis([0 0.003])
        xlim(ISIbins([1 end]));ylim(ISIbins([1 end]))
        
end
    
NiceSave('exISImap',savefolder,[]);



%% Rate correlation

dt = 1; %s, time window
r = 1; %cell rate
n = 0:20;

pN = ((r.*dt).^n .* exp(-r.*dt))./factorial(n);

meanN = (1./sum(pN)) .* sum(n.*pN);
%varN = 
CV = mean


