%Poisson Neuron Statistics
savefolder = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/Modeling/Figures';


%Generate global input flutuations

%A = 1;
theta = 1;
sigma = 1;
duration = 2000;
dt = 0.0001;
save_dt = dt;
numsignals = 1;

[drive,t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);


%%
A = logspace(-2,0.5,25);
sig = logspace(-2,0.5,25);
[MU,SIG] = meshgrid(A,sig);

clear spiketimes ISIs CV2 meanCV2 meanrate
%%
for ss = 1:length(sig)
    ss
    for mm = 1:length(A)
        exprate = MU(ss,mm).*exp(drive.*SIG(ss,mm));

        [ s ] = PoissonRateSpikeBins(exprate,dt);
        spiketimes{ss,mm} = t(s)';
        ISIs{ss,mm} = diff(spiketimes{ss,mm});
        CV2{ss,mm} = 2.*abs(ISIs{ss,mm}(2:end)-ISIs{ss,mm}(1:end-1))./...
            (ISIs{ss,mm}(2:end)+ISIs{ss,mm}(1:end-1));
        meanCV2(ss,mm) = mean(CV2{ss,mm});
        meanrate(ss,mm) = 1./mean(ISIs{ss,mm});
    end
end


%%
[~,sortCV2] = sort(meanCV2(:));
[~,sortrate] = sort(meanrate(:));

%%
viewwin = [0 50];
figure
subplot(4,1,2)
    plot(t,drive,'k')
    xlim(viewwin)
    xlabel('t (fluctuation time scale)')
subplot(4,1,1)
    plotSpikeRaster(spiketimes(sortCV2),'PlotType','scatter');
    xlim(viewwin)
    axis xy
    
subplot(3,3,7)
    plot(log10(meanrate(:)),(meanCV2(:)),'k.')
    LogScale('x',10)
    axis tight
    xlabel('Mean Rate');ylabel('<CV2>')
subplot(3,3,8)
    imagesc(log10(A),log10(sig),log10(meanrate))
    axis xy
    xlabel('A');ylabel('Sigma')
    LogScale('xy',10)
    colorbar
    LogScale('c',10)
subplot(3,3,9)
    imagesc(log10(A),log10(sig),meanCV2)
    axis xy
    xlabel('A');ylabel('Sigma')
    LogScale('xy',10)
    colorbar
NiceSave('popfluct',savefolder,'poisspop')
    %%
figure
subplot(2,2,1)
    plot(log10(meanrate(:)),(meanCV2(:)),'.')
    LogScale('x',10)
subplot(2,2,3)
    imagesc(log10(A),log10(sig),log10(meanrate))
    axis xy
    xlabel('A');ylabel('Sigma')
    LogScale('xy',10)
    colorbar
    LogScale('c',10)
subplot(2,2,4)
    imagesc(log10(A),log10(sig),meanCV2)
    axis xy
    xlabel('A');ylabel('Sigma')
    LogScale('xy',10)
    colorbar
NiceSave('popfluctstats',savefolder,'poisspop')