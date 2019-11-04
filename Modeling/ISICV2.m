savefolder = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/Modeling/Figures';


rate = 1;
poissISIs = exprnd(rate,5000000,1);
gaussISIs = 0.25*randn(5000000,1)+rate;
gaussISIs(gaussISIs <= 0) = [];

CV2_poiss = 2.*abs(poissISIs(2:end)-poissISIs(1:end-1))./(poissISIs(2:end)+poissISIs(1:end-1));
CV2_gausss = 2.*abs(gaussISIs(2:end)-gaussISIs(1:end-1))./(gaussISIs(2:end)+gaussISIs(1:end-1));


gaussSpiketimes = cumsum(gaussISIs);
poissSpiketimes = cumsum(poissISIs);

%%
    CV2bins = linspace(0,2,100);
    %CV2hist{ee} = hist(CV2,CV2bins);
    
    ISIbins = linspace(-3.5,1.5,100);
    %ISIhist_ex{ee} = hist(log10(ISIs_ex{ee}),CV2bins);
    
    JointHist_poiss = hist3([log10(poissISIs(2:end)),CV2_poiss],{ISIbins,CV2bins});
	JointHist_gauss = hist3([log10(gaussISIs(2:end)),CV2_gausss],{ISIbins,CV2bins});
    
    
    ISIhist_poiss = hist(log10(poissISIs(2:end)),ISIbins);
    ISIhist_gauss = hist(log10(gaussISIs(2:end)),ISIbins);
%%
figure

subplot(6,3,1)
plot(poissSpiketimes,ones(size(poissSpiketimes)),'.')
xlim([0 10])

subplot(6,3,2)
plot(gaussSpiketimes,ones(size(gaussSpiketimes)),'.')
xlim([0 10])


subplot(4,3,4)
hist((poissISIs))
%LogScale('x',10)
xlim([0 6])
xlabel('ISI')

subplot(4,3,5)
hist((gaussISIs))
xlim([0 6])
%LogScale('x',10)
xlabel('ISI')

subplot(4,3,7)
hist(log10(poissISIs))
xlim([-2.5 2.5])
LogScale('x',10)
xlabel('ISI')

subplot(4,3,8)
hist(log10(gaussISIs))
xlim([-2.5 2.5])
LogScale('x',10)
xlabel('ISI')


subplot(4,3,10)
hist(CV2_poiss)
xlabel('CV2')

subplot(4,3,11)
hist(CV2_gausss)
xlabel('CV2')

subplot(4,3,6)
    imagesc(ISIbins,CV2bins,JointHist_poiss')
    axis xy
    hold on
    plot(ISIbins,bz_NormToRange(ISIhist_poiss,0.5),'w','linewidth',1)
    %axis tight
    %xlim([-3 1.5]);ylim([0 2])
    box off
    LogScale('x',10,'exp',true,'nohalf',true)
    xlabel('ISI');ylabel('CV2')
    

subplot(4,3,12)
    imagesc(ISIbins,CV2bins,JointHist_gauss')
    axis xy
    hold on
    plot(ISIbins,bz_NormToRange(ISIhist_gauss,0.5),'w','linewidth',1)
    %axis tight
    %xlim([-3 1.5]);ylim([0 2])
    box off
    LogScale('x',10,'exp',true,'nohalf',true)
    xlabel('ISI');ylabel('CV2')

NiceSave('ISICV2',savefolder,'poisspop')
