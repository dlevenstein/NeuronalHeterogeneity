figfolder = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/Modeling/Figures/ToyPoisson';

%% Two models: continuous rate vs modal rate
ASrate = 40;
GSrate = 0.5;
pAS_peak = 0.75;
%width = (2/3).*pi;

k = 4;
%peaknorm = exp(k.*cos(0))./(2.*pi.*besseli(0,k)); %Build in.
%tuningcurve_continuous = @(th) exp(k.*cos(th))./(2.*pi.*besseli(0,k));
tuningcurve_continuous = @(th) GSrate + (ASrate-GSrate).*...
    exp(k.*cos(th))./(exp(k.*cos(0)));


tuningcurve_probabilistic = @(th) pAS_peak.*exp(k.*cos(th))./(exp(k.*cos(0)));
ratecurve_probabilistic = @(th) GSrate + (ASrate-GSrate).*...
    (rand(size(th)) < tuningcurve_probabilistic(th));

%%
theta_X = linspace(-pi,pi,100);

figure
subplot(2,2,1)
plot(theta_X,(tuningcurve_continuous(theta_X)),'k')
%LogScale('y',10)
box off; axis tight

subplot(2,2,2)
plot(theta_X,(tuningcurve_probabilistic(theta_X)),'k')
%LogScale('y',10)
box off; axis tight
%% "Head Direction"

simtime = 5000; %s
dt = 0.0001;
timescale = 50.^-1;
sigma = 4*pi;

[Xvar.theta,Xvar.timestamps] = OUNoise(timescale,sigma,simtime,dt,dt,1);


%%
%dt = 0.0001;
[ s_continuous ] = PoissonRateSpikeBins(tuningcurve_continuous(Xvar.theta),dt);
[ s_probabilistic ] = PoissonRateSpikeBins(ratecurve_probabilistic(Xvar.theta),dt);

spiketimes_continuous = {Xvar.timestamps(s_continuous)};
spiketimes_probabilistic = {Xvar.timestamps(s_probabilistic)};


%%
Xvar.data = mod(Xvar.theta,2*pi);
numXbins = 100;
[ISIbyHD_continuous] = bz_ConditionalISI(spiketimes_probabilistic,Xvar,...
    'showfig',true,'GammaFit',false,'numXbins',numXbins,'numISIbins',100,...
    'normtype','none','Xwin',[0 2.*pi],'Xbinoverlap',3);


%numXbins = 100;
[ISIbyHD_probabilstic] = bz_ConditionalISI(spiketimes_continuous,Xvar,...
    'showfig',true,'GammaFit',false,'numXbins',numXbins,'numISIbins',100,...
    'normtype','none','Xwin',[0 2.*pi],'Xbinoverlap',3);


%%
figure
subplot(3,2,1)
plot(theta_X,(tuningcurve_continuous(theta_X)),'k')
%LogScale('y',10)
box off; axis tight

subplot(3,2,2)
plot(theta_X,(tuningcurve_probabilistic(theta_X)),'k')
%LogScale('y',10)
box off; axis tight

subplot(3,2,3)
imagesc(ISIbyHD_continuous.Dist.Xbins,ISIbyHD_continuous.Dist.Ybins,ISIbyHD_continuous.Dist.pYX')
%bounds = ylim(gca);
hold on
imagesc(ISIbyHD_continuous.Dist.Xbins-2*pi,ISIbyHD_continuous.Dist.Ybins,ISIbyHD_continuous.Dist.pYX')
xlim([-1.5*pi 1.5*pi])
LogScale('y',10)
ylabel('ISI (s)')
%bz_AddRightRateAxis
%colorbar

subplot(3,2,4)
imagesc(ISIbyHD_probabilstic.Dist.Xbins,ISIbyHD_probabilstic.Dist.Ybins,ISIbyHD_probabilstic.Dist.pYX')
%bounds = ylim(gca);
hold on
imagesc(ISIbyHD_probabilstic.Dist.Xbins-2*pi,ISIbyHD_probabilstic.Dist.Ybins,ISIbyHD_probabilstic.Dist.pYX')
xlim([-1.5*pi 1.5*pi])
LogScale('y',10)
ylabel('ISI (s)')
%bz_AddRightRateAxis

subplot(6,1,5)
plot(Xvar.timestamps,mod(Xvar.theta,2*pi),'k')
hold on
plot(spiketimes_continuous{1},ones(size(spiketimes_continuous{1})).*2*pi +0.5,'.');
plot(spiketimes_probabilistic{1},ones(size(spiketimes_probabilistic{1})).*2*pi +1,'.');
xlim([0 4])

subplot(6,1,6)
plot(Xvar.timestamps,tuningcurve_continuous(Xvar.theta),'k')
hold on
plot(Xvar.timestamps,ratecurve_probabilistic(Xvar.theta),'k--')
box off
xlim([0 4])

NiceSave('ProbContTuning',figfolder,'ToyPoisson')