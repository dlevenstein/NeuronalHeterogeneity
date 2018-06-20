%FluctuatingPoissonVariability

mu = 1;
theta = 1;
sigma = 1;
duration = 2000;
dt = 0.0005;
save_dt = dt;
numsignals = 1;

numbins = 20;
thetas = logspace(-1.5,2.5,numbins);
sigmas = logspace(-1.5,2,numbins);

%Threshold-Linear Rate Fluctuations
for tt = 1:length(thetas)
    tt
    for ss = 1:length(sigmas)
        theta = thetas(tt);
        sigma = sigmas(ss);

%Generate OU noise drive with given std, time scale 
[drive,t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
exprate = exp(drive);%exprate = zscore(exprate);
linrate = drive+mu;linrate(linrate<0)=0;

[ s ] = PoissonRateSpikeBins(linrate,dt);
spiketimes = t(s);
ISIs = diff(spiketimes);
CV2 = 2.*abs(ISIs(2:end)-ISIs(1:end-1))./(ISIs(2:end)+ISIs(1:end-1));
meanCV2_lin(tt,ss) = mean(CV2);
    end
end

%% Figure and Examples
viewwin = [0 30];
t_idx = t>=viewwin(1) & t<=viewwin(2);
examplethetas = [0.3 0.3 10 200];
examplesigmas = [0.3 30  10 3];

for ee = 1:length(examplethetas)
    theta = examplethetas(ee);
    sigma = examplesigmas(ee);
    %Generate OU noise drive with given std, time scale 
    [drive,t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
    drive = drive+mu;
    linrate_ex{ee} = drive;linrate_ex{ee}(linrate_ex{ee}<0)=0;
    [ s ] = PoissonRateSpikeBins(linrate_ex{ee},dt);
    spiketimes_ex{ee} = t(s);
    ISIs_ex{ee} = diff(spiketimes_ex{ee});
    CV2 = 2.*abs(ISIs_ex{ee}(2:end)-ISIs_ex{ee}(1:end-1))./(ISIs_ex{ee}(2:end)+ISIs_ex{ee}(1:end-1));
    meanCV2_ex{ee} = mean(CV2);
end

figure
subplot(3,3,1)
    imagesc(log10(thetas),log10(sigmas),meanCV2_lin')
    hold on
    plot(log10(examplethetas),log10(examplesigmas),'ro')
    axis xy
    xlabel('Time Scale');ylabel('Magnitude')
    colorbar
    LogScale('x',10)
    LogScale('y',10)
    caxis([0.9 1.2])
    
for ee = 1:length(examplethetas)
    subplot(4,6,3.*(ee+1)+[1 2]+6)
        plot(t(t_idx),(linrate_ex{ee}(t_idx)),'k')
        hold on
        plot(spiketimes_ex{ee},max((linrate_ex{ee}(t_idx))).*ones(size(spiketimes_ex{ee})),'k.')
        xlabel('t');ylabel('Rate')
        axis tight
        xlim(viewwin)
        box off
        
    subplot(4,6,ee.*3+12)
        plot(log10(ISIs_ex{ee}(1:end-1)),log10(ISIs_ex{ee}(2:end)),'.')
        title(['<CV2>: ',num2str(meanCV2_ex{ee})])
        
end



%% Lognormal Rate Fluctuations

dt = 0.00002;
save_dt = dt;
numbins = 20;
thetas = logspace(-1.5,2.5,numbins);
sigmas = logspace(-1,1,numbins);

for tt = 1:length(thetas)
    tt
    for ss = 1:length(sigmas)
        theta = thetas(tt);
        sigma = sqrt(log(1+sigmas(ss).^2));
        mu = log(1./sqrt(sigmas(ss).^2+1));
        
    %Generate OU noise drive with given std, time scale 
    [drive,t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
    exprate = exp(drive+mu);

    [ s ] = PoissonRateSpikeBins(exprate,dt);
    spiketimes = t(s);
    ISIs = diff(spiketimes);
    CV2 = 2.*abs(ISIs(2:end)-ISIs(1:end-1))./(ISIs(2:end)+ISIs(1:end-1));
    meanCV2_log(tt,ss) = mean(CV2);
    end
end


%% Figure and Examples
viewwin = [0 30];
t_idx = t>=viewwin(1) & t<=viewwin(2);
examplethetas = [0.3 0.3 10 200];
examplesigmas = [0.3 30  10 3];
clear exprate

for ee = 1:length(examplethetas)
    theta = examplethetas(ee);
    sigma = sqrt(log(1+examplesigmas(ee).^2));
    mu = log(1./sqrt(examplesigmas(ee).^2+1));
        
    %Generate OU noise drive with given std, time scale 
    [drive,t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
    exprate{ee} = exp(drive+mu);
    [ s ] = PoissonRateSpikeBins(exprate{ee},dt);
    spiketimes_ex{ee} = t(s);
    ISIs_ex{ee} = diff(spiketimes_ex{ee});
    CV2 = 2.*abs(ISIs_ex{ee}(2:end)-ISIs_ex{ee}(1:end-1))./(ISIs_ex{ee}(2:end)+ISIs_ex{ee}(1:end-1));
    meanCV2_ex{ee} = mean(CV2);
end
%%
figure
subplot(3,3,1)
    imagesc(log10(thetas),log10(sigmas),meanCV2_log')
    hold on
    plot(log10(examplethetas),log10(examplesigmas),'ro')
    axis xy
    xlabel('Time Scale');ylabel('Magnitude')
    colorbar
    LogScale('x',10)
    LogScale('y',10)
    caxis([0.9 1.2])
    
for ee = 1:length(examplethetas)
    subplot(4,6,3.*(ee+1)+[1 2]+6)
        plot(t(t_idx),(exprate{ee}(t_idx)),'k')
        hold on
        plot(spiketimes_ex{ee},max((exprate{ee}(t_idx))).*ones(size(spiketimes_ex{ee})),'k.')
        xlabel('t');ylabel('Rate')
        axis tight
        xlim(viewwin)
        box off
        
    subplot(4,6,ee.*3+12)
        plot(log10(ISIs_ex{ee}(1:end-1)),log10(ISIs_ex{ee}(2:end)),'.')
        title(['<CV2>: ',num2str(meanCV2_ex{ee})])
        LogScale('xy',10)
        
end


%%
    mean(exprate)
    std(exprate)
    
%%
viewwin = [0 40];
t_idx = t>=viewwin(1) & t<=viewwin(2);

figure
subplot(4,4,1:3)
plot(t(t_idx),drive(t_idx),'k')
xlabel('t');ylabel('Drive')
xlim(viewwin)

subplot(4,4,5:7)
plot(t(t_idx),(exprate(t_idx)),'k')
hold on
plot(spiketimes,max((exprate(t_idx))).*ones(size(spiketimes)),'k.')
xlabel('t');ylabel('Rate (Exp)')
axis tight
xlim(viewwin)

subplot(4,4,8)
plot(log10(ISIs(1:end-1)),log10(ISIs(2:end)),'.')
title(['<CV2>: ',num2str(meanCV2)])
%meanCV2

subplot(4,4,9:11)
plot(t(t_idx),(linrate(t_idx)),'k')
xlabel('t');ylabel('Rate (ThLin)')
xlim(viewwin)