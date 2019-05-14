%FluctuatingPoissonVariability
figfolder ='/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/Modeling/Figures/PoissonFluctuations';

%% Threshold linear
% mu = 1;
% theta = 1;
% sigma = 1;
% duration = 2000;
% dt = 0.0005;
% save_dt = dt;
% numsignals = 1;
% 
% numbins = 20;
% thetas_lin = logspace(-1.5,2.5,numbins);
% sigmas_lin = logspace(-1.5,2,numbins);
% 
% %Threshold-Linear Rate Fluctuations
% for tt = 1:length(thetas_lin)
%     tt
%     for ss = 1:length(sigmas_lin)
%         theta = thetas_lin(tt);
%         sigma = sigmas_lin(ss);
% 
% %Generate OU noise drive with given std, time scale 
% [drive,t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
% exprate = exp(drive);%exprate = zscore(exprate);
% linrate = drive+mu;linrate(linrate<0)=0;
% 
% [ s ] = PoissonRateSpikeBins(linrate,dt);
% spiketimes = t(s);
% ISIs = diff(spiketimes);
% CV2 = 2.*abs(ISIs(2:end)-ISIs(1:end-1))./(ISIs(2:end)+ISIs(1:end-1));
% meanCV2_lin(tt,ss) = mean(CV2);
%     end
% end
% 
% %% Figure and Examples
% viewwin = [0 30];
% t_idx = t>=viewwin(1) & t<=viewwin(2);
% examplethetas = [0.3 0.3 10 200];
% examplesigmas = [0.3 30  10 3];
% 
% for ee = 1:length(examplethetas)
%     theta = examplethetas(ee);
%     sigma = examplesigmas(ee);
%     %Generate OU noise drive with given std, time scale 
%     [drive,t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
%     drive = drive+mu;
%     linrate_ex{ee} = drive;linrate_ex{ee}(linrate_ex{ee}<0)=0;
%     [ s ] = PoissonRateSpikeBins(linrate_ex{ee},dt);
%     spiketimes_ex{ee} = t(s);
%     ISIs_ex{ee} = diff(spiketimes_ex{ee});
%     CV2 = 2.*abs(ISIs_ex{ee}(2:end)-ISIs_ex{ee}(1:end-1))./(ISIs_ex{ee}(2:end)+ISIs_ex{ee}(1:end-1));
%     meanCV2_ex{ee} = mean(CV2);
% end
% 
% figure
% subplot(3,3,1)
%     imagesc(log10(thetas_lin),log10(sigmas_lin),meanCV2_lin')
%     hold on
%     plot(log10(examplethetas),log10(examplesigmas),'ro')
%     axis xy
%     xlabel('Time Scale');ylabel('Magnitude')
%     colorbar
%     LogScale('x',10)
%     LogScale('y',10)
%     caxis([0.9 1.2])
%     
% for ee = 1:length(examplethetas)
%     subplot(4,6,3.*(ee+1)+[1 2]+6)
%         plot(t(t_idx),(linrate_ex{ee}(t_idx)),'k')
%         hold on
%         plot(spiketimes_ex{ee},max((linrate_ex{ee}(t_idx))).*ones(size(spiketimes_ex{ee})),'k.')
%         xlabel('t');ylabel('Rate')
%         axis tight
%         xlim(viewwin)
%         box off
%         
%     subplot(4,6,ee.*3+12)
%         plot(log10(ISIs_ex{ee}(1:end-1)),log10(ISIs_ex{ee}(2:end)),'.')
%         title(['<CV2>: ',num2str(meanCV2_ex{ee})])
%         
% end
% 


%% Lognormal Rate Fluctuations: CV2
%Time units: 1 estimated spike/time bin
dt = 0.0002;
duration = 5000; 
save_dt = dt;
numbins = 20;
numsignals = 1;
thetas_log = logspace(-1.5,2.5,numbins);
sigmas_log = logspace(-1,1,numbins);
trialduration = 20; 
numtrials = 100;

for tt = 1:length(thetas_log)
    tt
    %Create OU with time scale theta, mean 0, std 1.
    theta = thetas_log(tt);
    %sigma=1;
    [drive,t] = OUNoise(theta,1,duration,dt,save_dt,numsignals);
    
    %Shift OU so that exp(OU) has mean 1 and std sigma
    for ss = 1:length(sigmas_log)
        sigma = sqrt(log(1+sigmas_log(ss).^2));
        mu = log(1./sqrt(sigmas_log(ss).^2+1));
        
        %Generate OU noise drive with given std, time scale 
        exprate = exp((drive.*sigma)+mu);
        %MR = mean(exprate) %VALIDATION
        %SR = std(exprate)
        [ s ] = PoissonRateSpikeBins(exprate,dt);
        spiketimes = t(s);
        ISIs = diff(spiketimes);
        CV2 = 2.*abs(ISIs(2:end)-ISIs(1:end-1))./(ISIs(2:end)+ISIs(1:end-1));
        meanCV2_log(tt,ss) = mean(CV2);
        
        %Get Rates and FF from 100 random "trials"
        binnedrate = bz_SpktToSpkmat({spiketimes},'binsize',trialduration);
        numspikes{tt,ss} = binnedrate.data(randsample(length(binnedrate.data),100));
        FF(tt,ss) = std(numspikes{tt,ss}).^2./mean(numspikes{tt,ss});
        
    end
end


%% Lognormal Rate Fluctuations: FanoFactor
% numtrials = 50;
% trialduration = 10; %Expected # spikes in trial... at mean rate
% 
% for tt = 1:length(thetas_log)
%     tt
%     %Create OU with time scale theta, mean 0, std 1.
%     theta = thetas_log(tt);
%     [drive,t] = OUNoise(theta,1,trialduration,dt,save_dt,numtrials);
% 
%     for ss = 1:length(sigmas_log)
%         sigma = sqrt(log(1+sigmas_log(ss).^2));
%         mu = log(1./sqrt(sigmas_log(ss).^2+1));
%         
%         %Generate OU noise drive with given std, time scale 
%         exprate = exp((drive*sigma)+mu);
% 
%         [ s ] = PoissonRateSpikeBins(exprate,dt);
%         numspikes{tt,ss} = sum(s,1);
%         FF(tt,ss) = std(numspikes{tt,ss}).^2./mean(numspikes{tt,ss});
%    
%     end
% end


%% Figure and Examples


examplethetas = [1   0.3  5 200];
examplesigmas = [0.3  3   3  3];
clear exprate
clear exprate_trials
numsignals = 1;
for ee = 1:length(examplethetas)
    ee
    theta = examplethetas(ee);
    sigma = sqrt(log(1+examplesigmas(ee).^2));
    mu = log(1./sqrt(examplesigmas(ee).^2+1));
        
    %Generate OU noise drive with given std, time scale 
    [drive,t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
    exprate{ee} = exp(drive+mu);
    [ s ] = PoissonRateSpikeBins(exprate{ee},dt);
    spiketimes_ex{ee} = t(s);
    
    %Get ISIs and CV2
    ISIs_ex{ee} = diff(spiketimes_ex{ee});
    CV2 = 2.*abs(ISIs_ex{ee}(2:end)-ISIs_ex{ee}(1:end-1))./(ISIs_ex{ee}(2:end)+ISIs_ex{ee}(1:end-1));
    meanCV2_ex{ee} = mean(CV2);
    
    %Get Rates and FF
    binnedrate = bz_SpktToSpkmat(spiketimes_ex(ee),'binsize',trialduration);
    numspikes_ex{ee} = binnedrate.data;
    FF_ex{ee} = std(numspikes_ex{ee}).^2./mean(numspikes_ex{ee});

    
end

%%

viewwin = [10 25];
t_idx = t>=viewwin(1) & t<=viewwin(2);

heatcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
figure
subplot(5,4,1)
%colormap(gca,heatcolormap)
    imagesc(log10(thetas_log),log10(sigmas_log),meanCV2_log')
    hold on
    plot(log10(examplethetas),log10(examplesigmas),'ro')
    plot(log10(thetas_log([1 end])),[0 0],'k--')
    plot([0 0],log10(sigmas_log([1 end])),'k--')
    axis xy
    xlabel('Time Scale');ylabel('Magnitude')
    colorbar
    LogScale('x',10)
    LogScale('y',10)
    caxis([0.9 1.3])
    title('<CV2>')
    
%FANO FACTOR PLOT    
subplot(5,4,2)
%colormap(gca,heatcolormap)
    imagesc(log10(thetas_log),log10(sigmas_log),log10(FF)')
    hold on
    plot(log10(examplethetas),log10(examplesigmas),'ro')
    plot(log10(thetas_log([1 end])),[0 0],'k--')
    plot([0 0],log10(sigmas_log([1 end])),'k--')
    axis xy
    xlabel('Time Scale');ylabel('Magnitude')
    colorbar
    LogScale('x',10)
    LogScale('y',10)
    caxis([-0.5 2])
    LogScale('c',10)

    title(['FF (binsize: ',num2str(trialduration),')'])
    
    
subplots = [4 7 10 13];
subplots2 = [9 15 21 27];
subplots3 = [7 11 15 19];
for ee = 1:length(examplethetas)
    subplot(5,3,subplots(ee))
        plot(t(t_idx),log10(exprate{ee}(t_idx)),'color',[0.5 0.5 0.5],'linewidth',1)
        hold on
        plot([spiketimes_ex{ee} spiketimes_ex{ee}]',...
            [2.*ones(size(spiketimes_ex{ee})) 2.3.*ones(size(spiketimes_ex{ee}))]',...
            'k')
        %xlabel('t');
        ylabel('Rate')
        axis tight
        
        xlim(viewwin)
        ylim([-2.5 2.5])
        

        
        bz_ScaleBar
        
        LogScale('y',10)
        box off
        if ee==4
            xlabel('Time (AU)')
        end
        
    subplot(5,6,subplots2(ee))
        plot(log10(ISIs_ex{ee}(1:end-1)),log10(ISIs_ex{ee}(2:end)),'k.','markersize',5)
        %title(['<CV2>: ',num2str(round(meanCV2_ex{ee},2))])
        xlim([-4 2]);ylim([-4 2])
        LogScale('xy',10)
        %xlim(
        if ee==4
            xlabel('ISI_n');ylabel('ISI_n_+_1')
        end
        
    subplot(5,4,subplots3(ee))
        histogram(numspikes_ex{ee},'facecolor','k','BinEdges', [linspace(0,3*trialduration,20),Inf])
        %title(['FF: ',num2str(round(FF_ex{ee},1))])
        box off
        axis tight
        if ee==4
            xlabel('Spike Count')
        end
        
    subplot(5,4,subplots3(ee)+1)
        histogram(log10(ISIs_ex{ee}),'facecolor','k','BinEdges', linspace(-3.5,1.5,25))
        %title(['FF: ',num2str(round(FF_ex{ee},1))])
        box off
        axis tight
        if ee==4
            xlabel('ISI')
        end
        
        LogScale('x',10)
        
        
end

NiceSave('CompareFFCV2',figfolder,'Poiss')


%% Reduction in variability with sinusoidal input

dt = 0.001;
%save_dt = dt;
numbins = 25;
duration = 5000; 
freqs = logspace(-2,2,numbins);
kappas = logspace(-1,3,numbins);
t = [0:dt:duration]';
for tt = 1:length(freqs)
    tt
    for ss = 1:length(kappas)
        freq = freqs(tt);
        amp = kappas(ss);
        
        phase = mod(t,1./freq).*2*pi.*(freq);
        rate = circ_vmpdf(phase, 0, amp);
        rate = rate./mean(rate);
    %Generate OU noise drive with given std, time scale 
    %exprate = exp(drive);
    %exprate = (exprate./mean(exprate));
    %exprate = exprate

    [ s ] = PoissonRateSpikeBins(rate,dt);
    spiketimes = t(s);
    ISIs = diff(spiketimes);
    CV2 = 2.*abs(ISIs(2:end)-ISIs(1:end-1))./(ISIs(2:end)+ISIs(1:end-1));
    meanCV2_sin(tt,ss) = mean(CV2);
    end
end
   

%% Examples


examplefreqs = [1 1 0.1 10];
exampleamps = [1  100  100  100];
clear exprate
clear exprate_trials
numsignals = 1;
for ee = 1:length(examplefreqs)
    ee
    freq = examplefreqs(ee);
     amp = exampleamps(ee);
        
     
    phase = mod(t,1./freq).*2*pi.*(freq);
    rate_ex{ee} = circ_vmpdf(phase, 0, amp);
    rate_ex{ee} = rate_ex{ee}./mean(rate_ex{ee});
    [ s ] = PoissonRateSpikeBins(rate_ex{ee},dt);
    spiketimes_ex{ee} = t(s);
    
    %Get ISIs and CV2
    ISIs_ex{ee} = diff(spiketimes_ex{ee});
    CV2_ex{ee} = 2.*abs(ISIs_ex{ee}(2:end)-ISIs_ex{ee}(1:end-1))./(ISIs_ex{ee}(2:end)+ISIs_ex{ee}(1:end-1));
    meanCV2_ex{ee} = mean(CV2_ex{ee});
    
    
end

%%
viewwin = [10 25];
t_idx = t>=viewwin(1) & t<=viewwin(2);

heatcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
figure
subplot(5,4,1)
%colormap(gca,heatcolormap)
    imagesc(log10(freqs),log10(kappas),meanCV2_sin')
    hold on
    plot(log10(examplefreqs),log10(exampleamps),'ro')
    plot(log10(freqs([1 end])),[0 0],'k--')
    plot([0 0],log10(kappas([1 end])),'k--')
    axis xy
    xlabel('Frequency');ylabel('Amplitude')
    colorbar
    LogScale('x',10)
    LogScale('y',10)
    caxis([0.9 1.3])
    title('<CV2>')
    
    
subplots = [4 7 10 13];
subplots2 = [9 15 21 27];
subplots3 = [7 11 15 19];
for ee = 1:length(examplefreqs)
    subplot(5,3,subplots(ee))
        plot(t(t_idx),(rate_ex{ee}(t_idx)),'color',[0.5 0.5 0.5],'linewidth',1)
        hold on
       % plot(t(t_idx),cos(2*pi.*(examplefreqs(ee)).*t(t_idx)),'color',[0.5 0.5 0.5],'linewidth',1)
        plot([spiketimes_ex{ee}(spiketimes_ex{ee}>=viewwin(1) & spiketimes_ex{ee}<=viewwin(2))...
            spiketimes_ex{ee}(spiketimes_ex{ee}>=viewwin(1) & spiketimes_ex{ee}<=viewwin(2))]',...
            [2.*ones(size(spiketimes_ex{ee}(spiketimes_ex{ee}>=viewwin(1) & spiketimes_ex{ee}<=viewwin(2))))...
            2.3.*ones(size(spiketimes_ex{ee}(spiketimes_ex{ee}>=viewwin(1) & spiketimes_ex{ee}<=viewwin(2))))]',...
            'k')
        %xlabel('t');
        ylabel('Rate')
        axis tight
        
        xlim(viewwin)
       % ylim([-2.5 2.5])
        

        
        bz_ScaleBar('AU')
        
        %LogScale('y',10)
        box off
        if ee==4
            xlabel('Time (AU)')
        end
        
    subplot(5,6,subplots2(ee))
        plot(log10(ISIs_ex{ee}(1:end-1)),log10(ISIs_ex{ee}(2:end)),'k.','markersize',5)
        %title(['<CV2>: ',num2str(round(meanCV2_ex{ee},2))])
        xlim([-4 2]);ylim([-4 2])
        LogScale('xy',10)
        %xlim(
        if ee==4
            xlabel('ISI_n');ylabel('ISI_n_+_1')
        end
        
    subplot(5,4,subplots3(ee))
        histogram(CV2_ex{ee},'facecolor','k','BinEdges', [linspace(0,2,20),Inf])
        %title(['FF: ',num2str(round(FF_ex{ee},1))])
        box off
        axis tight
        if ee==4
            xlabel('CV2')
        end
        
    subplot(5,4,subplots3(ee)+1)
        histogram(log10(ISIs_ex{ee}),'facecolor','k','BinEdges', linspace(-3.5,1.5,25))
        %title(['FF: ',num2str(round(FF_ex{ee},1))])
        box off
        axis tight
        if ee==4
            xlabel('ISI')
        end
        
        LogScale('x',10)
        
        
end

NiceSave('OscCV2',figfolder,'Poiss')
