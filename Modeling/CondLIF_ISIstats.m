

%TO DO: set DT, calculate probability of number of spikes
%% FOr Figure
% R = 0.3;
% [ S ] = PoissonRateSpikeBins(R,0.1,1000);
% figure
% stem(S,'k','linewidth',1)



%% Add the approprate folders to the path
%Path of the SOSpikingModel repository

repopath = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity'; 
%repopath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity'; 
addpath(genpath(repopath))

figfolder = [repopath,'/Modeling/Figures/EIBalance'];
savefolder = [repopath,'/Modeling/Simulation_Data'];

SAVESIM = false;

%%
%Neuron Parameters
cellparams.v_r = -55;
cellparams.v_th = -45;
cellparams.g_L = 182/18; %(nS)
cellparams.C = 182; %pF

cellparams.E_e = 0;
cellparams.E_i = -70;
cellparams.E_L = -65;
cellparams.E_h = 0;

%Input Parameters
synparams.w_e =1; %synaptic weight (nS)
%Set w_i such that IPSC~EPSC
w_i = synparams.w_e .* -(cellparams.v_r-cellparams.E_e)./(cellparams.v_r-cellparams.E_i);

synparams.w_i = w_i; %synaptic weight (nS)
synparams.tau_se = 5; %ms
synparams.tau_si = 5; %ms
synparams.K_e = 1;
synparams.K_i = 1;

rates.R_e = 1;
rates.R_i = 1;
rates.g_h = 0;


%% Simulations in Rate space
R_es = logspace(2,5,40);
R_is = logspace(2,5,40);

for ee = 1:length(R_es)
    ee
    rates.R_e = R_es(ee);
    for ii = 1:length(R_is)
        ii
        rates.R_i = R_is(ii);
    [ spikestats(ee,ii) ] = NoisyInputSims( cellparams,synparams,rates,'showfig',false );
    
    spkrate(ee,ii) = spikestats(ee,ii).rate;
    ISICV(ee,ii) = spikestats(ee,ii).ISI_CV;
    end
end
    
if SAVESIM
    save([savefolder,'/isistats.mat']);
end
    
%%
numspks = arrayfun(@(X) sum(~isnan(X.ISIs)),spikestats);

%% Calculate V_inf and Gamma
v_inf = -65:10:-35;
%v_inf = cellparams.v_th;

[R,V] = meshgrid(R_es./1000,v_inf);

D_e = cellparams.E_e - V;
D_i = cellparams.E_i - V;
D_L = cellparams.E_L - V;

%A_e = D_e.*synparams.w_e.*synparams.tau_se.*synparams.K_e;
%A_i = D_i.*synparams.w_i.*synparams.tau_si.*synparams.K_i;

g_e = synparams.w_e.*synparams.tau_se.*synparams.K_e.*R;
%g_i = synparams.w_i.*synparams.tau_si.*synparams.K_i.*R_is;

Ri_Vinf = -(g_e.*D_e+cellparams.g_L.*D_L)./(D_i.*synparams.w_i.*synparams.tau_si);
Ri_Vinf = Ri_Vinf.*1000;
Ri_Vinf(Ri_Vinf<0) = nan;

g_e = synparams.w_e.*synparams.tau_se.*synparams.K_e.*R_es./1000;
g_i = synparams.w_i.*synparams.tau_si.*synparams.K_i.*R_is./1000;
[G_e,G_i] = meshgrid(g_e,g_i);
gamma = 1+(1/cellparams.g_L).*(G_e+G_i);




%% Simulations in Vinf/Gamma space

Vinfs = linspace(-55,-35,30);
Gammas = logspace(0,2,30);

[V,G] = meshgrid(Vinfs,Gammas);
[ erate,irate ] = CondLIFReparm( V,G,cellparams,synparams );
for vv = 1:length(Vinfs)
    vv
    parfor ggg = 1:length(Gammas)
        %gg
    [ spikestats_VG(vv,ggg) ] = NoisyInputSims( cellparams,synparams,...
        [erate(vv,ggg) irate(vv,ggg) 0],'showfig',false );
    
    spkrate_VG(vv,ggg) = spikestats_VG(vv,ggg).rate;
    ISICV_VG(vv,ggg) = spikestats_VG(vv,ggg).ISI_CV;
    end
end
%%
if SAVESIM
    save([savefolder,'/isistats4.mat']);
end
%%
numspks_VG = arrayfun(@(X) sum(~isnan(X.ISIs)),spikestats_VG);

%% examples

subthreshex_lowG.Vinf = -50;
subthreshex_lowG.Gammma = 5;
[ subthreshex_lowG.R_e,subthreshex_lowG.R_i ] = ...
    CondLIFReparm( subthreshex_lowG.Vinf,subthreshex_lowG.Gammma,cellparams,synparams );
rates.R_e = subthreshex_lowG.R_e;
rates.R_i = subthreshex_lowG.R_i;
[subthreshex_lowG.spikestats,subthreshex_lowG.fig] = NoisyInputSims( cellparams,synparams,rates,...
    'showfig',true,'figfolder',figfolder ); 
%%
% subthreshex_highG.Vinf = -50;
% %subthreshex_highG.Gammma = 30;
% subthreshex_highG.Gammma = 30;
% [ subthreshex_highG.R_e,subthreshex_highG.R_i ] = ...
%     CondLIFReparm( subthreshex_highG.Vinf,subthreshex_highG.Gammma,cellparams,synparams );
% rates.R_e = subthreshex_highG.R_e;
% rates.R_i = subthreshex_highG.R_i;
% [subthreshex_highG.spikestats,subthreshex_highG.fig] = NoisyInputSims( cellparams,synparams,rates,...
%     'showfig',true,'figfolder',figfolder ); 

%%
supthreshex.Vinf = -40;
supthreshex.Gammma = 5;
[ supthreshex.R_e,supthreshex.R_i ] = ...
    CondLIFReparm( supthreshex.Vinf,supthreshex.Gammma,cellparams,synparams );
rates.R_e = supthreshex.R_e;
rates.R_i = supthreshex.R_i;
[supthreshex.spikestats,supthreshex.fig] = NoisyInputSims( cellparams,synparams,rates,...
    'showfig',true,'figfolder',figfolder );

%%
logISIbins = linspace(0,3.5,25);
supthreshex.ISIhist = hist(log10(supthreshex.spikestats.ISIs),logISIbins);
%subthreshex_highG.ISIhist = hist(log10(subthreshex_highG.spikestats.ISIs),logISIbins);
subthreshex_lowG.ISIhist = hist(log10(subthreshex_lowG.spikestats.ISIs),logISIbins);

supthreshex.returnISIhist = hist3([log10(supthreshex.spikestats.ISIs(1:end-1))',...
    log10(supthreshex.spikestats.ISIs(2:end))'],...
    {logISIbins,logISIbins});
% subthreshex_highG.returnISIhist = hist3([log10(subthreshex_highG.spikestats.ISIs(1:end-1))',...
%     log10(subthreshex_highG.spikestats.ISIs(2:end))'],...
%     {logISIbins,logISIbins});
subthreshex_lowG.returnISIhist = hist3([log10(subthreshex_lowG.spikestats.ISIs(1:end-1))',...
    log10(subthreshex_lowG.spikestats.ISIs(2:end))'],...
    {logISIbins,logISIbins});
%%
figure
histmap = flipud(colormap(gray));
colormap(gcf,histmap)
subplot(6,4,1)
    bar(logISIbins,supthreshex.ISIhist,'facecolor','k')
    
    box off
    xlim(logISIbins([1 end]))
    LogScale('x',10)
    xlabel('ISI (ms)')
subplot(6,6,4)
    imagesc(logISIbins,logISIbins,supthreshex.returnISIhist)
    axis xy
    LogScale('xy',10)
    xlabel('ISI_n');ylabel('ISI_n+1')
    
% subplot(4,3,7)
%     bar(logISIbins,subthreshex_highG.ISIhist)
%     LogScale('x',10)
%     box off
%     xlim(logISIbins([1 end]))
%     xlabel('ISI (ms)')
% subplot(4,4,11)
%     imagesc(logISIbins,logISIbins,subthreshex_highG.returnISIhist)
%     axis xy
%     LogScale('xy',10)
%     xlabel('ISI_n');ylabel('ISI_n+1')
    
subplot(6,4,9)
    bar(logISIbins,subthreshex_lowG.ISIhist,'facecolor','k')
    
    box off
    xlim(logISIbins([1 end]))
    LogScale('x',10)
    xlabel('ISI (ms)')
subplot(6,6,16)
    imagesc(logISIbins,logISIbins,subthreshex_lowG.returnISIhist)
    axis xy
    LogScale('xy',10)
    xlabel('ISI_n');ylabel('ISI_n+1')

NiceSave('ExSpikeStats',figfolder,'CondLIF')
%% Figure
spklim = 900;
figure
    subplot(3,3,1)
        h = imagesc(log10(R_es),log10(R_is),log10(spkrate'.*1000));
        hold on
        set(h, 'AlphaData', (numspks>spklim)') 
        plot(log10(R_es),log10(Ri_Vinf),'k--','linewidth',0.5)
        plot(log10(R_es),log10(Ri_Vinf(3,:)),'k','linewidth',1)
        %plot(log10(subthreshex_highG.R_e),log10(subthreshex_highG.R_i),'r*')
        plot(log10(subthreshex_lowG.R_e),log10(subthreshex_lowG.R_i),'r*')
        plot(log10(supthreshex.R_e),log10(supthreshex.R_i),'r*')
        colorbar
        axis xy
        caxis([-1 3.5])
        LogScale('xy',10)
        LogScale('c',10)
        title('Mean Rate');
        xlim([2 4.5]);ylim([2 4.5])
        xlabel('K_eR_e (Hz)');ylabel('K_iR_i (Hz)')
    subplot(3,3,2)
        h = imagesc(log10(R_es),log10(R_is),ISICV');
        hold on
        set(h, 'AlphaData', (numspks>spklim)') 
        plot(log10(R_es),log10(Ri_Vinf),'k--','linewidth',0.5)
        plot(log10(R_es),log10(Ri_Vinf(3,:)),'k','linewidth',1)
        %plot(log10(subthreshex_highG.R_e),log10(subthreshex_highG.R_i),'r*')
        plot(log10(subthreshex_lowG.R_e),log10(subthreshex_lowG.R_i),'r*')
        plot(log10(supthreshex.R_e),log10(supthreshex.R_i),'r*')
        colorbar
        caxis([0.5 1.5])
        xlim([2 4.5]);ylim([2 4.5])
        axis xy
        LogScale('xy',10)
        title('CV_I_S_I');
        xlabel('K_eR_e (Hz)');ylabel('K_iR_i (Hz)')
        

    subplot(3,3,7)
        h = imagesc(Vinfs,log10(Gammas),log10(spkrate_VG.*1000));
        set(h, 'AlphaData', (numspks_VG>spklim)) 
        hold on
        plot(cellparams.v_th.*[1 1],log10(Gammas([1 end])),'k','linewidth',1)
        %plot((subthreshex_highG.Vinf),log10(subthreshex_highG.Gammma),'r*')
        plot((subthreshex_lowG.Vinf),log10(subthreshex_lowG.Gammma),'r*')
        plot((supthreshex.Vinf),log10(supthreshex.Gammma),'r*')
        axis xy
        colorbar
        title('Mean Rate');
        caxis([-1 3.5])
        LogScale('y',10)
        LogScale('c',10)
        ylim([0 2])
        xlabel('V_i_n_f (mV)');ylabel('Gamma')
    subplot(3,3,8)
        h = imagesc(Vinfs,log10(Gammas),(ISICV_VG));
        set(h, 'AlphaData', (numspks_VG>spklim)) 
        hold on
        plot(cellparams.v_th.*[1 1],log10(Gammas([1 end])),'k','linewidth',1)
       %plot((subthreshex_highG.Vinf),log10(subthreshex_highG.Gammma),'r*')
        plot((subthreshex_lowG.Vinf),log10(subthreshex_lowG.Gammma),'r*')
        plot((supthreshex.Vinf),log10(supthreshex.Gammma),'r*')
        axis xy
        colorbar
        caxis([0.5 1.5])
        ylim([0 2])
        title('CV_I_S_I');
        LogScale('y',10)
        xlabel('V_i_n_f (mV)');ylabel('Gamma')

NiceSave('ISIbyInputRates',figfolder,'CondLIF')

%% Figure: converting Vinf<->Gamma
irateplot = irate;irateplot(irateplot<0)=nan;
%iheatmap = 

figure
subplot(3,3,1)
    imagesc(log10(R_es),log10(R_is),log10(gamma)')
    hold on
            plot(log10(R_es),log10(Ri_Vinf),'k--','linewidth',0.5)
            plot(log10(R_es),log10(Ri_Vinf(3,:)),'k','linewidth',1)
    LogScale('xy',10)
    xlim([2 4.5]);ylim([2 4.5])
    axis xy
    colorbar
    caxis([0 2])
    title('Gamma')
    LogScale('c',10)
    xlabel('K_eR_e (Hz)');ylabel('K_iR_i (Hz)')
    
subplot(4,4,3)
    h = imagesc(Vinfs,log10(Gammas),log10(erate));
        set(h, 'AlphaData', (irate>0)) 
    axis xy
    colorbar
    caxis([1 5])
    LogScale('c',10)
    LogScale('y',10)
    
    title('E Rate')
    xlabel('V_i_n_f');ylabel('Gamma')
    
subplot(4,4,4)
    h = imagesc(Vinfs,log10(Gammas),log10(irateplot));
    set(h, 'AlphaData', (irate>0)) 
    axis xy
    colorbar
    caxis([1 5])
    LogScale('c',10)
    title('I Rate')
    LogScale('y',10)
    
    xlabel('V_i_n_f');ylabel('Gamma')

NiceSave('VinfGamma',figfolder,'CondLIF')
        

        %%

pickgamma = [5 15 20];

numbins = 50;
clear ISIhist
ISIhist.logISIbins = linspace(-3.5,1.5,numbins);
ISIhist.FI = zeros(numbins,length(Vinfs),length(pickgamma));
for gg = 1:length(pickgamma)
    for vv = 1:length(Vinfs)
        ISIhist.FI(:,vv,gg) = hist(log10(spikestats_VG(pickgamma(gg),vv).ISIs./1000),ISIhist.logISIbins);
    end
end

pickRi = [5 15 20];

numbins = 50;
ISIhist.FI_EI = zeros(numbins,length(Vinfs),length(pickRi));
for ii = 1:length(pickRi)
    for ee = 1:length(R_es)
        ISIhist.FI_EI(:,ee,ii) = hist(log10(spikestats(ee,pickRi(ii)).ISIs./1000),ISIhist.logISIbins);
    end
end


%%
plotspkerate_VG = spkrate_VG;
plotspkerate_VG(numspks_VG<spklim) = nan;

plotISICV_VG = ISICV_VG;
plotISICV_VG(numspks_VG<spklim) = nan;

plotspkerate = spkrate;
plotspkerate(numspks<spklim) = nan;

plotISICV = ISICV;
plotISICV(numspks<spklim) = nan;


figure
    subplot(6,3,1)
        plot(Vinfs,log10(plotspkerate_VG(pickgamma,:).*1000),'linewidth',2)
        %legend({num2str(Gammas(pickgamma(1))),num2str(Gammas(pickgamma(2))),num2str(Gammas(pickgamma(3)))},'location','southeast')
        LogScale('y',10)
        ylim([-1 3.5]);xlim(Vinfs([1 end])) 
        xlabel('V_I_n_f');ylabel('Rate (Hz)')
        ylim([-1 3.5]);xlim(Vinfs([1 end])) 
        box off
        
    subplot(6,3,4)
        plot(Vinfs,(plotspkerate_VG(pickgamma,:).*1000),'linewidth',2)
        %legend({num2str(Gammas(pickgamma(1))),num2str(Gammas(pickgamma(2))),num2str(Gammas(pickgamma(3)))},'location','northwest')
        xlabel('V_I_n_f');ylabel('Rate (Hz)')
        box off
        ylim([0 1500]);xlim(Vinfs([1 end])) 
        
    subplot(6,3,5)
        plot(Vinfs,(plotspkerate_VG(pickgamma,:).*1000),'linewidth',2)
        %legend({num2str(Gammas(pickgamma(1))),num2str(Gammas(pickgamma(2))),num2str(Gammas(pickgamma(3)))},'location','northwest')
        xlabel('V_I_n_f');ylabel('Rate (Hz)')
        box off
        ylim([0 30]);xlim([-55 -45]) 
                legend({num2str(round(Gammas(pickgamma(1)))),...
            num2str(round(Gammas(pickgamma(2)))),num2str(round(Gammas(pickgamma(3))))},'location','eastoutside')
        
    subplot(6,3,7)
        plot(Vinfs,(plotISICV_VG(pickgamma,:)),'linewidth',2)
        %LogScale('y',10)
        xlabel('V_I_n_f');ylabel('CV_I_S_I')
        xlim(Vinfs([1 end])) ;ylim([0 1.5])
        box off
        
        for gg = 1:length(pickgamma)
    subplot(6,3,3*gg+7)
    imagesc(Vinfs,ISIhist.logISIbins,ISIhist.FI(:,:,gg))
    LogScale('y',10)
    axis xy
        end
        
        
    subplot(6,3,3)
        plot(log10(R_es),log10(plotspkerate(:,pickRi).*1000),'linewidth',2)
        %legend({num2str(R_is(pickRi(1))),num2str(R_is(pickRi(2))),num2str(R_is(pickRi(3)))},'location','southeast')
        axis tight
        xlim([2.7 4.5])
        LogScale('xy',10)
        xlabel('R_e');ylabel('Rate (Hz)')
        %axis tight
        box off
        
    subplot(6,3,6)
        plot(log10(R_es),(plotspkerate(:,pickRi).*1000),'linewidth',2)
        %legend({num2str(R_is(pickRi(1))),num2str(R_is(pickRi(2))),num2str(R_is(pickRi(3)))},'location','northwest')
        xlabel('R_e');ylabel('Rate (Hz)')
        %axis tight
        
        axis tight
        xlim([2.7 4.5])
        LogScale('x',10)
        box off
        ylim([0 1000])
        %NiceSave('FICurve',figfolder,'CondLIF')
        
    subplot(6,3,9)
        plot(log10(R_es),(plotISICV(:,pickRi)),'linewidth',2)
        %legend({num2str(round(R_is(pickRi(1)))),num2str(round(R_is(pickRi(2)))),num2str(round(R_is(pickRi(3))))},'location','eastoutside')
        axis tight
        xlim([2.7 4.5])
        LogScale('x',10)
        xlabel('R_e');ylabel('CV_I_S_I')
        %
        box off
        
        for gg = 1:length(pickRi)
    subplot(6,3,3*gg+9)
    imagesc(Vinfs,ISIhist.logISIbins,ISIhist.FI_EI(:,:,gg))
    LogScale('y',10)
    axis xy
        end
        
NiceSave('FICurve_VG',figfolder,'CondLIF')

%%
figure

NiceSave('FICurve_EI',figfolder,'CondLIF')


%% BELOW FOR LATER! MOVE TO OTHER SCRIPT
%% Fit exponential to the datas. 

expfittype = fittype( 'A.*k.^(V-V0)','dependent',{'R'},'independent',{'V'},'coefficients',{'A','k','V0'})
%%
options = fitoptions(expfittype);
options.StartPoint = [1,1,-50];
options.Lower = [0 1 -70];
options.Upper = [10 10 -30];
options.MaxIter = 1e15;
options.MaxFunEvals = 1e15;
expfit = fit(Vinfs',(spkrate_VG(pickgamma(2),:).*1000)',expfittype,options)

%%
PLfittype = fittype( '(exp(k)).*heaviside(V-V0).*(V-V0).^n','dependent',{'R'},'independent',{'V'},'coefficients',{'n','k','V0'})
%%
options = fitoptions(PLfittype);
options.StartPoint = [1,1,-50];
options.Lower = [1 -6 -70];
options.Upper = [10 1 -30];
options.MaxIter = 1e15;
options.MaxFunEvals = 1e15;

whichgamma = 1;
vmax = -43;
%figure
for gg = 2:length(Gammas)
PLfit = fit(Vinfs(Vinfs<vmax & ~isnan(spkrate_VG(gg,:)))',(spkrate_VG(gg,(Vinfs<vmax& ~isnan(spkrate_VG(gg,:)))).*1000)',PLfittype,options);
V0_eachfit(gg) = PLfit.V0;
n(gg) = PLfit.n;
k(gg) = PLfit.k;

%subplot(5,6,gg)
%plot(PLfit,Vinfs,(spkrate_VG(gg,:).*1000))
end

%%
figure
subplot(2,2,1)
plot(log(Gammas),V0_eachfit,'k.')
xlabel('G');ylabel('V0')
%LogScale('x',10)
subplot(2,2,2)
plot(log(Gammas),n,'k.')
%LogScale('x',10)
xlabel('ln(G)');ylabel('n')
subplot(2,2,3)
plot(log(Gammas),(k),'k.')
%LogScale('x',10)
xlabel('ln(G)');ylabel('ln(k)')

%%
linPLfittype = fittype( 'exp(mk.*log(gamma)+bk).*heaviside(V-V0).*(V-V0).^(mn.*log(gamma)+bn)',...
    'dependent',{'R'},'independent',{'V','gamma'},'coefficients',{'mn','bn','mk','bk','V0'});
%%

options = fitoptions(linPLfittype);
options.StartPoint = [0.5,1,-1,0,-50];
options.Upper = [3 3  0  2 -30];
options.Lower = [0 0 -5 -4 -70];
options.MaxIter = 1e15;
options.MaxFunEvals = 1e15;
%% Reshape everything
%condition: V<-43, gamma<100 ~isnan(rate)
V_fit = V(:);
G_fit = G(:);
rate_fit = spkrate_VG(:).*1000;

conditions = V_fit<=-44 & G_fit<50 &G_fit>1 & ~isnan(rate_fit);

V_fit = V_fit(conditions);
G_fit = G_fit(conditions);
rate_fit = rate_fit(conditions);

%%
linPLfit = fit([V_fit G_fit],rate_fit,linPLfittype,options)
%%
figure
plot(linPLfit,[V_fit G_fit],rate_fit)

%%
figure
mn = linPLfit.mn;
bn = linPLfit.bn;
mk = linPLfit.mk;
bk= linPLfit.bk;
V0= linPLfit.V0;
subplot(2,2,4)
    imagesc(Vinfs,log10(Gammas),log10(exp(mk.*log(G)+bk).*heaviside(V-V0).*(V-V0).^(mn.*log(G)+bn)))
    axis xy
    colorbar
    caxis([-1 5])
    ylim([0 2])
    LogScale('y',10)
    xlabel('V_I_n_f');ylabel('Gamma')

subplot(2,2,3)
    imagesc(Vinfs,log10(Gammas),log10(spkrate_VG.*1000))
    axis xy
    colorbar
    caxis([-1 5])
    ylim([0 2])
        LogScale('y',10)
    xlabel('V_I_n_f');ylabel('Gamma')    
subplot(3,3,1)
    plot(log(Gammas),V0_eachfit,'k.')
    hold on
    plot(log(Gammas),V0.*ones(size(Gammas)),'k')
    xlabel('G');ylabel('V0')
    %LogScale('x',10)
    xlim([0 log(100)])
subplot(3,3,2)
    plot(log(Gammas),n,'k.')
    hold on
    plot(log(Gammas),mn.*log(Gammas)+bn,'k')
    %LogScale('x',10)
    xlabel('ln(G)');ylabel('n')
    xlim([0 log(100)])
subplot(3,3,3)
    plot(log(Gammas),(k),'k.')
    hold on
    plot(log(Gammas),mk.*log(Gammas)+bk,'k')
    %LogScale('x',10)
    xlabel('ln(G)');ylabel('ln(k)')
    xlim([0 log(100)])
%%
xtest = linspace(-2,2,100);
figure
plot(xtest,0.002.*heaviside(xtest-1).*(xtest-1).^2,'k')


%% FI curve
R_es = logspace(2,3,20);
rates.R_i = 500;
for ee = 1:length(R_es)
    ee
    rates.R_e = R_es(ee);
    [ spikestats_FI(ee) ] = NoisyInputSims( cellparams,synparams,rates,'showfig',false );
    spkrate_FI(ee) = spikestats_FI(ee).rate;
    ISICV_FI(ee) = spikestats_FI(ee).ISI_CV;
end
   
%%
isihist.logbins = linspace(-3,log10(20),200);
allISIs = cat(1,spikestats_FI(:).ISIs);
isihist.hist = hist(log10(allISIs'),isihist.logbins);
%isihist.CV = nanmean(allISIs,2);
%ISICV2_FI
%%
figure
subplot(4,2,1)
    imagesc(log10(R_es),isihist.logbins,isihist.hist)
    hold on
    plot(log10(R_es),log10(1./spkrate_FI),'o-')
    axis xy
    LogScale('xy',10)
subplot(4,2,3)
    plot(log10(R_es),ISICV_FI,'o-')
    xlim(log10(R_es([1 end])))
    LogScale('x',10)



%%
synparams.w_e = 400; %synaptic weight (pS?)
synparams.w_i = 200; %synaptic weight (pS?)
synparams.tau_se = 0.01;
synparams.tau_si = 0.01;
synparams.K_e = 500;
synparams.K_i = 500;

rates.R_e = 1;
rates.R_i = 5;
rates.g_h = 0;

cellparams.v_r = -55;
cellparams.v_th = -45;
cellparams.g_L = 182/18;
cellparams.C = 182;


cellparams.E_e = 0;
cellparams.E_i = -70;
cellparams.E_L = -65;
cellparams.E_h = 0;


w_es = linspace(0,2000,21);
w_is = linspace(0,2000,21);

for ee = 1:length(w_es)
    ee
    synparams.w_e = w_es(ee);
    for ii = 1:length(w_is)
        ii
        synparams.w_i = w_is(ii);
    [ spikestats(ee,ii) ] = NoisyInputSims( cellparams,synparams,rates,'showfig',false );
    
    spkrate(ee,ii) = spikestats(ee,ii).rate;
    ISICV(ee,ii) = spikestats(ee,ii).ISI_CV;
    end
end


%% Calculate sub-superthrehsold
B_e = D_e.*synparams.tau_se.*synparams.K_e;
B_i = D_i.*synparams.tau_si.*synparams.K_i;

%% esitmate (interpolate) R = RE line




%%

figure
newmap = [[1 1 1];colormap(gcf)];
colormap(newmap)
    subplot(2,2,1)
        imagesc(w_es,w_is,log10(spkrate)')
        hold on
        plot(w_es,-(w_es.*B_e.*rates.R_e)./(B_i.*rates.R_i) - D_L./(B_i.*rates.R_i),...
            'r','LineWidth',2)
        colorbar
        axis xy
        %LogScale('xy',10)
        xlabel('w_e');ylabel('w_i')
    subplot(2,2,2)
        imagesc(w_es,w_is,ISICV')
        hold on
        plot(w_es,-(w_es.*B_e.*rates.R_e)./(B_i.*rates.R_i) - D_L./(B_i.*rates.R_i),...
            'r','LineWidth',2)
        colorbar
        caxis([0 1.5])
        axis xy
       % LogScale('xy',10)
        xlabel('w_e');ylabel('w_i')