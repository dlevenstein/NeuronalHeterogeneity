

%TO DO: set DT, calculate probability of number of spikes


%% Add the approprate folders to the path
%Path of the SOSpikingModel repository

repopath = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity'; 
%repopath = '/Users/jonathangornet/Documents/GitHub/SOSpikingModel'; 
addpath(genpath(repopath))

figfolder = [repopath,'Modeling/Figures/EIBalance'];
savefolder = [repopath,'Modeling/Simulation_Data'];


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
synparams.tau_se = 5; %mS
synparams.tau_si = 5; %mS
synparams.K_e = 1;
synparams.K_i = 1;

rates.R_e = 1;
rates.R_i = 1;
rates.g_h = 0;


R_es = logspace(1,5,25);
R_is = logspace(1,5,25);
%%
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

[R,V] = meshgrid(R_es,v_inf);

D_e = cellparams.E_e - V;
D_i = cellparams.E_i - V;
D_L = cellparams.E_L - V;

%A_e = D_e.*synparams.w_e.*synparams.tau_se.*synparams.K_e;
%A_i = D_i.*synparams.w_i.*synparams.tau_si.*synparams.K_i;

g_e = synparams.w_e.*synparams.tau_se.*synparams.K_e.*R;
%g_i = synparams.w_i.*synparams.tau_si.*synparams.K_i.*R_is;

Ri_Vinf = -(g_e.*D_e+cellparams.g_L.*D_L)./(D_i.*synparams.w_i.*synparams.tau_si);


g_e = synparams.w_e.*synparams.tau_se.*synparams.K_e.*R_es./1000;
g_i = synparams.w_i.*synparams.tau_si.*synparams.K_i.*R_is./1000;
[G_e,G_i] = meshgrid(g_e,g_i);
gamma = 1+(1/cellparams.g_L).*(G_e+G_i);

%%
figure
subplot(2,2,1)
imagesc(log10(R_es),log10(R_is),log10(gamma)')
hold on
        plot(log10(R_es),log10(Ri_Vinf),'k--','linewidth',0.5)
        plot(log10(R_es),log10(Ri_Vinf(3,:)),'k','linewidth',1)
LogScale('xy',10)
axis xy
colorbar
xlabel('K_eR_e (Hz)');ylabel('K_iR_i (Hz)')

NiceSave('VinfGamma',figfolder,'CondLIF')

%% examples

%Vinf = -50;
%Gammma = 2;
%g = 1;
subthreshex_lowG.R_e = 500;
subthreshex_lowG.R_i = 300;
rates.R_e = subthreshex_lowG.R_e;
rates.R_i = subthreshex_lowG.R_i;
[subthreshex_lowG.spikestats,subthreshex_lowG.fig] = NoisyInputSims( cellparams,synparams,rates,...
    'showfig',true,'figfolder',figfolder ); 
%%
subthreshex_highG.R_e = 1500;
subthreshex.R_i = 800;
rates.R_e = subthreshex_highG.R_e;
rates.R_i = subthreshex_highG.R_i;
[subthreshex_highG.spikestats,subthreshex_highG.fig] = NoisyInputSims( cellparams,synparams,rates,...
    'showfig',true,'figfolder',figfolder ); 

%%
supthreshex.R_e = 1500;
supthreshex.R_i = 500;
rates.R_e = supthreshex.R_e;
rates.R_i = supthreshex.R_i;
[supthreshex.spikestats,supthreshex.fig] = NoisyInputSims( cellparams,synparams,rates,...
    'showfig',true,'figfolder',figfolder );


%% test
test.R_e = 10000;
test.R_i = 5000;
rates.R_e = test.R_e;
rates.R_i = test.R_i;
[test.spikestats,test.fig] = NoisyInputSims( cellparams,synparams,rates,...
    'showfig',true);

%%
logISIbins = linspace(-1.5,1,25);
supthreshex.ISIhist = hist(log10(supthreshex.spikestats.ISIs),logISIbins);
subthreshex.ISIhist = hist(log10(subthreshex.spikestats.ISIs),logISIbins);

supthreshex.returnISIhist = hist3([log10(supthreshex.spikestats.ISIs(1:end-1))',...
    log10(supthreshex.spikestats.ISIs(2:end))'],...
    {logISIbins,logISIbins});
subthreshex.returnISIhist = hist3([log10(subthreshex.spikestats.ISIs(1:end-1))',...
    log10(subthreshex.spikestats.ISIs(2:end))'],...
    {logISIbins,logISIbins});
%%
figure
subplot(4,3,7)
    bar(logISIbins,supthreshex.ISIhist)
    LogScale('x',10)
    box off
    xlim(logISIbins([1 end]))
subplot(4,4,11)
    imagesc(logISIbins,logISIbins,supthreshex.returnISIhist)
    axis xy
subplot(4,3,10)
    bar(logISIbins,subthreshex.ISIhist)
    LogScale('x',10)
    box off
    xlim(logISIbins([1 end]))
subplot(4,4,15)
    imagesc(logISIbins,logISIbins,subthreshex.returnISIhist)
    axis xy
%NEED:
%ISI hist (log scale)
%subthreshex.ISIs
%ISI return map (log scale)
%%
spklim = 100;
figure
    subplot(2,2,1)
        h = imagesc(log10(R_es),log10(R_is),log10(spkrate'));
        hold on
        set(h, 'AlphaData', (numspks>spklim)') 
        plot(log10(R_es),log10(Ri_Vinf),'k--','linewidth',0.5)
        plot(log10(R_es),log10(Ri_Vinf(2,:)),'k','linewidth',1)
        plot(log10(subthreshex.R_e),log10(subthreshex.R_i),'r*')
        plot(log10(supthreshex.R_e),log10(supthreshex.R_i),'r*')
        colorbar
        axis xy
        LogScale('xy',10)
        title('Mean Rate');
        xlabel('K_eR_e (Hz)');ylabel('K_iR_i (Hz)')
    subplot(2,2,2)
        h = imagesc(log10(R_es),log10(R_is),ISICV');
        hold on
        set(h, 'AlphaData', (numspks>spklim)') 
        plot(log10(R_es),log10(Ri_Vinf),'k--','linewidth',0.5)
        plot(log10(R_es),log10(Ri_Vinf(2,:)),'k','linewidth',1)
        colorbar
        caxis([0 1])
        axis xy
        LogScale('xy',10)
        title('CV_I_S_I');
        xlabel('K_eR_e (Hz)');ylabel('K_iR_i (Hz)')
        
        

NiceSave('ISIbyInputRates',figfolder,'CondLIF')



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