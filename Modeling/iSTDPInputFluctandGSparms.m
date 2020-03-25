function iSTDPandGSparms(savepath)

%%
savepath = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/Modeling/Simulation_Data/iSTDPInputFluctandGSparms';

%%
TimeParams.dt = 0.1;
TimeParams.SimTime = 10000;

%Poisson Rate (add to Brunel sim)
g = 5;
%g = 4; %Initial strength of Inhibitoon (relative to excitation)

clear parms


parms.EPopNum = 1000;
parms.IPopNum = 250;
parms.u_0 = 0;

parms.V_rest = 0;
%parms.delay_s = 1.2;
parms.delay_s = 1.8.*rand(parms.EPopNum+parms.IPopNum,parms.EPopNum+parms.IPopNum)+1.2;
parms.g = g;

parms.V_th =20;
parms.tau_m = 20;
parms.V_reset = 10;
parms.t_ref = 1;

%Feedforward parameters
parms.N_FF = 1000;
parms.K_FF = 250;
parms.J_FF = 0.5;

%Conectivity: In degree
gamma = 0.25;
parms.Kee = 250;
parms.Kie = parms.Kee;
parms.Kei = parms.Kee.*gamma;
parms.Kii = parms.Kee.*gamma;


netname = 'weaklybalanced';
switch netname
    case 'weaklybalanced'
        parms.J = 0.4;
        parms.ex_rate = 20;
        parms.ex_rate = 40; %Not same as other script...
    case 'strongrecurrent'
        parms.J = 1.5; 
        parms.ex_rate = 20;
    case 'CA1like'
        parms.J = 0.4;
        parms.ex_rate = 30; 
        parms.Kee = 0; 
end


parms.LearningRate = 1e-2;
parms.TargetRate = [sort(exp(randn(parms.EPopNum,1)));nan(parms.IPopNum,1)]; %Target Rate for Excitatory cells (units of Hz)
parms.tauSTDP = 20;    %Time Constant for the STDP curve (Units of ms)


%v_th = th/(C_e.*J.*tau);
gammas_E = [0.1 0.3 0.5 0.7];
%gammas_E = [0.5];
clear SimValues
for gg = 1:length(gammas_E)
    loopparms = parms;
    loopparms.Kei = loopparms.Kee.*gammas_E(gg);
    loopparms.Kii = loopparms.Kee.*gammas_E(gg);
    loopparms.g = (6./4)./gammas_E(gg);
tic 
[SimValues{gg}] = Run_LIF_iSTDP(loopparms,TimeParams,'showprogress',true,...
    'cellout',true,'save_dt',100);
toc
end
%%
for gg = 1:length(gammas_E)
    PlotSimRaster(SimValues{gg},TimeParams.SimTime-[1000 0])
    NiceSave('iSTDPRaster',savepath,['gamma',num2str(gammas_E(gg))])
end
%% Save/load
filename = fullfile(savepath,['TrainedNet_',netname]);
save(filename,'SimValues','parms','TimeParams','netname','gammas_E')
%load(filename)


%% ISI 
for gg = 1:length(gammas_E) 
clear spikes
spikes.times = cellfun(@(X) X./1000,SimValues{gg}.spikesbycell,'UniformOutput',false);
spikes.UID = 1:length(SimValues{gg}.spikesbycell);
CellClass = cell(1,length(spikes.times));
CellClass(SimValues{gg}.EcellIDX) = {'E'};
CellClass(SimValues{gg}.IcellIDX) = {'I'};
timewindows.initialization = [0 20];
timewindows.equib = [20 TimeParams.SimTime./1000];
ISIstats(gg) = bz_ISIStats(spikes,'ints',timewindows,'showfig',true,'cellclass',CellClass);
end
%%
NiceSave('ISIStats',savepath,netname)
%%
figure
subplot(2,2,1)
plot(log10(parms.TargetRate),log10(ISIstats.summstats.equib.meanrate),'k.','markersize',2)
axis tight
box off
hold on
xlabel('Target Rate');ylabel('Simulated Rate')
UnityLine
LogScale('xy',10)
NiceSave('RateTarget',savepath,netname)

%% Fit the gamma model...
%Use high rate cells...
AScost = 0.0000; %0.13 for data
for gg = 1:length(gammas_E)
    logISIhist = ISIstats(gg).ISIhist.equib.log;
    usecells = randsample(SimValues{gg}.EcellIDX([end-500:end]),100);
    logISIhist = logISIhist(usecells,:)';
    logtimebins = ISIstats(gg).ISIhist.logbins;
    logISIhist = logISIhist./mode(diff(logtimebins));

GammaFit(gg) = bz_FitISISharedGammaModes(logISIhist,logtimebins,...
    'numAS',1,...
    'figfolder',savepath,'basePath',savepath,'figname',netname,...
    'AScost_lambda',AScost,'AScost_p',1/2,'ASguess',true,'MScost',3);

end

%%
figure %GS_CV values
subplot(3,2,1)
BoxAndScatterPlot({log10(GammaFit(1).sharedfit.GSCVs),...
    log10(GammaFit(2).sharedfit.GSCVs),log10(GammaFit(3).sharedfit.GSCVs),...
    log10(GammaFit(4).sharedfit.GSCVs)},'labels',parms.Kee.*gammas_E)
hold on
ylabel('Ground State CV');xlabel('K_E_I')
plot(xlim(gca),[0 0],'k--')
LogScale('y',10)
box off
NiceSave('CV',savepath,netname)
