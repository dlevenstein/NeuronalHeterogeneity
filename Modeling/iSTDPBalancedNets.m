function iSTDPBalancednets(savepath)

%%
savepath = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/Modeling/Simulation_Data';

%%
TimeParams.dt = 0.1;
TimeParams.SimTime = 100000;

%Poisson Rate (add to Brunel sim)
g = 5;
g = 4; %Initial strength of Inhibitoon (relative to excitation)
v_norm = 2;

clear parms

parms.EPopNum = 1000;
parms.IPopNum = 250;
parms.u_0 = 20.*v_norm;
parms.u_0 = 0;

parms.V_rest = 0;
parms.delay_s = 1.2;
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


netname = 'CA1like';
switch netname
    case 'weaklybalanced'
        parms.J = 0.4;
        parms.ex_rate = 20;
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

tic 
[SimValues] = Run_LIF_iSTDP(parms,TimeParams,'showprogress',true,...
    'cellout',true,'save_dt',100);
toc
%%

PlotSimRaster(SimValues,TimeParams.SimTime-[400 0])
NiceSave('iSTDPRaster',pwd,netname)
%% Save/load
filename = fullfile(savepath,['TrainedNet_',netname]);
save(filename,'SimValues','parms','TimeParams','netname')
%load(filename)
%%
%NiceSave('iSTDPRaster',pwd,netname)
%%
figure
subplot(2,2,1)
imagesc(SimValues.WeightMat_initial)
caxis([-2 0.2])
crameri('vik','pivot',0)
colorbar
subplot(2,2,2)
imagesc(SimValues.WeightMat)
colorbar

caxis([-2 0.2])
crameri('vik','pivot',0)

%% ISI 
clear spikes
spikes.times = cellfun(@(X) X./1000,SimValues.spikesbycell,'UniformOutput',false);
spikes.UID = 1:length(SimValues.spikesbycell);
CellClass = cell(1,length(spikes.times));
CellClass(SimValues.EcellIDX) = {'E'};
CellClass(SimValues.IcellIDX) = {'I'};
timewindows.initialization = [0 10];
timewindows.equib = [10 TimeParams.SimTime./1000];
ISIstats = bz_ISIStats(spikes,'ints',timewindows,'showfig',true,'cellclass',CellClass);

%%
NiceSave('ISIStats',pwd,netname)
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
NiceSave('RateTarget',pwd,netname)