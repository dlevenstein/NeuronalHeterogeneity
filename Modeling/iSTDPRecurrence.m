function iSTDPRecurrence(savepath)

%%
savepath = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/Modeling/Simulation_Data/Recurrence';

%%
TimeParams.dt = 0.1;
TimeParams.SimTime = 80000;

%Poisson Rate (add to Brunel sim)
g = 5;
%g = 4; %Initial strength of Inhibitoon (relative to excitation)
v_norm = 2;

clear parms

parms.EPopNum = 1000;
parms.IPopNum = 250;
parms.u_0 = 0;

parms.V_rest = 0;
parms.delay_s = 8.9.*rand(parms.EPopNum+parms.IPopNum,1)+1.1; %grid later
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


parms.J = 0.4;
parms.ex_rate = 20;
%

%Train with fluctuating rate
duration = TimeParams.SimTime;
dt = 0.1;
save_dt = 1;
numsignals = 1;

theta = 1./2000; %1s (1000ms) timescale
sigma = 10;

[ X,T ] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
ex_rate_fun = @(t) interp1(T,X,t,'nearest')+parms.ex_rate;
%
parms.ex_rate = ex_rate_fun;
%%
parms.LearningRate = 1e-2;
parms.TargetRate = [sort(exp(randn(parms.EPopNum,1)));nan(parms.IPopNum,1)]; %Target Rate for Excitatory cells (units of Hz)
parms.tauSTDP = 20;    %Time Constant for the STDP curve (Units of ms)

[SimValues] = Run_LIF_iSTDP(parms,TimeParams,'showprogress',true,...
    'cellout',true,'save_dt',1000);

%%
twin = [10000 TimeParams.SimTime];
spikemat = PlotSimRaster(SimValues,twin,'ratebin',100);

for ss = 1:4
    sextileidx = [(ss-1)*round((parms.EPopNum./6))+1 : (ss)*floor((parms.EPopNum./6))];
    spikemat.sextilerates(ss,:) = nanmean(1000*(spikemat.E.data(:,sextileidx)),2);
end

%%

% subplot(4,1,3)
% plot(spikemat.E.timestamps,bz_NormToRange(spikemat.sextilerates))

tvec = linspace(twin(1),twin(2),10000);
subplot(4,1,4)
hold off
    plot(tvec,ex_rate_fun(tvec),'linewidth',2)
    %xlim(TimeParams.SimTime-[2000 0])
    ylabel('Input Rate (Hz)')
    
 %%
 corrmat = corr(spikemat.E.data,'type','spearman');
 
 %%
 figure
 imagesc(corrmat)
 caxis([-0.1 0.1])
 crameri('berlin','pivot',0)
 ColorbarWithAxis([-0.1 0.1],'Correlation')
 xlabel('Cell (sorted by target FR');
 ylabel('Cell (sorted by target FR')
 axis xy
%%
TimeParams.SimTime = 20000;

%v_th = th/(C_e.*J.*tau);
inputrates = logspace(0.5,1.5,8);
for rr = 1:length(inputrates)
    loopparms = parms;
    loopparms.ex_rate = inputrates(rr);
tic 
[SimValues_inputs(rr)] = Run_LIF_iSTDP(loopparms,TimeParams,'showprogress',true,...
    'cellout',true,'save_dt',2,'J_mat',SimValues.WeightMat);
toc
end

%%

for rr = 1:length(inputrates)
PlotSimRaster(SimValues_inputs(rr),TimeParams.SimTime-[2000 0])
% subplot(4,1,4)
% hold off
%     plot(SimValues_fluct{rr}.t,ex_rate_fun{rr}(SimValues_fluct{rr}.t),'linewidth',2,'color',sigmacolors(rr,:))
%     xlim(TimeParams.SimTime-[2000 0])
end


%%

for ss = 1:length(inputrates) 
    %clear spikes
    spikes(ss).times = cellfun(@(X) X./1000,SimValues_inputs(ss).spikesbycell,'UniformOutput',false);
    spikes(ss).UID = 1:length(SimValues_inputs(ss).spikesbycell);
    CellClass = cell(1,length(spikes(ss).times));
    CellClass(SimValues_inputs(ss).EcellIDX) = {'E'};
    CellClass(SimValues_inputs(ss).IcellIDX) = {'I'};
    %timewindows.initialization = [0 20];
    %timewindows.equib = [0 TimeParams.SimTime];
    ISIstats(ss) = bz_ISIStats(spikes(ss),'showfig',true,'cellclass',CellClass);
    [popCCG(ss)] = PopCCG(spikes(ss),'showfig',true,'cellclass',CellClass);
    
    close all
end

%%
numratebins = 40;
numrvoltbins = 100;
FIStats.ratedist.bins = linspace(-1,1.5,numratebins);
FIStats.voltagedist.bins = linspace(-20,20,numrvoltbins);
FIStats.ratedist.pop = zeros(numratebins,length(inputrates));
FIStats.ratedist.cell = zeros(numratebins,length(inputrates));
FIStats.voltagedist.allcells = zeros(numrvoltbins,length(inputrates));
FIStats.ISIdist.bins = ISIstats(1).ISIhist.logbins;
FIStats.ISIdist.meanE = zeros(length(FIStats.ISIdist.bins),length(inputrates));

FIStats.popCCG.bins = popCCG(1).t_ccg;

for ss = 1:length(inputrates) 
    FIStats.ratedist.cell(:,ss) = hist(log10(ISIstats(ss).summstats.ALL.meanrate(SimValues_inputs(ss).EcellIDX)),FIStats.ratedist.bins);
    FIStats.ratedist.Icell(:,ss) = hist(log10(ISIstats(ss).summstats.ALL.meanrate(SimValues_inputs(ss).IcellIDX)),FIStats.ratedist.bins);
    
    allvoltage = SimValues_inputs(ss).V(SimValues_inputs(ss).EcellIDX,:);
    FIStats.voltagedist.allcells(:,ss) = hist(allvoltage(:),FIStats.voltagedist.bins);
    
    spikemat = PlotSimRaster(SimValues_inputs(ss),[0000 2000]);
    FIStats.ratedist.pop(:,ss) = hist(log10(spikemat.E.poprate),FIStats.ratedist.bins); 
    
    FIStats.ISIdist.meanE(:,ss) = mean(ISIstats(ss).ISIhist.ALL.log(SimValues_inputs(ss).EcellIDX,:),1);
    
    FIStats.popCCG.E(:,ss) = popCCG(ss).pop.E(:,1);
    FIStats.popCCG.I(:,ss) = popCCG(ss).pop.I(:,2);
end

%%
figure
subplot(3,3,1)
imagesc(log10(inputrates),FIStats.ratedist.bins,FIStats.ratedist.cell)
axis xy
LogScale('xy',10)

subplot(3,3,2)
imagesc(log10(inputrates),FIStats.ratedist.bins,FIStats.ratedist.Icell)
axis xy
LogScale('xy',10)


subplot(3,3,3)
imagesc(log10(inputrates),FIStats.voltagedist.bins,FIStats.voltagedist.allcells)
axis xy
LogScale('x',10)

subplot(3,3,4)
imagesc(log10(inputrates),FIStats.ratedist.bins,FIStats.ratedist.pop)
axis xy
LogScale('xy',10)


subplot(3,3,5)
imagesc(log10(inputrates),FIStats.ISIdist.bins,FIStats.ISIdist.meanE)
axis xy
LogScale('xy',10)

subplot(3,3,6)
imagesc(log10(inputrates),FIStats.popCCG.bins,FIStats.popCCG.I)
axis xy

subplot(3,3,7)
imagesc(log10(inputrates),FIStats.popCCG.bins,FIStats.popCCG.E)
axis xy
%LogScale('xy',10)

%%
%save('-v7.3')

%%

PlotSimRaster(SimValues,TimeParams.SimTime-[400 0])
NiceSave('iSTDPRaster',pwd,netname)
%% Save/load
filename = fullfile(savepath,['TrainedNet_',netname]);
%save(filename,'SimValues','parms','TimeParams','netname')
load(filename)

%% Effective Indegree
K_net = sum(SimValues.WeightMat,2);
K_E = sum(SimValues.WeightMat(:,SimValues.EcellIDX),2);
K_I = sum(SimValues.WeightMat(:,SimValues.IcellIDX),2);
K_EI = K_E./K_I;

figure
subplot(2,2,1)
plot(log10(parms.TargetRate),K_net,'.')
subplot(2,2,2)
plot(log10(parms.TargetRate),K_EI,'.')
subplot(2,2,3)
plot(log10(parms.TargetRate),K_E,'.')
subplot(2,2,4)
plot(log10(parms.TargetRate),K_I,'.')
%%
PlotSimRaster(SimValues,TimeParams.SimTime-[2000 0])
NiceSave('iSTDPRaster_late',pwd,netname)

PlotSimRaster(SimValues,[0000 2000])
NiceSave('iSTDPRaster_early',pwd,netname)
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