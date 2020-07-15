function iSTDPRingnet(savepath)

%%
savepath = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/Modeling/Simulation_Data';



%%
TimeParams.dt = 0.1;
TimeParams.SimTime = 250000;

%%

clear parms


%% Tuning Curve as a function of input sparseness*, weight, N

%%


parms.EPopNum = 1200;
parms.IPopNum = 300;
parms.u_0 = 0;

%Conectivity: In degree
gamma = 0.5; %initally 0.25 to match 4x less inhibitory cells
gammaI = 2; %relative E->I connectivity
parms.Kee = 300;
parms.Kie = parms.Kee.*gammaI;
parms.Kei = parms.Kee.*gamma;
parms.Kii = parms.Kee.*gamma;

parms.V_rest = 0;
parms.delay_s =  4.*rand(parms.EPopNum+parms.IPopNum,1)+1;;
parms.g = 2;


parms.V_th =20;
parms.tau_m = 20;
parms.V_reset = 10;
parms.t_ref = 1;


alpha = 0.5; %root K scaling
parms.J = (parms.V_th-parms.V_rest)./(parms.Kee.^alpha);

parms.LearningRate = 1e-2;
parms.LearningRate = 1e-2.*parms.J;%Learning rate is O(1/100) of synaptic weight
parms.TargetRate = [sort(exp(randn(parms.EPopNum,1)));nan(parms.IPopNum,1)]; %Target Rate for Excitatory cells (units of Hz)
parms.tauSTDP = 20;    %Time Constant for the STDP curve (Units of ms)





%% Feedforward parameters
% parms.N_FF = 1000;
% parms.K_FF = [250 750];
% parms.J_FF = [0.4 0.1];
%Feedforward parameters
parms.N_FF = 2400;
parms.K_FF = 600;
%Root K scaling for FF
%parms.J_FF = 0.1;
parms.J_FF = (parms.V_th-parms.V_rest)./(parms.K_FF.^0.5); %1/RootK scaling
%parms.J_FF = (parms.V_th-parms.V_rest)./(parms.K_FF); %1/K scaling
%% Ring input

[theta,T] = OUNoise(100000.^-1,4*pi,TimeParams.SimTime,TimeParams.dt,TimeParams.dt.*5,1);
neuronIDX = 2.*pi.*[1:parms.N_FF]./parms.N_FF;

meanrate = 10;

k = 10;
bess = besseli(0,4);
tol = 0.25;
expmean = mean(exp(k.*cos(neuronIDX-theta(1)))./(2.*pi.*bess));

%Rate of the input cells given theta
R_th = @(th) meanrate.*exp(k.*cos(neuronIDX-th))./(2.*pi.*bess)./expmean;
%Rate of the input cells given time (and timeseries theta)
R = @(t) R_th(theta(abs(T-t)<tol));
parms.ex_rate = R;
%%

figure
imagesc(T,neuronIDX,R(T)')
colorbar
hold on
plot(T,mod(theta,2*pi),'k.','markersize',1)
%%
%v_th = th/(C_e.*J.*tau);

tic 
[SimValues] = Run_LIF_iSTDP(parms,TimeParams,'showprogress',true,...
    'cellout',true,'save_dt',1000,'estrate',15);
toc
%%
save('RingSim')
%% Get Sorting by max(Input(theta))


th = linspace(0,2*pi,100);
inputtuning = zeros(length(th),parms.EPopNum+parms.IPopNum);
W = SimValues.FF_mat;

for tt = 1:length(th)
    inputtuning(tt,:) = R_th(th(tt))*W;
end
[~,peakth] = max(inputtuning);
peakth = th(peakth);
[~,sortpeak] = sort(peakth(1:parms.EPopNum));
sortpeak = [sortpeak ,(parms.EPopNum+[1:parms.IPopNum])];
%%
figure
imagesc(inputtuning(:,sortpeak))

%%
netname='ring';
overlay_HD = parms.EPopNum.*mod(theta,2*pi)./(2*pi);
overlay_HD(abs(diff(overlay_HD))>100) = nan;
PlotSimRaster(SimValues,TimeParams.SimTime-[3000 0],...
    'cellsort',sortpeak,'overlay',[T overlay_HD])
NiceSave('iSTDPRaster_late',pwd,netname)

PlotSimRaster(SimValues,[0000 5000],...
    'cellsort',sortpeak,'overlay',[T overlay_HD])
NiceSave('iSTDPRaster_early',pwd,netname)

% PlotSimRaster(SimValues,[0000 5000],...
%     'cellsort',sortpeak,'overlay',[T overlay_HD])
% NiceSave('iSTDPRaster_early',pwd,netname)







