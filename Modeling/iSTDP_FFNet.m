function iSTDPRingnet(savepath)

%%
savepath = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/Modeling/Simulation_Data/iSTDP_FFNet/';
if ~exist(savepath,'dir')
    mkdir(savepath)
end
%%
display(['Will save to ',savepath])

clusterpar = false;
if clusterpar
    pc = parcluster('local');
    % % store temporary files in the 'scratch' drive on the cluster, labeled by job ID
    pc.JobStorageLocation = strcat(getenv('SCRATCH'), '/', getenv('SLURM_JOB_ID'));
    % % enable MATLAB to utilize the multiple cores allocated in the job script
    % % SLURM_NTASKS_PER_NODE is a variable set in the job script by the flag --tasks-per-node
    % % we use SLURM_NTASKS_PER_NODE - 1, because one of these tasks is the original MATLAB script itself
    parpool(pc, str2num(getenv('SLURM_NTASKS_PER_NODE'))-1);
end
%%
TimeParams.dt = 0.1;
TimeParams.SimTime = 100000;
%TimeParams.SimTime = 1000;

%Poisson Rate (add to Brunel sim)
%g = 5;
%g = 4; %Initial strength of Inhibitoon (relative to excitation)

clear parms


parms.EPopNum = 100;
parms.IPopNum = 100;
parms.u_0 = 0;

%Conectivity: In degree
%gamma = 0.5; %initally 0.25 to match 4x less inhibitory cells
%gammaI = 2; %relative E->I connectivity
parms.Kee = 0;
parms.Kie = 0;
parms.Kei = 50;
parms.Kii = 50;


parms.V_rest = 0;
%parms.delay_s = 1.2;
%parms.delay_s = 1.8.*rand(parms.EPopNum+parms.IPopNum,parms.EPopNum+parms.IPopNum)+1.2;
parms.delay_s = 4.*rand(parms.EPopNum+parms.IPopNum,1)+1; %grid later
parms.g = 1;

parms.V_th =20;
parms.tau_m = 20;
parms.V_reset = 10;
parms.t_ref = 1;

%Feedforward parameters
parms.N_FF = 1200;
parms.K_FF = 300;
parms.J_FF = (parms.V_th-parms.V_rest)./(parms.K_FF.^0.5); %1/RootK scaling




alpha = 1; %Root K scaling
parms.J = (parms.V_th-parms.V_rest)./(parms.Kee.^alpha);
parms.J = (parms.V_th-parms.V_rest)./(parms.Kei.^alpha); 

parms.LearningRate = parms.g.*parms.J.*2e-2; %O(1/100) initial synaptic weight
%parms.TargetRate = [sort(exp(randn(parms.EPopNum,1)));nan(parms.IPopNum,1)]; %Target Rate for Excitatory cells (units of Hz)
parms.TargetRate = [sort(exp(randn(parms.EPopNum,1)));20.*ones(parms.IPopNum,1)]; %Target Rate for Excitatory cells (units of Hz)
parms.tauSTDP = 20;    %Time Constant for the STDP curve (Units of ms)

%%
%v_th = th/(C_e.*J.*tau);
v_th = 1000*(parms.V_th-parms.V_rest)/(parms.K_FF.*parms.J_FF.*parms.tau_m);
v_rel = 2; %3 times the threshold rate
meanrate = v_rel.*v_th;

%%

% %Initial Training with no fluctuating inputs
% %meanrate = 10;
% 
% %Train with fluctuating rate
% meanrate = v_th.*2;
% duration = TimeParams.SimTime;
% OU_simdt = 0.1;
% OU_savedt = 1;
% numsignals = 1;
% 
% theta = 1./5000; %5s (1000ms) timescale
% sigma = 1.*v_th;
% %sigma = 1;
% 
% %disp('Making OU noise...')
% [ X,T ] = OUNoise(theta,sigma,duration,OU_simdt,OU_savedt,numsignals);
% parms.ex_rate = @(t) interp1(T,X,t,'nearest')+meanrate;
    

%%
parms.ex_rate = meanrate;
[SimValues] = Run_LIF_iSTDP(parms,TimeParams,'showprogress',true,...
    'cellout',true,'save_dt',100,'plotEIweight',true);
%NiceSave('TrainingFigure',savepath,[])

%%
figure
plot(SimValues.t,SimValues.V(1,:),'k')
end