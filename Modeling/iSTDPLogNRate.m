function iSTDPLogNRate(repopath,whichnet)

repopath = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity'; 
%whichnet = 'Uniform';
addpath(genpath(repopath))

modelrepopath = '/Users/dlevenstein/Project Repos/SOSpikingModel';
addpath(genpath(modelrepopath))


figfolder = fullfile(repopath,'Modeling','Figures','iSTDP');
savedatafolder = fullfile(repopath,'Modeling','Simulation_Data','LogN/');
%% Set Cluster Number
%numClusters = 5; %This is if you have 5 clusters, so cluster_number inputs are 1-5 

%% Load Scripts Here 

netfolder = fullfile(repopath,'Modeling','StuffForDan','ExperimentRateLongTrained/');


        load([netfolder 'Lognormal_m_1_s_10_EE_LognormalRates_Noise_10ms_50pA_K_IE_250_Spikes.mat']);
        %savefilename = 'AdaptationVCurrentSpikes.mat';

%% Network Parameters
PopParams = SimValues.PopParams;
% clear PopParams
% 
% PopParams.EPopNum = 2000;
% PopParams.IPopNum = 500;
% 
% %Noise Properties 
% PopParams.I_e  = 0;         %External input
% PopParams.sigma = 0;        %noise magnitude: variance
% PopParams.theta = 0.1;      %noise time scale (1/ms)
% 
%Neuron properties
PopParams.E_L     = [-65 -67];          %rev potential: leak (mV)
PopParams.g_L     = [182/18 119/8];     %leak conductance (nS)
PopParams.C       = [182 119];          %capacitance (pF)
PopParams.V_th    = [-45 -47];          %spike threshold (mV)
PopParams.V_reset = [-55 -55];          %reset potential (mV)
PopParams.t_ref   = 0.5;                %refractory period (ms)

%Synaptic Properties
PopParams.E_e     = 0;              %rev potential: E (mV)
PopParams.E_i     = -80;            %rev potential: I (mV)
PopParams.tau_s   = [5 5];          %synaptic decay timescale (1/ms)
% 
%Adaptation Properties (loop through a or b in loop module)
PopParams.E_w     = -70;            %rev potential: adaptation (mV)
PopParams.a       = 0;              %adaptation decay timescale (1/ms)
PopParams.b       = 0;              %adaptation activation rate (1/ms)
PopParams.tau_w   = 500;            %adaptation decay (ms)
PopParams.gwnorm  = 0;              %magnitude of adaptation

PopParams.t_syn = rand(PopParams.EPopNum+PopParams.IPopNum,1)*0.4+0.05;                %Synaptic Delay
% 
 PopParams.LearningRate = 0;
% PopParams.TargetRateI = nan; %Target E rate nan Hz (Turns off plasticity)
% PopParams.TargetRateE = nan; %Target I rate nan Hz (Turns off plasticity)
% PopParams.tauSTDP = 20;

%% Adaptation and Noise Parameters 

% PopParamsAnalysis              = PopParams;
% PopParamsAnalysis.LearningRate = 0;                      %Learning rate
% PopParamsAnalysis.tau_w        = 500;    %300                %adaptation decay (ms)
% PopParamsAnalysis.sigma        = 100;     %100               %Noise variance (pA) (Set to Covariance Matrix to add covariance
 PopParams.W            = SimValues.WeightMat(:,:,end); %Synaptic Weights
% PopParamsAnalysis.gwnorm       = 1;                      %Adaptation norm
% PopParamsAnalysis.t_syn        = rand(PopParams.EPopNum+PopParams.IPopNum,1)*0.4+0.05;                      %Synaptic Delay (ms)

%% Note Adaptation Equation
%dwdt =  (- w + a.*(V - E_w))./tau_w; (Line 529)
%w(spikeneurons) = w(spikeneurons) + b(spikeneurons); at time of spike (Line 559)
%g_w = gwnorm.*w; (Line 648)
%% Starting Values

%PopParamsAnalysis.V0 = PopParams_in.V0;
%PopParamsAnalysis.w0 = PopParams_in.w0;
%PopParamsAnalysis.s0 = PopParams_in.s0;


PopParams.p0spike = 0.02; %start ON

TimeParams.dt      = 0.05;
TimeParams.SimTime = 5000;


%%

SimValues_new = AdLIFfunction_iSTDP(PopParams,TimeParams,'cellout',true,...
    'showprogress',true,'showfig',true,'save_weights',TimeParams.SimTime,...
    'save_dt',TimeParams.SimTime,'useGPU',false,'defaultNeuronParams',false);%,...
    %'recordInterval',[(0:RecordTime:RecordTime) + (TimeParams.SimTime - RecordTime)]');
%% Set up for parallel in cluster
% pc = parcluster('local');
% % store temporary files in the 'scratch' drive on the cluster, labeled by job ID
% pc.JobStorageLocation = strcat(getenv('SCRATCH'), '/', getenv('SLURM_JOB_ID'));
% % enable MATLAB to utilize the multiple cores allocated in the job script
% % SLURM_NTASKS_PER_NODE is a variable set in the job script by the flag --tasks-per-node
% % we use SLURM_NTASKS_PER_NODE - 1, because one of these tasks is the original MATLAB script itself
% parpool(pc, str2num(getenv('SLURM_NTASKS_PER_NODE'))-1);
%%
spikes_In.times = cellfun(@(X) X./1000,spikesbycell,'UniformOutput',false);
spikes_In.UID = 1:(PopParams.EPopNum+PopParams.IPopNum);
celltype = cell(size(spikes_In.UID));
celltype(SimValues.EcellIDX) = {'E'};
celltype(SimValues.IcellIDX) = {'I'};
ISIStats = bz_ISIStats(spikes_In,'cellclass',celltype,'showfig',true,...
    'savecellinfo',true,'basePath',savedatafolder,'figfolder',figfolder,...
    'forceRedetect',true);

%%
baseName = 'LogNRate'
histcolors = flipud(gray);
figure
subplot(3,3,1)

imagesc(ISIStats.ISIhist.logbins,[1 length(ISIStats.sorts.ALL.rateE)],...
    ISIStats.ISIhist.ALL.log(ISIStats.sorts.ALL.rateE,:))
hold on
%plot(ISIstats.ISIhist.logbins,bz_NormToRange(ISIstats.meanISIhist.Sup.dark.,0.3),'k','linewidth',1)
plot(log10(1./ISIStats.summstats.ALL.meanrate(ISIStats.sorts.ALL.rateE)),...
    [1:length(ISIStats.sorts.ALL.rateE)],'k.','markersize',2)
colormap(gca,histcolors)
caxis([0 0.15])
LogScale('x',10,'exp',true)
ylabel({'Sup Cells','Sort by Rate'})
xlabel('ISI (s)')

NiceSave('ISIdists',figfolder,baseName)
