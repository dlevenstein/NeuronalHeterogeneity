function iSTDPRecurrence(savepath)

%%
%savepath = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/Modeling/Simulation_Data/Recurrence';
if ~exist(savepath,'dir')
    mkdir(savepath)
end
%%
display(['Will save to ',savepath])

pc = parcluster('local');

    % sto
% % store temporary files in the 'scratch' drive on the cluster, labeled by job ID
pc.JobStorageLocation = strcat(getenv('SCRATCH'), '/', getenv('SLURM_JOB_ID'));
% % enable MATLAB to utilize the multiple cores allocated in the job script
% % SLURM_NTASKS_PER_NODE is a variable set in the job script by the flag --tasks-per-node
% % we use SLURM_NTASKS_PER_NODE - 1, because one of these tasks is the original MATLAB script itself
parpool(pc, str2num(getenv('SLURM_NTASKS_PER_NODE'))-1);
%%
TimeParams.dt = 0.1;

clear parms

parms.EPopNum = 1000;
parms.IPopNum = 250;
parms.u_0 = 0;

parms.V_rest = 0;
parms.delay_s = 8.9.*rand(parms.EPopNum+parms.IPopNum,1)+1.1; %grid later
parms.g = 6; %Initial strength of Inhibitoon (relative to excitation)

parms.V_th =20;
parms.tau_m = 20;
parms.V_reset = 10;
parms.t_ref = 1;

%Feedforward parameters
parms.N_FF = 1000;
parms.K_FF = 250;
%Root K scaling for FF
%parms.J_FF = 0.1;
parms.J_FF = (parms.V_th-parms.V_rest)./(parms.K_FF.^0.5);

%Conectivity: In degree
gamma = 0.25; %4x less inhibitory cells
parms.Kee = 250;
parms.Kie = parms.Kee;
parms.Kei = parms.Kee.*gamma;
parms.Kii = parms.Kee.*gamma;


parms.LearningRate = 1e-2;
parms.TargetRate = [sort(exp(randn(parms.EPopNum,1)));nan(parms.IPopNum,1)]; %Target Rate for Excitatory cells (units of Hz)
parms.tauSTDP = 20;    %Time Constant for the STDP curve (Units of ms)


%%
numJs = 11;
alphas = linspace(0,1,numJs)
%alphas = 0.5;
Js = (parms.V_th-parms.V_rest)./(parms.Kee.^alphas)
%Js = logspace(log10(rootKscale),1,numJs);

%
%%
numInputs = 10;
%v_th = th/(C_e.*J.*tau);
v_th = 1000*(parms.V_th-parms.V_rest)/(parms.K_FF.*parms.J_FF.*parms.tau_m);
inputrates = logspace(-0.5,1,numInputs).*v_th;


%%
parfor jj = 1:numJs
    
    TimeParams_Jloop = TimeParams;
    TimeParams_Jloop.SimTime = 120000;
    TimeParams_Jloop.SimTime = 100;

    parms_Jloop = parms;
    parms_Jloop.J = Js(jj);

    %Train with fluctuating rate
    meanrate = v_th.*3;
    duration = TimeParams_Jloop.SimTime;
    OU_simdt = 0.1;
    OU_savedt = 1;
    numsignals = 1;

    theta = 1./2000; %1s (1000ms) timescale
    sigma = v_th;

    %disp('Making OU noise...')
    [ X,T ] = OUNoise(theta,sigma,duration,OU_simdt,OU_savedt,numsignals);
    parms_Jloop.ex_rate = @(t) interp1(T,X,t,'nearest')+meanrate;
    %disp('DONE!')
    %
    %%

    [tempj] = Run_LIF_iSTDP(parms_Jloop,TimeParams_Jloop,'showprogress',false,...
        'cellout',true,'save_dt',1000,'estrate',50);
    SimValues_train{jj} = tempj;
    
    NiceSave('TrainingFigure',savepath,['alpha',num2str(round(alphas(jj),1))])

    disp('J sim done')
    %% Different inputs
    TimeParams_Iloop = TimeParams;
    TimeParams_Iloop.SimTime = 30000;
    TimeParams_Iloop.SimTime = 30;
    for rr = 1:numInputs
        
        parms_Iloop = parms;
        parms_Iloop.ex_rate = inputrates(rr);
        %tic 
        disp('Starting Input Sim')
        [temp] = Run_LIF_iSTDP(parms_Iloop,TimeParams_Iloop,'showprogress',true,...
            'cellout',true,'save_dt',2,'J_mat',SimValues_train{jj}.WeightMat);
        %toc
        disp('Input sim done')
        NiceSave('SimFig',savepath,['alpha',num2str(round(alphas(jj),1)),'input',num2str(round(inputrates(rr),1))])
        disp('savedfig')
        SimValues_inputs{jj,rr} = temp;
        disp('tempassigned')
    end

end
%%
if ~exist(savepath,'dir')
    mkdir(savepath)
end
savefilename = fullfile(savepath,'simresults.mat');

save(savefilename,'-v7.3')

%% Plot Rasters


% 
% 
% %%
% 
% for ss = 1:length(inputrates) 
%     %clear spikes
%     spikes(ss).times = cellfun(@(X) X./1000,SimValues_inputs(ss).spikesbycell,'UniformOutput',false);
%     spikes(ss).UID = 1:length(SimValues_inputs(ss).spikesbycell);
%     CellClass = cell(1,length(spikes(ss).times));
%     CellClass(SimValues_inputs(ss).EcellIDX) = {'E'};
%     CellClass(SimValues_inputs(ss).IcellIDX) = {'I'};
%     %timewindows.initialization = [0 20];
%     %timewindows.equib = [0 TimeParams.SimTime];
%     ISIstats(ss) = bz_ISIStats(spikes(ss),'showfig',true,'cellclass',CellClass);
%     [popCCG(ss)] = PopCCG(spikes(ss),'showfig',true,'cellclass',CellClass);
%     
%     close all
% end
% 
% %%
% numratebins = 40;
% numrvoltbins = 100;
% FIStats.ratedist.bins = linspace(-1,1.5,numratebins);
% FIStats.voltagedist.bins = linspace(-20,20,numrvoltbins);
% FIStats.ratedist.pop = zeros(numratebins,length(inputrates));
% FIStats.ratedist.cell = zeros(numratebins,length(inputrates));
% FIStats.voltagedist.allcells = zeros(numrvoltbins,length(inputrates));
% FIStats.ISIdist.bins = ISIstats(1).ISIhist.logbins;
% FIStats.ISIdist.meanE = zeros(length(FIStats.ISIdist.bins),length(inputrates));
% 
% FIStats.popCCG.bins = popCCG(1).t_ccg;
% 
% for ss = 1:length(inputrates) 
%     FIStats.ratedist.cell(:,ss) = hist(log10(ISIstats(ss).summstats.ALL.meanrate(SimValues_inputs(ss).EcellIDX)),FIStats.ratedist.bins);
%     FIStats.ratedist.Icell(:,ss) = hist(log10(ISIstats(ss).summstats.ALL.meanrate(SimValues_inputs(ss).IcellIDX)),FIStats.ratedist.bins);
%     
%     allvoltage = SimValues_inputs(ss).V(SimValues_inputs(ss).EcellIDX,:);
%     FIStats.voltagedist.allcells(:,ss) = hist(allvoltage(:),FIStats.voltagedist.bins);
%     
%     spikemat = PlotSimRaster(SimValues_inputs(ss),[0000 2000]);
%     FIStats.ratedist.pop(:,ss) = hist(log10(spikemat.E.poprate),FIStats.ratedist.bins); 
%     
%     FIStats.ISIdist.meanE(:,ss) = mean(ISIstats(ss).ISIhist.ALL.log(SimValues_inputs(ss).EcellIDX,:),1);
%     
%     FIStats.popCCG.E(:,ss) = popCCG(ss).pop.E(:,1);
%     FIStats.popCCG.I(:,ss) = popCCG(ss).pop.I(:,2);
% end
% 
% %%
% figure
% subplot(3,3,1)
% imagesc(log10(inputrates),FIStats.ratedist.bins,FIStats.ratedist.cell)
% axis xy
% LogScale('xy',10)
% 
% subplot(3,3,2)
% imagesc(log10(inputrates),FIStats.ratedist.bins,FIStats.ratedist.Icell)
% axis xy
% LogScale('xy',10)
% 
% 
% subplot(3,3,3)
% imagesc(log10(inputrates),FIStats.voltagedist.bins,FIStats.voltagedist.allcells)
% axis xy
% LogScale('x',10)
% 
% subplot(3,3,4)
% imagesc(log10(inputrates),FIStats.ratedist.bins,FIStats.ratedist.pop)
% axis xy
% LogScale('xy',10)
% 
% 
% subplot(3,3,5)
% imagesc(log10(inputrates),FIStats.ISIdist.bins,FIStats.ISIdist.meanE)
% axis xy
% LogScale('xy',10)
% 
% subplot(3,3,6)
% imagesc(log10(inputrates),FIStats.popCCG.bins,FIStats.popCCG.I)
% axis xy
% 
% subplot(3,3,7)
% imagesc(log10(inputrates),FIStats.popCCG.bins,FIStats.popCCG.E)
% axis xy
% %LogScale('xy',10)
% NiceSave('FIStats',pwd,'HiJ')
% 
% %%
% save('-v7.3')
% 
% %%
% 
% PlotSimRaster(SimValues,TimeParams.SimTime-[400 0])
% NiceSave('iSTDPRaster',pwd,netname)
% %% Save/load
% filename = fullfile(savepath,['TrainedNet_',netname]);
% %save(filename,'SimValues','parms','TimeParams','netname')
% load(filename)
% 
% %% Effective Indegree
% K_net = sum(SimValues.WeightMat,2);
% K_E = sum(SimValues.WeightMat(:,SimValues.EcellIDX),2);
% K_I = sum(SimValues.WeightMat(:,SimValues.IcellIDX),2);
% K_EI = K_E./K_I;
% 
% figure
% subplot(2,2,1)
% plot(log10(parms.TargetRate),K_net,'.')
% subplot(2,2,2)
% plot(log10(parms.TargetRate),K_EI,'.')
% subplot(2,2,3)
% plot(log10(parms.TargetRate),K_E,'.')
% subplot(2,2,4)
% plot(log10(parms.TargetRate),K_I,'.')
% %%
% PlotSimRaster(SimValues,TimeParams.SimTime-[2000 0])
% NiceSave('iSTDPRaster_late',pwd,netname)
% 
% PlotSimRaster(SimValues,[0000 2000])
% NiceSave('iSTDPRaster_early',pwd,netname)
% %%
% figure
% subplot(2,2,1)
% imagesc(SimValues.WeightMat_initial)
% caxis([-2 0.2])
% crameri('vik','pivot',0)
% colorbar
% subplot(2,2,2)
% imagesc(SimValues.WeightMat)
% colorbar
% 
% caxis([-2 0.2])
% crameri('vik','pivot',0)
% 
% %% ISI 
% clear spikes
% spikes.times = cellfun(@(X) X./1000,SimValues.spikesbycell,'UniformOutput',false);
% spikes.UID = 1:length(SimValues.spikesbycell);
% CellClass = cell(1,length(spikes.times));
% CellClass(SimValues.EcellIDX) = {'E'};
% CellClass(SimValues.IcellIDX) = {'I'};
% timewindows.initialization = [0 10];
% timewindows.equib = [10 TimeParams.SimTime./1000];
% ISIstats = bz_ISIStats(spikes,'ints',timewindows,'showfig',true,'cellclass',CellClass);
% 
% %%
% NiceSave('ISIStats',pwd,netname)
% %%
% figure
% subplot(2,2,1)
% plot(log10(parms.TargetRate),log10(ISIstats.summstats.equib.meanrate),'k.','markersize',2)
% axis tight
% box off
% hold on
% xlabel('Target Rate');ylabel('Simulated Rate')
% UnityLine
% LogScale('xy',10)
% NiceSave('RateTarget',pwd,netname)