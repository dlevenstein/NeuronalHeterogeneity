function iSTDPInputFluctandGSparms(savepath)

%%
%savepath = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/Modeling/Simulation_Data/iSTDPInputFluctandGSparms/';
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
TimeParams.SimTime = 100000;
%TimeParams.SimTime = 1000;

%Poisson Rate (add to Brunel sim)
%g = 5;
%g = 4; %Initial strength of Inhibitoon (relative to excitation)

clear parms


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
%parms.delay_s = 1.2;
%parms.delay_s = 1.8.*rand(parms.EPopNum+parms.IPopNum,parms.EPopNum+parms.IPopNum)+1.2;
parms.delay_s = 4.*rand(parms.EPopNum+parms.IPopNum,1)+1; %grid later
parms.g = 2;

parms.V_th =20;
parms.tau_m = 20;
parms.V_reset = 10;
parms.t_ref = 1;

%Feedforward parameters
parms.N_FF = 2400;
parms.K_FF = 600;
parms.J_FF = (parms.V_th-parms.V_rest)./(parms.K_FF.^0.5); %1/RootK scaling




alpha = 0.5; %Root K scaling
parms.J = (parms.V_th-parms.V_rest)./(parms.Kee.^alpha);

parms.LearningRate = 1e-2;
parms.TargetRate = [sort(exp(randn(parms.EPopNum,1)));nan(parms.IPopNum,1)]; %Target Rate for Excitatory cells (units of Hz)
parms.tauSTDP = 20;    %Time Constant for the STDP curve (Units of ms)

%%
numInputs = 1;
%v_th = th/(C_e.*J.*tau);
v_th = 1000*(parms.V_th-parms.V_rest)/(parms.K_FF.*parms.J_FF.*parms.tau_m);
v_rel = 3; %3 times the threshold rate
meanrate = v_rel.*v_th;

%%

%Initial Training with no fluctuating inputs
%meanrate = 10;
parms.ex_rate = meanrate;
[SimValues] = Run_LIF_iSTDP(parms,TimeParams,'showprogress',true,...
    'cellout',true,'save_dt',1000,'plotEIweight',true);
NiceSave('TrainingFigure',savepath,[])

%% Then add fluctuations


TimeParams.SimTime = 100000;
%TimeParams.SimTime = 100;

duration = TimeParams.SimTime;
dt = 0.1;
save_dt = 1;
numsignals = 1;

sigmas = [logspace(-1,0.5,5)].*v_th;
thetas = 1./logspace(1,3,5);
[Si,Th] = meshgrid(sigmas,thetas);
numsims = length(sigmas).*length(thetas);
%%
parfor ss = 1:numsims
%for ss = 1:length(sigmas)
    %sig = S(ss)
    
    [ X,T ] = OUNoise(Th(ss),Si(ss),duration,dt,save_dt,numsignals);
    ex_rate_fun{ss} = @(t) interp1(T,X,t,'nearest')+meanrate;
    
    
    loopparms = parms;
    loopparms.ex_rate = ex_rate_fun{ss};
    [SimValues_fluct{ss}] = Run_LIF_iSTDP(loopparms,TimeParams,'showprogress','parloop',...
        'cellout',true,'save_dt',2,'J_mat',SimValues.WeightMat);
    NiceSave('NoiseFigure',savepath,['S',num2str(round(Si(ss),1)),'T',num2str(round(1./Th(ss),1))])
end
%%
clear pc
savefilename = fullfile(savepath,'simresults.mat');

save(savefilename,'-v7.3')
disp('mat file saved')


%%
% sigmacolors = crameri('imola',length(sigmas)+1);
% sigmacolors(end,:) = [];
% %%
% 
% for ss = 1:length(sigmas)
% PlotSimRaster(SimValues_fluct{ss},TimeParams.SimTime-[2000 0])
% subplot(4,1,4)
% hold off
%     plot(SimValues_fluct{ss}.t,ex_rate_fun{ss}(SimValues_fluct{ss}.t),'linewidth',2,'color',sigmacolors(ss,:))
%     xlim(TimeParams.SimTime-[2000 0])
% end
% 
%     %%
%     for ss = 1:length(sigmas) 
% %clear spikes
% spikes(ss).times = cellfun(@(X) X./1000,SimValues_fluct{ss}.spikesbycell,'UniformOutput',false);
% spikes(ss).UID = 1:length(SimValues_fluct{ss}.spikesbycell);
% CellClass = cell(1,length(spikes(ss).times));
% CellClass(SimValues_fluct{ss}.EcellIDX) = {'E'};
% CellClass(SimValues_fluct{ss}.IcellIDX) = {'I'};
% %timewindows.initialization = [0 20];
% %timewindows.equib = [0 TimeParams.SimTime];
% ISIstats(ss) = bz_ISIStats(spikes(ss),'showfig',true,'cellclass',CellClass);
%     end
% 
%     %% CCG
% for ss = 1:length(sigmas)
%     [popCCG(ss)] = PopCCG(spikes(ss),'showfig',true,'cellclass',CellClass);
% end
% %%
% AScost = 0.0001; %0.13 for data
% for ss = 1:length(sigmas)
%     logISIhist = ISIstats(ss).ISIhist.ALL.log;
%     usecells = randsample(SimValues_fluct{ss}.EcellIDX([100:end]),250);
%     logISIhist = logISIhist(usecells,:)';
%     logtimebins = ISIstats(ss).ISIhist.logbins;
%     logISIhist = logISIhist./mode(diff(logtimebins));
% 
% GammaFit(ss) = bz_FitISISharedGammaModes(logISIhist,logtimebins,...
%     'numAS',1,...
%     'figfolder',savepath,'basePath',savepath,'figname',['sigma',num2str(sigmas(ss))],...
%     'AScost_lambda',AScost,'AScost_p',1/2,'ASguess',true,'MScost',3);
% 
% end
% 
% 
% %%
% classnames = {'E','I'};
% figure
% for ss = 1:length(sigmas)
%     
%         PlotSimRaster(SimValues_fluct{ss},TimeParams.SimTime-[2000 0])
%         subplot(2,1,1)
%             plot(SimValues_fluct{ss}.t,bz_NormToRange(ex_rate_fun{ss}(SimValues_fluct{ss}.t),1,[-50 100]),'linewidth',2,'color',sigmacolors(ss,:))
%         
%         subplot(4,3,10)
%             plot(-GammaFit(ss).sharedfit.GSlogrates,log10(GammaFit(ss).sharedfit.GSCVs),'.')
%             hold on
%             plot(ISIstats(ss).ISIhist.logbins([1 end]),[0 0],'k--')
%             plot(-GammaFit(ss).sharedfit.ASlogrates,log10(GammaFit(ss).sharedfit.ASCVs),'o')
%             xlim(ISIstats(ss).ISIhist.logbins([1 end]))
%             box off
%             
%             
%         subplot(4,3,7)
%             imagesc(ISIstats(ss).ISIhist.logbins,SimValues.EcellIDX,...
%                 ISIstats(ss).ISIhist.ALL.log(ISIstats(ss).sorts.ALL.ratebyclass,:))
%             %colorbar
%             caxis([0 0.15])
%             
% %         for tt = 1:2
% %             subplot(3,3,tt+7)
% %                 imagesc(popCCG(ss).t_ccg,[0 length(spikes(ss).times)],log10(popCCG(ss).cells.(classnames{tt}))')
% %                 hold on
% % 
% %                 plot(popCCG(ss).t_ccg,bz_NormToRange(-popCCG(ss).pop.(classnames{tt})(:,1),[1 length(SimValues_fluct{ss}.EcellIDX)-1]),'linewidth',2);%,'color',cellcolor{tt})
% % 
% %                 for tt2 = 2:numclasses
% %                     plot(xlim(gca),sum(inclasscells{tt2-1}).*[1 1],'w')
% %                     plot(popCCG(ss).t_ccg,bz_NormToRange(-popccg(:,tt2,tt),sum(inclasscells{tt2-1})+[1 sum(inclasscells{tt2})-1]),'linewidth',2);%,'color',cellcolor{tt})
% %                 end
% %                 %plot(t_ccg,bz_NormToRange(-popccg(:,tt,2),sum(CellClass.pE)+[1 sum(CellClass.pI)]),'color',cellcolor{tt})
% %                 %plot(xlim(gca),sum(CellClass.pE).*[1 1],'w')
% %                 colorbar
% %                 title(classnames{tt})
% %                 xlabel('t lag')
% %         end
%             
% 	NiceSave('RasterEtc',savepath,[netname,'_Sigma',num2str(round(sigmas(ss)))],'includeDate',true)
% 
% end
% 
% 
% 
% %%
% figure %GS_CV values
% subplot(3,3,1)
% BoxAndScatterPlot({log10(GammaFit(1).sharedfit.GSCVs),...
%     log10(GammaFit(2).sharedfit.GSCVs),log10(GammaFit(3).sharedfit.GSCVs),...
%     log10(GammaFit(4).sharedfit.GSCVs)},'labels',round(sigmas,1),'colors',sigmacolors)
% hold on
% ylim([-0.5 0.5])
% ylabel('Ground State CV');xlabel('Fluctuation Amplitude (Hz)')
% plot(xlim(gca),[0 0],'k--')
% LogScale('y',10)
% box off
% 
% 
% for ss = 1:length(sigmas)
%     subplot(6,3,2)
%     hold on
%     plot(popCCG(ss).t_ccg,popCCG(ss).pop.E(:,1),'color',sigmacolors(ss,:),'linewidth',1)
%     plot([0 0],ylim(gca),'k')
%     xlim([-0.1 0.1])
%     ylabel('E Rate');xlabel('t lag (E Spikes)')
%     
%     subplot(6,3,5)
%     hold on
%     plot(popCCG(ss).t_ccg,popCCG(ss).pop.I(:,1),'color',sigmacolors(ss,:),'linewidth',1)
%     plot([0 0],ylim(gca),'k')
%     xlim([-0.1 0.1])
%     ylabel('I Rate');xlabel('t lag (E Spikes)')
%     
%     subplot(6,3,3)
%     hold on
%     plot(popCCG(ss).t_ccg,popCCG(ss).pop.E(:,2),'color',sigmacolors(ss,:),'linewidth',1)
%     plot([0 0],ylim(gca),'r')
%     xlim([-0.1 0.1])
%     ylabel('E Rate');xlabel('t lag (I Spikes)')
%     
%     subplot(6,3,6)
%     hold on
%     plot(popCCG(ss).t_ccg,popCCG(ss).pop.I(:,2),'color',sigmacolors(ss,:),'linewidth',1)
%     plot([0 0],ylim(gca),'r')
%     xlim([-0.1 0.1])
%     ylabel('I Rate');xlabel('t lag (I Spikes)')
% end
% 
% NiceSave('CV_CCG',savepath,netname)
% 

