function iSTDPandGSparms(savepath)

%%
savepath = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/Modeling/Simulation_Data/iSTDPInputFluctandGSparms/';

%%
TimeParams.dt = 0.1;
TimeParams.SimTime = 75000;
%TimeParams.SimTime = 100;

%Poisson Rate (add to Brunel sim)
g = 5;
%g = 4; %Initial strength of Inhibitoon (relative to excitation)

clear parms


parms.EPopNum = 1000;
parms.IPopNum = 250;
parms.u_0 = 0;

parms.V_rest = 0;
%parms.delay_s = 1.2;
%parms.delay_s = 1.8.*rand(parms.EPopNum+parms.IPopNum,parms.EPopNum+parms.IPopNum)+1.2;
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



netname = 'weaklybalanced';
parms.J = 0.5;
% switch netname
%     case 'weaklybalanced'
%         parms.J = 0.4;
%         %parms.ex_rate = 20;
%         %parms.ex_rate = 40; %Not same as other script...
%     case 'strongrecurrent'
%         parms.J = 1.5; 
%         parms.ex_rate = 20;
%     case 'CA1like'
%         parms.J = 0.4;
%         parms.ex_rate = 30; 
%         parms.Kee = 0; 
% end


parms.LearningRate = 0.5e-2;
parms.TargetRate = [sort(exp(randn(parms.EPopNum,1)));nan(parms.IPopNum,1)]; %Target Rate for Excitatory cells (units of Hz)
parms.tauSTDP = 20;    %Time Constant for the STDP curve (Units of ms)


%Initial Training with no fluctuating inputs
meanrate = 10;
parms.ex_rate = meanrate;
[SimValues] = Run_LIF_iSTDP(parms,TimeParams,'showprogress',true,...
    'cellout',true,'save_dt',1000);

%% Then add fluctuations


TimeParams.SimTime = 150000;

theta = 1./100; %100ms timescale
duration = TimeParams.SimTime;
dt = 0.1;
save_dt = 1;
numsignals = 1;

sigmas = [logspace(0,1,4)];
for ss = 1:length(sigmas)
    [ X,T ] = OUNoise(theta,sigmas(ss),duration,dt,save_dt,numsignals);
    
    ex_rate_fun{ss} = @(t) interp1(T,X,t,'nearest')+meanrate;
end
%parfor (ss = 1:length(sigmas),2)
for ss = 1:length(sigmas)
    ss
    loopparms = parms;
    loopparms.ex_rate = ex_rate_fun{ss};
    [SimValues_fluct{ss}] = Run_LIF_iSTDP(loopparms,TimeParams,'showprogress',true,...
        'cellout',true,'save_dt',2,'J_mat',SimValues.WeightMat);

end
%%
sigmacolors = crameri('imola',length(sigmas)+1);
sigmacolors(end,:) = [];
%%

for ss = 1:length(sigmas)
PlotSimRaster(SimValues_fluct{ss},TimeParams.SimTime-[2000 0])
subplot(4,1,4)
hold off
    plot(SimValues_fluct{ss}.t,ex_rate_fun{ss}(SimValues_fluct{ss}.t),'linewidth',2,'color',sigmacolors(ss,:))
    xlim(TimeParams.SimTime-[2000 0])
end

    %%
    for ss = 1:length(sigmas) 
%clear spikes
spikes(ss).times = cellfun(@(X) X./1000,SimValues_fluct{ss}.spikesbycell,'UniformOutput',false);
spikes(ss).UID = 1:length(SimValues_fluct{ss}.spikesbycell);
CellClass = cell(1,length(spikes(ss).times));
CellClass(SimValues_fluct{ss}.EcellIDX) = {'E'};
CellClass(SimValues_fluct{ss}.IcellIDX) = {'I'};
%timewindows.initialization = [0 20];
%timewindows.equib = [0 TimeParams.SimTime];
ISIstats(ss) = bz_ISIStats(spikes(ss),'showfig',true,'cellclass',CellClass);
    end

    
%%
AScost = 0.0001; %0.13 for data
for ss = 1:length(sigmas)
    logISIhist = ISIstats(ss).ISIhist.ALL.log;
    usecells = randsample(SimValues_fluct{ss}.EcellIDX([100:end]),250);
    logISIhist = logISIhist(usecells,:)';
    logtimebins = ISIstats(ss).ISIhist.logbins;
    logISIhist = logISIhist./mode(diff(logtimebins));

GammaFit(ss) = bz_FitISISharedGammaModes(logISIhist,logtimebins,...
    'numAS',1,...
    'figfolder',savepath,'basePath',savepath,'figname',['sigma',num2str(sigmas(ss))],...
    'AScost_lambda',AScost,'AScost_p',1/2,'ASguess',true,'MScost',3);

end


%%
figure
for ss = 1:length(sigmas)
    
        PlotSimRaster(SimValues_fluct{ss},TimeParams.SimTime-[2000 0])
        subplot(2,1,1)
            plot(SimValues_fluct{ss}.t,bz_NormToRange(ex_rate_fun{ss}(SimValues_fluct{ss}.t),1,[-50 100]),'linewidth',2,'color',sigmacolors(ss,:))
        
        subplot(4,2,7)
            plot(-GammaFit(ss).sharedfit.GSlogrates,log10(GammaFit(ss).sharedfit.GSCVs),'.')
            hold on
            plot(ISIstats(ss).ISIhist.logbins([1 end]),[0 0],'k--')
            plot(-GammaFit(ss).sharedfit.ASlogrates,log10(GammaFit(ss).sharedfit.ASCVs),'o')
            xlim(ISIstats(ss).ISIhist.logbins([1 end]))
            
            
        subplot(4,2,5)
            imagesc(ISIstats(ss).ISIhist.logbins,SimValues.EcellIDX,...
                ISIstats(ss).ISIhist.ALL.log(ISIstats(ss).sorts.ALL.ratebyclass,:))
            %colorbar
            caxis([0 0.15])
	NiceSave('RasterEtc',savepath,[netname,'_Sigma',num2str(round(sigmas(ss)))],'includeDate',true)

end

%%
%ss = 1;
for ss = 1:length(sigmas)
    [popCCG(ss)] = PopCCG(spikes(ss),'showfig',true,'cellclass',CellClass);
end


%%
figure %GS_CV values
subplot(3,3,1)
BoxAndScatterPlot({log10(GammaFit(1).sharedfit.GSCVs),...
    log10(GammaFit(2).sharedfit.GSCVs),log10(GammaFit(3).sharedfit.GSCVs),...
    log10(GammaFit(4).sharedfit.GSCVs)},'labels',round(sigmas,1),'colors',sigmacolors)
hold on
ylim([-0.5 0.5])
ylabel('Ground State CV');xlabel('Fluctuation Amplitude (Hz)')
plot(xlim(gca),[0 0],'k--')
LogScale('y',10)
box off


for ss = 1:length(sigmas)
    subplot(6,3,2)
    hold on
    plot(popCCG(ss).t_ccg,popCCG(ss).pop.E(:,1),'color',sigmacolors(ss,:),'linewidth',1)
    plot([0 0],ylim(gca),'k')
    xlim([-0.1 0.1])
    ylabel('E Rate');xlabel('t lag (E Spikes)')
    
    subplot(6,3,5)
    hold on
    plot(popCCG(ss).t_ccg,popCCG(ss).pop.I(:,1),'color',sigmacolors(ss,:),'linewidth',1)
    plot([0 0],ylim(gca),'k')
    xlim([-0.1 0.1])
    ylabel('I Rate');xlabel('t lag (E Spikes)')
    
    subplot(6,3,3)
    hold on
    plot(popCCG(ss).t_ccg,popCCG(ss).pop.E(:,2),'color',sigmacolors(ss,:),'linewidth',1)
    plot([0 0],ylim(gca),'r')
    xlim([-0.1 0.1])
    ylabel('E Rate');xlabel('t lag (I Spikes)')
    
    subplot(6,3,6)
    hold on
    plot(popCCG(ss).t_ccg,popCCG(ss).pop.I(:,2),'color',sigmacolors(ss,:),'linewidth',1)
    plot([0 0],ylim(gca),'r')
    xlim([-0.1 0.1])
    ylabel('I Rate');xlabel('t lag (I Spikes)')
end

NiceSave('CV_CCG',savepath,netname)
%% Pop autocorrelation
ss = 4;

allspikes = SimValues_fluct{ss}.spikesbycell;
%allspikes{end+1} = SimValues_fluct{ss}.spikes(:,1); 

[allccg,t_ccg] = CCG(allspikes,[],'binSize',2,'duration',250,'norm','rate');  

%%
Espikes = ismember(SimValues_fluct{ss}.spikes(:,2),SimValues_fluct{ss}.EcellIDX);
Ispikes = ismember(SimValues_fluct{ss}.spikes(:,2),SimValues_fluct{ss}.IcellIDX);
popSpikes{1} = SimValues_fluct{ss}.spikes(Espikes,1);
popSpikes{2} = SimValues_fluct{ss}.spikes(Ispikes,1);
[popccg,t_ccg] = CCG(popSpikes,[],'binSize',2,'duration',250,'norm','rate'); 
%%
clear cellpopccg
clear cellspikes
%cellpopccg = struct('t');
numcellsccg = SimValues_fluct{ss}.PopParams.EPopNum+SimValues_fluct{ss}.PopParams.IPopNum;
numcellsccg = 10;
for cc = 1:numcellsccg
    bz_Counter(cc,numcellsccg,'Calculating PopCCG: Cell');
    thiscellspikes = ismember(SimValues_fluct{ss}.spikes(:,2),cc);
    Espikes = ismember(SimValues_fluct{ss}.spikes(:,2),SimValues_fluct{ss}.EcellIDX);
    Ispikes = ismember(SimValues_fluct{ss}.spikes(:,2),SimValues_fluct{ss}.IcellIDX);
    cellspikes = SimValues_fluct{ss}.spikes(thiscellspikes,1);
    otherEspikes = SimValues_fluct{ss}.spikes(Espikes&~thiscellspikes,1);
    otherIspikes = SimValues_fluct{ss}.spikes(Ispikes&~thiscellspikes,1);
    [temp,t] = CCG({cellspikes,otherEspikes,otherIspikes},[],'binSize',2,'duration',400,'norm','rate'); 
    Epop(:,cc) = temp(:,1,2).*1000./sum(~ismember(SimValues_fluct{ss}.EcellIDX,cc));
    Ipop(:,cc) = temp(:,1,3).*1000./sum(~ismember(SimValues_fluct{ss}.IcellIDX,cc));
end

%%
figure
subplot(2,2,1)
imagesc(log10(cellpopccg.Epop)')
colorbar


subplot(2,2,2)
imagesc(log10(cellpopccg.Ipop'))
%%
% figure
% plot(temp(:,1,2),'k')
% hold on
% plot(temp(:,1,3),'r')
%%
figure
plot(t_ccg,popccg(:,1,1),'k')
hold on
plot(t_ccg,popccg(:,2,2),'r')
plot(t_ccg,popccg(:,1,2),'r--')
plot(t_ccg,popccg(:,2,1),'k--')

%%
for cc = 1:length(allspikes)
    
end
%Time units: ms
%%
figure
plot(t_ccg,allccg(:,end,end))
%%
imagesc(t_ccg,[0 1],log10(squeeze(allccg(:,900,:))'.*1000))
colorbar
%caxis([0 100])
%%
% %v_th = th/(C_e.*J.*tau);
% gammas_E = [0.1 0.3 0.5 0.7];
% %gammas_E = [0.5];
% clear SimValues
% for gg = 1:length(gammas_E)
%     loopparms = parms;
%     loopparms.Kei = loopparms.Kee.*gammas_E(gg);
%     loopparms.Kii = loopparms.Kee.*gammas_E(gg);
%     loopparms.g = (6./4)./gammas_E(gg);
% tic 
% [SimValues{gg}] = Run_LIF_iSTDP(loopparms,TimeParams,'showprogress',true,...
%     'cellout',true,'save_dt',100);
% toc
% end
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
