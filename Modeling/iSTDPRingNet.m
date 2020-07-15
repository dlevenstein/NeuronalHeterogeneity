function iSTDPRingnet(savepath)

%%
savepath = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/Modeling/Simulation_Data';



%%
TimeParams.dt = 0.1;
TimeParams.SimTime = 10000;

%%
%Feedforward parameters
clear parms
parms.N_FF = 1000;
parms.K_FF = [250 750];
parms.J_FF = [0.4 0.1];
%% Ring input

[theta,T] = OUNoise(40000.^-1,4*pi,TimeParams.SimTime,TimeParams.dt,TimeParams.dt.*5,1);
neuronIDX = 2.*pi.*[1:parms.N_FF]./parms.N_FF;

meanrate = 30;

k = 2;
bess = besseli(0,4);
tol = 0.25;
expmean = mean(exp(k.*cos(neuronIDX-theta(1)))./(2.*pi.*bess));

%Rate of the input cells given theta
R_th = @(th) meanrate.*exp(k.*cos(neuronIDX-th))./(2.*pi.*bess)./expmean;
%Rate of the input cells given time (and timeseries theta)
R = @(t) R_th(theta(abs(T-t)<tol));

%%

figure
imagesc(T,neuronIDX,R(T)')
colorbar
hold on
plot(T,mod(theta,2*pi),'k.','markersize',1)

%% Tuning Curve as a function of input sparseness*, weight, N

%%

%Poisson Rate (add to Brunel sim)
%g = 5;
%g = 4; %Initial strength of Inhibitoon (relative to excitation) and I->I strength
v_norm = 2;



parms.EPopNum = 1000;
parms.IPopNum = 250;
parms.u_0 = 20.*v_norm;
parms.u_0 = 0;

parms.V_rest = 0;
%parms.delay_s = 1.8.*rand(parms.EPopNum+parms.IPopNum,parms.EPopNum+parms.IPopNum)+1.2;
parms.delay_s = 1.8.*rand(parms.EPopNum+parms.IPopNum,1)+1.2;

%parms.delay_s = 1.2;
parms.g = g;

parms.V_th =20;
parms.tau_m = 20;
parms.V_reset = 10;
parms.t_ref = 1;

%Conectivity: In degree
gamma = 0.25;
parms.Kee = 250;
parms.Kie = parms.Kee;
parms.Kei = parms.Kee.*gamma;
parms.Kii = parms.Kee.*gamma;



parms.LearningRate = 1e-2;
parms.TargetRate = [sort(exp(1.1.*randn(parms.EPopNum,1)-1.1));nan(parms.IPopNum,1)]; %Target Rate for Excitatory cells (units of Hz)
parms.TargetRate = [sort(exp(randn(parms.EPopNum,1)));nan(parms.IPopNum,1)]; %Target Rate for Excitatory cells (units of Hz)

parms.tauSTDP = 20;    %Time Constant for the STDP curve (Units of ms)

figure
hist(log10(parms.TargetRate))
LogScale('x',10)

netname = 'ring';
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
        %Increase K?...
    case 'ring'
        g = 10;
        parms.g = g;
        parms.J = 0.05; %Weight Exc->Inh %note: 0.05 gives tuning
        parms.ex_rate = R; 
        
        parms.Kee = 500;
        parms.Kie = parms.Kee;
        parms.Kei = parms.Kee.*gamma;
        parms.Kii = parms.Kee.*gamma;
        parms.Kee = 0; 
        
        %Slow down iSTDP?
        parms.LearningRate = 1e-3;
        parms.LearningRate = 2.5e-3;
end

%%
%v_th = th/(C_e.*J.*tau);

tic 
[SimValues] = Run_LIF_iSTDP(parms,TimeParams,'showprogress',true,...
    'cellout',true,'save_dt',1000,'estrate',15);
toc

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
overlay_HD = parms.EPopNum.*mod(theta,2*pi)./(2*pi);
overlay_HD(abs(diff(overlay_HD))>100) = nan;
PlotSimRaster(SimValues,TimeParams.SimTime-[5000 0],...
    'cellsort',sortpeak,'overlay',[T overlay_HD])
NiceSave('iSTDPRaster_late',pwd,netname)

PlotSimRaster(SimValues,[0000 5000],...
    'cellsort',sortpeak,'overlay',[T overlay_HD])
NiceSave('iSTDPRaster_early',pwd,netname)

%% Calculate ISI tuning curve
headdir.timestamps = T./1000;
headdir.data = theta;
spikes.times = cellfun(@(X) X./1000,SimValues.spikesbycell,'UniformOutput',false);
%%
[ISIbyHD] = bz_ConditionalISI(spikes.times,headdir,...
    'showfig',false,'GammaFit',false,'numXbins',10,'numISIbins',100,...
    'normtype','none','Xwin',[0 2.*pi]);



%%
figure
plot(log10(parms.TargetRate),squeeze(ISIbyHD.MutInf),'.')

%%
[~,ISIbyHD.fieldpeak] = max(ISIbyHD.Dist.SpikeRate,[],2);
ISIbyHD.fieldpeak = ISIbyHD.Dist.Xbins(ISIbyHD.fieldpeak);

%%
spikes.numcells = length(spikes.times);
for cc = 1:spikes.numcells
    position_norm = headdir;
    position_norm.data = headdir.data-ISIbyHD.fieldpeak(:,:,cc);
    position_norm.data = mod(position_norm.data+pi,2.*pi)-pi;
[ISIbyHD_align(cc)] = bz_ConditionalISI(spikes.times{cc},position_norm,...
    'showfig',false,'GammaFit',false,'numXbins',15,'numISIbins',100,...
    'normtype','none','Xwin',[-pi pi],'minX',10);
end



ISIbyHD_align_mean = bz_CollapseStruct( ISIbyHD_align(SimValues.EcellIDX),3,'mean',true);


ISIbyHD_align = bz_CollapseStruct( ISIbyHD_align,3,'justcat',true);

%%
figure
subplot(2,2,1)
    imagesc(ISIbyHD_align_mean.Dist.Xbins,ISIbyHD_align_mean.Dist.Ybins,ISIbyHD_align_mean.Dist.pYX')
    hold on
    imagesc(ISIbyHD_align_mean.Dist.Xbins+2.*pi,ISIbyHD_align_mean.Dist.Ybins,ISIbyHD_align_mean.Dist.pYX')
    plot(ISIbyHD_align_mean.Dist.Xbins,-log10(ISIbyHD_align_mean.Dist.SpikeRate),'r')
    plot(ISIbyHD_align_mean.Dist.Xbins+2.*pi,-log10(ISIbyHD_align_mean.Dist.SpikeRate),'r')
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to PF Peak (m)')
    xlim([-pi 3.*pi])
subplot(2,2,2)
    imagesc(ISIbyHD_align_mean.Dist.Xbins,ISIbyHD_align_mean.Dist.Ybins,ISIbyHD_align.Dist.pYX(:,:,800)')
    hold on
    imagesc(ISIbyHD_align_mean.Dist.Xbins+2.*pi,ISIbyHD_align_mean.Dist.Ybins,ISIbyHD_align.Dist.pYX(:,:,800)')
    plot(ISIbyHD_align_mean.Dist.Xbins,-log10(ISIbyHD_align.Dist.SpikeRate(:,:,800)),'r')
    plot(ISIbyHD_align_mean.Dist.Xbins+2.*pi,-log10(ISIbyHD_align.Dist.SpikeRate(:,:,800)),'r')
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to PF Peak (m)')
    xlim([-pi 3.*pi])

%% Save/load
filename = fullfile(savepath,['TrainedNet_',netname]);
save(filename,'SimValues','parms','TimeParams','netname')
%load(filename)

%%
figure
subplot(2,2,1)
imagesc(SimValues.WeightMat_initial(sortpeak,sortpeak))
caxis([-2 0.2])
crameri('vik','pivot',0)
colorbar
subplot(2,2,2)
imagesc(SimValues.WeightMat(sortpeak,sortpeak))
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
timewindows.initialization = [0 100];
timewindows.equib = [100 TimeParams.SimTime./1000];
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

%% Tuning Curve

ISIstats.allspikes.position = cellfun(@(X) interp1(T./1000,mod(theta,2*pi),X,'nearest'),...
    ISIstats.allspikes.times,'UniformOutput',false);
ISIstats.allspikes.position_norm = cellfun(@(X,Y) interp1(T./1000,mod(theta-Y+pi,2*pi),X,'nearest'),...
    ISIstats.allspikes.times,num2cell(peakth),'UniformOutput',false);
ISIstats.allspikes.ISInp1 = cellfun(@(X) [X(2:end);nan],...
    ISIstats.allspikes.ISIs,'UniformOutput',false);

%%
%Conditioning by theta occupancy...
% [ ISIbyPos ] = cellfun(@(X,Y,Z,W) ConditionalHist( [Z;Z],log10([X;Y]),...
%     'Xbounds',[0 2*pi],'numXbins',25,'Ybounds',[-3 2],'numYbins',125,'minX',15,...
%     'conditionby',mod(theta-W+pi,2*pi)),...
%     ISIstats.allspikes.ISIs,ISIstats.allspikes.ISInp1,...
%     ISIstats.allspikes.position_norm,num2cell(peakth),...
%     'UniformOutput',false);

[ ISIbyPos ] = cellfun(@(X,Y,Z,Q) ConditionalHist( [Z(Q);Z(Q)],log10([X(Q);Y(Q)]),...
    'Xbounds',[0 2*pi],'numXbins',20,'Ybounds',[-2.5 1],'numYbins',125,'minX',20),...
    ISIstats.allspikes.ISIs,ISIstats.allspikes.ISInp1,...
    ISIstats.allspikes.position_norm,ISIstats.allspikes.instate.equib,...
    'UniformOutput',false);

ISIbyPos = cat(1,ISIbyPos{:});
ISIbyPos = CollapseStruct( ISIbyPos,3);
ISIbyPos.pYX = ISIbyPos.pYX./mode(diff(T));
ISIbyPos.rate = sum(ISIbyPos.pYX,2);
ISIbyPos.meanpYX.E = nanmean(ISIbyPos.pYX(:,:,SimValues.EcellIDX),3);
ISIbyPos.meanpYX.I = nanmean(ISIbyPos.pYX(:,:,SimValues.IcellIDX),3);

sextgroups = discretize(SimValues.EcellIDX,6);
for ss = 1:6
    ISIbyPos.meanpYX.sextiles(:,:,ss) = nanmean(ISIbyPos.pYX(:,:,sextgroups==ss),3);
    ISIstats.Jointhist.sextiles(ss,:,:) = nanmean(ISIstats.Jointhist.equib.log(sextgroups==ss,:,:),1);

end
%%
classes = {'E','I'};
excell = [650 750 850 950];
figure
for cc = 1:2
subplot(2,2,cc)
imagesc(ISIbyPos.Xbins(1,:,1)-pi,ISIbyPos.Ybins(1,:,1),ISIbyPos.meanpYX.(classes{cc})')
hold on
bz_piTickLabel('x')
%plot(th,bz_NormToRange(-inputtuning(:,excell)),'r')
LogScale('y',10)
title(classes{cc})
end

for ss = 1:5
    subplot(4,5,10+ss)
imagesc(ISIbyPos.Xbins(1,:,1)-pi,ISIbyPos.Ybins(1,:,1),ISIbyPos.meanpYX.sextiles(:,:,ss)')
hold on
%plot(th,bz_NormToRange(-inputtuning(:,excell)),'r')
LogScale('y',10,'exp',true)
bz_piTickLabel('x')
if ss>1
    set(gca,'yticklabel',[])
else
    ylabel('ISI (s)')
end
title(ss)

%     subplot(4,4,12+ss)
% imagesc(ISIbyPos.Xbins(1,:,1)-pi,ISIbyPos.Ybins(1,:,1),ISIbyPos.pYX(:,:,excell(ss))')
% hold on
% %plot(th,bz_NormToRange(-inputtuning(:,excell)),'r')
% LogScale('y',10)
% bz_piTickLabel('x')
% if ss>1
%     set(gca,'yticklabel',[])
% end
end
NiceSave('TuningCurves',pwd,netname)

%%
figure
for ss = 1:6
    subplot(6,4,(ss-1)*4+1)
   imagesc(squeeze(ISIstats.Jointhist.sextiles(ss,:,:) )')
   axis xy
end
