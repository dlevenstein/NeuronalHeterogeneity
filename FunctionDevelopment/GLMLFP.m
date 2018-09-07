function [ GLMFP ] = GLMLFP( spiketimes,specgram,varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Specgram: structure with
%   specgram.data (complex valued spectrogram)
%   specgram.timestamps
%   specgram.freqs
%   specgram.samplingRate (Hz)
% Spiketimes
%
%
%specgram = wavespec;
%
%%
% parse args
p = inputParser;
addParameter(p,'intervals',[0 Inf],@isnumeric)

parse(p,varargin{:})
intervals = p.Results.intervals;

%%
dt = 1/specgram.samplingRate;
spkmat = bz_SpktToSpkmat(spiketimes,'binsize',dt);
spkmat.specgram = interp1(specgram.timestamps,specgram.data,spkmat.timestamps,'nearest');

status = InIntervals(spkmat.timestamps,intervals);
spkmat.specgram = spkmat.specgram(status,:);
spkmat.data = spkmat.data(status);
spkmat.timestamps = spkmat.timestamps(status);

spkmat_in = double(spkmat.data);

abX = zscore(log10(abs(spkmat.specgram)));
%abX = NormToInt(log10(abs(specgram.data)),'modZ');
angX = angle(spkmat.specgram);
%% Binn power for nonparametric power dependence
%nphasebins = 10;
npowerbins = 10;
nfreqs = length(specgram.freqs);

poweredges = linspace(-1.75,1.75,npowerbins+1);
powercenters = poweredges(1:end-1)+0.5.*diff(poweredges(1:2));
poweredges(1) = -Inf; poweredges(end) = Inf;
[~,~,powerBIN] = histcounts(abX,poweredges);

binnedpowers = zeros(size(powerBIN,1),npowerbins);
for ff = 1:nfreqs
    for tt = 1:length(powerBIN)
    binnedpowers(tt,powerBIN(tt,ff),ff) = 1;
    end
end

%X = reshape(binnedphasepowers,length(phaseBIN),[]);
%%
% figure
% subplot(1,3,1)
% imagesc(X(1:5000,:))
% 
% subplot(1,6,3)
% imagesc(spkmat(1:5000,:))
% 
% subplot(1,3,3)
% imagesc(powercenters,1:5000,binnedpowers(1:5000,:))
% hold on
% plot(abX(1:5000),1:5000,'r')
% xlabel('Power')
%% 
%pGLMwts1 = glmfit(X,spkmat,'poisson');
%%
% rateterm =pGLMwts1(1);
% phasekernel = pGLMwts1(2:end);
% 
% phasekernel = reshape(phasekernel,nphasebins,npowerbins);

%%
% figure
% imagesc(phasekernel)

%% Making the explanatory matrix

TH = angX;
POW = binnedpowers;

%%
%Make a variable to hold the predicted kernels:
%nfreq = 1;
kernelPredict = zeros(1+nfreqs.*npowerbins+nfreqs,1)';
powidx = 2:(nfreqs.*npowerbins+1);
phaseidx = (nfreqs.*npowerbins+2):(1+nfreqs.*npowerbins+nfreqs);
%Kernels: [rate, powerdependence, phase]

%Make a function for log(rate) (i.e. input!, pre-nonlinearity) as f'n of kernels (a) and design matrix (units spk/s)
%A = @(a,S,R) (S*a(2:ntfilt+1) + R*a(ntfilt+2:ntfilt+nthist+1) + log(a(1)));
%A = @(a,POW,TH) (POW*a(2:npowerbins+1)'.*cos(TH+a(npowerbins+2)) + log(a(1)));
A = @(a,POW,TH) (LFPCouplingKernel(POW,TH,a(phaseidx),a(powidx)) + log(a(1)));
%Make a function for loglikelihood as f'n of kernel
%Recall....
%        Poisson likelihood:      P(s|r) = (r*dt)^s/s! exp(-(r.*dt))  
%     (probability of observing spike count s in window of duration dt given rate r)
%     giving log-likelihood:  log P(s|r) =  s log (r*dt) - (r*dt)   
%nlogL = @(k) -(spkmat'*(A(k,S,R).*dt) - sum(exp(A(k,S,R)).*dt));
%nlogL = @(k) -(spkmat_in'*(A(k,POW,TH).*dt) - sum(exp(A(k,POW,TH)).*dt));
nlogL = @(k) -(spkmat_in'*(A(k,POW,TH).*dt) - sum(exp(A(k,POW,TH)).*dt));

%%
% powerkerneltest = linspace(0,1,npowerbins);
% 
% test = POW*powerkerneltest'.*sin(TH+0);
% spkmat.*test;

%%
% allkerneltest = [0.1,powerkerneltest,0];
% 
% -(spkmat_in'*(A(allkerneltest,POW,TH).*dt) - sum(exp(A(allkerneltest,POW,TH)).*dt));
% 
% nlogL(allkerneltest)
%%
%%minimize nloglik!
%Constrain: power kernel>=0
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');%,'UseParallel',true);
%try also: 'Algorithm','active-set'
options.MaxFunctionEvaluations = 2e4;
%%
%Add constraint: no <0 coupling
lb = -inf(size(kernelPredict));
ub = inf(size(kernelPredict));
lb(powidx) = 0;
ub(powidx) = 5;
lb(phaseidx) = -3*pi;
ub(phaseidx) = 3*pi;
lb(1) = 1e-4; ub(1) = 100;
initial = zeros(size(kernelPredict));
initial(1) = 1;
kernelPredict = fmincon(nlogL,initial,[],[],[],[],lb,ub,[],options);

%%
R0 = kernelPredict(1);
Rpower = kernelPredict(powidx);
Rpower = reshape(Rpower,npowerbins,nfreqs);
Rphase = kernelPredict(phaseidx);
figure
plot(powercenters,Rpower)
%% Check the result
predictedinput =A(kernelPredict,POW,TH);
predictedrate = exp(predictedinput);
figure
subplot(2,2,1)
    plot(spkmat.timestamps,predictedinput)
    hold on
    %plot(t,input)
    xlim([1 5])
    %legend('Predicted Drive','Drive')
subplot(2,2,2)
    plot(spkmat.timestamps,predictedrate)
    hold on
    %plot(t,rate)
    xlim([1 5])
    legend('Predicted Rate','Rate')

% %Make a variable to hold the predicted kernels
% kernelPredict = zeros(1+ntfilt+nthist,1);
% 
% %Make a function for rate (input!, pre-nonlinearity) as f'n of kernels (a) and design matrix (units spk/s)
% A = @(a,S,R) (S*a(2:ntfilt+1) + R*a(ntfilt+2:ntfilt+nthist+1) + log(a(1)));
% 
% %Make a function for loglikelihood as f'n of kernel
% %Recall....
% %        Poisson likelihood:      P(s|r) = (r*dt)^s/s! exp(-(r.*dt))  
% %     (probability of observing spike count s in window of duration dt given rate r)
% %     giving log-likelihood:  log P(s|r) =  s log (r*dt) - (r*dt)   
% nlogL = @(k) -(spkmat'*(A(k,S,R).*dt) - sum(exp(A(k,S,R)).*dt));
% 
% %%minimize nloglik!
% initial = 0.1.*ones(size(kernelPredict));
% kernelPredict = fminunc(nlogL,initial);
%%
GLMFP.R0 = R0;
GLMFP.Rpower = Rpower;
GLMFP.Rphase = Rphase;
%GLMFP.freqs 
GLMFP.powerbins = powercenters;
GLMFP.predDrive =predictedinput;
GLMFP.predRate = predictedrate;
GLMFP.timestamps = spkmat.timestamps;
%GLMFP.MSE?

end

