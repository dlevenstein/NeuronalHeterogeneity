function [ GLMFP ] = GLMLFP_raw( spiketimes,lfp,varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% lfp: structure with
%   specgram.data (lfp signal)
%   specgram.timestamps
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
dt = 1/lfp.samplingRate;
spkmat = bz_SpktToSpkmat(spiketimes,'binsize',dt);
spkmat.lfp = interp1(lfp.timestamps,double(lfp.data),spkmat.timestamps,'nearest');

status = InIntervals(spkmat.timestamps,intervals);
spkmat.lfp = spkmat.lfp(status,:);
spkmat.data = spkmat.data(status);
spkmat.timestamps = spkmat.timestamps(status);

spkmat_in = double(spkmat.data);

%%
numchans = size(lfp.data,2);

%%
pGLMwts1 = glmfit(spkmat.lfp,spkmat_in,'poisson','constant','on');

%%
ratepred_pGLM = exp(pGLMwts1(1) + spkmat.lfp*pGLMwts1(2:numchans+1));

%%
% xwin = [500 505];
% figure
% subplot(2,1,1)
% plot(spkmat.timestamps,spkmat.lfp,'k')
% xlim(xwin)
% subplot(2,1,2)
% plot(spkmat.timestamps,ratepred_pGLM./dt,'k')
% xlim(xwin)

%%
GLMFP.R0 = pGLMwts1(1);
GLMFP.Rlfp = pGLMwts1(2:numchans+1);
%GLMFP.predDrive =predictedinput;
GLMFP.predRate = ratepred_pGLM;
GLMFP.timestamps = spkmat.timestamps;
GLMFP.dt = dt;
%GLMFP.ratefun = A
%GLMFP.MSE?
%%
% % rateterm =pGLMwts1(1);
% % phasekernel = pGLMwts1(2:end);
% % 
% % phasekernel = reshape(phasekernel,nphasebins,npowerbins);
% 
% %%
% % figure
% % imagesc(phasekernel)
% 
% %% Making the explanatory matrix
% 
% %abX = abX;
% 
% %%
% %Make a variable to hold the predicted kernels:
% %nfreq = 1;
% kernelPredict = zeros(1+numchans,1)';
% lfpidx = 2:(numchans+1);
% 
% %Kernels: [rate, powerdependence, phase]
% 
% %Make a function for log(rate) (i.e. input!, pre-nonlinearity) as f'n of kernels (a) and design matrix (units spk/s)
% %A = @(a,S,R) (S*a(2:ntfilt+1) + R*a(ntfilt+2:ntfilt+nthist+1) + log(a(1)));
% %A = @(a,POW,TH) (POW*a(2:npowerbins+1)'.*cos(TH+a(npowerbins+2)) + log(a(1)));
% A = @(a,L) (L*a(lfpidx) + log(a(1)));
% %Make a function for loglikelihood as f'n of kernel
% %Recall....
% %        Poisson likelihood:      P(s|r) = (r*dt)^s/s! exp(-(r.*dt))  
% %     (probability of observing spike count s in window of duration dt given rate r)
% %     giving log-likelihood:  log P(s|r) =  s log (r*dt) - (r*dt)   
% %nlogL = @(k) -(spkmat'*(A(k,S,R).*dt) - sum(exp(A(k,S,R)).*dt));
% %nlogL = @(k) -(spkmat_in'*(A(k,POW,TH).*dt) - sum(exp(A(k,POW,TH)).*dt));
% nlogL = @(k) -(spkmat_in'*(A(k,spkmat.lfp).*dt) - sum(exp(A(k,spkmat.lfp)).*dt));
% 
% %%
% % powerkerneltest = linspace(0,1,npowerbins);
% % 
% % test = POW*powerkerneltest'.*sin(TH+0);
% % spkmat.*test;
% 
% %%
% % allkerneltest = [0.1,powerkerneltest,0];
% % 
% % -(spkmat_in'*(A(allkerneltest,POW,TH).*dt) - sum(exp(A(allkerneltest,POW,TH)).*dt));
% % 
% % nlogL(allkerneltest)
% %%
% %%minimize nloglik!
% %Constrain: power kernel>=0
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');%,'UseParallel',true);
% %try also: 'Algorithm','active-set','sqp'
% %Decrease tolerance.....
% options.MaxFunctionEvaluations = 2e4;
% %%
% %Add constraint: no <0 coupling
% lb = -inf(size(kernelPredict));
% ub = inf(size(kernelPredict));
% lb(1) = 1e-5; ub(1) = 100;
% 
% %Aeq = zeros(size(kernelPredict));
% %Aeq(ratepowidx) = 1;
% %beq = 0;
% 
% initial = zeros(size(kernelPredict));
% initial(1) = 1;
% kernelPredict = fmincon(nlogL,initial,[],[],[],[],lb,ub,[],options);
% %Aeq: zeros(size(kernelPredict) Aeq(ratepowidx)=1  beq: 0
% %%
% R0 = kernelPredict(1);
% Rlfp = kernelPredict(lfpidx);
% 
% 
% 
% %figure
% %plot(powercenters,Rpower)
% %% Check the result
% predictedrate = exp(A(kernelPredict,spkmat.lfp));
% figure
% % subplot(2,2,1)
% %     plot(spkmat.timestamps,predictedinput)
% %     hold on
% %     %plot(t,input)
% %     xlim([1 5])
%     %legend('Predicted Drive','Drive')
% subplot(2,2,2)
%     plot(spkmat.timestamps,predictedrate)
%     hold on
%     %plot(t,rate)
%     xlim([1 5])
%     legend('Predicted Rate','Rate')
% 
% % %Make a variable to hold the predicted kernels
% % kernelPredict = zeros(1+ntfilt+nthist,1);
% % 
% % %Make a function for rate (input!, pre-nonlinearity) as f'n of kernels (a) and design matrix (units spk/s)
% % A = @(a,S,R) (S*a(2:ntfilt+1) + R*a(ntfilt+2:ntfilt+nthist+1) + log(a(1)));
% % 
% % %Make a function for loglikelihood as f'n of kernel
% % %Recall....
% % %        Poisson likelihood:      P(s|r) = (r*dt)^s/s! exp(-(r.*dt))  
% % %     (probability of observing spike count s in window of duration dt given rate r)
% % %     giving log-likelihood:  log P(s|r) =  s log (r*dt) - (r*dt)   
% % nlogL = @(k) -(spkmat'*(A(k,S,R).*dt) - sum(exp(A(k,S,R)).*dt));
% % 
% % %%minimize nloglik!
% % initial = 0.1.*ones(size(kernelPredict));
% % kernelPredict = fminunc(nlogL,initial);
%%
% GLMFP.R0 = R0;
% GLMFP.Rlfp = Rlfp;
% %GLMFP.predDrive =predictedinput;
% GLMFP.predRate = predictedrate;
% GLMFP.timestamps = spkmat.timestamps;
% GLMFP.dt = dt;
% %GLMFP.ratefun = A
% %GLMFP.MSE?

end

