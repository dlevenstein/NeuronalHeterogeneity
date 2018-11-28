function [ EIGLM ] = GLMEI( spikes,CellClass,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%%
% parse args
p = inputParser;
addParameter(p,'intervals',[0 Inf],@isnumeric)
addParameter(p,'refcell','All')
addParameter(p,'smoothwin',0.05)

parse(p,varargin{:})
intervals = p.Results.intervals;
refcell = p.Results.refcell;
smoothwin = p.Results.smoothwin;

%%
%intervals = SleepState.ints.NREMstate;
%refcell = 9;
%smoothwin = 0.05;
%%

dt = 0.001;
spkmat = bz_SpktToSpkmat(spikes,'binsize',dt);

status = InIntervals(spkmat.timestamps,intervals);
spkmat.data = spkmat.data(status,:);
spkmat.timestamps = spkmat.timestamps(status);

%% Break spikes into this (reference) cell and e/iMUA

spkmat_in = double(spkmat.data(:,refcell));

%Get MUA from other cells
Edx = CellClass.pE; Edx(refcell)=0;
Idx = CellClass.pI; Idx(refcell)=0;
EMUA = sum(spkmat.data(:,Edx),2);
IMUA = sum(spkmat.data(:,Idx),2);
EMUA = movsum(EMUA,smoothwin./dt)./sum(Edx)./smoothwin;
IMUA = movsum(IMUA,smoothwin./dt)./sum(Idx)./smoothwin;

EMUA = double(EMUA);
IMUA = double(IMUA);

%% Figure of MUA
% viewints = bz_RandomWindowInIntervals(intervals,10);
% figure
% subplot(3,1,1)
% plot(spkmat.timestamps,EMUA,'k')
% hold on
% plot(spkmat.timestamps,IMUA,'r')
% xlim(viewints)
% subplot(3,1,2)
% plot(spkmat.timestamps,spkmat_in,'k')
% xlim(viewints)

%% Set up the kernel and indices for parameters
kernelPredict = zeros(3,1)';
ridx = 1;
Eidx = 2;
Iidx = 3;

%% Functions.... clean this up

%Function for log(rate) from E, I and kernels (a)
A = @(a,E,I) (a(Eidx)*E + a(Iidx)*I + log(a(ridx)));% + ...
               % + H*a(histidx));
% V = @(a,POW,TH,abX) (... %(abX*a(ratepowidx) + ...
%     LFPCouplingKernel(POW,TH,a(phaseidx),a(powidx)) + ...
%     log(a(1)));

%Make function for rate as f'n of voltage
%R = @(k) k(kidx).*max(V(k,Erate,Irate)-Vth,0).^k(alphaidx) + k(1);
%R = @(k) k(alphaidx).*exp(k(kidx).*(V(k,Erate,Irate)-Vth));
%Make a function for loglikelihood as f'n of Voltage
% rate = exp(V+Vth)? (but what sets 0... etc?) 
%Recall....
%        Poisson likelihood:      P(s|r) = (r*dt)^s/s! exp(-(r.*dt))  
%     (probability of observing spike count s in window of duration dt given rate r)
%     giving log-likelihood:  log P(s|r) =  s log (r*dt) - (r*dt)   
%nlogL = @(k) -(spkmat'*(A(k,S,R).*dt) - sum(exp(A(k,S,R)).*dt));
%nlogL = @(k) -(spkmat_in'*(A(k,POW,TH).*dt) - sum(exp(A(k,POW,TH)).*dt));
%nlogL = @(k) -(spkmat_in'*(A(k,POW,TH,abX).*dt) - sum(exp(A(k,POW,TH,abX)).*dt));
nlogL = @(k) -(spkmat_in'*(A(k,EMUA,IMUA).*dt) - sum(exp(A(k,EMUA,IMUA)).*dt));

%% Constraints and fitting parameters
%%minimize nloglik!
%Constrain: power kernel>=0
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');%,'UseParallel',true);
%try also: 'Algorithm','active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 2e4;

%Add constraint: no <0 coupling
lb = -inf(size(kernelPredict));
ub = inf(size(kernelPredict));
% lb(Eidx) = 0;
% lb(Iidx) = 0;
% lb(kidx) = 1e-6;
% lb(alphaidx) = 0.5;
% lb(ELidx) = -70;
% ub(ELidx) = -20;
lb(ridx) = 1e-5; ub(ridx) = 100; 
initial = zeros(size(kernelPredict));

initial(ridx) = 1e-3;
initial(Eidx) = 0;
initial(Iidx) = 0;



%%
kernelPredict = fmincon(nlogL,initial,[],[],[],[],lb,ub,[],options);

%%
%NlogL = nlogL(kernelPredict);
R0 = kernelPredict(ridx);
RE = kernelPredict(Eidx);
RI = kernelPredict(Iidx);
predictedrate = exp(A(kernelPredict,EMUA,IMUA));

EIGLM.R0 = R0;
EIGLM.RE = RE;
EIGLM.RI = RI;
EIGLM.predRate = predictedrate;
EIGLM.timestamps = spkmat.timestamps;
EIGLM.dt = dt;
EIGLM.nlogL = nlogL(kernelPredict);

%%
% A(kernelPredict)
% %%
% viewints = bz_RandomWindowInIntervals(intervals,5);
% 
% figure
% subplot(4,1,1)
% plot(spkmat.timestamps,Erate,'k')
% hold on
% plot(spkmat.timestamps,Irate,'r')
% xlim(viewints)
% subplot(4,1,2)
% plot(spkmat.timestamps,spkmat_in,'k')
% xlim(viewints)
% subplot(4,1,3)
% plot(spkmat.timestamps,A(initial,Erate,Irate))
% xlim(viewints)
% 
% subplot(4,1,4)
% plot(spkmat.timestamps,exp(A(initial,Erate,Irate)))
% xlim(viewints)
% 
% %%
% figure
% subplot(4,1,1)
% plot(spkmat.timestamps,EMUA,'k')
% hold on
% plot(spkmat.timestamps,IMUA,'r')
% xlim(viewints)
% subplot(4,1,2)
% plot(spkmat.timestamps,spkmat_in,'k')
% xlim(viewints)
% subplot(4,1,3)
% plot(spkmat.timestamps,A(kernelPredict,EMUA,IMUA))
% xlim(viewints)
% 
% subplot(4,1,4)
% plot(spkmat.timestamps,exp(A(kernelPredict,EMUA,IMUA)))
% xlim(viewints)
% 
% %%
% Vplot = linspace(E_I,E_E,100);
% figure
% subplot(2,2,1)
% %plot(Vplot,initial(kidx).*max(Vplot-Vth,0).^initial(alphaidx) + initial(1));
% plot(Vplot,initial(alphaidx).*exp(initial(kidx).*(Vplot-Vth)))
% 
% %k(alphaidx).*exp(k(alphaidx).*(V(k,Erate,Irate)-Vth))
% subplot(2,2,2)
% plot(Vplot,kernelPredict(alphaidx).*exp(kernelPredict(kidx).*(Vplot-Vth)))
% 
% %plot(Vplot,kernelPredict(kidx).*max(Vplot-Vth,0).^kernelPredict(alphaidx) + kernelPredict(1));
end

