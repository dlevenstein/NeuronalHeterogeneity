function [ output_args ] = cGLMEI( spikes,CellClass,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%%
% parse args
p = inputParser;
addParameter(p,'intervals',[0 Inf],@isnumeric)

parse(p,varargin{:})
intervals = p.Results.intervals;

%%
intervals = SleepState.ints.NREMstate;
refcell = 9;
smoothwin = 0.5;
%%

dt = 0.001;
spkmat = bz_SpktToSpkmat(spikes,'binsize',dt);

status = InIntervals(spkmat.timestamps,intervals);
spkmat.data = spkmat.data(status,:);
spkmat.timestamps = spkmat.timestamps(status);

%%
spkmat_in = double(spkmat.data(:,refcell));

%%
Edx = CellClass.pE; Eidx(refcell)=0;
Idx = CellClass.pI; Iidx(refcell)=0;
Erate = sum(spkmat.data(:,Edx),2);
Irate = sum(spkmat.data(:,Idx),2);
Erate = movsum(Erate,smoothwin./dt)./sum(Edx)./smoothwin;
Irate = movsum(Irate,smoothwin./dt)./sum(Idx)./smoothwin;

%%

spkmat_in = double(spkmat.data(:,refcell));
Erate = double(Erate);
Irate = double(Irate);
%%
viewints = bz_RandomWindowInIntervals(intervals,10);
figure
subplot(3,1,1)
plot(spkmat.timestamps,Erate,'k')
hold on
plot(spkmat.timestamps,Irate,'r')
xlim(viewints)
subplot(3,1,2)
plot(spkmat.timestamps,spkmat_in,'k')
xlim(viewints)

%%
kernelPredict = zeros(1+2+3,1)';
Eidx = 2;
Iidx = 3;
kidx = 4;
alphaidx = 5;
ELidx = 6;


%%
E_E = 0;
E_I = -70;
%E_L = -50; %Should be parameter.... sets baseline voltage in absense of input
Vth = -45; %Voltage threshold

%Function for voltage from E, I, spike history and kernels (a)
V = @(a,E,I) ((a(Eidx)*E*E_E + a(Iidx)*I*E_I + 1*a(ELidx)) ./ ...
                (a(Eidx)*E + a(Iidx)*I + 1));% + ...
               % + H*a(histidx));
% V = @(a,POW,TH,abX) (... %(abX*a(ratepowidx) + ...
%     LFPCouplingKernel(POW,TH,a(phaseidx),a(powidx)) + ...
%     log(a(1)));

%Make function for rate as f'n of voltage
%R = @(k) k(kidx).*max(V(k,Erate,Irate)-Vth,0).^k(alphaidx) + k(1);
R = @(k) k(alphaidx).*exp(k(kidx).*(V(k,Erate,Irate)-Vth));
%Make a function for loglikelihood as f'n of Voltage
% rate = exp(V+Vth)? (but what sets 0... etc?) 
%Recall....
%        Poisson likelihood:      P(s|r) = (r*dt)^s/s! exp(-(r.*dt))  
%     (probability of observing spike count s in window of duration dt given rate r)
%     giving log-likelihood:  log P(s|r) =  s log (r*dt) - (r*dt)   
%nlogL = @(k) -(spkmat'*(A(k,S,R).*dt) - sum(exp(A(k,S,R)).*dt));
%nlogL = @(k) -(spkmat_in'*(A(k,POW,TH).*dt) - sum(exp(A(k,POW,TH)).*dt));
%nlogL = @(k) -(spkmat_in'*(A(k,POW,TH,abX).*dt) - sum(exp(A(k,POW,TH,abX)).*dt));
nlogL = @(k) -(spkmat_in'*(log(R(k).*dt)) - sum(R(k).*dt));

%%
%%minimize nloglik!
%Constrain: power kernel>=0
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');%,'UseParallel',true);
%try also: 'Algorithm','active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 2e4;
%%
%Add constraint: no <0 coupling
lb = -inf(size(kernelPredict));
ub = inf(size(kernelPredict));
lb(Eidx) = 0;
lb(Iidx) = 0;
lb(kidx) = 1e-6;
lb(alphaidx) = 0.5;
lb(ELidx) = -70;
ub(ELidx) = -20;
lb(1) = 1e-10; ub(1) = 1e-4; %should just get rid of this...
initial = zeros(size(kernelPredict));

initial(1) = 1e-10;
initial(Eidx) = 2;
initial(Iidx) = 0.5;
initial(kidx) = 0.005;
initial(alphaidx) = 2;
initial(ELidx) = -50;


%%
kernelPredict = fmincon(nlogL,initial,[],[],[],[],lb,ub,[],options);

%%
nlogL(initial)

%%
R(kernelPredict)
%%
viewints = bz_RandomWindowInIntervals(intervals,10);

figure
subplot(4,1,1)
plot(spkmat.timestamps,Erate,'k')
hold on
plot(spkmat.timestamps,Irate,'r')
xlim(viewints)
subplot(4,1,2)
plot(spkmat.timestamps,spkmat_in,'k')
xlim(viewints)
subplot(4,1,3)
plot(spkmat.timestamps,V(initial,Erate,Irate))
xlim(viewints)

subplot(4,1,4)
plot(spkmat.timestamps,R(initial))
xlim(viewints)

%%
figure
subplot(4,1,1)
plot(spkmat.timestamps,Erate,'k')
hold on
plot(spkmat.timestamps,Irate,'r')
xlim(viewints)
subplot(4,1,2)
plot(spkmat.timestamps,spkmat_in,'k')
xlim(viewints)
subplot(4,1,3)
plot(spkmat.timestamps,V(kernelPredict,Erate,Irate))
xlim(viewints)

subplot(4,1,4)
plot(spkmat.timestamps,R(kernelPredict))
xlim(viewints)

%%
Vplot = linspace(E_I,E_E,100);
figure
subplot(2,2,1)
%plot(Vplot,initial(kidx).*max(Vplot-Vth,0).^initial(alphaidx) + initial(1));
plot(Vplot,initial(alphaidx).*exp(initial(kidx).*(Vplot-Vth)))

%k(alphaidx).*exp(k(alphaidx).*(V(k,Erate,Irate)-Vth))
subplot(2,2,2)
plot(Vplot,kernelPredict(alphaidx).*exp(kernelPredict(kidx).*(Vplot-Vth)))

%plot(Vplot,kernelPredict(kidx).*max(Vplot-Vth,0).^kernelPredict(alphaidx) + kernelPredict(1));
end

