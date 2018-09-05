function [ spiketimes,rate,input ] = PhaseModulatedSpiking( freq,amp1,phase )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

figfolder = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/Modeling/Figures';
%%
freq1 = 3;
freq2 = 30;
%phase = pi;

T = 5000; %s
dt = 0.001; %s
%t = dt:dt:T;

%amp = 1;

[ amp1,t] = OUNoise(0.5,1.2,T,dt,dt,1);
[ amp2,t] = OUNoise(0.5,1.2,T,dt,dt,1);

freqinput1 = amp1.*cos(2.*pi.*freq1.*t);
freqinput2 = amp2.*cos(2.*pi.*freq2.*t);

rate_0 = 10; %Hz

input = freqinput1 + freqinput2 + log(rate_0);

spkmat = zeros(size(input));
rate = zeros(size(input));

for tt = 1:length(t)
    rate(tt) = exp(input(tt));
    spkmat(tt) = rand(1)<=rate(tt).*dt;
end
spktimes = find(spkmat).*dt;

%%

figure
plot(t,input,'k')
xlim([1 10])

%% Get the LFP explanatory matrix
lfp.data = input+0.3.*randn(size(input));
lfp.timestamps = t;
lfp.samplingRate = 1./dt;
wavespec = bz_WaveSpec(lfp,'nfreqs',10,'frange',[1 100]);


%%
tic
[ GLMFP ] = GLMLFP( {spktimes},wavespec );
toc
%% Simulate ISIs
for tt = 1:length(GLMFP.timestamps)
    GLMFP.spkmat(tt) = rand(1)<=GLMFP.predRate(tt);
end
GLMFP.spktimes = find(GLMFP.spkmat).*dt;
%% ISI
isis = diff(spktimes);
GLMFP.isis = diff(GLMFP.spktimes);

isidist.bins = linspace(-3,1,40);
isidist.hist = hist(log10(isis),isidist.bins);

isidist.histsim = hist(log10(GLMFP.isis),isidist.bins);
%%

xwin = [0 5];
figure
subplot(2,2,3)
%plot(GLMFP.powerbins,GLMFP.Rpower,'k','linewidth',2)
imagesc(log10(wavespec.freqs),GLMFP.powerbins,GLMFP.Rpower)
colorbar
axis xy
axis tight
xlabel('Freq (Hz)');ylabel('Power')
LogScale('x',10)
subplot(4,1,2)
plot(GLMFP.timestamps,GLMFP.predRate,'k')
hold on
plot(t,rate./1000,'r')
legend('Predicted','Actual')
ylabel('Rate (spk/ms)')
xlabel('t (s)')
xlim(xwin)

subplot(4,1,1)
plot(lfp.timestamps,lfp.data,'k')
hold on
plot(spktimes,1.*ones(size(spktimes)),'k.')
ylabel('"LFP"')
xlabel('t (s)')
xlim(xwin)

subplot(2,2,4)
plot(isidist.bins,isidist.hist,'r','linewidth',2)
hold on
plot(isidist.bins,isidist.histsim,'k','linewidth',2)
legend('Actual','Simulated','location','northwest')
LogScale('x',10)
xlabel('ISI (s)')

NiceSave('PhaseCouplingISIToy',figfolder,'PhasePoisson')
%%
xwin = [0 5];
figure
subplot(2,1,1)
imagesc(wavespec.timestamps,log10(wavespec.freqs),abs(wavespec.data'))
xlim(xwin)
axis xy
LogScale('y',10)
subplot(2,1,2)
plot(t,lfp.data,'k')
%plot(t,amp)
hold on
plot(spktimes,5.*ones(size(spktimes)),'k.')
xlim(xwin)
box off

NiceSave
%% Make the explanatory matrix

%With the "raw values"
reX = real(wavespec.data);
imX = imag(wavespec.data);
abX = log10(abs(wavespec.data));
angX = angle(wavespec.data);

X = [reX,imX,abX];


%With bins
numPbins = 20;
phaseedges = linspace(-pi,pi,numPbins+1);
phasecenters = phaseedges(1:end-1)+diff(phaseedges(1:2));
[~,~,phaseBIN] = histcounts(angX,phaseedges);
[~,~,powerBIN] = histcounts(abX);

%binnesphases = interp1(phasebins,1:numPbins,angX,'nearest');
%binnedphases = zeros(size(binnesphases,1),size(binnesphases,2),numPbins);
%binnedphases(:,:,binnesphases) = 1;

binnedphases = zeros(size(phaseBIN,1),numPbins);
for tt = 1:length(phaseBIN)
binnedphases(tt,phaseBIN(tt)) = 1;
end



%%
X = binnedphases;
%%
figure
subplot(1,3,1)
imagesc(binnedphases(1:500,:))




subplot(1,3,2)
imagesc(imX(1:500,:))
subplot(1,3,3)
imagesc(abX(1:500,:))
%% Use glmfit
pGLMwts1 = glmfit(X,spkmat,'poisson');
rateterm =pGLMwts1(1);
phasekernel = pGLMwts1(2:numPbins+1);
%%
figure
plot(phasecenters,phasekernel)
%%
pGLMconst1 = pGLMwts1(1);
pGLMre1 = pGLMwts1(2:1+length(wavespec.freqs));
pGLMim1 = pGLMwts1(length(wavespec.freqs)+2:2*length(wavespec.freqs)+1);
pGLMab1 = pGLMwts1(2*length(wavespec.freqs)+2:end);

%%
figure
subplot(2,2,1)
plot(wavespec.freqs,pGLMre1)
title('RealPart')
xlabel('f (Hz)')

subplot(2,2,2)
plot(wavespec.freqs,pGLMim1)
title('ImagPart')
xlabel('f (Hz)')

subplot(2,2,3)
plot(wavespec.freqs,pGLMab1)
title('Abs')
xlabel('f (Hz)')



%%

%% 
%Make a variable to hold the predicted kernels
kernelPredict = zeros(1+ntfilt+nthist,1);

%Make a function for rate as f'n of kernels (a) and design matrix (units spk/s)
A = @(a,S,R) (S*a(2:ntfilt+1) + R*a(ntfilt+2:ntfilt+nthist+1) + log(a(1)));

%Make a function for loglikelihood as f'n of kernel
%Recall....
%        Poisson likelihood:      P(s|r) = (r*dt)^s/s! exp(-(r.*dt))  
%     (probability of observing spike count s in window of duration dt given rate r)
%     giving log-likelihood:  log P(s|r) =  s log (r*dt) - (r*dt)   
nlogL = @(k) -(spkmat'*(A(k,S,R).*dt) - sum(exp(A(k,S,R)).*dt));

%%minimize nloglik!
initial = 0.1.*ones(size(kernelPredict));
kernelPredict = fminunc(nlogL,initial);

end

