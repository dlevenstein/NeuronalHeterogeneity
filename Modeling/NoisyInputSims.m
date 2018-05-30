function [ spikestats,isifig ] = NoisyInputSims( cellparams,synparams,rates,varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%%
p = inputParser;
addParameter(p,'showfig',true,@islogical)
addParameter(p,'figfolder',false)
addParameter(p,'ntrials',1000,@isnumeric)
parse(p,varargin{:})
SHOWFIG = p.Results.showfig;
figfolder = p.Results.figfolder;
ntrials = p.Results.ntrials; %number of "neurons" or trials

%%

onsettransient = 0.1; %for conductance

%Parameters: DEV
% w_e = 200; %synaptic weight (pS?)
% w_i = 100; %synaptic weight (pS?)
% tau_se = 0.01;
% tau_si = 0.01;
% R_e = 1;
% R_i = 4;
% K_e = 500;
% K_i = 500;
% 
% v_r = -55;
% v_th = -45;
% g_L = 182/18;
% C = 182;
% g_h = 0;
% 
% E_e = 0;
% E_i = -70;
% E_L = -65;
% E_h = 0;

w_e = synparams.w_e; %synaptic weight (pS?)
w_i = synparams.w_i; %synaptic weight (pS?)
tau_se = synparams.tau_se;
tau_si = synparams.tau_si;
K_e = synparams.K_e;
K_i = synparams.K_i;

R_e = rates.R_e;
R_i = rates.R_i;
g_h = rates.g_h;

v_r = cellparams.v_r;
v_th = cellparams.v_th;
g_L = cellparams.g_L;
C = cellparams.C;

E_e = cellparams.E_e;
E_i = cellparams.E_i;
E_L = cellparams.E_L;
E_h = cellparams.E_h;

%Calculate poisson spike times
T = 25;% max simulation time, s
%dt = T./max([(R_e.*K_e),(R_i.*K_i)])./20;
dt = 0.0005;

%RK_e = repmat(R_e.*K_e,ntrials,1);
%RK_i = repmat(R_i.*K_i,ntrials,1);
RK_e = R_e.*K_e;
RK_i = R_i.*K_i;
%% Simulate
TimeStamps = [-onsettransient:dt:T]';

%[ s_e ] = PoissonRateSpikeBins(RK_e,dt,length(TimeStamps));
%[ s_i ] = PoissonRateSpikeBins(RK_i,dt,length(TimeStamps));
lambda_e = RK_e.*dt;
lambda_i = RK_i.*dt;
s_e = poissrnd(lambda_e,length(TimeStamps),ntrials);
s_i = poissrnd(lambda_i,length(TimeStamps),ntrials);

%Initialize conductance, voltage matrices
g_e = nan(length(TimeStamps),ntrials);
g_i = nan(length(TimeStamps),ntrials);
v = nan(length(TimeStamps),ntrials);

%Initial Conditions
ISIs = nan(1,ntrials);
g_e(1,:)=zeros(size(g_e(1,:)));
g_i(1,:)=zeros(size(g_i(1,:)));

for tt = 2:length(TimeStamps)
    if TimeStamps(tt)>0 && all(isnan(v(tt-1,:)))
        v(tt-1,:)=v_r.*ones(size(v(tt-1,:)));
    end
    %Conductance state and Driving potential
    gamma = g_e(tt-1,:) + g_i(tt-1,:) + g_L + g_h;
    H = (g_e(tt-1,:).*E_e + g_i(tt-1,:).*E_i + g_L.*E_L + g_h.*E_h)./gamma;
    
    %Dynamical equations
    dg_edt = -g_e(tt-1,:)./tau_se;
    dg_idt = -g_i(tt-1,:)./tau_si;
    dvdt = (-v(tt-1,:) + H)./(C./gamma);
    
    dg_e = dg_edt.*dt + w_e.*s_e(tt,:);
    dg_i = dg_idt.*dt + w_i.*s_i(tt,:);
    dv = dvdt.*dt;
   
    %effects
    g_e(tt,:) = g_e(tt-1,:)+dg_e;
    g_i(tt,:) = g_i(tt-1,:)+dg_i;
    v(tt,:) = v(tt-1,:)+dv;
    
    spikeneurons = v(tt,:)>v_th;
    v(tt,spikeneurons) = nan;
    ISIs(spikeneurons) = TimeStamps(tt);
    
    if (all(isnan(v(tt,:))) | all(v(tt,:)<(v_r))) && TimeStamps(tt)>0
        %Can put a check - if no spikes by certain time?
       break
    end
    
    
end
   

%% Outputs
%put in bins.... (input?)
[ISIdist.counts,ISIdist.bins] = hist(ISIs,30);
meanV = nanmean(v,2);
stdV = nanstd(v,[],2);

%Calculate expected H? tau (C/gamma)?

%Mean/std of condutance
%Mean/CV(std?) of ISI distribution

%%
spikestats.ISIs = ISIs;
spikestats.rate = 1./nanmean(ISIs);
spikestats.ISI_CV = nanstd(ISIs)./nanmean(ISIs);
%%
timewin = TimeStamps([1 end]);
%timewin = [0 3];
if SHOWFIG
    isifig = figure;
    subplot(2,1,1)
    plot(TimeStamps,v,'color',0.5.*[1 1 1])
    hold on
    plot(TimeStamps,v(:,1),'color',0.2.*[1 1 1])
    plot(TimeStamps,meanV,'k')
    plot(TimeStamps,meanV+stdV,'k--')
    plot(TimeStamps,meanV-stdV,'k--')
    plot(TimeStamps([1 end]),v_th.*[1 1],'r--')
    bar(ISIdist.bins,ISIdist.counts./(0.2.*max(ISIdist.counts))+v_th,...
        'BaseValue',v_th,'facecolor','k');
    %b(1).BaseValue = v_th;
    box off
    xlabel('t (s)')
    ylabel('V')
    xlim(timewin)

    subplot(4,1,3)
    plot(TimeStamps,g_e(:,1),'g')
    xlabel('t (s)')
    ylabel('g_e')
    box off
    xlim(timewin)

    subplot(4,1,4)
    plot(TimeStamps,g_i(:,1),'r')
    xlabel('t (s)')
    ylabel('g_i')
    box off
    xlim(timewin)
    
    if figfolder
        NiceSave('SimISI',figfolder,'CondLIF')
    end
end

end

