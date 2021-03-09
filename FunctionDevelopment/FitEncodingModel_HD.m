function [] = FitEncodingModel_HD(s,x,dt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Possible here: make a guess that's close! (Use continuous model fit?....)



init = zeros(6,1);
init(1) = pi;    %x0
init(2) = 1;    %k
init(3) = 0.8;  %pAS_0
init(4) = 0.1;  %pAS_pi
init(5) = 5;   %rAS
init(6) = 1;    %rGS

lb(1) = -5*pi;    %x0
lb(2) = 0.1;    %k
lb(3) = 0;   %pAS_0
lb(4) = 0;  %pAS_pi
lb(5) = 1;   %rAS
lb(6) = 1e-2;    %rGS

ub(1) = 5*pi;    %x0
ub(2) = 10;    %k
ub(3) = 1;   %pAS_0
ub(4) = 1;  %pAS_pi
ub(5) = 200;   %rAS
ub(6) = 20;    %rGS

%likelihood = @(a,x,s) double(EncodingModal_VMPoisson(s,x,a(1),a(2),a(3),a(4),a(5),a(6),dt));
likelihood = @(a,x,s) (EncodingContinuous_VMPoisson(double(s),x,a(1),a(2),a(5),a(6),dt));
nlogL = @(k) -(nansum(log(likelihood(k,x,s))));

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');%,'UseParallel',true);
%try also: 'Algorithm','active-set','sqp'
options.MaxFunctionEvaluations = 2e4;

%%
kernelPredict = fmincon(nlogL,init,[],[],[],[],lb,ub,[],options);

%%
nlogL(init)
%%
likelihood(init,x,s)

%% Plot: fitted vs observed tuning curve



end

