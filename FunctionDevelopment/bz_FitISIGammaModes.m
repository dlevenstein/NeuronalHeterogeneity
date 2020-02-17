function [lambdas,ks,weights,fiterror] = bz_FitISIGammaModes(logbins,logISIhist,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%   INPUTS
%       logbins
%       logISIhist
%
%   Options
%       'logbase'       (default: 10)
%       'numpad'        (number of bins below/above to pad)
%       'maxNmodes'  
%       'returnNmodes'
%       'showfig'
%       'lambdabounds'
%
%   OUTPUTS
%       lambdas: the beta parameter (mean ISI is ks/lambdas)
%       ks: the alpha parameter     (CV is 1/ks)
%% Input Parser

% parse args
p = inputParser;
addParameter(p,'returnNmodes',3)
addParameter(p,'showfig',true)
addParameter(p,'logbase',10)
addParameter(p,'maxNmodes',10)
addParameter(p,'lambdabounds',[-4 10])
addParameter(p,'numpad',15)

parse(p,varargin{:})
numpad = p.Results.numpad;
logbase = p.Results.logbase;
maxNmodes = p.Results.maxNmodes;
returnNmodes = p.Results.returnNmodes;
SHOWFIG = p.Results.showfig;
lambdabounds = p.Results.lambdabounds; %units: loglambda, e


%% DEV
%excell =5;

%logbins = ISIStats.ISIhist.logbins;
%logISIhist = ISIStats.ISIhist.NREMstate.log(excell,:);

%%
taubins = logbins.*log(logbase); %Convert to log base e
%Pad the ISI distirbution with zeros
taubins = [linspace(min(taubins)-5,min(taubins),numpad),taubins,...
    linspace(max(taubins),max(taubins)+5,numpad)];
logISIhist = [zeros(1,numpad),logISIhist,zeros(1,numpad)];

%% The multigammafunction of all parameters
multigamfun = @(lambkweit) sum(...
    lambkweit(end/(3/2)+1:end).*...
    (exp(lambkweit(1:end/3)).^(1./lambkweit(end/3+1:end/(3/2)))).*...
    (exp((1./lambkweit(end/3+1:end/(3/2)))*taubins)) ./ ...
    (gamma((1./lambkweit(end/3+1:end/(3/2)))).*...
    exp(exp(lambkweit(1:end/3))*exp(taubins))),1);


%% Try log rates and 1/k

trymodes = [1:maxNmodes];
fiterror = zeros(size(trymodes));

for nummodes = trymodes

%Initialize parms
init = [linspace(-1.5,5.5,nummodes)';...    %Lambda 
    0.8.*ones(nummodes,1);         %K 
    ones(nummodes,1)./(3.*nummodes)];             %Weights (normalize later)


difffun = @(lambkweit) sum((logISIhist-multigamfun(lambkweit)).^2);


ub = [lambdabounds(2).*ones(nummodes,1);...    %Lambda
    inf(nummodes,1);         %K 
    5.*ones(nummodes,1)];             %Weights (normalize later)
lb =  [lambdabounds(1).*ones(nummodes,1);...    %Lambda
    zeros(nummodes,1);         %K 
    zeros(nummodes,1)];             %Weights (normalize later)
options = optimoptions('fmincon','Algorithm','sqp','Display','off');%,'UseParallel',true);
%try also: 'Algorithm','active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 1e5;
options.MaxIterations = 1000; 

%fitparms = fminsearch(difffun,init)
fitparms{nummodes} = fmincon(difffun,init,[],[],[],[],lb,ub,[],options);
fiterror(nummodes) = difffun(fitparms{nummodes});

end

%TO DO: Pad the edges with zeros....
%%
lambdas = exp(fitparms{returnNmodes}(1:returnNmodes));  %log lambda
ks  = 1./fitparms{returnNmodes}(returnNmodes+1:end-returnNmodes); %1/k
weights  = fitparms{returnNmodes}(end-returnNmodes+1:end);
%%
if SHOWFIG
   
    
timebins = taubins./log(logbase);
meanISI = ks./lambdas;
%Initialize parms
init = [linspace(-1.5,5.5,returnNmodes)';...    %Lambda 
    0.8.*ones(returnNmodes,1);         %K 
    ones(returnNmodes,1)./(3.*returnNmodes)];             %Weights (normalize later)
    
figure
subplot(4,2,1)
plot(timebins,logISIhist,'k','linewidth',2)
hold on
plot(timebins,multigamfun(fitparms{returnNmodes}),'r','linewidth',2)
for mm = 1:returnNmodes
    plot(timebins,multigamfun(fitparms{returnNmodes}(returnNmodes.*[0;1;2]+mm)),'r')
end
LogScale('x',logbase)

subplot(2,2,2)
stem(log10(1./meanISI),weights)
LogScale('x',10)
xlabel('Rate (Hz)')
ylabel('weight')

subplot(2,2,4)
scatter(log10(1./meanISI),1./ks,10)
%crameri('bilbao')
LogScale('x',10)
xlabel('Lambda (Timescale)');ylabel('1/k (CV)')

subplot(2,2,3)
plot(trymodes,log10(fiterror),'o-')
xlabel('Number of Modes')
ylabel('Error')

subplot(4,2,3)
plot(timebins,logISIhist,'k','linewidth',2)
hold on
plot(timebins,multigamfun(init),'r','linewidth',2)
for mm = 1:returnNmodes
    plot(timebins,multigamfun(init(returnNmodes.*[0;1;2]+mm)),'r')
end
title('Initalization')


end

end

