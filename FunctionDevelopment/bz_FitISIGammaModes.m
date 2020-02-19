function [lambdas,ks,weights,fiterror] = bz_FitISIGammaModes(ISIs,varargin)
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
%       'ISIs'          for measuing AIC (maybe just make hist from this?)
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
addParameter(p,'lambdabounds',[-4 7])
addParameter(p,'numpad',15)
addParameter(p,'minISIs',200)

parse(p,varargin{:})
numpad = p.Results.numpad;
logbase = p.Results.logbase;
maxNmodes = p.Results.maxNmodes;
returnNmodes = p.Results.returnNmodes;
SHOWFIG = p.Results.showfig;
lambdabounds = p.Results.lambdabounds; %units: loglambda, e
minISIs = p.Results.minISIs;
%% 
if length(ISIs)<minISIs
    lambdas = nan(returnNmodes,1);
    ks = nan(returnNmodes,1);
    weights = nan(returnNmodes,1);
    fiterror = nan(1,maxNmodes);
    return
end
%% DEV
%excell =5;

%logbins = ISIStats.ISIhist.logbins;
%logISIhist = ISIStats.ISIhist.NREMstate.log(excell,:);

%%
% taubins = logbins.*log(logbase); %Convert to log base e
% %Pad the ISI distirbution with zeros
% taubins = [linspace(min(taubins)-5,min(taubins),numpad),taubins,...
%     linspace(max(taubins),max(taubins)+5,numpad)];
% logISIhist = [zeros(1,numpad),logISIhist,zeros(1,numpad)];

taubins = linspace(-10,8,500);
logISIhist = hist(log(ISIs),taubins);
logISIhist = logISIhist./(sum(logISIhist).*mode(diff(taubins)));
%% The multigammafunction of all parameters

%Sum of loggammas
% multigamfun = @(lambkweit,tau) sum(...
%     lambkweit(end/(3/2)+1:end).*...
%     (exp(lambkweit(1:end/3)).^(1./lambkweit(end/3+1:end/(3/2)))).*...
%     (exp((1./lambkweit(end/3+1:end/(3/2)))*tau)) ./ ...
%     (gamma((1./lambkweit(end/3+1:end/(3/2)))).*...
%     exp(exp(lambkweit(1:end/3))*exp(tau))),1);

    %Sum of loggammas
    function plogt = multigamfun(lambkweit,tau)
        lambda = exp(lambkweit(1:end/3));  %beta
        k = 1./lambkweit(end/3+1:end/(3/2)); %alpha (log transform for fitting?)
        weight = lambkweit(end/(3/2)+1:end);
        
        plogt = sum(...
            weight.*...
            (lambda.^(k).*exp(k*tau)) ./ ...
            (gamma(k).*exp(lambda*exp(tau))),1);
          
    end


%% Try log rates and 1/k

trymodes = [1:maxNmodes];
fiterror = zeros(size(trymodes));
initweightfactor = 10;

for nummodes = trymodes

%Initialize parms
init = [linspace(-1.5,5,nummodes)';...    %Lambda 
    0.8.*ones(nummodes,1);         %K 
    ones(nummodes,1)./(nummodes)];             %Weights (normalize later)


difffun = @(lambkweit) sum((logISIhist-multigamfun(lambkweit,taubins)).^2);


ub = [lambdabounds(2).*ones(nummodes,1);...    %Lambda
    6.*ones(nummodes,1);         %K 
    ones(nummodes,1)];             %Weights (normalize later)
lb =  [lambdabounds(1).*ones(nummodes,1);...    %Lambda
    zeros(nummodes,1);         %K 
    zeros(nummodes,1)];             %Weights (normalize later)

%Constraint: weights sum to 1
Aeq = [zeros(1,nummodes) zeros(1,nummodes) ones(1,nummodes)];
beq = 1;

options = optimoptions('fmincon','Algorithm','sqp','Display','off');%,'UseParallel',true);
%try also: 'Algorithm','active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 1e5;
options.MaxIterations = 1000; 

fitparms{nummodes} = fmincon(difffun,init,[],[],Aeq,beq,lb,ub,[],options);
fiterror(nummodes) = difffun(fitparms{nummodes});

end

%% AIC/BIC
%liklihood = sum(log(multigamfun(fitparms{returnNmodes},ISIs')));
%%
lambdas = exp(fitparms{returnNmodes}(1:returnNmodes));  %log lambda
ks  = 1./fitparms{returnNmodes}(returnNmodes+1:end-returnNmodes); %1/k
weights  = fitparms{returnNmodes}(end-returnNmodes+1:end);
%%
if SHOWFIG
   
    
timebins = taubins./log(logbase);
meanISI = ks./lambdas;
%Initialize parms
init = [linspace(-1.5,5,returnNmodes)';...    %Lambda 
    0.8.*ones(returnNmodes,1);         %K 
    ones(returnNmodes,1)./(returnNmodes)];             %Weights (normalize later)
    
figure
subplot(4,2,1)
plot(timebins,logISIhist,'k','linewidth',2)
hold on
plot(timebins,multigamfun(fitparms{returnNmodes},taubins),'r','linewidth',2)
for mm = 1:returnNmodes
    plot(timebins,multigamfun(fitparms{returnNmodes}(returnNmodes.*[0;1;2]+mm),taubins),'r')
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
xlabel('Rate (Hz)');ylabel('1/k (CV)')

subplot(2,2,3)
plot(trymodes,log10(fiterror),'o-')
xlabel('Number of Modes')
ylabel('Error')

subplot(4,2,3)
plot(timebins,logISIhist,'k','linewidth',2)
hold on
plot(timebins,multigamfun(init,taubins),'r','linewidth',2)
for mm = 1:returnNmodes
    plot(timebins,multigamfun(init(returnNmodes.*[0;1;2]+mm),taubins),'r')
end
title('Initalization')


end

end

