function [lambdas,ks,weights,fiterror,returnNmodes] = bz_FitISIGammaModes(ISIs,varargin)
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
%       'returnNmodes'  'auto' (selects based on drop in error)
%       'promthresh'    for auto selection (default: 0.01))
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
addParameter(p,'returnNmodes',6)
addParameter(p,'autoNmodes',true)
addParameter(p,'showfig',true)
addParameter(p,'logbase',10)
addParameter(p,'maxNmodes',10)
%addParameter(p,'lambdabounds',[-5 8])
addParameter(p,'logratebounds',[-3 3])
addParameter(p,'numpad',15)
addParameter(p,'minISIs',250)
addParameter(p,'promthresh',0.01)
addParameter(p,'sequentialreduce',false)
%addParameter(p,'lasso',0)


parse(p,varargin{:})
numpad = p.Results.numpad;
logbase = p.Results.logbase;
maxNmodes = p.Results.maxNmodes;
returnNmodes = p.Results.returnNmodes;
SHOWFIG = p.Results.showfig;
logratebounds = p.Results.logratebounds; %units: loglambda, e
minISIs = p.Results.minISIs;
sequentialreduce = p.Results.sequentialreduce;
promthresh = p.Results.promthresh;
autoNmodes = p.Results.autoNmodes;
%lasso = p.Results.lasso; %lasso doesn't work because weights sum to 1...
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
    function plogt = multigamfun(lambkweit,tau)
        %Parameters: 1) log10 rate  2) log10CV . 3) weights
        k = 1./(10.^lambkweit(end/3+1:end/(3/2))); %alpha (log transform)
        %lambda = exp(lambkweit(1:end/3));  %beta
        lambda = (10.^lambkweit(1:end/3)).*k;
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

for nummodes = maxNmodes:-1:1

%Initialize parms
if nummodes ==maxNmodes || ~sequentialreduce
%     init = [linspace(-1.5,5.5,nummodes)';...    %Lambda 
%         -0.3.*ones(nummodes,1);         %K  (used to be CV=0.8...)
%         ones(nummodes,1)./(nummodes)];             %Weights (normalize later)
    init = [linspace(-0.5,1.5,nummodes)';...    %Lambda 
        -0.3.*ones(nummodes,1);         %K  (used to be CV=0.8...)
        ones(nummodes,1)./(nummodes)];             %Weights (normalize later)
else  %Remove the lowest weight, keep rate/cv, renormalize the weights, refit
    init = fitparms{nummodes+1};
    [~,lowestweightmode] = min(init(end/(3/2)+1:end));
    init(lowestweightmode + (nummodes+1).*[0 1 2]) = [];
    init(end/(3/2)+1:end) = ones(nummodes,1)./(nummodes);
end


difffun = @(lambkweit) sum((logISIhist-multigamfun(lambkweit,taubins)).^2);


ub = [logratebounds(2).*ones(nummodes,1);...    %Lambda
    1.5.*ones(nummodes,1);         %K (optimization parameter is log(CV)
    ones(nummodes,1)];             %Weights (normalize later)
lb =  [logratebounds(1).*ones(nummodes,1);...    %Lambda
    -2.5*ones(nummodes,1);         %K (optimization parameter is log(CV)
    zeros(nummodes,1)];             %Weights (normalize later)

%Constraint: weights sum to 1
Aeq = [zeros(1,nummodes) zeros(1,nummodes) ones(1,nummodes)];
beq = 1;

options = optimoptions('fmincon','Algorithm','sqp','Display','off','UseParallel',true);
%try also: 'Algorithm','active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 1e5;
options.MaxIterations = 1000; 

fitparms{nummodes} = fmincon(difffun,init,[],[],Aeq,beq,lb,ub,[],options);
fiterror(nummodes) = difffun(fitparms{nummodes});

end


%% Pick the number of modes to return based on drop in error
errordrop = [0 diff(log10(fiterror))];
[~,putNmodes,~,P] = findpeaks(-errordrop);
putNmodes(P<promthresh) = [];
savereturnNmodes = returnNmodes;
if autoNmodes
    returnNmodes = max(putNmodes(putNmodes<=savereturnNmodes));
    if isempty(returnNmodes)
        [~,putNmodes,~,P] = findpeaks(-errordrop);
        returnNmodes = max(putNmodes(putNmodes<=savereturnNmodes));
        %Need a better way to do this.... prominence is dependent on error
        %which is dependent on number of spikes
        if isempty(returnNmodes)
            returnNmodes = savereturnNmodes;
        end
    end
end


% AIC/BIC using liklihood...
%%

ks  = 1./(10.^fitparms{returnNmodes}(returnNmodes+1:end-returnNmodes)); %1/k
lambdas = (10.^fitparms{returnNmodes}(1:returnNmodes)).*ks;  %log lambda
weights  = fitparms{returnNmodes}(end-returnNmodes+1:end);

ks(returnNmodes+1:savereturnNmodes) = nan;
lambdas(returnNmodes+1:savereturnNmodes) = nan;
weights(returnNmodes+1:savereturnNmodes) = nan;

%%
if SHOWFIG
   
    
timebins = taubins./log(logbase);
meanISI = ks./lambdas;
%Initialize parms
init = [linspace(-1.5,5.5,returnNmodes)';...    %Lambda 
    -0.3.*ones(returnNmodes,1);         %K  (used to be CV=0.8...)
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

subplot(4,2,4)
stem(log10(1./meanISI),weights)
LogScale('x',10)
xlabel('Rate (Hz)')
ylabel('weight')

subplot(2,2,4)
scatter(log10(1./meanISI),1./ks,10)
%crameri('bilbao')
LogScale('x',10)
xlabel('Rate (Hz)');ylabel('1/k (CV)')

subplot(4,2,5)
plot(trymodes,errordrop,'o-')
hold on
plot(trymodes(putNmodes),errordrop(putNmodes),'r^')
xlabel('Number of Modes')
ylabel('Total Squared Error')
%LogScale('y',10)

subplot(4,2,7)
plot(trymodes,log10(fiterror),'o-')
xlabel('Number of Modes')
ylabel('Total Squared Error')
LogScale('y',10)

if ~sequentialreduce
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

end

