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
%       'Nestimatemethod'  'ascending' or 'descending'
%       'autoNmodes'    'LargeInflection','TSEthresh'
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
addParameter(p,'Nestimatemethod','descending')
addParameter(p,'showfig',true)
addParameter(p,'logbase',10)
addParameter(p,'maxNmodes',10)
%addParameter(p,'lambdabounds',[-5 8])
addParameter(p,'logratebounds',[-3 3])
addParameter(p,'numpad',15)
addParameter(p,'minISIs',300)
addParameter(p,'promthresh',0.05)
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
Nestimatemethod = p.Results.Nestimatemethod;
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
%numbins = min(max(round(length(ISIs)./30),150),350);
numbins = 250;
taubins = linspace(-10,8,numbins);
timebins = taubins./log(logbase);
logISIhist = hist(log(ISIs),taubins);
logISIhist = logISIhist./(sum(logISIhist).*mode(diff(taubins)));
%% The multigammafunction of all parameters

    %Sum of loggammas
    function plogt = multigamfun(lambkweit,tau)
        %Parameters: 1) log10 rate  2) log10CV . 3) weights
        k = 1./(lambkweit(end/3+1:end/(3/2))); %alpha
        %lambda = exp(lambkweit(1:end/3));  %beta
        lambda = (10.^lambkweit(1:end/3)).*k;
        weight = lambkweit(end/(3/2)+1:end);
        
        plogt = sum(...
            weight.*...
            (lambda.^(k).*exp(k*tau)) ./ ...
            (gamma(k).*exp(lambda*exp(tau))),1);
        %plogt(isnan(plogt)) = 0;
    end


%% Try log rates and 1/k

trymodes = [1:maxNmodes];
fiterror = zeros(size(trymodes));
initweightfactor = 10;

switch Nestimatemethod
    case 'descending'
        trymodes = [maxNmodes:-1:1];
    case 'ascending'
        trymodes = [1:maxNmodes];
end

for nummodes = trymodes
    
    sigmodeweight = 0.02;
    switch Nestimatemethod
        case 'descending'  
            %Initialize parms
            if nummodes ==maxNmodes || ~sequentialreduce
                init = [linspace(-1,2,nummodes)';...    %log mean rate  
                    0.8.*ones(nummodes,1);         % CV
                    ones(nummodes,1)./(nummodes)];             %Weights
            else  
                init = fitparms{nummodes+1};
                [~,lowestweightmode] = min(init(end/(3/2)+1:end));
                %Remove the lowest weight
                init(lowestweightmode + (nummodes+1).*[0 1 2]) = [];
                sigmodes = (init(end/(3/2)+1:end)>sigmodeweight);
                %renormalize the weights
                init(end/(3/2)+1:end) = init(end/(3/2)+1:end)+ones(nummodes,1)./(nummodes);
                init(end/(3/2)+1:end) = init(end/(3/2)+1:end)./sum(init(end/(3/2)+1:end));
                %redistirbute the insignificant rates between the significant modes
                %sigrates = init(sigmodes);
                %init(~sigmodes) = linspace(-0.5,2,sum(~sigmodes));
                %redistirbute the insignificant rates to residual peaks
                if ~all(sigmodes)
                    resid = movmean(logISIhist-multigamfun(init,taubins),10);
                    [~,peakresid] = findpeaks(resid,'NPeaks',sum(~sigmodes),'SortStr','descend');
                    peakresid = (timebins(peakresid)); %units: log10(ISI)
                    if length(peakresid)<(sum(~sigmodes))
                        peakresid(end+1:(sum(~sigmodes))) = linspace(-2.5,2,sum(~sigmodes)-length(peakresid));
                    end
                    init(find(~sigmodes))=-peakresid';
                end
                %set insignificant mode CV back to 0.8
                init(find(~sigmodes)+nummodes) = 0.8;
            end

        case 'ascending'
            if nummodes == 1
                init = [0;0.8;1];    %log mean rate ; CV ; Weights
            else
                init = fitparms{nummodes-1};
                sigmodes = (init(end/(3/2)+1:end)>sigmodeweight);
                %Remove the insignificant modes
                if any(~sigmodes)
                    init(find(~sigmodes) + (nummodes-1).*[0 1 2]) = [];
                end
                %add a mode at the peaks of the (smoothed) residual (1+insigmodes)
                resid = movmean(logISIhist-multigamfun(init,taubins),10);
%                 if nummodes == 2 %running into issues with low rate
%                 peaks...
%                     %[~,peakresid] = findpeaks(resid,'NPeaks',1+sum(~sigmodes));
%                     peakresid = -1;
%                 else
                    [~,peakresid] = findpeaks(resid,'NPeaks',1+sum(~sigmodes),'SortStr','descend');
                    peakresid = (timebins(peakresid)); %units: log10(ISI)
%                end
                
                if length(peakresid)<(1+sum(~sigmodes))
                    peakresid(end+1:(1+sum(~sigmodes))) = linspace(-2.5,2,1+sum(~sigmodes)-length(peakresid));
                end
                init = [init(1:end/3) ; -peakresid';...  %convert to log10(1/ISI)
                    init(end/3+1:end/(3/2)); 0.8.*ones(size(peakresid'));...
                    init(end/(3/2)+1:end)+(1./nummodes); ones(size(peakresid'))./(nummodes)]; %keep the old weights, add the new weights at fraction
                %set CVs of the new modes
                %renormalize the weights
                %init(end/(3/2)+1:end) = ones(nummodes,1)./(nummodes);
                init(end/(3/2)+1:end) = init(end/(3/2)+1:end)./sum(init(end/(3/2)+1:end));
%                 %% Debug figure
%                 figure
%                 subplot(2,2,1)
%                 plot(timebins,resid)
%                 subplot(2,2,2)
%                     plot(timebins,logISIhist,'k','linewidth',2)
%                     hold on
%                     plot(timebins,multigamfun(fitparms{nummodes-1},taubins),'r','linewidth',2)
%                     plot(timebins,multigamfun(init,taubins),'r--','linewidth',2)
%                     plot(peakresid,zeros(size(peakresid)),'r+');
%                      %for mm = 1:returnNmodes
%                      %    plot(timebins,multigamfun(fitparms{nummodes-1}(returnNmodes.*[0;1;2]+mm),taubins),'r')
%                      %end
            end
    end


difffun = @(lambkweit) sum((logISIhist-multigamfun(lambkweit,taubins)).^2);

lb =  [logratebounds(1).*ones(nummodes,1);...    %Lambda (optimization parameter is log rate)
    0*ones(nummodes,1);         %K (optimization parameter is log(CV)
    zeros(nummodes,1)];             %Weights
ub = [logratebounds(2).*ones(nummodes,1);...    %Lambda
    10.*ones(nummodes,1);         %K (optimization parameter is log(CV)
    ones(nummodes,1)];             %Weights

%Constraint: weights sum to 1
Aeq = [zeros(1,nummodes) zeros(1,nummodes) ones(1,nummodes)];
beq = 1;

options = optimoptions('fmincon','Algorithm','sqp','Display','off');%,'UseParallel',true);
%try also: 'Algorithm','active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 1e5;
options.MaxIterations = 1000; 

ICs{nummodes} = init;
fitparms{nummodes} = fmincon(difffun,init,[],[],Aeq,beq,lb,ub,[],options);
fiterror(nummodes) = difffun(fitparms{nummodes});

end

switch Nestimatemethod
    case 'descending' 
        trymodes = fliplr(trymodes);
end

%% Pick the number of modes to return based on drop in error
fiterror_raw = fiterror;
fiterror_norm = fiterror./min(fiterror);
%fiterror_norm = zscore((fiterror));
errordrop = [0 diff(log10(fiterror))];
inflection = diff(errordrop);
%[~,putNmodes,~,P] = findpeaks(-errordrop);

%%
% figure
% plot(trymodes(2:end),diff(errordrop))
%%
 %Take the largest N peak over some prominence
%putNmodes(P<promthresh) = [];
savereturnNmodes = returnNmodes;
if any(autoNmodes) %If no peaks over that prominence, take the largest prominence peak
switch autoNmodes
    case 'LargeInflection'    
    
    [~,putNmodes,~,P] = findpeaks(inflection);
    
    %Largest N peak with prominence over threshold (for inflection)
    [~,whichmode] = max(putNmodes(putNmodes<=savereturnNmodes & P>promthresh));
    returnNmodes = putNmodes(whichmode);
    
    if isempty(whichmode)
        %Largest prominence peak
        [~,whichmode] = max(P(putNmodes<=savereturnNmodes));
        returnNmodes = putNmodes(whichmode);
    end
    
    case 'dTSEDrop'
        %Largest prominence peak      (for error drop)
        [~,whichmode] = max(P(putNmodes<=savereturnNmodes));
        returnNmodes = putNmodes(whichmode);
        if isempty(whichmode)
            %Largest prominence peak
            [~,whichmode] = max(P(putNmodes<=savereturnNmodes));
            returnNmodes = putNmodes(whichmode);
        end
        
        
    case 'TSEthresh'
        TSEthresh = 0.1;
        returnNmodes = find(log(fiterror_norm)<TSEthresh,1,'first');
        putNmodes = returnNmodes;
end
    

    if isempty(returnNmodes)  %if no peaks
        returnNmodes = savereturnNmodes;
    end
end


% AIC/BIC using liklihood...
%%

ks  = 1./(fitparms{returnNmodes}(returnNmodes+1:end-returnNmodes)); %1/k
lambdas = (10.^fitparms{returnNmodes}(1:returnNmodes)).*ks;  %log lambda
weights  = fitparms{returnNmodes}(end-returnNmodes+1:end);

ks(returnNmodes+1:savereturnNmodes) = nan;
lambdas(returnNmodes+1:savereturnNmodes) = nan;
weights(returnNmodes+1:savereturnNmodes) = nan;

%%
if SHOWFIG
   
   
meanISI = ks./lambdas;
%Initialize parms
init = [linspace(-1.5,5.5,returnNmodes)';...    %Lambda 
    -0.3.*ones(returnNmodes,1);         %K  (used to be CV=0.8...)
    ones(returnNmodes,1)./(returnNmodes)];             %Weights (normalize later)
    
figure

subplot(4,2,1)
plot(timebins,logISIhist,'color',[0.5 0.5 0.5],'linewidth',2)
hold on
plot(timebins,multigamfun(fitparms{returnNmodes},taubins),'k','linewidth',2)
for mm = 1:returnNmodes
    plot(timebins,multigamfun(fitparms{returnNmodes}(returnNmodes.*[0;1;2]+mm),taubins),'r')
end
axis tight
box off
LogScale('x',logbase)

subplot(4,2,3)
stem(log10(1./meanISI),weights)
LogScale('x',10)
xlabel('Rate (Hz)')
ylabel('weight')

subplot(2,2,2)
scatter(log10(1./meanISI),1./ks,10)
%crameri('bilbao')
LogScale('x',10)
xlabel('Rate (Hz)');ylabel('1/k (CV)')

subplot(6,3,11)
plot(trymodes,errordrop,'o-')
hold on
plot(trymodes(putNmodes),errordrop(putNmodes),'r^')
xlabel('Number of Modes')
ylabel('dTSE')

subplot(6,3,12)
plot(trymodes(1:end-1),inflection,'o-')
hold on
plot([0 trymodes(end)],[0 0],'r--')
plot(trymodes(putNmodes),errordrop(putNmodes),'r^')
xlabel('Number of Modes')
ylabel('inflect.')
%LogScale('y',10)

subplot(6,3,10)
plot(trymodes,log10(fiterror),'o-')
xlabel('Number of Modes')
ylabel('logTSE')
%LogScale('y',10)

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

for nn = 1:6
subplot(6,3,nn+12)
plot(timebins,logISIhist,'k','linewidth',2)
hold on
plot(timebins,multigamfun(ICs{nn},taubins),'r--','linewidth',1)
plot(timebins,multigamfun(fitparms{nn},taubins),'r','linewidth',2)
for mm = 1:nn
    plot(timebins,multigamfun(fitparms{nn}(nn.*[0;1;2]+mm),taubins),'r')
end

axis tight
box off
LogScale('x',logbase)

end

end

end
