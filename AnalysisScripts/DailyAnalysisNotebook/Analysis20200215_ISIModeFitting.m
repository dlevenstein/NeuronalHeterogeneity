function [ ] = AnalysisXXXXXXXX(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = pwd;
basePath = '/Users/dlevenstein/Dropbox/research/Datasets/Cicero_09102014';
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
%spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
%CellClass = bz_LoadCellinfo(basePath,'CellClass');
%SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');



%% Make a fake distribution: logGamma
taus = linspace(-7,3,100);
lambda = [0.5;500]; %Multiple rates
weights = [1;0.1];
ks =[1;20];

total = sum(weights.*(lambda.^ks).*(exp(ks*taus)) ./ (gamma(ks).*exp(lambda*exp(taus))),1);

figure
plot(taus,total,'k')
NiceSave('TwoGammas',figfolder,baseName,'includeDate',true)

%% Fit the fake distribution
nummodes = length(lambda);
%The multigammafunction of all parameters
multigamfun = @(lambkweit) sum(...
    lambkweit(end/(3/2)+1:end).*...
    (lambkweit(1:end/3).^lambkweit(end/3+1:end/(3/2))).*...
    (exp(lambkweit(end/3+1:end/(3/2))*taus)) ./ ...
    (gamma(lambkweit(end/3+1:end/(3/2))).*...
    exp(lambkweit(1:end/3)*exp(taus))),1);

%Initialize parms
init = [logspace(-2,2,nummodes)';...    %Lambda (make log scale later... following logPoiss
    ones(nummodes,1);         %K (make inverse later - reparamaterize to CV?
    ones(nummodes,1)./nummodes];             %Weights (normalize later)

difffun = @(lambkweit) sum((total-multigamfun(lambkweit)).^2);

ub = inf(size(init));
lb = zeros(size(init));
options = optimoptions('fmincon','Algorithm','sqp');%,'UseParallel',true);
%try also: 'Algorithm','active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 1e5;
options.MaxIterations = 1000; 

%fitparms = fminsearch(difffun,init)
fitparms = fmincon(difffun,init,[],[],[],[],lb,ub,[],options);
%%
figure
plot(taus,total,'k','linewidth',2)
hold on
plot(taus,multigamfun([lambda;ks;weights]),'g--','linewidth',2)
plot(taus,multigamfun(fitparms),'r')


%% Now let's try with a real ISI dist...
excell =5;
logbase = 10;
logISIbins = ISIStats.ISIhist.logbins;
taubins = logISIbins.*log(logbase); %Convert to log base e
testdist = ISIStats.ISIhist.NREMstate.log(excell,:);

numpad = 15;
taubins = [linspace(min(taubins)-5,min(taubins),numpad),taubins,...
    linspace(max(taubins),max(taubins)+5,numpad)];
testdist = [zeros(1,numpad),testdist,zeros(1,numpad)];

%%

trymodes = 1:3;
fiterror = zeros(size(trymodes));
for nummodes = trymodes
nummodes
%Initialize parms
init = [logspace(-2,2,nummodes)';...    %Lambda (make log scale later... following logPoiss
    ones(nummodes,1);         %K (make inverse later - reparamaterize to CV?
    ones(nummodes,1)./nummodes];             %Weights (normalize later)

%The multigammafunction of all parameters
multigamfun = @(lambkweit) sum(...
    lambkweit(end/(3/2)+1:end).*...
    (lambkweit(1:end/3).^lambkweit(end/3+1:end/(3/2))).*...
    (exp(lambkweit(end/3+1:end/(3/2))*taubins)) ./ ...
    (gamma(lambkweit(end/3+1:end/(3/2))).*...
    exp(lambkweit(1:end/3)*exp(taubins))),1);

difffun = @(lambkweit) sum((testdist-multigamfun(lambkweit)).^2);


ub = inf(size(init));
lb = zeros(size(init));
options = optimoptions('fmincon','Algorithm','sqp');%,'UseParallel',true);
%try also: 'Algorithm','active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 1e5;
options.MaxIterations = 1000; 

%fitparms = fminsearch(difffun,init)
fitparms = fmincon(difffun,init,[],[],[],[],lb,ub,[],options);
fiterror(nummodes) = difffun(fitparms);

end



%%
figure
subplot(2,2,1)
plot(taubins,testdist,'k','linewidth',2)
hold on
plot(taubins,multigamfun(fitparms),'r')
subplot(2,2,2)
plot(trymodes,log10(fiterror),'o-')



%% Try log rates and 1/k

trymodes = [1:10,3];
fiterror = zeros(size(trymodes));
for nummodes = trymodes
nummodes
%Initialize parms
init = [linspace(-1.5,5.5,nummodes)';...    %Lambda 
    0.8.*ones(nummodes,1);         %K 
    ones(nummodes,1)./(3.*nummodes)];             %Weights (normalize later)

%The multigammafunction of all parameters
multigamfun = @(lambkweit) sum(...
    lambkweit(end/(3/2)+1:end).*...
    (exp(lambkweit(1:end/3)).^(1./lambkweit(end/3+1:end/(3/2)))).*...
    (exp((1./lambkweit(end/3+1:end/(3/2)))*taubins)) ./ ...
    (gamma((1./lambkweit(end/3+1:end/(3/2)))).*...
    exp(exp(lambkweit(1:end/3))*exp(taubins))),1);

difffun = @(lambkweit) sum((testdist-multigamfun(lambkweit)).^2);


ub = [10.*ones(nummodes,1);...    %Lambda
    inf(nummodes,1);         %K 
    5.*ones(nummodes,1)];             %Weights (normalize later)
lb =  [-4.*ones(nummodes,1);...    %Lambda
    zeros(nummodes,1);         %K 
    zeros(nummodes,1)];             %Weights (normalize later)
options = optimoptions('fmincon','Algorithm','sqp');%,'UseParallel',true);
%try also: 'Algorithm','active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 1e5;
options.MaxIterations = 1000; 

%fitparms = fminsearch(difffun,init)
fitparms = fmincon(difffun,init,[],[],[],[],lb,ub,[],options);
fiterror(nummodes) = difffun(fitparms);

end

%TO DO: Pad the edges with zeros....
%%
rates = fitparms(1:nummodes);  %log lambda
ks  = fitparms(nummodes+1:end-nummodes); %1/k
weights  = fitparms(end-nummodes+1:end);
%%
figure
subplot(4,2,1)
plot(taubins,testdist,'k','linewidth',2)
hold on
plot(taubins,multigamfun(fitparms),'r','linewidth',2)
for mm = 1:nummodes
    plot(taubins,multigamfun(fitparms(nummodes.*[0;1;2]+mm)),'r')
end

subplot(2,2,2)
stem(rates,weights)
%LogScale('x',exp(1))
xlabel('Lambda (Timescale)')
ylabel('weight')

subplot(2,2,4)
scatter(rates,ks,10)
%crameri('bilbao')
xlabel('Lambda (Timescale)');ylabel('1/k (CV)')

subplot(2,2,3)
plot(trymodes,log10(fiterror),'o-')
xlabel('Number of Modes')
ylabel('Error')

subplot(4,2,3)
plot(taubins,testdist,'k','linewidth',2)
hold on
plot(taubins,multigamfun(init),'r','linewidth',2)
for mm = 1:nummodes
    plot(taubins,multigamfun(init(nummodes.*[0;1;2]+mm)),'r')
end
title('Initalization')

