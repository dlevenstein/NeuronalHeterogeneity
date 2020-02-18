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
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = pwd;
basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
%spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};



%% Make a fake distribution

taus = linspace(-7,3,100);
lambda = [1;10;50]; %Multiple rates
weights = [1;1;1];

total = sum(weights.*lambda*exp(taus) ./ exp(lambda*exp(taus)),1);

%%
nummodes = length(lambda);
multiexpfun = @(lambweit) sum(lambweit(end-nummodes+1:end).*lambweit(1:nummodes)*exp(taus) ./ exp(lambweit(1:nummodes)*exp(taus)),1);

init = [0.1;1;10;1;1;1];
difffun = @(lambweit) sum((total-multiexpfun(lambweit)).^2);

ub = inf(size(init));
lb = zeros(size(init));
%options = optimoptions('fmincon','Display','iter','Algorithm','sqp');%,'UseParallel',true);

%fitparms = fminsearch(difffun,init)
fitparms = fmincon(difffun,init,[],[],[],[],lb,ub);
%%
figure
plot(taus,total,'k','linewidth',2)
hold on
plot(taus,multiexpfun([lambda;weights]),'g--','linewidth',2)
plot(taus,multiexpfun(fitparms),'r')

%% Now let's try with a real ISI dist...
excell = 3;
logbase = 10;
logISIbins = ISIStats.ISIhist.logbins;
taubins = logISIbins.*log(logbase); %Convert to log base e
testdist = ISIStats.ISIhist.NREMstate.log(excell,:);


%%

trymodes = 1:20;
fiterror = zeros(size(trymodes));
for nummodes = trymodes

init = [logspace(-2,2,nummodes)';0.1.*ones(nummodes,1)];
multiexpfun = @(lambweit) sum(lambweit(end-nummodes+1:end).*lambweit(1:nummodes)*exp(taubins) ./ exp(lambweit(1:nummodes)*exp(taubins)),1);


difffun = @(lambweit) sum((testdist-multiexpfun(lambweit)).^2);


ub = [1000.*ones(nummodes,1);ones(nummodes,1)];
lb = [zeros(nummodes,1);-ones(nummodes,1)];
options = optimoptions('fmincon','Algorithm','sqp');%,'UseParallel',true);
%try also: 'Algorithm','active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 2e4;
options.MaxIterations = 1000; 

%fitparms = fminsearch(difffun,init)
fitparms = fmincon(difffun,init,[],[],[],[],lb,ub,[],options);
fiterror(nummodes) = difffun(fitparms);
end
%%
figure
plot(trymodes,log10(fiterror),'o-')
%%
figure
plot(taubins,testdist,'k','linewidth',2)
hold on
plot(taubins,multiexpfun(fitparms),'r')

%% Try log rates

maxmodes = 6;
trymodes = 1:maxmodes;
fiterror = zeros(size(trymodes));
for nummodes = trymodes
nummodes
init = [linspace(-5,5,nummodes)';ones(nummodes,1)./nummodes];
multiexpfun = @(lambweit) sum(lambweit(end/2+1:end).*exp(lambweit(1:end/2))*exp(taubins) ./ exp(exp(lambweit(1:end/2))*exp(taubins)),1);


difffun = @(lambweit) sum((testdist-multiexpfun(lambweit)).^2);


ub = [10.*ones(nummodes,1);ones(nummodes,1)];
lb = [-10.*ones(nummodes,1);zeros(nummodes,1)];
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

% nummodes = 4;
% init = [linspace(-5,5,nummodes)';ones(nummodes,1)./nummodes];
% ub = [10.*ones(nummodes,1);ones(nummodes,1)];
% lb = [-10.*ones(nummodes,1);zeros(nummodes,1)];
% fitparms = fmincon(difffun,init,[],[],[],[],lb,ub,[],options);
%%
rates = fitparms(1:nummodes);
weights  = fitparms(end-nummodes+1:end);
%%
figure
subplot(2,2,1)
plot(taubins,testdist,'k','linewidth',2)
hold on
plot(taubins,multiexpfun(fitparms),'r','linewidth',2)
for mm = 1:nummodes
    plot(taubins,multiexpfun(fitparms(nummodes.*[0;1]+mm)),'r')
end

subplot(2,2,2)
stem(rates,weights)
LogScale('x',exp(1))
xlabel('Rate (Hz)')

subplot(2,2,3)
plot(trymodes,log10(fiterror),'o-')
xlabel('Number of Modes')
ylabel('Error')


%%
taus = linspace(-7,3,100);
lambda = [0.5;400]; %Multiple rates
weights = [1;0.2];
ks =[1;20];

total = sum(weights.*(lambda.^ks).*(exp(ks*taus)) ./ (gamma(ks).*exp(lambda*exp(taus))),1);

figure
plot(taus,total,'k')
NiceSave('TwoGammas',figfolder,baseName,'includeDate',true)

