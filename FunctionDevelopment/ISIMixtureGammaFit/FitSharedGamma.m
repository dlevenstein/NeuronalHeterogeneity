function [sharedfit,costval,costval_full] = FitSharedGamma(logISIhist,taubins,varargin)
%FITSHAREDGAMMA Summary of this function goes here
% 
%   INPUT
%       logISIhist      [numtimebins x numcells]  probability density (N/(sum*dbin))
%       taubins         [numtimebins x 1] base e
%   (options)
%       'init_struct'
%       'numAS'         (only needed if no initial guess provided)
%       'AScost_lambda'
%       'AScost_p'
%       'MScost'
%       'MSthresh'
%%
p = inputParser;
addParameter(p,'numAS',3)
addParameter(p,'init_struct',[])
addParameter(p,'AScost_lambda',0)
addParameter(p,'AScost_p',2/3)
addParameter(p,'MScost',0)
addParameter(p,'MSthresh',0.002)
addParameter(p,'display','iter')

parse(p,varargin{:})
numAS = p.Results.numAS;
init_struct = p.Results.init_struct;

AScost_lambda = p.Results.AScost_lambda;
AScost_p = p.Results.AScost_p;
MScost = p.Results.MScost;
MSthresh = p.Results.MSthresh;
display = p.Results.display;
%%
numcells = size(logISIhist,2);
%% If there's no initial guess

if isempty(init_struct)
    init_struct.GSlogrates = -0.5.*ones(1,numcells);
    init_struct.GSCVs = 1.5.*ones(1,numcells);
    init_struct.GSweights = 0.5.*ones(1,numcells);

    init_struct.ASlogrates = linspace(1,2.5,numAS);
    init_struct.ASCVs = 0.3.*ones(1,numAS);
    
    init_struct.ASweights  = 0.5.*ones(numcells,numAS)./(numAS);
else
    numAS = length(init_struct.ASlogrates);
end

%%
init = convertGSASparms(init_struct);

%Upper/Lower Bounds
clear lb ub
lb.GSlogrates = -2.*ones(1,numcells);
lb.GSCVs =      zeros(1,numcells);
lb.GSweights =  zeros(1,numcells);
lb.ASlogrates = 0.*ones(1,numAS); %formerly 0.3
lb.ASCVs =      zeros(1,numAS);
lb.ASweights  = zeros(numcells,numAS);
lb = convertGSASparms(lb);

ub.GSlogrates = 2.*ones(1,numcells);
ub.GSCVs =      4.*ones(1,numcells);
ub.GSweights =  ones(1,numcells);
ub.ASlogrates = 3.*ones(1,numAS);
ub.ASCVs =      2.*ones(1,numAS);
ub.ASweights  = ones(numcells,numAS);
ub = convertGSASparms(ub);

%Make the constraint matrix for all weights to add to 1
Aeq = zeros(numcells,length(ub));
Aeq_ASonly = zeros(numcells,length(ub));
Beq = ones(numcells,1);
for cc = 1:numcells
    thiscell.GSlogrates = zeros(1,numcells);
    thiscell.GSCVs =      zeros(1,numcells);
    thiscell.GSweights =  zeros(1,numcells);
    thiscell.ASlogrates = zeros(1,numAS);
    thiscell.ASCVs =      zeros(1,numAS);
    thiscell.ASweights  = zeros(numcells,numAS);
    thiscell.ASweights(cc,:) = 1;
    Aeq_ASonly(cc,:) = convertGSASparms(thiscell);
    thiscell.GSweights(cc) = 1;
    Aeq(cc,:) = convertGSASparms(thiscell);
end
Aeq_ASonly(Aeq_ASonly~=1)=0;
Aeq(Aeq~=1)=0;

options = optimoptions('fmincon','Algorithm','sqp' ,'UseParallel',false,'Display',display);%
%try also: 'Algorithm','interior-point''active-set'
%Decrease tolerance.....
options.MaxFunctionEvaluations = 1e8;
options.MaxIterations = 1000; 

%% Fit all the distributions together

%The Loss Function for each cell
cellloss = @(GSASparm_vect) sum((logISIhist-GSASmodel(GSASparm_vect,taubins,numcells,numAS)).^2);

%Loss function for only refreactory spikes
%Only for high spike density (positive)
smallISI = MSthresh; %<2ms ISIs penalize (make parameter)
sub1msbins = taubins<=log(smallISI); %Which bins are small enough for small-time cost
% cellloss_ref = @(GSASparm_vect) sum(...
%     ((logISIhist(sub1msbins,:)-GSASmodel(GSASparm_vect,taubins(sub1msbins),numcells,numAS))...
%     .*((logISIhist(sub1msbins,:)-GSASmodel(GSASparm_vect,taubins(sub1msbins),numcells,numAS))>0)).^2);
cellloss_ref = @(GSASparm_vect) sum(...
    (logISIhist(sub1msbins,:)-GSASmodel(GSASparm_vect,taubins(sub1msbins),numcells,numAS)).^2);

%Total loss function with regularization etc
costfun = @(GSASparm_vect) sum(cellloss(GSASparm_vect) ...
    + AScost_lambda.*(abs(Aeq_ASonly*GSASparm_vect)').^(AScost_p)...; %L1/2 norm on AS weights to promote sparseness
    + MScost.*cellloss_ref(GSASparm_vect)); 

%Fitting
fitparms = fmincon(costfun,init,[],[],Aeq,Beq,lb,ub,[],options);
costval_full = costfun(fitparms); %Get the total loss
costval = cellloss(fitparms); %Get the loss for each cell

%Convert back into structure for output
sharedfit = convertGSASparms(fitparms,numcells,numAS);

%Sort AS modes by mean rate from high to low
[~,ASratesort] = sort(sharedfit.ASlogrates,'descend');
sharedfit.ASlogrates = sharedfit.ASlogrates(ASratesort);
sharedfit.ASCVs = sharedfit.ASCVs(ASratesort);
sharedfit.ASweights = sharedfit.ASweights(:,ASratesort);

end

