function [ConditionalISI,ConditionalISIModes] = bz_ConditionalISI(spikes,conditionalvariable,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%
%
%Future: multiple conditional variables... (fit together)
%%
p = inputParser;
addParameter(p,'ints',[0 inf])
addParameter(p,'normtype','percentile')
addParameter(p,'minX',50)
addParameter(p,'numISIbins',120)
addParameter(p,'ISIbounds',[0.001 100])
addParameter(p,'numXbins',10)
addParameter(p,'Xwin',[0 1])
addParameter(p,'GammaFitParms',[])
addParameter(p,'GammaFit',false)
addParameter(p,'MutInf',true)
addParameter(p,'ISIDist',true)
addParameter(p,'showfig',true)
addParameter(p,'figname',[])
addParameter(p,'basePath',pwd,@isstr)
addParameter(p,'figfolder',false)
parse(p,varargin{:})
ints = p.Results.ints;
normtype = p.Results.normtype;
minX = p.Results.minX;
numISIbins = p.Results.numISIbins;
logISIbounds = log10(p.Results.ISIbounds);
numXbins = p.Results.numXbins;
Xwin = p.Results.Xwin;
GFParms = p.Results.GammaFitParms;
DO_GammaFit = p.Results.GammaFit;
SHOWFIG = p.Results.showfig;
figname = p.Results.figname;
basePath = p.Results.basePath;
figfolder = p.Results.figfolder;
DO_ISIDist = p.Results.ISIDist;
DO_MutInf = p.Results.MutInf;

%%
baseName = bz_BasenameFromBasepath(basePath);

%%
ISIs.n = diff(spikes);
ISIs.times = spikes(2:end-1);
ISIs.np1 = ISIs.n(2:end);
ISIs.n = ISIs.n(1:end-1);

%% Restrict ISIs and conditional variable to intervals
conditionalvariable.dt = mode(diff(conditionalvariable.timestamps));
inints = InIntervals(conditionalvariable.timestamps,ints);
conditionalvariable.data = conditionalvariable.data(inints);
conditionalvariable.timestamps = conditionalvariable.timestamps(inints);


inints = InIntervals(ISIs.times,ints);
ISIs.n = ISIs.n(inints);
ISIs.np1 = ISIs.np1(inints);
ISIs.times = ISIs.times(inints);

if ~strcmp(normtype,'none')
    [conditionalvariable.data] = NormToInt(conditionalvariable.data,normtype);
end


ISIs.x = interp1(conditionalvariable.timestamps,conditionalvariable.data,ISIs.times,'nearest');
%% Conditional Distribution/Rate, and MutInfo
if DO_ISIDist
    [ ConditionalISIDist ] = ConditionalHist( [ISIs.x;ISIs.x],log10([ISIs.n;ISIs.np1]),...
        'Xbounds',Xwin,'numXbins',numXbins,'Ybounds',logISIbounds,'numYbins',numISIbins,'minX',minX);

    ConditionalISIDist.Xocc = hist(conditionalvariable.data,ConditionalISIDist.Xbins);
    ConditionalISIDist.SpikeRate = ConditionalISIDist.Xhist./(ConditionalISIDist.Xocc.*conditionalvariable.dt.*2);

    %Convert to prob density
    ConditionalISIDist.pYX = ConditionalISIDist.pYX./mode(diff(ConditionalISIDist.Xbins));
end

if DO_MutInf
    ConditionalISIDist.MutInf = mutualinfo(log10([ISIs.n;ISIs.np1]),[ISIs.x;ISIs.x]);
end

if DO_GammaFit
    %% Set up everything for fitting


    %Initial Conditions
    %init_struct.GSlogrates = -log10(meanISI)-0.5;
    init_struct.GSlogrates = GFParms.GSlogrates.*ones(1,numXbins);
    init_struct.GSCVs = GFParms.GSCVs.*ones(1,numXbins);
    init_struct.GSweights = GFParms.GSweights.*ones(1,numXbins);

    % if ASguess
    %     init_struct.ASlogrates = educatedGuess.logrates(1:numAS);
    %     init_struct.ASCVs = educatedGuess.CVs(1:numAS);
    % else
        init_struct.ASlogrates = GFParms.ASlogrates;
        init_struct.ASCVs = GFParms.ASCVs;
    % end
    init_struct.ASweights  = repmat(GFParms.ASweights,numXbins,1);
    init = convertGSASparms(init_struct);



    %%

    taubins = ConditionalISIDist.Ybins./log10(exp(1));
    logISIhist = ConditionalISIDist.pYX'.* mode(diff(ConditionalISIDist.Xbins))./mode(diff(taubins)); %convert to dtau
    numXbins = numXbins;
    numAS = length(GFParms.ASweights);

    %%


    %Upper/Lower Bounds
    clear lb ub
    lb.GSlogrates = -2.*ones(1,numXbins);
    lb.GSCVs =      zeros(1,numXbins);
    lb.GSweights =  zeros(1,numXbins);
    lb.ASlogrates = 0.3.*ones(1,numAS);
    lb.ASCVs =      zeros(1,numAS);
    lb.ASweights  = zeros(numXbins,numAS);
    lb = convertGSASparms(lb);

    ub.GSlogrates = 2.*ones(1,numXbins);
    ub.GSCVs =      4.*ones(1,numXbins);
    ub.GSweights =  ones(1,numXbins);
    ub.ASlogrates = 3.*ones(1,numAS);
    ub.ASCVs =      2.*ones(1,numAS);
    ub.ASweights  = ones(numXbins,numAS);
    ub = convertGSASparms(ub);

    %Make the constraint matrix for all weights to add to 1
    Aeq = zeros(numXbins,length(ub));
    Aeq_ASonly = zeros(numXbins,length(ub));
    Beq = ones(numXbins,1);
    for cc = 1:numXbins
        thiscell.GSlogrates = zeros(1,numXbins);
        thiscell.GSCVs =      zeros(1,numXbins);
        thiscell.GSweights =  zeros(1,numXbins);
        thiscell.ASlogrates = zeros(1,numAS);
        thiscell.ASCVs =      zeros(1,numAS);
        thiscell.ASweights  = zeros(numXbins,numAS);
        thiscell.ASweights(cc,:) = 1;
        Aeq_ASonly(cc,:) = convertGSASparms(thiscell);
        thiscell.GSweights(cc) = 1;
        Aeq(cc,:) = convertGSASparms(thiscell);
    end
    Aeq_ASonly(Aeq_ASonly~=1)=0;
    Aeq(Aeq~=1)=0;

    options = optimoptions('fmincon','Algorithm','sqp' ,'UseParallel',false,'Display','none');%
    %try also: 'Algorithm','interior-point''active-set'
    %Decrease tolerance.....
    options.MaxFunctionEvaluations = 1e8;
    options.MaxIterations = 1000; 

    %% Fit all the distributions together
    AScost_lambda = 0;
    AScost_p = 1/2;
    MScost = 3;
    sub1msbins = ConditionalISIDist.Xbins<=-2.7;

    costfun = @(GSASparm_vect) sum(sum((logISIhist-GSASmodel(GSASparm_vect,taubins,numXbins,numAS)).^2)) ...
        + AScost_lambda.*sum((abs(Aeq_ASonly*GSASparm_vect)).^(AScost_p))...; %L1/2 norm on AS weights to promote sparseness
        + MScost.*sum(sum((logISIhist(sub1msbins,:)-GSASmodel(GSASparm_vect,taubins(sub1msbins),numXbins,numAS)).^2)); 

    fitparms = fmincon(costfun,init,[],[],Aeq,Beq,lb,ub,[],options);
    ConditionalISIModes = convertGSASparms(fitparms,numXbins,numAS);

    %% Mode Correlations
    [ConditionalISIModes.GSCorr,ConditionalISIModes.GScorr_p] = corr(ConditionalISIModes.GSweights',ConditionalISIDist.Xbins','type','Pearson');
    [ConditionalISIModes.ASCorr,ConditionalISIModes.AScorr_p] = corr(ConditionalISIModes.ASweights,ConditionalISIDist.Xbins','type','Pearson');
    ConditionalISIModes.GSCorr = ConditionalISIModes.GSCorr';
    ConditionalISIModes.GScorr_p = ConditionalISIModes.GScorr_p';
    ConditionalISIModes.ASCorr = ConditionalISIModes.ASCorr';
    ConditionalISIModes.AScorr_p = ConditionalISIModes.AScorr_p';

end
%%
if SHOWFIG
testmodel = GSASmodel(ConditionalISIModes,taubins,numXbins,numAS);
figure
subplot(2,2,1)
imagesc(ConditionalISIDist.Xbins,ConditionalISIDist.Ybins,testmodel)
hold on
plot(ConditionalISIDist.Xbins,-ConditionalISIModes.GSlogrates,'ro')
plot(ConditionalISIDist.Xbins([1 end]),-GFParms.GSlogrates.*[1 1],'r--')
%colorbar
caxis([0 0.4])


subplot(2,2,2)
yyaxis left
imagesc(ConditionalISIDist.Xbins,ConditionalISIDist.Ybins,ConditionalISIDist.pYX')
bounds = ylim(gca);
hold on
LogScale('y',10)
ylabel('ISI (s)')
yyaxis right
plot(ConditionalISIDist.Xbins,log10(ConditionalISIDist.SpikeRate),'r','linewidth',2)
ylim(-fliplr(bounds))
LogScale('y',10,'nohalf',true)
ylabel('Rate (Hz)')
%colorbar

subplot(2,2,3)
%plot(sharedfit.ASweights
plot(ConditionalISIModes.ASlogrates(ConditionalISIModes.AScorr_p<0.05),ConditionalISIModes.ASCorr(ConditionalISIModes.AScorr_p<0.05),'o')
hold on
plot(ConditionalISIModes.ASlogrates(ConditionalISIModes.AScorr_p>0.05),ConditionalISIModes.ASCorr(ConditionalISIModes.AScorr_p>0.05),'.')
%hold on
plot(GFParms.GSlogrates,ConditionalISIModes.GSCorr,'o')
plot(xlim(gca),[0 0],'k--')
LogScale('x',10)
xlabel('Rate (Hz)')
ylabel('Weight Correlation')
box off

if figfolder
    NiceSave(['ConditionalISI',figname],figfolder,baseName);
end

end
end

