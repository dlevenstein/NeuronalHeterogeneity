function [GammaFit] = bz_FitISISharedGammaModes_new(spikes,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%   INPUTS
%       spikes          a buzcode spikes structure
%                       -or-
%                       [numtimebins x numcells]  probability density (N/(sum*dbin))
%                       dbin must be in units of base e (sorry.... o_O) 
%
%   Options
%       'logtimebins'   put in the time bins, if using a probabilty density
%                       input (log10)
%       'logbase'       (default: 10)
%       'numAS'         number of activated states
%       'figfolder'     a folder to save the figure in
%       'basePath'
%       'AScost_lambda'
%       'AScost_p'
%       'ASguess'       A L1/2 norm cost on Activated states that tries to
%                       push weights to 0 and avoid overlapping AS modes
%       'showfig'
%       'display_results'  fitting display parameter: 'iter' or 'final'

%
%   OUTPUTS
%       GammaFit
%% Input Parser

% parse args
p = inputParser;
addParameter(p,'logbase',10)
addParameter(p,'numAS',6)
addParameter(p,'maxAS',6)
addParameter(p,'showfig',true)
addParameter(p,'figfolder',false)
addParameter(p,'figname',[])
addParameter(p,'basePath',pwd,@isstr)
addParameter(p,'AScost_lambda',0)
addParameter(p,'AScost_p',1)
addParameter(p,'MScost',0)
addParameter(p,'MSthresh',0.002)
addParameter(p,'ASguess',[])
addParameter(p,'meanFR',[])
addParameter(p,'usecells',[])
addParameter(p,'holdweights',false)

addParameter(p,'singlefit',false,@islogical)

addParameter(p,'init_struct',[])
addParameter(p,'ints',[-Inf Inf])

addParameter(p,'logtimebins',[])

addParameter(p,'spkthresh',250)
addParameter(p,'forceRedetect',false,@islogical);
addParameter(p,'savecellinfo',false,@islogical)
addParameter(p,'savenumAS','all')
addParameter(p,'display_results','iter')
addParameter(p,'UseParallel',false)

parse(p,varargin{:})
logbase = p.Results.logbase;
numAS = p.Results.numAS;
maxAS = p.Results.maxAS;
SHOWFIG = p.Results.showfig;
figfolder = p.Results.figfolder;
figname = p.Results.figname;
basePath = p.Results.basePath;
AScost_lambda = p.Results.AScost_lambda;
AScost_p = p.Results.AScost_p;
ASguess = p.Results.ASguess;
MScost = p.Results.MScost;
MSthresh = p.Results.MSthresh;
meanFR = p.Results.meanFR;
SINGLEFIT = p.Results.singlefit;
init_struct = p.Results.init_struct;
ints = p.Results.ints;
logtimebins = p.Results.logtimebins;
spkthresh = p.Results.spkthresh;
usecells = p.Results.usecells;
holdweights = p.Results.holdweights;
forceRedetect = p.Results.forceRedetect;
SAVECELLINFO = p.Results.savecellinfo;
display_results = p.Results.display_results;
UseParallel = p.Results.UseParallel;



educatedGuess.logrates =  [2.0 2.25 0.9  1.5 1.0 0.5]; %HiGamma Burst Theta LowGamma Irregular10Hz ThetaDouble
educatedGuess.CVs =       [0.2 0.05 0.05 0.3 0.6 0.05];

%%
% File naming
baseName = bz_BasenameFromBasepath(basePath);
cellinfofilename = fullfile(basePath,[baseName,'.GammaFit.cellinfo.mat']); %Update in a bit and below
if strcmp(figfolder,'detectionfigures')
    figfolder = [basePath,filesep,'DetectionFigures'];
end

if exist(cellinfofilename,'file') && ~forceRedetect
    GammaFit = bz_LoadCellinfo(basePath,'GammaFit');%Update in a bit
    return
end
% if length(ISIs)<minISIs
%     return
% end
%% Get the ISI stats from the spikes
if isempty(logtimebins)
    %Note: don't do this if the histogram is given...
    useints.ints = ints;
    [ ISIstats ] = bz_ISIStats( spikes,'ints',useints,...
        'shuffleCV2',false,...
        'savecellinfo',false,'basePath',basePath,'forceRedetect',true,...
        'numISIbins',150,'logISIbounds',[0.0001 500]);

    %Curate the histograms
    numspks = cellfun(@sum,ISIstats.allspikes.instate.ints);
    logtimebins = ISIstats.ISIhist.logbins;
    logISIhist = ISIstats.ISIhist.ints.log;
    usecells = find(usecells & numspks>spkthresh);
    logISIhist = logISIhist(usecells,:)';
    logISIhist = logISIhist./mode(diff(logtimebins));
    meanFR = ISIstats.summstats.ints.meanrate(usecells);
    
    taubins = logtimebins./log10(exp(1));
    %Put the logISIhist in probabilty density with bins of size e
    logISIhist = logISIhist.* mode(diff(logtimebins))./mode(diff(taubins)); %convert to dtau
else
    logISIhist = spikes;
    
    %Make the tau (base e) bins
    taubins = logtimebins./log10(exp(1));
end

%% DEV

sub1msbins = logtimebins<=-2.7;

numcells = size(logISIhist,2);

%% First guess: no AS, GS rate = mean Rate

%If don't want to loop, set maxAS to 0. Put initialization inside loop
if ~isempty(init_struct)
    INITGIVEN = true;
    maxAS = 0;
else
    INITGIVEN = false;
end

for aa = 1:(maxAS+1)
    %aa =5;
    %% Making the initialization structure
    %If init_struct is given, skip all this
    if ~INITGIVEN    
        init_struct(aa).GSlogrates = log10(meanFR)-0.25; %Formerly -0.5.... why?
        init_struct(aa).GSCVs = ones(1,numcells);  %Formerly 1.5.... why?
        init_struct(aa).GSweights = 0.5.*ones(1,numcells);

        if aa == 1 %Start: no AS, GS rate = mean Rate
            init_struct(aa).ASlogrates = [];
            init_struct(aa).ASCVs = [];
            init_struct(aa).ASweights = [];

            if ASguess %Former initialization - educated guess
                init_struct(aa).ASlogrates = educatedGuess.logrates(1:numAS);
                init_struct(aa).ASCVs = educatedGuess.CVs(1:numAS);
                init_struct(aa).ASweights  = 0.5.*ones(numcells,numAS)./(numAS);
                %init_struct.ASlogrates = linspace(1,2.5,numAS);
                %init_struct.ASCVs = 0.3.*ones(1,numAS);
            end
        else
            init_struct(aa).ASlogrates = [sharedfit(aa-1).ASlogrates -missingpeak(aa-1)];
            init_struct(aa).ASCVs = [sharedfit(aa-1).ASCVs 0.05];
            init_struct(aa).ASweights  = 0.5.*ones(numcells,aa-1)./(aa-1);
            
            if holdweights && aa>3 %Keep the old weights proportional, but add new mode
                renormweights = (aa-2).*0.5./(aa-1).*...
                    abs(sharedfit(aa-1).ASweights./(sum(sharedfit(aa-1).ASweights,2)));
                init_struct(aa).ASweights(:,1:aa-2) = renormweights;
            end
        end
    end
    %% Fit
    tic
    [sharedfit(aa),costval(aa,:)] = FitSharedGamma(logISIhist,taubins,...
        'MScost',MScost,'MSthresh',MSthresh,'AScost_p',AScost_p,'AScost_lambda',AScost_lambda,...
        'init_struct',init_struct(aa),'display_results',display_results,...
        'UseParallel',true); 
    computetime(aa) = toc

    %%
    weightthresh = 0.01;
    allweights = sharedfit(aa).ASweights;
    allweights(log10(allweights)<-5) = 0.00001;
    numsigAS = sum(allweights>weightthresh,2);



    %% Shared Fit Figure (PUt into FitSharedGamma?)
    lowthreshcolor = [0.95 0.95 0.95];
    numrepeats = 3;
    %excell = excells;
    histcolors = [repmat([1 1 1],numrepeats,1);makeColorMap(lowthreshcolor,[0 0 0])];


    fitISI = GSASmodel(sharedfit(aa),taubins,numcells,aa-1);
    initISI = GSASmodel(init_struct(aa),taubins,numcells,aa-1);
    meanISIdist = mean(logISIhist,2);
    meanISIdist_fit = mean(fitISI,2);
    GScolor = [0.6 0.4 0];
    [~,sortGSrate] = sort(sharedfit(aa).GSlogrates);

    AS_lb = 0; %Lower bound of AS mode rate (note:inverse. log time units)
    missingdensity_mean = meanISIdist-meanISIdist_fit;
    [~,missingpeak(aa)] = max(missingdensity_mean(logtimebins<AS_lb));
    missingpeak(aa) = logtimebins(missingpeak(aa));

    missingdensity_all = movmean(logISIhist,3)-fitISI;
    [~,missingpeak_all] = max(missingdensity_all);
    missingpeak_all = logtimebins(missingpeak_all);
    peakdensity = hist(missingpeak_all,logtimebins);
    peakdensity = movsum(peakdensity,3);

    % figure
    % imagesc(missingdensity_all)

    %%
    figure

        subplot(4,3,6)
            plot(logtimebins,missingdensity_mean)
            hold on
            plot(AS_lb.*[1 1],ylim(gca),'r')
            box off
            plot(xlim(gca),[0 0],'k--')
            xlim([-3 2])
            LogScale('x',10,'exp',true)

        subplot(4,3,9)
            hist(missingpeak_all,logtimebins)
            hold on;box off
            plot(logtimebins,peakdensity)
            xlim([-3 2])
            LogScale('x',10,'exp',true)

    %     subplot(3,3,7)
    %         plot(log10(meanFR),sharedfit(aa).GSlogrates,'.','color',GScolor)
    %         hold on
    %         UnityLine
    %         xlabel('Mean FR (Hz)');ylabel('GS Rate (Hz)')
    %         LogScale('xy',10)

        subplot(3,3,1)

            imagesc(logtimebins,[1 numcells],logISIhist(:,sortGSrate)')
            hold on
            plot(log10(MSthresh).*[1 1],ylim(gca),'r')
            plot(logtimebins,-bz_NormToRange(meanISIdist,0.3)+numcells,'k','linewidth',2)
            plot(-log10(meanFR(sortGSrate)),[1:numcells],'k.')
            plot(-sharedfit(aa).GSlogrates(sortGSrate),[1:numcells],'.','color',GScolor)
            %colorbar
            colormap(gca,histcolors)
            title(figname)
            xlim([-3 2])
            ylabel(['Cells (',num2str(numcells),')'])
            set(gca,'yticklabel',[])
            LogScale('x',10,'exp',true)

        subplot(3,3,7)
            imagesc(logtimebins,[1 numcells],fitISI(:,sortGSrate)')
            hold on
            plot(log10(MSthresh).*[1 1],ylim(gca),'r')
            plot(logtimebins,-bz_NormToRange(meanISIdist_fit,0.3)+numcells,'k','linewidth',2)
            %colorbar
            colormap(gca,histcolors)
            xlim([-3 2])
            LogScale('x',10,'exp',true)

        subplot(3,3,2)
            imagesc(logtimebins,[1 numcells],initISI(:,sortGSrate)')
            hold on
           % plot(logtimebins,-bz_NormToRange(meanISIdist_fit,0.3)+numcells,'k','linewidth',2)
            %colorbar
            plot(xlim(gca),[0 0],'k--')
            colormap(gca,histcolors)
            xlim([-3 2])
            title('Initialization')
            LogScale('x',10,'exp',true)

        subplot(4,3,3)
            plot(logtimebins,meanISIdist_fit,'k--','linewidth',1)
            hold on;box off
            plot(logtimebins,meanISIdist,'k','linewidth',2)
            xlim([-3 2])
            LogScale('x',10,'exp',true)


        subplot(3,3,4)
            plot(-sharedfit(aa).ASlogrates,log10(sharedfit(aa).ASCVs),'ko')
            hold on
            scatter(-sharedfit(aa).GSlogrates,log10(sharedfit(aa).GSCVs),...
                10.*sharedfit(aa).GSweights+0.00001,GScolor,...
                'filled')
            xlim([-3 2])
             LogScale('xy',10)
            plot(xlim(gca),[0 0],'k--')
            box off

        subplot(6,3,17)
            hist(log10(allweights(:)))
            hold on;box off
            plot(log10(weightthresh).*[1 1],ylim(gca),'k--')
            LogScale('x',10)
            xlabel('Weight');ylabel('# Modes (All Cells)')

        subplot(6,3,14)
            hist(numsigAS,linspace(0,aa-1,aa))
            hold on
            box off
            xlabel('# Modes');ylabel('# Cells')

        subplot(3,3,5)
            plot(linspace(0,size(costval,1)-1,size(costval,1)),mean(log10(costval),2),'ko-')
            hold on
            plot(linspace(0,size(costval,1)-1,size(costval,1)),log10(costval),'k.')
            box off
            axis tight
            LogScale('y',10,'nohalf',true)
            xlabel('# Modes');ylabel('Errror')


    if figfolder
        NiceSave(['SharedGammaModes',num2str(aa-1),'AS_',figname],figfolder,baseName);
    end

end

%% Here: Pick the best number of AS modes (temporary)
%savenumAS
% if strcmp(savenumAS,'all')
% else
%     sharedfit = sharedfit(savenumAS);
% end

%% The single cell fit (after picking...)
if SINGLEFIT
    %% Fit Each cell distribution, starting from the group dists
    display('Group fit Complete! Now Fitting each cell independently')
    for cc = 1:numcells
        bz_Counter(cc,numcells,'Cell')
        %Initial Conditions
        cinit_struct.GSlogrates = sharedfit(numAS+1).GSlogrates(cc);
        cinit_struct.GSCVs = sharedfit(numAS+1).GSCVs(cc);
        cinit_struct.GSweights = sharedfit(numAS+1).GSweights(cc);

        cinit_struct.ASlogrates = sharedfit(numAS+1).ASlogrates;
        cinit_struct.ASCVs = sharedfit(numAS+1).ASCVs;
        cinit_struct.ASweights  = sharedfit(numAS+1).ASweights(cc,:);
        cinit = convertGSASparms(cinit_struct);

        %Upper/Lower Bounds
        clear clb cub
        clb.GSlogrates = -2.*ones(1,1);
        clb.GSCVs =      zeros(1,1);
        clb.GSweights =  zeros(1,1);
        clb.ASlogrates = 0.*ones(1,numAS); %Was 0.3
        clb.ASCVs =      zeros(1,numAS);
        clb.ASweights  = zeros(1,numAS);
        clb = convertGSASparms(clb);

        cub.GSlogrates = 2.*ones(1,1);
        cub.GSCVs =      4.*ones(1,1);
        cub.GSweights =  ones(1,1);
        cub.ASlogrates = 2.5.*ones(1,numAS);
        cub.ASCVs =      2.*ones(1,numAS);
        cub.ASweights  = ones(1,numAS);
        cub = convertGSASparms(cub);

        %Make the constraint matrix for all weights to add to 1
        cAeq = zeros(1,length(cub));
        cAeq_ASonly = zeros(1,length(cub));
        cBeq = ones(1,1);
        thiscell.GSlogrates = zeros(1,1);
        thiscell.GSCVs =      zeros(1,1);
        thiscell.GSweights =  zeros(1,1);
        thiscell.ASlogrates = zeros(1,numAS);
        thiscell.ASCVs =      zeros(1,numAS);
        thiscell.ASweights  = ones(1,numAS);
        cAeq_ASonly = convertGSASparms(thiscell)';
        thiscell.GSweights =  ones(1,1);
        cAeq = convertGSASparms(thiscell)';
        cAeq_ASonly(cAeq_ASonly~=1)=0;
        cAeq(cAeq~=1)=0;

        options = optimoptions('fmincon','Algorithm','sqp','UseParallel',false,'Display','off');%
        %try also: 'Algorithm','active-set', 'sqp'
        %Decrease tolerance.....
        options.MaxFunctionEvaluations = 1e8;
        options.MaxIterations = 1000; 

        %% Fit the single-cell distribution 
        thisdist = logISIhist(:,cc);

        cdifffun = @(GSASparm_vect) sum(sum((thisdist-GSASmodel(GSASparm_vect,taubins,1,numAS)).^2)) ...
        + AScost_lambda.*sum((abs(cAeq_ASonly*GSASparm_vect)).^(AScost_p))... ; %L2/3 norm on AS weights to promote sparseness
        + MScost.*sum(sum((thisdist(sub1msbins)-GSASmodel(GSASparm_vect,taubins(sub1msbins),1,numAS)).^2)); 

        fitparms_singlecell = fmincon(cdifffun,cinit,[],[],cAeq,cBeq,clb,cub,[],options);
    %     if(zerospkcells(cc))
    %         fitparms_singlecell = nan(size(fitparms_singlecell));
    %     end
        GammaFit.singlecell(cc) = convertGSASparms(fitparms_singlecell,1,numAS);


    end
    
    singlecell_all = bz_CollapseStruct(GammaFit.singlecell,1);    
end
%%
%GammaFit.singlecell(zerospkcells) = [];
GammaFit.taubins = taubins;
GammaFit.ISIdists = logISIhist;
GammaFit.sharedfit = sharedfit;
GammaFit.logtimebins = logtimebins;
GammaFit.numcells = numcells;
GammaFit.initialconditions = init_struct;
GammaFit.costval = costval;
GammaFit.computetime = computetime;

GammaFit.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');
GammaFit.detectorinfo.detectorname = 'bz_FitISISharedGammaModes';
GammaFit.detectorinfo.detectionparms = p.Results;
%Collapse the structure

%%
%Move this into function (note: might need subfunctions)
GammaFit.cellstats.meanrate = meanFR;

try
GammaFit.cellstats.UID = spikes.UID(usecells);
catch
end
if isfield(spikes,'region')
    GammaFit.cellstats.region = spikes.region(usecells);
end

    

% if length(GammaFit.(statenames{ss}).cellstats.meanrate) ~= ...
%         length(GammaFit.(statenames{ss}).singlecell)
%     error('bad number of cells')
% end

%% Saving
% numstates = 1;
% if numstates >1
%     %joint indexing across different sets of intervals
%     GammaFit.(statenames{1}).cellstats.NW = false(size(GammaFit.(statenames{1}).cellstats.meanrate));
%     GammaFit.(statenames{2}).cellstats.NW = false(size(GammaFit.(statenames{2}).cellstats.meanrate));
%     [~,NW1,NW2]=intersect(usecells{1},usecells{2});
%     GammaFit.(statenames{1}).cellstats.NW(NW1)=true;
%     GammaFit.(statenames{2}).cellstats.NW(NW2)=true;
% end
% 
% if SAVECELLINFO
%     save(cellinfofilename,'GammaFit')
% end


%% SKIP FOR NOW
if SHOWFIG | figfolder
    
GScolor = [0.6 0.4 0];

%figure
%hist(singlecell_all.GSweights)
%%
[~,sortGSrate] = sort(sharedfit(numAS+1).GSlogrates);
for aa = 1:numAS
    [~,sortASweight{aa}] = sort(sharedfit(numAS+1).ASweights(:,aa));
end

%% ISI dist
meanISIdist = mean(logISIhist,2);
%%
%sortrate = sort(ISIStats.summstats.WAKEstate.meanrate);
fitISI = GSASmodel(sharedfit(numAS+1),taubins,numcells,numAS);
figure
subplot(2,2,1)
imagesc(logtimebins,[1 numcells],logISIhist(:,sortGSrate)')
hold on
plot(logtimebins,-bz_NormToRange(meanISIdist,0.3)+numcells,'k','linewidth',2)
%colorbar
title(figname)
xlim([-3 2])

subplot(2,2,2)
imagesc(logtimebins,[1 numcells],fitISI(:,sortGSrate)')
%colorbar
xlim([-3 2])





subplot(2,2,3)
plot(-sharedfit(numAS+1).ASlogrates,log10(sharedfit(numAS+1).ASCVs),'o')
hold on
% scatter(-sharedfit.GSlogrates,sharedfit.GSCVs,...
%     5.*singlecell_all.GSweights+0.00001,log10(ISIStats.summstats.WAKEstate.meanrate(ISIStats.sorts.WAKEstate.ratepE)),...
%     'filled')
scatter(-sharedfit(numAS+1).GSlogrates,log10(sharedfit(numAS+1).GSCVs),...
    10.*singlecell_all.GSweights+0.00001,...
    'filled')
for aa = 1:numAS
scatter(-singlecell_all.ASlogrates(:,aa),log10(singlecell_all.ASCVs(:,aa)),25.*singlecell_all.ASweights(:,aa)+0.00001,...
    'filled')
end
%axis tight
box off
plot(logtimebins([1 end]),[0 0],'k--')
xlabel('Mean ISI (s)');ylabel('CV')
xlim([-3 2])
LogScale('x',10)

subplot(2,2,4)
plot(singlecell_all.GSlogrates,sharedfit(numAS+1).GSlogrates,'.')
xlabel('Single-cell GS rate');ylabel('Group Fit GS rate')

if figfolder
    NiceSave(['GammaModes',figname],figfolder,baseName);
end

%%
figure
subplot(3,3,1)
imagesc(taubins,[1 numcells],logISIhist(:,sortGSrate)')
%colorbar
for aa = 1:numAS
    subplot(3,3,3+aa)
    imagesc(taubins,[1 numcells],logISIhist(:,sortASweight{aa})')
end


%% Single cell example(s)
numex = 3;
excells = randi(numcells,numex);
figure
for ee = 1:numex
    excell = excells(ee);
subplot(3,3,ee)
    plot(logtimebins,logISIhist(:,excell),'color',[0.5 0.5 0.5],'linewidth',2)
    hold on
    plot(logtimebins,fitISI(:,excell),'k--','linewidth',1)
    plot(GammaFit.logtimebins,...
        GSASmodel(GammaFit.singlecell(excell),...
        GammaFit.taubins),...
        'k','linewidth',2)
    hold on
    plot(logtimebins,...
        LogGamma(singlecell_all.GSlogrates(excell),...
        singlecell_all.GSCVs(excell),...
        singlecell_all.GSweights(excell)',taubins'),'color',GScolor);
    for aa = 1:numAS
        plot(logtimebins,...
            LogGamma(singlecell_all.ASlogrates(excell,aa),...
            singlecell_all.ASCVs(excell,aa),...
            singlecell_all.ASweights(excell,aa)',taubins'),'k');
    end
    box off
    axis tight
    %title(

subplot(3,3,3+ee)
scatter(-singlecell_all.ASlogrates(excell,:),log10(singlecell_all.ASCVs(excell,:)),...
    100*singlecell_all.ASweights(excell,:)+0.00001,'k','filled')
hold on
scatter(-singlecell_all.GSlogrates(excell),log10(singlecell_all.GSCVs(excell)),...
    100*singlecell_all.GSweights(excell)+0.00001,GScolor,'filled')
plot(logtimebins([1 end]),[0 0],'k--')
ylabel('CV');xlabel('mean ISI (s)')
xlim(logtimebins([1 end]));ylim([-2 1])
LogScale('x',10,'exp',true)

end
if figfolder
    NiceSave(['CellExample',figname],figfolder,baseName);
end


%% Example figures for process

if ~isfield(GammaFit.cellstats,'UID')
    exUID = randsample(1:length(GammaFit.cellstats.meanrate),1);
else
    exUID = randsample(GammaFit.cellstats.UID,1);
end

% bz_PlotISIDistModes(GammaFit.(statenames{ss}),GammaFit.(statenames{ss}).cellstats.UID,'whichShare',pp)
 lowthreshcolor = [0.95 0.95 0.95];
numrepeats = 3;
%excell = excells;
histcolors = [repmat([1 1 1],numrepeats,1);makeColorMap(lowthreshcolor,[0 0 0])];   

figure
for pp = 1:numAS+1
    fitISI = GSASmodel(GammaFit.sharedfit(pp),...
        GammaFit.taubins,GammaFit.numcells,pp-1);
    [~,sortGSrate] = sort(GammaFit.sharedfit(pp).GSlogrates);

    subplot(3,7,pp)
        imagesc(GammaFit.logtimebins,[1 GammaFit.numcells],fitISI(:,sortGSrate)')
        hold on
        %plot(log10(MSthresh).*[1 1],ylim(gca),'r')
        %plot(logtimebins,-bz_NormToRange(meanISIdist_fit,0.3)+numcells,'k','linewidth',2)
        %colorbar
        colormap(gca,histcolors)
        xlim([-3 2])
        set(gca,'xtick',[]);set(gca,'ytick',[])
        title([num2str(pp-1),' AS Modes'])

    subplot(3,7,pp+7)
        bz_PlotISIDistModes(GammaFit,'all','showSingleFits',true,...
            'whichShare',pp,'dotscale',10,'dotscaleAS',150)
        ylim([-1.5 1.6])
        LogScale('y',10,'nohalf',true)
        if pp>1
            set(gca,'yticklabels',[])
            ylabel('')
        end
        box off

    subplot(3,7,pp+14)
        bz_PlotISIDistModes(GammaFit,exUID,'whichShare',pp)
        ylim([-1.5 1.6])
        LogScale('y',10,'nohalf',true)
        if pp>1
            set(gca,'yticklabels',[])
            ylabel('')
        end
        box off
end

if figfolder
    NiceSave(['CompareNAS_',figname],figfolder,baseName);
end


end
end
