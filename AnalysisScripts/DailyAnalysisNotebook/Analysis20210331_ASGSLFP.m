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
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
%basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Cicero_09102014');
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
%sessionInfo = bz_getSessionInfo(basePath);
%spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
%CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
%statecolors = {'k','b','r',[0.6 0.6 0.6]};

% try
%     celltypes = CellClass.celltypes;
% catch
%     celltypes = unique(CellClass.label);
% end
% cellcolor = {'k','r'};


%%
load([basePath,'/GammaProcessed1/hmm_out.mat'])
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');
%%
ModeHMM.WAKEstate = WAKEall;
ModeHMM.NREMstate = NREMall;

numcells = length(WAKEall);
spkthresh = 50;
MeanReturn.logbins = linspace(-3,2,50);
%get next ISI (nan for last one in the state)
%Cat all the cells
for ss = 1:2
for cc = 1:numcells

ModeHMM.(states{ss})(cc).next_isi = cellfun(@(prevISI) [prevISI(2:end) nan],ModeHMM.(states{ss})(cc).state_isi,'UniformOutput',false);
ModeHMM.(states{ss})(cc).next_state = cellfun(@(prevState) [prevState(2:end) nan],ModeHMM.(states{ss})(cc).decoded_mode,'UniformOutput',false);

ModeHMM.(states{ss})(cc).prev_isi = cellfun(@(prevISI) prevISI(1:end),ModeHMM.(states{ss})(cc).state_isi,'UniformOutput',false);
ModeHMM.(states{ss})(cc).prev_state = cellfun(@(prevState) prevState(1:end),ModeHMM.(states{ss})(cc).decoded_mode,'UniformOutput',false);

ModeHMM.(states{ss})(cc).prev_isi = cat(2,ModeHMM.(states{ss})(cc).prev_isi{:});
ModeHMM.(states{ss})(cc).next_isi = cat(2,ModeHMM.(states{ss})(cc).next_isi{:});

ModeHMM.(states{ss})(cc).prev_state = cat(2,ModeHMM.(states{ss})(cc).prev_state{:});
ModeHMM.(states{ss})(cc).next_state = cat(2,ModeHMM.(states{ss})(cc).next_state{:});
ModeHMM.(states{ss})(cc).state_spk = cat(2,ModeHMM.(states{ss})(cc).state_spk{:});


    
end

end

%%
%ignorepairs = false(numcells.*7);
ss = 1;
for sm = 1:6
    CellClass.celltypes{sm} = ['Mode',num2str(sm)];
    CellClass.(CellClass.celltypes{sm}) = false(numcells.*6);
    for cc = 1:numcells
        cellidx = cc + (sm-1).*numcells;
        spikes_modes.UID(cellidx) = cellidx;
        if sm == 6
            inmode = ModeHMM.(states{ss})(cc).next_state == sm & ModeHMM.(states{ss})(cc).prev_state==sm;
%         elseif sm == 7 %All AS spikes
%              inmode = ~(ModeHMM.(states{ss})(cc).next_state == 6 & ModeHMM.(states{ss})(cc).prev_state==6);
        else
            inmode = ModeHMM.(states{ss})(cc).next_state == sm | ModeHMM.(states{ss})(cc).prev_state==sm;
         end
        spikes_modes.times{cellidx} = ModeHMM.(states{ss})(cc).state_spk(inmode)';
        CellClass.label{cellidx} = ['Mode',num2str(sm)];
        CellClass.(CellClass.celltypes{sm})(cellidx) = true;
%         for sm2 = 1:6
%             cellidx2 = cc + (sm2-1).*numcells;
%             ignorepairs(cellidx,cellidx2)=true;
%         end
    end
    
end
spikes_modes.numcells = numcells.*6;

%% Load the LFP in this spike group
downsamplefactor = 2;
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);




%% Spike coupling with wavelets

%Take only subset of time (random intervals) so wavelets doesn't break
%computer
usetime = 4000;
winsize = 25;
if sum(diff(SleepState.ints.WAKEstate,1,2))>usetime
    nwin = round(usetime./winsize);
    %winsize = 30; %s
    try
        windows = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,winsize,nwin );
    catch
        windows = SleepState.ints.WAKEstate;
    end
else
    windows = SleepState.ints.WAKEstate;
end
%%
subpop = CellClass.Mode1 + CellClass.Mode2.*2 + CellClass.Mode3.*3 + ...
    CellClass.Mode4.*4 + CellClass.Mode5.*5 + CellClass.Mode6.*6;
%Calculate pop-phase coupling for all channels
[SpikeLFPCoupling] = ...
    bz_GenSpikeLFPCoupling(spikes_modes,lfp,...
    'int',windows,'frange',[1 312],'ncyc',12,...
    'cellclass',CellClass.label,'synchwin',0.002,'synchdt',0.002,...
    'nfreqs',150,'ISIpower',false,'spikeLim',5000);
    %close all

%%
modenames = {'AS1','AS2','AS3','AS4','AS5','GS'};

morder = [2 4 1 5 3 6];
GScolor = [0.6 0.4 0];
modecolors = crameri('bamako',5);
%modecolors = [modecolors;GScolor];
modecolors = {modecolors(1,:),modecolors(2,:),modecolors(3,:),modecolors(4,:),modecolors(5,:),...
    GScolor};

figure
subplot(2,2,1)
    for sm = 1:6
        hold on
    plot(log10(SpikeLFPCoupling.freqs),SpikeLFPCoupling.pop.(CellClass.celltypes{morder(sm)}).phasemag,...
        'color',modecolors{morder(sm)},'linewidth',2)

    end
    LogScale('x',10)
    xlabel('f (hz)');ylabel('Spike-Phase Coupling')
    
subplot(2,2,2)
    for sm = 1:6
        hold on
    plot(log10(SpikeLFPCoupling.freqs),SpikeLFPCoupling.pop.(CellClass.celltypes{morder(sm)}).powercorr  ,...
        'color',modecolors{morder(sm)},'linewidth',2)

    end
    LogScale('x',10)
    xlabel('f (hz)');ylabel('Spike-Phase Coupling')

NiceSave('GSLFPCoupling',figfolder,baseName);

end
