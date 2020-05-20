function [PopConditional,PopConditional_phase] = ISIModulationAnalysis(basePath,figfolder)
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
%reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
% basePath = '/Users/dl2820/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
% %basePath = pwd;
% %basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {[0 0 0],[0 0 1],[1 0 0],[0.6 0.6 0.6]};

try
celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};



    
%% Load the LFP

downsamplefactor = 2;
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%%

% ss = 1;
% state = states{ss};
% %Take only subset of time (random intervals) so wavelets doesn't break
% %computer (total 625s)
% usetime = 7200;%2500
% winsize = 25;
% if sum(diff(SleepState.ints.(state),1,2))>usetime
%     nwin = round(usetime./winsize);
%     %winsize = 30; %s
%     try
%         windows = bz_RandomWindowInIntervals( SleepState.ints.(state),winsize,nwin );
%     catch
%         windows = SleepState.ints.(state);
%     end
% else
%     windows = SleepState.ints.(state);
% end
        
%% Calculate Residuals
dt = 0.01;
winsize = 10;
frange = [1 312];
nfreqs = 150;
[specslope] = bz_PowerSpectrumSlope(lfp,winsize,dt,'spectype','wavelet',...
    'nfreqs',nfreqs,'showfig',true,'showprogress',true,'frange',frange,...
    'saveMat',basePath,'saveName',['Chan',num2str(lfp.channels)],...
    'saveFolder','WavPSS');
    

%%
for tt = 1:length(celltypes)
    popspikes.(celltypes{tt}) = cat(1,spikes.times{CellClass.(celltypes{tt})});
end

%%
clear fPower


fPower.timestamps = specslope.timestamps;
%fPower(length(specslope.freqs)+1).data = [];
clear PopConditionalISIDist PopConditionalISIDist_phase
for ff = 1:length(specslope.freqs)
	bz_Counter(ff,length(specslope.freqs),'Freq')
    fPower.data = specslope.resid(:,ff);
    for tt = 1:length(celltypes)
        for ss = 1:3
        [PopConditionalISIDist.(states{ss}).(celltypes{tt})(ff)] = ConditionalISI(popspikes.(celltypes{tt}),fPower,...
            'ints',SleepState.ints.(states{ss}),...
            'showfig',false,'ISIDist',false);
        end
    end
    fPower.data = specslope.phase(:,ff);
    for tt = 1:length(celltypes)
        for ss = 1:3
        [PopConditionalISIDist_phase.(states{ss}).(celltypes{tt})(ff)] = ConditionalISI(popspikes.(celltypes{tt}),fPower,...
            'ints',SleepState.ints.(states{ss}),...
            'showfig',false,'normtype','none','Xwin',[-pi pi],'ISIDist',false);
        end
    end
end



%%
PopConditional = bz_CollapseStruct(PopConditionalISIDist,3,'justcat',true);
PopConditional_phase = bz_CollapseStruct(PopConditionalISIDist_phase,3,'justcat',true);

%%

figure
for tt = 1:length(celltypes)
subplot(2,2,1+(tt-1)*2)
hold on
for ss = 1:2
plot(log10(specslope.freqs),squeeze(PopConditional.(states{ss}).(celltypes{tt}).MutInf),'color',statecolors{ss},'linewidth',1)
end
LogScale('x',10)
xlabel('f (Hz)')
ylabel({(celltypes{tt}),'MI[ISI;Power]'})
title('ISI-Power Modulation')


subplot(2,2,2+(tt-1)*2)
hold on
for ss = 1:2
plot(log10(specslope.freqs),squeeze(PopConditional_phase.(states{ss}).(celltypes{tt}).MutInf),'color',statecolors{ss},'linewidth',1)
end
LogScale('x',10)
xlabel('f (Hz)')
ylabel('MI[ISI;Phase]')
title('ISI-Phase Modulation')
end

NiceSave('PopISIMod',figfolder,baseName)

end
