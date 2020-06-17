function [MIMap] = ISILFPMap(basePath,figfolder)
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



%%
%Regions = unique(sessionInfo.region(~cellfun(@isempty,sessionInfo.region)));
Regions = unique(spikes.region(~cellfun(@isempty,spikes.region)));


%% Load the LFP
for rr = 1:length(Regions)
inregionchanIDX = ismember(sessionInfo.region,Regions{rr});
inregionchan = sessionInfo.channels(inregionchanIDX);
downsamplefactor = 2;
lfp = bz_GetLFP(inregionchan,... %note: may have to load separately here for RAM...
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
dt = 1;
winsize = 1;
frange = [1 312];
nfreqs = 150;
[specslope] = bz_PowerSpectrumSlope(lfp,winsize,dt,'spectype','fft',...
    'nfreqs',nfreqs,'showfig',false,'showprogress',false,'frange',frange);%,...
    %'saveMat',basePath,'saveName',['Chan',num2str(lfpchannel)],...
    %'saveFolder','WavPSS');
    

%%
inregioncellIDX = ismember(spikes.region,Regions{rr});
MIMap.(Regions{rr}).UIDs = spikes.UID(inregioncellIDX);
MIMap.(Regions{rr}).ChanID = inregionchan;
MIMap.freqs = specslope.freqs;
%%
for tt = 1:length(celltypes)
    popspikes.(Regions{rr}).(celltypes{tt}) = cat(1,spikes.times{CellClass.(celltypes{tt})&inregioncellIDX});
end

%%
clear fPower

%cc = 1; %channel

fPower.timestamps = specslope.timestamps;
%fPower(length(specslope.freqs)+1).data = [];
clear PopConditionalISIDist PopConditionalISIDist_phase PopConditionalISIDist_power
for ff = 1:length(specslope.freqs)
	bz_Counter(ff,length(specslope.freqs),'Freq')
    for cc = 1:length(inregionchan)
        fPower.data = specslope.resid(:,ff,cc);
        for tt = 1:length(celltypes)
            for ss = 1:3
            [PopConditionalISIDist] = bz_ConditionalISI(popspikes.(Regions{rr}).(celltypes{tt}),fPower,...
                'ints',SleepState.ints.(states{ss}),...
                'showfig',false,'ISIDist',false);
            %test(rr,ss,tt,ff,cc) = PopConditionalISIDist.MutInf;
            MIMap.(Regions{rr}).(states{ss}).(celltypes{tt})(ff,cc) = PopConditionalISIDist.MutInf;
            end
        end
    end
    
end

%%
%rr = 2;
[~,SGsIDX] = cellfun(@(SGAll) (ismember(SGAll,inregionchan)),...
                   sessionInfo.spikeGroups.groups,'UniformOutput',false); %Just cells in the right region
SGsThisRegion = cellfun(@(SGAll) SGAll(ismember(SGAll,inregionchan)),...
                   sessionInfo.spikeGroups.groups,'UniformOutput',false); %Just cells in the right region
SGLength = cellfun(@(SG) length(SG),SGsThisRegion);
SGnum = find(SGLength>0);
SGorder = [SGsIDX{:}]; %The spikegroup order of all cells 
SGorder = SGorder(SGorder~=0);
%SGsort = 
%Sort by group
%save MI map - make sure it has channel numbers and UID of cells in each
%region/population. and freqs. and spike group ID/sortings.
%Save figures in detection figures
MIMap.(Regions{rr}).SGorder;
MIMap.(Regions{rr}).SGLength;
MIMap.(Regions{rr}).SGnum;
%%
figure
    for ss = 1:3
        for tt = 1:length(celltypes)
            subplot(2,3,ss+(tt-1).*3)
                imagesc(log2(MIMap.freqs),[1 length(MIMap.(Regions{rr}).ChanID)],...
                    MIMap.(Regions{rr}).(states{ss}).(celltypes{tt})(:,SGorder)')
                hold on
                plot(log2(MIMap.freqs([1 end])),cumsum(SGLength)'*[1 1]+0.5,'w')
                
                set(gca,'YTickLabel',[]);
                if ss == 1
                    ylabel({(celltypes{tt}),'Spike Group'})
                    set(gca,'YTick',cumsum(SGLength(SGnum))-0.5.*SGLength(SGnum)+0.5)
                    set(gca,'YTickLabel',SGnum)
                end
                if tt == 1
                    title(states{ss})
                end
                if tt == length(celltypes)
                    xlabel('f (Hz)');
                    LogScale('x',2)
                else
                    set(gca,'XTickLabel',[]);
                end
        end
    end
    
    
NiceSave(['ISIMod_',Regions{rr}],figfolder,baseName)
end

end