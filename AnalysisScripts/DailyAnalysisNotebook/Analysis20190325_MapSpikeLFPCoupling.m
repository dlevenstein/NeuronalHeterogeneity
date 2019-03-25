function [ ] = Analysis20190325(basePath,figfolder)
% Date XX/XX/20XX
%
%Goal: create a function to map pop-LFP coupling for every channel in a
%recording
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC/Achilles_10252013';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};

%% The necessary inputs
states = fieldnames(SleepState.ints);
numstates = length(states);
%% 
% loop spike groups, do each spike group in a row
%Get the region of each spike group. Only couple to cells in the region
spikeGroups = sessionInfo.spikeGroups;

%% get the region of the spikegroup
for gg = 1:length(spikeGroups.groups) 
    spikeGroups.region(gg) = unique(sessionInfo.region(ismember(sessionInfo.channels,spikeGroups.groups{gg})));
end



%%
regions = unique(spikeGroups.region);

AllCellClass = cell(1,spikes.numcells);
for rr = 1:length(regions)
    for cc = 1:length(celltypes)
        classname = [regions{rr},'_',celltypes{cc}];
        AllCellClass(CellClass.(celltypes{cc}) & strcmp(spikes.region,regions{rr})) = {classname};
    end
end
allclasses = unique(AllCellClass);
classcolors = {'k','r','k','r'};
classline = {'-','-','--','--'};
%%

for gg = 1:length(spikeGroups.groups)
    display(['Mapping Spike Group ',num2str(gg)])
    %% Load the LFP in this spike group
    downsamplefactor = 2;
    lfp = bz_GetLFP(spikeGroups.groups{gg},...
        'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

    %% Restrict to state
    for ss = 1:numstates
        state = states{ss};


        %% Spike coupling with wavelets

        %Take only subset of time (random intervals) so wavelets doesn't break
        %computer (total 625s)
        usetime = 3000;%2500
        winsize = 25;
        if sum(diff(SleepState.ints.(state),1,2))>usetime
            nwin = round(usetime./winsize);
            %winsize = 30; %s
            try
                windows = bz_RandomWindowInIntervals( SleepState.ints.(state),winsize,nwin );
            catch
                windows = SleepState.ints.(state);
            end
        else
            windows = SleepState.ints.(state);
        end

        %Calculate pop-phase coupling for all channels
        [SpikeLFPCoupling(gg).(state)] = ...
            bz_GenSpikeLFPCoupling(spikes,lfp,'channel',spikeGroups.groups{gg},...
            'int',windows,'frange',[1 312],'ncyc',15,...
            'cellclass',AllCellClass,'synchwin',0.002,'synchdt',0.002,...
            'nfreqs',150,'ISIpower',true,'spikeLim',20000);
            close all

    end
    
end

%% Saving...
saveswitch = false;

savename = fullfile(basePath,[baseName,'.SpikeLFPCoupling.chaninfo.mat']);

if saveswitch
    save(savename,'SpikeLFPCoupling');
end

%%
%Need: region for each spike group. Split populations by region and cell
%type
%Histogram of peak channel for cells (by population)

%Histogram of proportion of cells coupled to each frequency? (need
%significance....)

%Use this to pick gamma(/ripple?) channel for each cell
for gg = 1:length(spikeGroups.groups)
	spikeGroups.ingroupcellIDX{gg} = ismember(spikes.maxWaveformCh,spikeGroups.groups{gg});
    spikeGroups.inregioncellIDX{gg} = ismember(spikes.region,spikeGroups.region{gg});    
    
    for ss = 1:numstates
        state = states{ss};
       SpikeLFPCoupling(gg).(state).ingroupmean.ratepowercorr = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ratepowercorr(spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outgroupmean.ratepowercorr = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ratepowercorr(...
           spikeGroups.inregioncellIDX{gg}&~spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outregionmean.ratepowercorr = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ratepowercorr(...
           ~spikeGroups.inregioncellIDX{gg},:,:),1));   

       SpikeLFPCoupling(gg).(state).ingroupmean.spikephasemag = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.spikephasemag(spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outgroupmean.spikephasemag = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.spikephasemag(...
           spikeGroups.inregioncellIDX{gg}&~spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outregionmean.spikephasemag = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.spikephasemag(...
           ~spikeGroups.inregioncellIDX{gg},:,:),1));

       SpikeLFPCoupling(gg).(state).ingroupmean.ISIpowermodulation = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ISIpowermodulation(spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outgroupmean.ISIpowermodulation = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ISIpowermodulation(...
           spikeGroups.inregioncellIDX{gg}&~spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outregionmean.ISIpowermodulation = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ISIpowermodulation(...
           ~spikeGroups.inregioncellIDX{gg},:,:),1));

       for tt = 1:2
           SpikeLFPCoupling(gg).(state).pop.(celltypes{tt}).ISIpowermodulation = ...
               squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ISIpowermodulation(CellClass.(celltypes{tt}),:,:),1));
       end
       
       for tt = 1:length(allclasses)
           SpikeLFPCoupling(gg).(state).pop.(allclasses{tt}).ISIpowermodulation = ...
               squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ISIpowermodulation(...
               ismember(AllCellClass,allclasses{tt}),:,:),1));
       end
    end
end


%% Try Collapsing

%SpikLFPCollapsed = bz_CollapseStruct(SpikeLFPCoupling,'match','justcat',true);


%% R vs L

order = [1 3 2 4];
for gg = 1:length(spikeGroups.groups)
figure
suptitle(['SpikeGroup ',num2str(gg),': ',spikeGroups.region{gg}])
for ss = 1:length(states)
for tt = 1:length(allclasses)
    subplot(4,3,ss+(order(tt)*3)-3)
        imagesc(log2(SpikeLFPCoupling(gg).(states{ss}).freqs),[0 1],...
            SpikeLFPCoupling(gg).(states{ss}).pop.(allclasses{tt}).ISIpowermodulation')
        if ss == 1
            ylabel(allclasses{tt})
        end
        LogScale('x',2)
        if tt ==4
            xlabel('f (Hz)')
        end
        colorbar
        if tt==1
            title(states{ss})
        end
end
end

NiceSave(['ISIPowerModulation_',num2str(gg)],figfolder,baseName,'includeDate',true)


end

%%
state = 'WAKEstate';
gg =4 ;
figure

for ss = 1:length(states)
for tt = 1:length(CellClass.celltypes)
    subplot(3,2,tt+(ss*2)-2)
        imagesc(log2(SpikeLFPCoupling(gg).(states{ss}).freqs),[0 1],...
            SpikeLFPCoupling(gg).(states{ss}).pop.(celltypes{tt}).ISIpowermodulation')
        title(celltypes{tt})
        LogScale('x',2)
        colorbar
        if tt==1
            ylabel(states{ss})
        end
end
end




end
