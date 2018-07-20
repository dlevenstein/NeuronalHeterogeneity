function [ LogACG ] = bz_LogACG( spiketimes,varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%
%   INPUTS
%       spikes   a cell array: each cell is the spike times of a neuron  
%
%   (options)
%       'ints'        A structure with one or more intervals in which to 
%                       calculate the LogACG.
%                       states.stateNAME = [start stop]
%                       Will calculate LogACG separately for each state
%                       (Can also 'load' from SleepState.states.mat)
%       'savecellinfo'  logical (default: false) save a cellinfo file?
%       'basePath'
%       'numbins'       default: 40
%
%DLevenstein 2018
%% Parse the inputs
defaultstates.ALL = [-Inf Inf];

% parse args
p = inputParser;
addParameter(p,'ints',defaultstates)
addParameter(p,'savecellinfo',false,@islogical)
addParameter(p,'basePath',pwd,@isstr)
addParameter(p,'numbins',40)

parse(p,varargin{:})
ints = p.Results.ints;
SAVECELLINFO = p.Results.savecellinfo;
basePath = p.Results.basePath;
nbins = p.Results.numbins;


%% Load the stuff
baseName = bz_BasenameFromBasepath(basePath);
cellinfofilename = fullfile(basePath,[baseName,'.LogACG.cellinfo.mat']);

if exist(cellinfofilename,'file') && ~forceRedetect
    LogACG = bz_LoadCellinfo(basePath,'LogACG');
    return
end

if strcmp(ints,'load')
    ints = bz_LoadStates(basePath,'SleepState');
    ints = ints.ints;
end
statenames = fieldnames(ints);
numstates = length(statenames);
%% Calculate the LogACG
numcells = length(spiketimes);

%Bins and stuff
logrange = [-2.5 1.5];
spiketimeresolution = 5e-05;  %s, Assuming 20khz wideband.... make this automatic.
tedges_log = logspace(logrange(1),logrange(2),nbins);
tedges_log = [0 tedges_log];
t_log = log10(tedges_log(2:end))-0.5.*diff(log10(tedges_log(2:3)));

for ss = 1:numstates
    %Get only the spikes in the current state, shift spike times to merge state intervals
    instatespiketimes = cellfun(@(X) Restrict(X,ints.(statenames{ss}),'shift','on'),...
    spiketimes,'UniformOutput',false);
    
    %Calculate the acg for each cell
    for cc = 1:numcells
    [linacg(:,cc),t] = CCG(instatespiketimes(cc),...
        [],'binSize',spiketimeresolution,'duration',70); %should calculate duration from logrange...
    end

    %Put things in the log bins
    posbins = t>0;
    counts_lin = linacg(posbins,:);
    t_lin = t(posbins);
    %Find the bins everyone should go in
    for tt = 1:nbins   
        inbins = t_lin>tedges_log(tt)&t_lin<=tedges_log(tt+1);
        counts_log(tt,:) = sum(counts_lin(inbins,:),1);
        binSize_log(tt) = sum(inbins).*spiketimeresolution;
    end
     
    %Normalize by duration of each bin and number of reference spikes
    for cc = 1:numcells
        numREFspikes = length(instatespiketimes{cc});%number of reference events for group
        counts_log(:,cc) = counts_log(:,cc)./numREFspikes./binSize_log';
    end

    LogACG.(statenames{ss}).acg = counts_log;
end

%% The output
if numstates==1
    LogACG = LogACG.(statenames{ss});
end
LogACG.t = t_log;

if SAVECELLINFO
    save(cellinfofilename,'LogACG')
end
%% Figure
% figure
% subplot(2,2,1)
% imagesc(t_log,[1 spikes.numcells], log10(counts_log(:,ISIStats.sorts.ALL.ratebyclass))')
% colorbar
% LogScale('x',10)
% 
% subplot(2,2,2)
% plot(t_log,counts_log(:,5))
% LogScale('x',10)
end

