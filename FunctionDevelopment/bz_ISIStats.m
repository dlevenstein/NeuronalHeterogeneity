function [ ISIstats ] = bz_ISIStats( spikes,varargin )
%ISIstats = bz_ISIStats(spikes,varargin) calculates the statistics 
%inter-spike intervals for the spiketimes in spikes.
%
%   INPUTS
%       spikes
%
%       (options)
%       'states'
%       'cellclasses'
%       'savecellinfo'
%       'basePath'
%       'figfolder'
%
%   OUTPUTS
%       
%
%DLevenstein
%% Parse the inputs

%% Dev
states = SleepState;
%%
statenames = {'NREMstate','REMstate','WAKEstate'};
numstates = length(statenames);


%% ISI and CV2 statistics
numcells = length(spikes.UID);

%Calculate ISI and CV2 for allspikes
spikes.ISIs = cellfun(@diff,spikes.times,'UniformOutput',false);
spikes.CV2 = cellfun(@(X) 2.*abs(X(2:end)-X(1:end-1))./(X(2:end)+X(1:end-1)),spikes.ISIs ,'UniformOutput',false);

%%


end

