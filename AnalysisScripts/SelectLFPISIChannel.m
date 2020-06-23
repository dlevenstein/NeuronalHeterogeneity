function [  ] = SelectLFPISIChannel( baseName,regionName_in,regionName_out,whichstate )

% baseName = 'Rat11-20150330';
% regionName_in = 'bla';
% regionName_out = 'bla';
%%
%Load from '/Users/dl2820/Project Repos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/ISILFPMap'
if ~exist('regionName_out','var')
    regionName_out = regionName_in;
end
if ~exist('whichstate','var')
    whichstate = 'WAKEstate';
end

%%
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISILFPMap'];

[ISILFPMap] = GetMatResults(figfolder,'ISILFPMap','baseNames',baseName);
%%
sorting = ISILFPMap.MIMap.(regionName_in).SGorder;
heatmap = ISILFPMap.MIMap.(regionName_in).(whichstate).pE(:,sorting);
ChanID = ISILFPMap.MIMap.(regionName_in).ChanID(sorting);

%Use bz_ChannelSelect to pick a channel
[selectedchannel] = bz_ChannelSelect(heatmap,ChanID)

%save.... where? re-save in AnalysisResults. THEN, in analysis script, save
%as channel tag
%Ideally it gets to sessionInfo.mat on cluster in a
%channeltag
%Question: what do to with PIR and BLA?
%%


%%
filename = fullfile(figfolder,[baseName,'.Analysisresults.ISILFPMap.mat']);
load(filename)
MIMap.selectedchans.(regionName_out).channel = selectedchannel;
MIMap.selectedchans.(regionName_out).regname = regionName_in;
save(filename,'MIMap');
end