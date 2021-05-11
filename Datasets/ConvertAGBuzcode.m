function [ checkspikes ] = ConvertAGBuzcode( basePath,figfolder )
%This function converts the AG dataset from crcns to buzcode format


%% Get the cell type cellinfo
baseName = bz_BasenameFromBasepath(basePath);
load(fullfile(basePath,[baseName,'_sessInfo.mat']))

%%
allIDs = [sessInfo.Spikes.PyrIDs;sessInfo.Spikes.IntIDs];
shankID = floor(allIDs/100);
cluID = mod(allIDs, 100);

%% Spikes and sessionInfo
%Load/save spikes, in the process will prompt user to put in regions and
%bad channels - get the spike groups for regions from file 
%Channel_Orderings.pdf and channels for badchannel (red) 
%NOTE: must subtract 1 from the channel in the pdf, because neuroscope 0-indexing
%instead of 1-indexing!!!!!!
spikes = bz_GetSpikes('basePath',basePath,'saveMat',true,...
    'getWaveforms','force','forceReload',true,'onlyLoad',[shankID cluID]);


%%
newCellClass.UID = spikes.UID;
for cc = 1:spikes.numcells
   AGcellID(cc) = 100.*spikes.shankID(cc) + spikes.cluID(cc); 
   ispE(cc) = any(sessInfo.Spikes.PyrIDs==AGcellID(cc));
   ispI(cc) = any(sessInfo.Spikes.IntIDs==AGcellID(cc));
   isnone(cc) = ~any(ispE | ispI);
end
%%
% use bz_cellclassification to make the CellClass structure and create
% detectionfigure, but use Andres' classification
bz_CellClassification(basePath,'keepKnown',true,...
    'knownE',spikes.UID(ispE),'knownI',spikes.UID(ispI),'forceReload',true)

%%
%Check if number of spikes are equal to andres' structure
checkspikes = length(spikes.spindices) == length(sessInfo.Spikes.SpikeTimes)


%% behavior
behaviorfielname = fullfile(basePath,[baseName,'.position.behavior.mat']);

position.timestamps = sessInfo.Position.TimeStamps';
position.samplingRate = 1./mean(diff(position.timestamps));
position.position.x = sessInfo.Position.TwoDLocation(:,1);
position.position.y = sessInfo.Position.TwoDLocation(:,2);
position.position.lin = sessInfo.Position.OneDLocation;
position.units = 'm';
position.behaviorinfo.MazeType = sessInfo.Position.MazeType;
position.behaviorinfo.substructnames = 'position';
position.behaviorinfo.processingfunction = 'ConvertAGBuzcode';

position.Epochs.PREEpoch = sessInfo.Epochs.PREEpoch;
position.Epochs.MazeEpoch = sessInfo.Epochs.MazeEpoch;
position.Epochs.POSTEpoch = sessInfo.Epochs.POSTEpoch;

%%
figure
subplot(2,2,1)
    scatter(position.position.x,position.position.y,1,position.timestamps)
    xlabel('x');ylabel('y');
    title('Position')
    
subplot(2,2,2)
    scatter(position.timestamps,position.position.lin,1,position.timestamps)
    xlabel('time (s)');ylabel('linear position)');
    title('Linearized Position')
    axis tight
NiceSave('Position',figfolder,baseName)
%%
save(behaviorfielname,'position')


end

