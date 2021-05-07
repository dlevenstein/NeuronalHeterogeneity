function [  ] = GGChansAndCells( ratPath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%ratPath = '/mnt/proraidDL/Database/GGData/Rat08';
[~,ratName] = fileparts(ratPath);

%SessionIndexingFile
% load('/mnt/proraidDL/Database/GGData/_extra/AllRats-Variables/sessionindexing.mat')
ratnumber = str2num(ratName(end-1:end));
% sessionindex = ratsessionindex(ratsessionindex(:,1)==ratnumber,2);

%SpikeParameters.  1:session, 2:shank, 3:cell,.... 10:MonoSynID
spikeparmsfile = fullfile(ratPath,[ratName,'-SpikeParameters.mat']);
load(spikeparmsfile);
try
    typefile = fullfile(ratPath,[ratName,'-finalType.mat']);
    load(typefile);
catch
    load('/mnt/proraidDL/Database/GGData/_extra/AllRats-Variables/AllRats-FinalType.mat')
    rat11 = finalType(:,1)==11;
    finalType = finalType(rat11,5);
end
shankIDs = SpikeParameters(:,2);
cluID = SpikeParameters(:,3);
sessionIDs = SpikeParameters(:,1);
%% Build the recording locations list
load([ratPath,'/RecordingLocations.mat'])
locationmap = [];
regionmap = [];
for rr = 2:5
    for reg = 1:11
        shanksetc = shankInfo{rr,reg};
        thislocationmap = [shankInfo{rr,12}.*ones(size(shanksetc,1),1),...
            shanksetc];%, repmat(shankInfo(1,reg),size(shanksetc,1),1)];
        locationmap = [locationmap ; thislocationmap];
        thisregionmap = repmat(shankInfo(1,reg),size(shanksetc,1),1);
        regionmap = [regionmap ; thisregionmap];
    end
end 
clear regionIDs
for cc = 1:length(sessionIDs)
    thiscell = locationmap(:,1)==ratnumber & ...
        locationmap(:,2)==sessionIDs(cc) & locationmap(:,3)==shankIDs(cc);
    if ~any(thiscell)
        %cc
        regionIDs{cc,1} = {''};
    else
        regionIDs{cc,1} = regionmap(thiscell);
    end
end
%%
baseNames = dir([ratPath,'/',ratName,'-*']);
baseNames = {baseNames([baseNames.isdir]).name};
basePaths = cellfun(@(X) fullfile(ratPath,X),baseNames,'UniformOutput',false);
%%

for bb = 1:length(basePaths)
    try
    display(['Loading ',baseNames{bb}])
    spikes = bz_GetSpikes('basepath',basePaths{bb});
    
    thissession = SpikeParameters(:,1) == (bb);
    clear ispE ispI isNone regions
    for cc = 1:spikes.numcells
        thiscell = (spikes.shankID(cc) == shankIDs & ...
            spikes.cluID(cc) == cluID & thissession);
        ispE(cc) = finalType(thiscell)==1;
        ispI(cc) = finalType(thiscell)==2;
        isNone(cc) = finalType(thiscell)==0;
        
        regions(cc) = regionIDs(thiscell);
%         if ~strcmp(thisregion{cc},spikes.region{cc})
%             keyboard
%         end
    end
        if ~any(ispE | ispI | isNone)
            error('cellsmissing')
        end
    %%
    %Pass through Cell classification for each cell
    bz_CellClassification(basePaths{bb},'keepKnown',true,...
    'knownE',spikes.UID(ispE),'knownI',spikes.UID(ispI),...
    'ignorecells',spikes.UID(isNone),'forceReload',true);

    %Save with the correct regins
    spikes.region = regions;
    save(spikes.filename,'spikes')
    catch
       display('ISSUE!') 
    end
    
    
end
%%


end

