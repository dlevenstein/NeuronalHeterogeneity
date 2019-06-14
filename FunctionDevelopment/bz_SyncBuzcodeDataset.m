function [ success ] = bz_SyncBuzcodeDataset( fromFolder,toFolder )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

rsynccommand = ['rsync -auvzP --exclude=*dat --exclude=*.clu* --exclude=*.fet.* ',...
    '--exclude=*.res.* --exclude=*.spk.* --exclude=.phy --exclude=*.npy'];

end

