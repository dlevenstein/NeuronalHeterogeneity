function [ success ] = bz_SyncBuzcodeDataset( fromFolder,toFolder )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

rsynccommand = ['rsync -auvzP --exclude=*dat --exclude=*.clu* --exclude=*.fet.* ',...
    '--exclude=*.res.* --exclude=*.spk.* --exclude=.phy --exclude=*.npy'];


% YS from server
 rsync -auvzP dl2820@bigpurple.nyumc.org:/gpfs/data/buzsakilab/DL/Database/YSData /mnt/proraidDL/Database/

%YS to server
rsync -auvzP --exclude=UnitSummary --exclude=revisions_cell_metrics --exclude=rez.mat /mnt/proraidDL/Database/YSData dl2820@bigpurple.nyumc.org:/gpfs/data/buzsakilab/DL/Database/


%AG to server 
rsync -auvzP --exclude=Old --exclude=Zip /mnt/proraidDL/Database/AGData dl2820@bigpurple.nyumc.org:/gpfs/data/buzsakilab/DL/Database/

%BW to server
rsync -auvzP --exclude=c3po --exclude=*_ACC /mnt/proraidDL/Database/BWData dl2820@bigpurple.nyumc.org:/gpfs/data/buzsakilab/DL/Database/


%AP from NYUshare (L/K to keep softlinks on send/recieve)
rsync -auvzPLK --exclude=*raw* --exclude=CluSAV --exclude=Analysis --exclude=*dat --exclude=*.clu* --exclude=*.fet.* --exclude=*.res.* --exclude=*.spk.* --exclude=*.m1m2.* --exclude=*.alg.* --exclude=*.mm.* /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onNYUShare/AP_THAL/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL/
%/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL/

%AP to cluster
rsync -auvzP /mnt/proraidDL/Database/APData /mnt/BigPurple/Database/

end

