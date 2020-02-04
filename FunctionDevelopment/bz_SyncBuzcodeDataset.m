function [ success ] = bz_SyncBuzcodeDataset( fromFolder,toFolder )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

rsynccommand = ['rsync -auvzP --exclude=*dat --exclude=*.clu* --exclude=*.fet.* ',...
    '--exclude=*.res.* --exclude=*.spk.* --exclude=.phy --exclude=*.npy'];


% YS from cluster
 rsync -auvzP dl2820@bigpurple.nyumc.org:/gpfs/data/buzsakilab/DL/Database/YSData /mnt/proraidDL/Database/

%YS to server
rsync -auvzP --exclude=UnitSummary --exclude=revisions_cell_metrics --exclude=rez.mat /mnt/proraidDL/Database/YSData dl2820@bigpurple.nyumc.org:/gpfs/data/buzsakilab/DL/Database/


%AG to server 
rsync -auvzP --exclude=Old --exclude=Zip /mnt/proraidDL/Database/AGData dl2820@bigpurple.nyumc.org:/gpfs/data/buzsakilab/DL/Database/

%BW to server
rsync -auvzP --exclude=c3po --exclude=*_ACC /mnt/proraidDL/Database/BWData dl2820@bigpurple.nyumc.org:/gpfs/data/buzsakilab/DL/Database/



%AP from NYUshare (L/K to keep softlinks on send/recieve)
rsync -auvzPLK --exclude=*raw* --exclude=CluSAV --exclude=Analysis --exclude=*dat --exclude=*.clu* --exclude=*.fet.* --exclude=*.res.* --exclude=*.spk.* --exclude=*.m1m2.* --exclude=*.alg.* --exclude=*.mm.* /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onNYUShare/AP_THAL/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL/


%GG first time from NYUshare
rsync -auvzPLK --exclude=Histology --exclude=Figures --exclude=Rat07 --exclude=*GLMoutput* --exclude=*.fil --exclude=*dat --exclude=*.clu* --exclude=*.fet.* --exclude=*.klg.* --exclude=*.res.* --exclude=*.spk.* --exclude=RawKK /mnt/NyuShare/Buzsakilabspace/Datasets/GirardeauG/ /mnt/proraidDL/Database/GGData/
%to cluster
rsync -auvzPLK --exclude=Histology --exclude=Figures --exclude=Rat07 --exclude=*GLMoutput* --exclude=*.fil --exclude=*dat --exclude=*.clu* --exclude=*.fet.* --exclude=*.klg.* --exclude=*.res.* --exclude=*.spk.* --exclude=RawKK /mnt/proraidDL/Database/GGData/ /mnt/BigPurple/Database/GGData/
%from cluster
rsync -auvzPLK --exclude=Histology --exclude=Figures --exclude=Rat07 --exclude=*GLMoutput* --exclude=*.fil --exclude=*dat --exclude=*.clu* --exclude=*.fet.* --exclude=*.klg.* --exclude=*.res.* --exclude=*.spk.* --exclude=RawKK /mnt/BigPurple/Database/GGData/ /mnt/proraidDL/Database/GGData/ 

%AP from cluster
rsync -auvzPLK /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onCluster/AP_THAL/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL/

%Everythign to cluster
rsync -auvzPLK /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onCluster/AP_THAL/
rsync -auvzPLK /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onCluster/BW_CTX/
rsync -auvzPLK /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onCluster/AG_HPC/
rsync -auvzPLK --exclude=UnitSummary --exclude=revisions_cell_metrics --exclude=rez.mat /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onCluster/YS_CTX/

%Everythign to NYUshare
rsync -auvzPLK --exclude=*.AnalysisResults.* /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onNYUShare/AP_THAL/
rsync -auvzPLK --exclude=*.AnalysisResults.* /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onNYUShare/BW_CTX/
rsync -auvzPLK --exclude=*.AnalysisResults.* /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onNYUShare/AG_HPC/
rsync -auvzPLK --exclude=*.AnalysisResults.* /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onNYUShare/YS_CTX/

%Everythign from cluster
rsync -auvzPLK /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onCluster/AP_THAL/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL/
rsync -auvzPLK /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onCluster/BW_CTX/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX/
rsync -auvzPLK /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onCluster/AG_HPC/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC/
rsync -auvzPLK /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onCluster/YS_CTX/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX/



%AP to NYUShare
rsync -auvzPLK /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL/ /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onNYUShare/AP_THAL/ 

%WM from NYUShare (ACH)
rsync -auvzP --exclude=*dat --exclude=Atropine --exclude=Kilosort* /mnt/NyuShare/dl2820/WMDataset/ /mnt/proraidDL/Database/WMData/AChPupil/


end

%FOr Roman
%rsync -auvzPLK --exclude=*.AnalysisResults.* /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX/ /mnt/NyuShare/Buzsakilabspace/LabShare/RomanHuszar/DATA/YSData_DL/
%rsync -auvzPLK --exclude=*.AnalysisResults.* /home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX/ /mnt/NyuShare/Buzsakilabspace/LabShare/RomanHuszar/DATA/BWData_DL/
