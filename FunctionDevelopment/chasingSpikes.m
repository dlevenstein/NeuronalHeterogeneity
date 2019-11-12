
function [spikes, spikesGroupAll] = chasingSpikes(d, varargin)
% Detect action potential in intracell traces
%
% USAGE
%
%   [spikes] = chasingSpikes(d, varargin)
%
% INPUTS
%
% d              - R x C matrix with intracell data. C is data, and R should
%               be greater than 1, where: R(1) is intracell signal and R(2)
%               is current signal. Following rows will be ignored.
% fs             - Sampling frequency, default 20000
% win            - Window size before/after spikes in seconds, default [0.002 0.02]
% peakProminence - Prominence for spike detection in mV, default 10.
% cellPosition   - R x C matrix indicating segment of d (in s) for N intracell
%               recording, with the format: [start_cell_1, end_cell_1; 
%               start_cell2, end_cell2; ...].
%
% OUTPUTS
%
% spikes         - timstamp in seconds
% spikesGroup    - R x C matrix containing all spikes.
%
%   Manu Valero 2018

% Parse options
p = inputParser;
addParameter(p,'fs',20000,@isnumeric);
addParameter(p,'peakProminence',10,@isnumeric);
addParameter(p,'cellPosition',[],@ismatrix);
addParameter(p,'win',[0.002 0.002],@isvector);

parse(p,varargin{:});

win = p.Results.win;
fs = p.Results.fs;
peakProminence = p.Results.peakProminence;
cellPos = p.Results.cellPosition;
if isempty(cellPos)
    cellPos = [1 numel(d(1,:))/fs];
end
hp = 10;
fprintf('%iHz high pass filtering... \n',hp);
hpFilt = designfilt('highpassiir','FilterOrder',8, 'PassbandFrequency',hp,'PassbandRipple',0.1, 'SampleRate',fs);
d_filt = filtfilt(hpFilt,d(1,:));
[peaks,locs]=findpeaks(d_filt,'MinPeakProminence',peakProminence);
win = win * fs;
cellPos = cellPos * fs;

ii = 1;
while ii <= size(cellPos,1)
    locsTemp = locs(locs > cellPos(ii,1) & locs < cellPos(ii,2));
    clear spikesGroup max_diff
    for jj = 1:length(locsTemp)
        if locsTemp(jj)-win(1)>0 && locsTemp(jj)+win(2)<length(d_filt)
            spikesGroup(:,jj)=d_filt(int32((locsTemp(jj)-win(1)):(locsTemp(jj)+win(2))));
        end
    end
    clear max_diff
    for jj = 1:size(spikesGroup,2)
        max_diff(jj) = max(diff(spikesGroup(:,jj)));
    end
    
    % removing artifacts
    fig = figure;
    histogram(max_diff);
    hold on
    ax =axis;
    p1 = plot([10 10],[ax(3) ax(4)]);
    draggable(p1,'constraint','h');
    btn = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
                            'Units','normalize','Position', [.89 .92 .10 .06],...
                            'Callback', 'uiresume(gcbf)');
    disp('Press continue when done...');
    uiwait(gcf);
    pos=mean(get(p1,'XData'));
    close(fig);
    spikesGroup(:,find(max_diff>pos)) = []; locsTemp(find(max_diff>pos))=[];

    % removing fake spikes
    fig = figure;
    hold on
    plot(spikesGroup);
    ylim([-20 100]) %axis tight;
    ax=axis;
    l1=plot([ax(1) ax(2)],50*ones(2,1),'--','LineWidth',1.5);
    l2=plot([ax(1) ax(2)],20*ones(2,1),'--','LineWidth',1.5);
    l3=plot([35 35],[ax(3) ax(4)],'--','LineWidth',1.5);
    l4=plot([45 45],[ax(3) ax(4)],'--','LineWidth',1.5);
    draggable(l1,'constraint','v');
    draggable(l2,'constraint','v');
    draggable(l3,'constraint','h');
    draggable(l4,'constraint','h');
    btn = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
                            'Units','normalize','Position', [.89 .92 .10 .06],...
                            'Callback', 'uiresume(gcbf)');
    disp('Press continue when done...');
    uiwait(gcf);
    yl(1)=mean(get(l1,'YData'));
    yl(2)=mean(get(l2,'YData'));
    xl(1)=mean(get(l3,'XData'));
    xl(2)=mean(get(l4,'XData'));
    fill([ax(1) ax(2) ax(2) ax(1)],[yl(1) yl(1) ax(4) ax(4)],[1 0 0], ...
        'EdgeColor','none','FaceAlpha',.5);
    fill([xl(1) xl(2) xl(2) xl(1)],[ax(3) ax(3) yl(2) yl(2)],[1 0 0], ...
        'EdgeColor','none','FaceAlpha',.5);
    pause(2);
    close(fig);
    
    [~,er]=find(spikesGroup>yl(1)); % exclusion zone 1
    spikesGroup(:,unique(er))=[]; locsTemp(unique(er))=[];
    [~,er]=find(spikesGroup(int32(xl(1):xl(2)),:)<yl(2)); % exclusion zone 2
    spikesGroup(:,unique(er))=[]; locsTemp(unique(er))=[];
    xtspk=linspace(-win(1)/fs * 1000,win(2)/fs * 1000,size(spikesGroup,1));
    
    h = figure;
    hold on
    plot(xtspk,spikesGroup,'color',[.8 .8 .8]);
    plot(xtspk,mean(spikesGroup,2),'k');
    xlabel('ms');
    ylabel('mV');
    
    opt = input('Press ''r'' to repeat clustering, otherwise press any key... ','s');
    if ~strcmpi(opt, 'r')
        spikes{ii} =  locsTemp;
        spikesGroupAll{ii} = spikesGroup;
        ii = ii + 1;
    end
    try close(h); end
end


end