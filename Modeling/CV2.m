    %%
    [X Y] = meshgrid(ISIstats.ISIhist.logbins,ISIstats.ISIhist.logbins);
    CV2test = 2.*abs(10.^(X)-10.^(Y))./(10.^(X)+10.^(Y));
    
bwcolormap = [makeColorMap([0 0 0.2],[0 0 0.9],[1 1 1]);makeColorMap([1 1 1],[0.9 0 0],[0.2 0 0])];
    
    figure
    subplot(3,3,1)
    imagesc(ISIstats.ISIhist.logbins,ISIstats.ISIhist.logbins,CV2test')
    axis xy
        colorbar
        colormap(bwcolormap)
        caxis([0 2])
        title('CV2')
    caxis([0 2])
    xlabel('ISI_n');ylabel('ISI_n_+_1')
    LogScale('xy',10)
    
    NiceSave('CV2map',figfolder,[]);