    %%
    logbins = linspace(-3,2,50);
    [X Y] = meshgrid(logbins,logbins);
    CV2test = 2.*abs(10.^(X)-10.^(Y))./(10.^(X)+10.^(Y));
    logCV2test = 2.*abs((X)-(Y))./((X)+(Y));
    
bwcolormap = [makeColorMap([0 0 0.2],[0 0 0.9],[1 1 1]);makeColorMap([1 1 1],[0.9 0 0],[0.2 0 0])];
    
    figure
    subplot(3,3,1)
        imagesc(logbins,logbins,CV2test')
        axis xy
            colorbar
            colormap(bwcolormap)
            caxis([0 2])
            title('CV2')
        caxis([0 2])
        xlabel('ISI_n');ylabel('ISI_n_+_1')
        LogScale('xy',10)
    

    subplot(3,3,2)
        imagesc(logbins,logbins,logCV2test')
        axis xy
            colorbar
            colormap(bwcolormap)
            caxis([0 2])
            title('logCV2')
        caxis([0 2])
        xlabel('ISI_n');ylabel('ISI_n_+_1')
        LogScale('xy',10)
    
    %NiceSave('CV2map',figfolder,[]);