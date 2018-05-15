    %%
    [X Y] = meshgrid(ISIstats.ISIhist.logbins,ISIstats.ISIhist.logbins);
    CV2test = 2.*abs(10.^(X)-10.^(Y))./(10.^(X)+10.^(Y));
    
    figure
    imagesc(ISIstats.ISIhist.logbins,ISIstats.ISIhist.logbins,CV2test')
    axis xy
    colorbar
    title('CV2')
    xlabel('ISI_n');ylabel('ISI_n_+_1')
    LogScale('xy',10)