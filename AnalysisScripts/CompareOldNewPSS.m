basePath = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/Dataset/OnCluster/171209_WT_EM1M3';
baseName = bz_BasenameFromBasepath(basePath);

figfolder = '/home/dlevenstein/ProjectRepos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/OptimizePSS';
    load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.mat']));
    load(fullfile(basePath,[baseName,'.PowerSpectrumSlope.lfp.mat']));
    
    %%
    specslope.resid(specslope.resid<0) = 0;
    PSpecSlope.Shortwin.OSCI(PSpecSlope.Shortwin.OSCI<0) = 0;
%%
wind = bz_RandomWindowInIntervals(specslope.timestamps([1 end]),20);
exchan = 20;
figure
subplot(4,1,1)
plot(specslope.timestamps,specslope.data(:,exchan),'k')
hold on
plot(PSpecSlope.Shortwin.timestamps,PSpecSlope.Shortwin.PSS(:,exchan),'r')
xlim(wind)


subplot(4,1,2)
plot(specslope.timestamps,bz_NormToRange(specslope.data(:,exchan),[0 1]),'k')
hold on
plot(PSpecSlope.Shortwin.timestamps,bz_NormToRange(PSpecSlope.Shortwin.PSS(:,exchan),[0 1]),'r')
xlim(wind)

subplot(4,1,3)
imagesc(specslope.timestamps,log10(specslope.freqs),specslope.specgram(:,:,exchan)')
hold on
plot(specslope.timestamps,bz_NormToRange(specslope.data(:,exchan),log10(specslope.freqs([1 end]))),'k')
xlim(wind)
%colorbar
axis xy
SpecColorRange( specslope.specgram(:,:,exchan) )

subplot(4,1,4)
imagesc(PSpecSlope.Shortwin.timestamps,(PSpecSlope.Shortwin.freqs),log10(PSpecSlope.Shortwin.SPEC(:,:,exchan)))
hold on
plot(PSpecSlope.Shortwin.timestamps,bz_NormToRange(PSpecSlope.Shortwin.PSS(:,exchan),(PSpecSlope.Shortwin.freqs([1 end]))),'r')
xlim(wind)
%colorbar
axis xy
NiceSave('OldNew',figfolder,baseName)

%subplot(PSpecSlope.Shortwin.movingwin

%%
wind = bz_RandomWindowInIntervals(specslope.timestamps([1 end]),20);

figure
subplot(4,1,1)
imagesc(specslope.timestamps,log10(specslope.freqs),specslope.specgram(:,:,exchan)')
hold on
plot(specslope.timestamps,bz_NormToRange(specslope.data(:,exchan),log10(specslope.freqs([1 end]))),'k')
xlim(wind)
%colorbar
SpecColorRange( specslope.specgram(:,:,exchan) )
axis xy

subplot(4,1,2)
imagesc(PSpecSlope.Shortwin.timestamps,(PSpecSlope.Shortwin.freqs),log10(PSpecSlope.Shortwin.SPEC(:,:,exchan)))
hold on
plot(PSpecSlope.Shortwin.timestamps,bz_NormToRange(PSpecSlope.Shortwin.PSS(:,exchan),(PSpecSlope.Shortwin.freqs([1 end]))),'r')
xlim(wind)
%colorbar
axis xy


subplot(4,1,3)
imagesc(specslope.timestamps,log10(specslope.freqs),specslope.resid(:,:,exchan)')
hold on
%plot(specslope.timestamps,bz_NormToRange(specslope.data(:,exchan),log10(specslope.freqs([1 end]))),'k')
xlim(wind)
caxis([0 1])
%colorbar
axis xy

subplot(4,1,4)
imagesc(PSpecSlope.Shortwin.timestamps,(PSpecSlope.Shortwin.freqs),log10(PSpecSlope.Shortwin.OSCI(:,:,exchan)))
hold on
%plot(PSpecSlope.Shortwin.timestamps,bz_NormToRange(PSpecSlope.Shortwin.PSS(:,exchan),(PSpecSlope.Shortwin.freqs([1 end]))),'r')
xlim(wind)
%colorbar
caxis([0 7])
axis xy
NiceSave('OldNew_OSCI',figfolder,baseName)


%%

wind = bz_RandomWindowInIntervals(specslope.timestamps([1 end]),25);

fftidx_new = find(specslope.timestamps>wind(1),1);
fftidx_old = find(PSpecSlope.Shortwin.timestamps>wind(1),1);


figure
subplot(4,1,1)
imagesc(specslope.timestamps,log10(specslope.freqs),specslope.specgram(:,:,exchan)')
hold on
plot(specslope.timestamps,bz_NormToRange(specslope.data(:,exchan),log10(specslope.freqs([1 end]))),'k')
xlim(wind)
colorbar
axis xy
SpecColorRange( specslope.specgram(:,:,exchan) )

subplot(4,1,2)
imagesc(PSpecSlope.Shortwin.timestamps,(PSpecSlope.Shortwin.freqs),log10(PSpecSlope.Shortwin.SPEC(:,:,exchan)))
hold on
plot(PSpecSlope.Shortwin.timestamps,bz_NormToRange(PSpecSlope.Shortwin.PSS(:,exchan),(PSpecSlope.Shortwin.freqs([1 end]))),'r')
xlim(wind)
colorbar
axis xy


%INTERCEPT + CONS?!

subplot(3,4,10)
plot(log2(specslope.freqs),(specslope.specgram(fftidx_new,:,exchan)),'k')
hold on
plot(log2(specslope.freqs),...
    log10(specslope.freqs).*specslope.data(fftidx_new)+specslope.intercept(fftidx_new,exchan))
%plot(log10(StatePlotMaterials.swFFTfreqs),(StatePlotMaterials.IRASAsmooth(:,fftidx)))
axis tight
box off
LogScale('x',2)
%ylim([3 5])
xlabel('f (Hz)');ylabel('log(power)')

subplot(3,4,11)
plot(log2(PSpecSlope.Shortwin.freqs),log10(PSpecSlope.Shortwin.SPEC(:,fftidx_old,exchan)),'k')
hold on
plot(log2(PSpecSlope.Shortwin.freqs),log10(PSpecSlope.Shortwin.FRAC(:,fftidx_old,exchan)),'r')
%plot(log2(PSpecSlope.Shortwin.freqs),...
%    log10(PSpecSlope.Shortwin.freqs).*PSpecSlope.Shortwin.PSS(fftidx_new)+PSpecSlope.Shortwin.intercept(fftidx_new,exchan))
%plot(log10(StatePlotMaterials.swFFTfreqs),(StatePlotMaterials.IRASAsmooth(:,fftidx)))
axis tight
box off
LogScale('x',2)
%ylim([3 5])
xlabel('f (Hz)');ylabel('log(power)')

NiceSave('OldNew_Fit',figfolder,baseName)
