reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DecodePositionAnalysis'];


[DecodeALL,baseNames] = GetMatResults(figfolder,'DecodePositionAnalysis');
%%
for rr = 1:length(DecodeALL)
    baseName = char(DecodeALL(rr).name)
    isTHAL(rr) = strcmp(baseName(1:5),'Mouse');
end
%%
regions = {'THAL','CA1'};
Decode.THAL = bz_CollapseStruct(DecodeALL(isTHAL));
Decode.CA1 = bz_CollapseStruct(DecodeALL(~isTHAL));

%%
for rr = 1:2
ISIbyDecPOS_norm.(regions{rr}) = bz_CollapseStruct(Decode.(regions{rr}).ISIbyDecPOS_norm,3,'justcat',true);
ISIbyDecPOS_normNREM.(regions{rr}) = bz_CollapseStruct(Decode.(regions{rr}).ISIbyDecPOS_normNREM,3,'justcat',true);
CellInfo.(regions{rr}) = bz_CollapseStruct(Decode.(regions{rr}).CellInfo,'match','justcat',true);

ISIbyDecPOS_norm_mean.(regions{rr}) = bz_CollapseStruct( Decode.(regions{rr}).ISIbyDecPOS_norm(CellInfo.(regions{rr}).keep),3,'mean',true);
ISIbyDecPOS_norm_meanNREM.(regions{rr}) = bz_CollapseStruct( Decode.(regions{rr}).ISIbyDecPOS_normNREM(CellInfo.(regions{rr}).keep),3,'mean',true);
end

%%
figure
for rr = 1:2
subplot(2,2,1+(rr-1)*2)
    imagesc(ISIbyDecPOS_norm_mean.(regions{rr}).Dist.Xbins,...
        ISIbyDecPOS_norm_mean.(regions{rr}).Dist.Ybins,...
        ISIbyDecPOS_norm_mean.(regions{rr}).Dist.pYX')
    hold on
    plot(ISIbyDecPOS_norm_mean.(regions{rr}).Dist.Xbins,-log10(ISIbyDecPOS_norm_mean.(regions{rr}).Dist.SpikeRate),'r')
    if rr == 1
        imagesc(ISIbyDecPOS_norm_mean.(regions{rr}).Dist.Xbins+2*pi,...
            ISIbyDecPOS_norm_mean.(regions{rr}).Dist.Ybins,...
            ISIbyDecPOS_norm_mean.(regions{rr}).Dist.pYX')
        imagesc(ISIbyDecPOS_norm_mean.(regions{rr}).Dist.Xbins-2*pi,...
            ISIbyDecPOS_norm_mean.(regions{rr}).Dist.Ybins,...
            ISIbyDecPOS_norm_mean.(regions{rr}).Dist.pYX')
        plot(ISIbyDecPOS_norm_mean.(regions{rr}).Dist.Xbins+2*pi,-log10(ISIbyDecPOS_norm_mean.(regions{rr}).Dist.SpikeRate),'r')
        plot(ISIbyDecPOS_norm_mean.(regions{rr}).Dist.Xbins-2*pi,-log10(ISIbyDecPOS_norm_mean.(regions{rr}).Dist.SpikeRate),'r')
    end
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    if rr==2
            xlabel('Decoded Position relative to PF Peak (m)')
            title('WAKE')
            xlim([-0.9 0.9]);
    elseif rr == 1
            xlabel('Decoded HD relative to Tuning Peak (m)')
            title('WAKE')
            bz_piTickLabel('x')
            xlim([-1.5*pi 1.5.*pi])
    end
    
subplot(2,2,2+(rr-1)*2)
    imagesc(ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.Xbins,...
        ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.Ybins,ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.pYX')
    hold on
    plot(ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.Xbins,-log10(ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.SpikeRate),'r')
    if rr == 1
        imagesc(ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.Xbins+2*pi,...
                ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.Ybins,ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.pYX')
        imagesc(ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.Xbins-2*pi,...
                ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.Ybins,ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.pYX')
        plot(ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.Xbins+2*pi,-log10(ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.SpikeRate),'r')
        plot(ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.Xbins-2*pi,-log10(ISIbyDecPOS_norm_meanNREM.(regions{rr}).Dist.SpikeRate),'r')

    end
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    if rr==2
            xlabel('Decoded Position relative to PF Peak (m)')
            title('NREM SWRs');xlim([-0.9 0.9]);
    elseif rr == 1
            xlabel('Decoded HD relative to Tuning Peak (m)')
            title('NREM')
            bz_piTickLabel('x')
            xlim([-1.5*pi 1.5.*pi])
    end

    
    
end
NiceSave('ISIByDecodedPos',figfolder,[])
