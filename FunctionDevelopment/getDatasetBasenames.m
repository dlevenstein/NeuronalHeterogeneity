function [baseNames] = getDatasetBasenames(varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
p = inputParser;
addParameter(p,'noPIfCTX',false)

parse(p,varargin{:})
noPIfCTX = p.Results.noPIfCTX;

%%
baseNames.CA1 = {'Achilles_10252013';'Achilles_11012013';'Buddy_06272013';'Cicero_09012014';'Cicero_09102014';'Cicero_09172014';'Gatsby_08022013';'Gatsby_08282013'};

switch noPIfCTX
    case false
        baseNames.fCTX = {'20140526_277um';'20140527_421um';'20140528_565um';...
            'BWRat17_121712';'BWRat17_121912';'BWRat18_020513';'BWRat19_032413';...
            'BWRat19_032513';'BWRat20_101013';'BWRat20_101513';'BWRat21_121113';...
            'BWRat21_121613';'BWRat21_121813';'Bogey_012615';'Dino_061814';...
            'Dino_061914';'Dino_062014';'Dino_072114';'Dino_072314';'Dino_072414';...
            'Rizzo_022615';'Rizzo_022715';'Splinter_020515';'Splinter_020915';'Templeton_032415'};
    case true
        baseNames.fCTX = {'20140526_277um';'20140527_421um';'20140528_565um';...
            'BWRat17_121712';'BWRat17_121912';'BWRat18_020513';'BWRat19_032413';...
            'BWRat19_032513';'BWRat20_101013';'BWRat20_101513';'BWRat21_121113';...
            'BWRat21_121613';'BWRat21_121813';'Bogey_012615';'Dino_061814';...
            'Dino_061914';'Dino_062014';'Dino_072114';'Dino_072314';'Dino_072414';...
            'Rizzo_022615';'Rizzo_022715';'Splinter_020515';'Splinter_020915';'Templeton_032415'};
end
baseNames.THAL = {'Mouse12-120806';'Mouse12-120807';'Mouse12-120808';'Mouse12-120809';'Mouse12-120810';'Mouse17-130125';'Mouse17-130128';'Mouse17-130129';'Mouse17-130130';'Mouse17-130131';'Mouse17-130201';'Mouse17-130202';'Mouse17-130203';'Mouse17-130204';'Mouse20-130514';'Mouse20-130515';'Mouse20-130516';'Mouse20-130517';'Mouse20-130520';'Mouse24-131213';'Mouse24-131216';'Mouse24-131217';'Mouse24-131218';'Mouse25-140123';'Mouse25-140124';'Mouse25-140128';'Mouse25-140129';'Mouse25-140130';'Mouse25-140131';'Mouse25-140203';'Mouse25-140204';'Mouse25-140205';'Mouse25-140206';'Mouse28-140310';'Mouse28-140311';'Mouse28-140312';'Mouse28-140313';'Mouse28-140314';'Mouse28-140317';'Mouse28-140318';'Mouse32-140822'};
baseNames.vCTX = {'YMV01_170818';'YMV02_170815';'YMV03_170818';'YMV04_170907';'YMV06_170913';'YMV07_170914';'YMV08_170922';'YMV09_171204';'YMV10_171213';'YMV11_171208';'YMV12_171211';'YMV13_180127';'YMV14_180128';'YMV15_180205';'YMV16_180206';'YMV17_180207';'YMV18_180208';'YMV19_180209'};
baseNames.BLA = {'Rat08-20130708';'Rat08-20130709';'Rat08-20130710';'Rat08-20130711';'Rat08-20130712';'Rat08-20130713';'Rat08-20130715';'Rat08-20130716';'Rat08-20130717';'Rat08-20130719';'Rat09-20140324';'Rat09-20140325';'Rat09-20140326';'Rat09-20140327';'Rat09-20140328';'Rat09-20140329';'Rat10-20140619';'Rat10-20140620';'Rat10-20140622';'Rat10-20140624';'Rat10-20140626';'Rat10-20140627';'Rat10-20140628';'Rat10-20140629';'Rat10-20140701';'Rat10-20140702';'Rat10-20140703';'Rat10-20140704';'Rat10-20140707';'Rat10-20140708';'Rat11-20150321';'Rat11-20150323';'Rat11-20150325';'Rat11-20150326';'Rat11-20150327';'Rat11-20150328';'Rat11-20150330'};
end

