#!/bin/bash
sbatch ./FullSharedFitjob.bash  'vCTX'  \
{'/gpfs/data/buzsakilab/DL/Database/YSData/YMV01/YMV01_170818',\
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV02/YMV02_170815',\
#/gpfs/data/buzsakilab/DL/Database/YSData/YMV03/YMV03_170818
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV04/YMV04_170907',\
#/gpfs/data/buzsakilab/DL/Database/YSData/YMV05/YMV05_170912
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV06/YMV06_170913',\
#/gpfs/data/buzsakilab/DL/Database/YSData/YMV07/YMV07_170914
#/gpfs/data/buzsakilab/DL/Database/YSData/YMV08/YMV08_170922
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV09/YMV09_171204',\
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV10/YMV10_171213',\
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV11/YMV11_171208',\
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV12/YMV12_171211',\
#/gpfs/data/buzsakilab/DL/Database/YSData/YMV13/YMV13_180127
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV14/YMV14_180128',\
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV15/YMV15_180205',\
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV16/YMV16_180206',\
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV17/YMV17_180207',\
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV18/YMV18_180208',\
'/gpfs/data/buzsakilab/DL/Database/YSData/YMV19/YMV19_180209'}  \
[]

sbatch ./FullSharedFitjob.bash  'fCTX'  \
{'/gpfs/data/buzsakilab/DL/Database/BWData/Bogey/Bogey_012615',\
'/gpfs/data/buzsakilab/DL/Database/BWData/BWRat17/BWRat17_121712',\
'/gpfs/data/buzsakilab/DL/Database/BWData/BWRat17/BWRat17_121912',\
'/gpfs/data/buzsakilab/DL/Database/BWData/BWRat18/BWRat18_020513',\
'/gpfs/data/buzsakilab/DL/Database/BWData/BWRat19/BWRat19_032413',\
'/gpfs/data/buzsakilab/DL/Database/BWData/BWRat19/BWRat19_032513',\
'/gpfs/data/buzsakilab/DL/Database/BWData/BWRat20/BWRat20_101013',\
'/gpfs/data/buzsakilab/DL/Database/BWData/BWRat20/BWRat20_101513',\
'/gpfs/data/buzsakilab/DL/Database/BWData/BWRat21/BWRat21_121113',\
'/gpfs/data/buzsakilab/DL/Database/BWData/BWRat21/BWRat21_121613',\
'/gpfs/data/buzsakilab/DL/Database/BWData/BWRat21/BWRat21_121813',\
'/gpfs/data/buzsakilab/DL/Database/BWData/Dino/Dino_061814',\
'/gpfs/data/buzsakilab/DL/Database/BWData/Dino/Dino_061914',\
'/gpfs/data/buzsakilab/DL/Database/BWData/Dino/Dino_062014',\
'/gpfs/data/buzsakilab/DL/Database/BWData/Dino/Dino_072114',\
'/gpfs/data/buzsakilab/DL/Database/BWData/Dino/Dino_072314',\
'/gpfs/data/buzsakilab/DL/Database/BWData/Dino/Dino_072414',\
'/gpfs/data/buzsakilab/DL/Database/BWData/JennBuzsaki22/20140526_277um',\
'/gpfs/data/buzsakilab/DL/Database/BWData/JennBuzsaki22/20140527_421um',\
'/gpfs/data/buzsakilab/DL/Database/BWData/JennBuzsaki22/20140528_565um',\
'/gpfs/data/buzsakilab/DL/Database/BWData/Rizzo/Rizzo_022615',\
'/gpfs/data/buzsakilab/DL/Database/BWData/Rizzo/Rizzo_022715',\
'/gpfs/data/buzsakilab/DL/Database/BWData/Splinter/Splinter_020515',\
'/gpfs/data/buzsakilab/DL/Database/BWData/Splinter/Splinter_020915',\
'/gpfs/data/buzsakilab/DL/Database/BWData/Templeton/Templeton_032415'}  \
[]

sbatch ./FullSharedFitjob.bash  'fCTX'  \
{'/gpfs/data/buzsakilab/DL/Database/AGData/Achilles/Achilles_11012013',\
'/gpfs/data/buzsakilab/DL/Database/AGData/Achilles/Achilles_10252013',\
'/gpfs/data/buzsakilab/DL/Database/AGData/Buddy/Buddy_06272013',\
'/gpfs/data/buzsakilab/DL/Database/AGData/Cicero/Cicero_09102014',\
'/gpfs/data/buzsakilab/DL/Database/AGData/Cicero/Cicero_09012014',\
'/gpfs/data/buzsakilab/DL/Database/AGData/Cicero/Cicero_09172014',\
'/gpfs/data/buzsakilab/DL/Database/AGData/Gatsby/Gatsby_08282013',\
'/gpfs/data/buzsakilab/DL/Database/AGData/Gatsby/Gatsby_08022013'}  \
[]



sbatch ./FullSharedFitjob.bash  'THAL'  \
{'/gpfs/data/buzsakilab/DL/Database/APData/Mouse12/Mouse12-120806',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse12/Mouse12-120807',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse12/Mouse12-120808',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse12/Mouse12-120809',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse12/Mouse12-120810',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse17/Mouse17-130125',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse17/Mouse17-130128',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse17/Mouse17-130129',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse17/Mouse17-130130',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse17/Mouse17-130131',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse17/Mouse17-130201',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse17/Mouse17-130202',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse17/Mouse17-130203',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse17/Mouse17-130204',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse20/Mouse20-130514',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse20/Mouse20-130515',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse20/Mouse20-130516',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse20/Mouse20-130517',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse20/Mouse20-130520',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse24/Mouse24-131213',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse24/Mouse24-131216',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse24/Mouse24-131217',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse24/Mouse24-131218',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse25/Mouse25-140123',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse25/Mouse25-140124',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse25/Mouse25-140128',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse25/Mouse25-140129',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse25/Mouse25-140130',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse25/Mouse25-140131',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse25/Mouse25-140203',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse25/Mouse25-140204',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse25/Mouse25-140205',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse25/Mouse25-140206',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse28/Mouse28-140310',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse28/Mouse28-140311',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse28/Mouse28-140312',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse28/Mouse28-140313',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse28/Mouse28-140314',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse28/Mouse28-140317',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse28/Mouse28-140318',\
'/gpfs/data/buzsakilab/DL/Database/APData/Mouse32/Mouse32-140822'}  \
[]


sbatch ./FullSharedFitjob.bash  'PIR'  \
{'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130708/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130709/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130710/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130711/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130712/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130713/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130715/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130716/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130717/',\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130718/,\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130719/',\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130720/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130722/,\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140324/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140325/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140326/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140327/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140328/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140329/',\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140331/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140401/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140402/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140403/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140404/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140405/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140407/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140408/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140409/,\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140619/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140620/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140622/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140624/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140626/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140627/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140628/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140629/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140701/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140702/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140703/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140704/',\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140705/,\  #double region on a shank...
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140707/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140708/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150321/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150323/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150325/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150326/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150327/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150328/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150330/'}  \
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150331/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150401/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150402/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150403/,\
'pir'

sbatch ./FullSharedFitjob.bash  'BLA'  \
{'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130708/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130709/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130710/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130711/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130712/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130713/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130715/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130716/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130717/',\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130718/,\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130719/',\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130720/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat08/Rat08-20130722/,\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140324/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140325/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140326/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140327/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140328/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140329/',\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140331/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140401/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140402/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140403/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140404/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140405/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140407/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140408/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat09/Rat09-20140409/,\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140619/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140620/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140622/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140624/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140626/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140627/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140628/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140629/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140701/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140702/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140703/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140704/',\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140705/,\  #double region on a shank...
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140707/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat10/Rat10-20140708/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150321/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150323/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150325/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150326/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150327/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150328/',\
'/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150330/'}  \
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150331/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150401/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150402/,\
#/gpfs/data/buzsakilab/DL/Database/GGData/Rat11/Rat11-20150403/,\
'bla'
