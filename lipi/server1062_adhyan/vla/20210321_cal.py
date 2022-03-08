##Also set this as general script

#### SCRIPT for calibrating the 20130423 data at L Band ####
import os

mspath = '/nas08-data02/rohit/20210322/cal_all/attempt4/'
msvis = '/nas08-data02/rohit/20210322/cal_all/attempt4/21A-118.sb39452573.eb39520716.59295.66536_cal_all.1s.ms'
calhs='20210322_L_100ms_cal.ms'
calprefix='/nas08-data02/rohit/20210322/cal_all/attempt4/calSUN_20210322-L'

os.system('rm -rf '+calprefix+'*')
os.system('rm -rf '+calhs)
os.system('rm -rf '+calhs+'.flagversions')

# change directory to mspath
os.chdir(mspath)

# calibration controls
doinfo = 1
doflag = 0
precal=1
docalib = 1
plotipcal=0
plotbp=0
docal20db = 0
docalsun = 0
doimage = 0
doplot = 0

refant='ea07'

## information step##
if doinfo:
###### MS information  ######
#### listobs ####
    listobs(vis=msvis,listfile=msvis+'.listobs',overwrite=True)
#### plotants ####
    figfile=msvis+'.ants.png'
    plotants(vis=msvis,figfile=figfile)

## flag step ##
if doflag:
    ##### Data Flagging ######### flagcmd: apply telescope flags #####flagcmd(vis=msvis, inpmode='list', inpfile="20130423_L_50ms.ms.onlineflags")#flagdata(vis=msvis, scan='1~6, 8, 28, 30, 50, 52')
    flagdata(vis=msvis, mode='quack', quackmode='beg', quackinterval=6.0)
# flagcmd(vis=msfile,inpmode='list',inpfile='manualflgs_20120310A.flgcmd')

if precal:
##### Calibrations ######
#### hanningsmooth ####
    hanningsmooth(vis=msvis,datacolumn='data',field='',outputvis=calhs)
#### check for antenna position ####
    gencal(vis=calhs,caltable=calprefix+'.antpos',caltype='antpos')
#### gain curve calibration ####
    gencal(vis=calhs,caltable=calprefix+'.gaincurve',caltype='gc')
    listcal(vis=calhs,caltable=calprefix+'.gaincurve')
# All solutions are unity, NOT going to use it
    gencal(vis=calhs,caltable=calprefix+'.rq',caltype='rq')

# setting flux density scale #
setjy(vis=calhs, field='3C48', model='3C48_L.im')

## try the CASA flag algorithm 
#flagmanager(vis=calhs, mode='save', versionname='before-auto-flag')
#flagmanager(vis=calhs, mode='restore', versionname='before-auto-flag')
#flagdata(vis=calhs, mode='tfcrop', spw='0~7',field='',scan='',
#         datacolumn='data', action='apply',
#         display='both', flagbackup=False)


#flaginfo
#flagInfo = flagdata(vis=calhs, mode='summary')

## Refant Antenna ##
#antennas='ea01,ea02,ea03,ea05,ea07,ea08,ea09,ea10,ea11,ea12,ea13,ea14,ea16,ea17,ea18,ea19,ea20,ea21,ea22,ea25,ea26,ea27,ea28'

## SPW 0
#antenna='ea01,ea24,ea28';spw0='0:32~37,0:44~46,0:66~70,0:90~110'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw0,antenna=antenna)
#antenna='ea02,ea10,ea19,ea20';spw0='0:30~35,0:42~48,0:65~70,0:75~80,0:90~115'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw0,antenna=antenna)
#antenna='ea11';spw0='0:30~35,0:42~48,0:53~57,0:75~80,0:90~115'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw0,antenna=antenna)
#antenna='ea03,ea18,ea21';spw0='0:30~35,0:90~115'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw0,antenna=antenna)
#antenna='ea04,ea05,ea07,ea09,ea14,ea16,ea26,ea27';spw0='0:30~35,0:65~70,0:90~115'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw0,antenna=antenna)
#antenna='ea08,ea12,ea17,ea22';spw0='0:65~70,0:90~115'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw0,antenna=antenna)
#antenna='ea13,ea25';spw0='0:90~115'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw0,antenna=antenna)
## SPW 1
#antenna='ea01,ea02,ea05,ea08,ea09,ea18,ea22';spw1='1:45~60,1:80~90,1:100~110'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw1,antenna=antenna)
#antenna='ea03,ea04,ea10,ea11,ea12,ea13,ea14,ea16,ea17,ea20,ea25,ea26,ea27';spw1='1:45~60,1:100~110'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw1,antenna=antenna)
#antenna='ea07,ea21';spw1='1:45~60,1:100~110,1:120~127'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw1,antenna=antenna)
#antenna='ea19,ea24,ea28';spw1='1:0~20,1:45~60,1:80~90,1:100~110'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw1,antenna=antenna)
## SPW 2
#spw2='2:10~40'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw2,antenna=antenna)
## SPW 3
## SPW 4
#antenna='ea01';spw4='4:25~38,4:60~64,4:75~85,4:92~100'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw4,antenna=antenna)
#antenna='ea02,ea03,ea04,ea05,ea07,ea08,ea09,ea10,ea11,ea12,ea13,ea17,ea18,ea19,ea20,ea21,ea22,ea24,ea25,ea26,ea27,ea28';spw4='4:25~65,4:75~85,4:92~100'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw4,antenna=antenna)
#antenna='ea14,ea16';spw4='4:30~38,4:47~55,4:58~65,4:75~85,4:92~100'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw4,antenna=antenna)
## SPW 5
#antenna='';spw5='5:0~5,5:50~75'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spw5,antenna=antenna)
####################################################################################

flagant0='ea01&ea02; ea01&ea04; ea01&ea05;ea01&ea09;ea01&ea13;ea02&ea09;ea02&ea13; ea02&ea26;ea03&ea08; ea03&ea28; ea04&ea05; ea04&ea09;ea04&ea13;ea04&ea18;ea04&ea19;ea04&ea20;ea05&ea09;ea05&ea13;ea05&ea18;ea05&ea19; ea05&ea20;ea05&ea22;ea05&ea27;ea07&ea10;ea08&ea26;ea09&ea12; ea09&ea13; ea09&ea18; ea09&ea19; ea09&ea20; ea09&ea21; ea09&ea22; ea09&ea26;  ea09&ea27; ea10&ea14; ea11&ea20; ea11&ea25; ea12&ea18;ea12&ea19; ea12&ea20; ea13&ea18; ea13&ea19; ea13&ea22; ea13&ea26;ea14&ea16; ea14&ea21; ea16&ea21; ea17&ea25; ea18&ea19;   ea18&ea20; ea18&ea22;ea19&ea22; ea19&ea27;  ea20&ea25; ea21&ea27; ea22&ea27; ea24&ea28'
flagant1='!ea01&ea07;!ea01&ea08;!ea01&ea10;!ea01&ea11;!ea01&ea14;!ea01&ea17;!ea01&ea24;!ea01&ea25;!ea02&ea07;!ea02&ea10;!ea02&ea11;!ea02&ea14;!ea02&ea16;!ea02&ea17;!ea02&ea20;!ea02&ea21;!ea02&ea24;!ea02&ea25;!ea02&ea27;!ea02&ea28;!ea03&ea04;!ea03&ea05;!ea03&ea07;!ea03&ea10;!ea03&ea11;!ea03&ea12;!ea03&ea13;!ea03&ea14;!ea03&ea16;!ea03&ea17;!ea03&ea18;!ea03&ea19;!ea03&ea20;!ea03&ea21;!ea03&ea22;!ea03&ea24;!ea03&ea25;!ea03&ea27;!ea04&ea07;!ea04&ea08;!ea04&ea10;!ea04&ea11;!ea04&ea14;!ea04&ea16;!ea04&ea17;!ea04&ea24;!ea04&ea25;!ea04&ea26;!ea04&ea28;!ea05&ea07;!ea05&ea08;!ea05&ea10;!ea05&ea11;!ea05&ea14;!ea05&ea16;!ea05&ea17;!ea05&ea24;!ea05&ea25;!ea05&ea28;!ea07&ea08;!ea07&ea09;!ea07&ea11;!ea07&ea12;!ea07&ea13;!ea07&ea14;!ea07&ea16;!ea07&ea17;!ea07&ea18;!ea07&ea19;!ea07&ea20;!ea07&ea22;!ea07& ea24;!ea07&ea25;!ea07&ea26;!ea07&ea27;!ea07&ea28;!ea07&ea10;!ea07&ea11;!ea07&ea12;!ea08&ea14;!ea08& ea16;!ea08&ea17;!ea08&ea18;!ea08&ea19;!ea08&ea20;!ea08&ea21;!ea08&ea24;!ea08&ea27;!ea09&ea10;!ea09&ea11;!ea09&ea17;!ea09&ea25;!ea10&ea11;!ea10&ea12;!ea10&ea13;!ea10&ea16;!ea10&ea17;!ea10&ea18;!ea10&ea19;!ea10&ea20;!ea10&ea21;!ea10&ea22;!ea10&ea24;!ea10&ea25;!ea10&ea26;!ea10&ea27;!ea10&ea28;!ea11&ea14;!ea11&ea16;!ea11&ea17;!ea11&ea18;!ea11&ea19;!ea11&ea21;!ea11&ea22;!ea11&ea24;!ea11&ea26;!ea11&ea27;!ea11&ea28;!ea12&ea14;!ea12&ea17;!ea12&ea22;!ea12&ea24;!ea12&ea25;!ea12&ea26;!ea12&ea28;!ea13&ea14;!ea13&ea16;!ea13&ea17;!ea13&ea25;!ea13&ea28;!ea14&ea17;!ea14&ea18;!ea14&ea19;!ea14&ea22;!ea14&ea24;!ea14&ea25;!ea14&ea26;!ea14&ea28;!ea16&ea17;!ea16&ea18;!ea16&ea22;!ea16&ea24;!ea16&ea25;!ea16&ea26;!ea16&ea28;!ea17&ea18;!ea17&ea19;!ea17&ea20;!ea17&ea21;!ea17&ea22;!ea17&ea24;!ea17&ea26;!ea17&ea27;!ea17&ea28;!ea18&ea24;!ea18&ea25;!ea18&ea26;!ea18&ea28;!ea19&ea24;!ea19&ea25;!ea19&ea26;!ea19&ea28;!ea20&ea24;!ea20&ea26;!ea21&ea24;!ea21&ea25;!ea21&ea26;!ea22&ea24;!ea22&ea25;!ea22&ea28;!ea24&ea25;!ea24&ea26;!ea24&ea27;!ea25&ea26;!ea25&ea27;!ea25&ea28;!ea26&ea28;!ea27&ea28'
flagant2='ea01&ea02; ea01&ea04; ea01&ea05; ea01&ea09; ea01&ea13; ea01&ea18; ea01&ea19; ea01&ea20; ea01&ea21; ea01&ea22; ea01&ea26; ea01&ea27;  ea02&ea04; ea02&ea05; ea02&ea09; ea02&ea13; ea02&ea26; ea02&ea28; ea03&ea08; ea03&ea09; ea03&ea13; ea03&ea26; ea03&ea28; ea04&ea05; ea04&ea09;  ea04&ea12; ea04&ea13; ea04&ea18; ea04&ea19; ea04&ea20; ea04&ea21; ea04&ea22; ea05&ea09; ea05&ea10; ea05&ea12; ea05&ea13; ea05&ea18; ea05&ea19; ea05&ea20; ea05&ea21; ea05&ea22; ea05&ea27; ea07&ea10; ea08&ea09; ea08&ea13; ea08&ea26; ea08&ea27; ea08&ea28; ea09&ea11; ea09&ea12; ea09&ea13; ea09&ea16; ea09&ea18; ea09&ea19; ea09&ea20; ea09&ea21; ea09&ea22; ea09&ea24; ea09&ea25; ea09&ea26; ea09&ea27; ea09&ea28; ea10&ea14; ea11&ea12;  ea11&ea20; ea12&ea16; ea12&ea18; ea12&ea19; ea12&ea20; ea13&ea19; ea13&ea20; ea13&ea21; ea13&ea22; ea13&ea26; ea13&ea27; ea14&ea16; ea16&ea18; ea16&ea21; ea16&ea27; ea17&ea25; ea18&ea19; ea18&ea20; ea18&ea21; ea18&ea27; ea19&ea22; ea19&ea27; ea19&ea28; ea20&ea25; ea21&ea27; ea22&ea27; ea24&ea26; ea24&ea28'
flagant3='ea02&ea13; ea04&ea09; ea05&ea09; ea09&ea22; ea12&ea18; ea12&ea20; ea14&ea16; ea16&ea27; ea19&ea22; ea19&ea27; ea21&ea27'
flagant4='ea01&ea07;!ea01&ea10;!ea01&ea14;!ea01&ea16;!ea01&ea17;!ea01&ea21;!ea01&ea24;!ea01&ea25;!ea01&ea28;!ea02&ea07;!ea02&ea10;!ea02&ea14;!ea02&ea16;!ea02&ea21;!ea02&ea24;!ea03&ea04;!ea03&ea05;!ea03&ea07;!ea03&ea10;!ea03&ea12;!ea03&ea14;!ea03&ea16;!ea03&ea27;!ea04&ea07;!ea04&ea10;!ea04&ea14;!ea04&ea16;!ea04&ea21;!ea04&ea24;!ea05&ea07;!ea05&ea10;!ea05&ea14;!ea05&ea28;!ea07&ea08;!ea07&ea09;!ea07&ea11;!ea07&ea12;!ea07&ea13;!ea07&ea17;!ea07&ea18;!ea07&ea19;!ea07&ea20;!ea07&ea22;!ea07&ea24;!ea07&ea25;!ea07&ea26;!ea07&ea27;!ea07&ea28;!ea08&ea10;!ea08&ea14;!ea08&ea16;!ea09&ea10;!ea09&ea24;!ea09&ea25;!ea10&ea11;!ea10&ea12;!ea10&ea13;!ea10&ea17;!ea10&ea18;!ea10&ea20;!ea10&ea22;!ea10&ea24;!ea10&ea25;!ea10&ea26;!ea10&ea27;!ea10&ea28;!ea11&ea14;!ea11&ea16;!ea11&ea25;!ea12&ea16;!ea12&ea24;!ea13&ea14;!ea13&ea16;!ea13&ea17;!ea13&ea21;!ea13&ea24;!ea13&ea28;!ea14&ea17;!ea14&ea24;!ea14&ea25;!ea14&ea28;!ea16&ea17;!ea16&ea24;!ea16&ea25;!ea16&ea26;!ea16&ea28;!ea17&ea18;!ea17&ea22;!ea17&ea27;!ea18&ea28;!ea19&ea24;!ea19&ea24;!ea20&ea25;!ea20&ea28;!ea21&ea24;!ea21&ea25;!ea21&ea26;!ea21&ea28;!ea22&ea28;!ea25&ea27;!ea26&ea27;!ea26&ea28;!ea27&ea28'
flagant5=' ea01&ea02; ea01&ea04; ea01&ea05; ea01&ea09; ea01&ea13; ea01&ea18; ea01&ea19; ea02&ea13; ea02&ea26; ea04&ea05; ea04&ea09; ea04&ea13; ea04&ea19; ea05&ea09;ea05&ea19; ea07&ea10; ea08&ea26; ea09&ea13; ea09&ea19; ea09&ea21; ea09&ea22; ea10&ea14; ea12&ea18; ea12&ea19; ea12&ea20; ea14&ea16; ea16&ea21; ea18&ea19; ea19&ea20; ea19&ea21; ea19&ea22; ea19&ea27;  ea21&ea27'
flagant6='ea01&ea02; ea01&ea04; ea01&ea05; ea01&ea09; ea01&ea12;  ea01&ea13;  ea01&ea26; ea01&ea27; ea02&ea13; ea02&ea26; ea04&ea05; ea04&ea09; ea04&ea18; ea04&ea19; ea05&ea09; ea02&ea18;ea02&ea19; ea07&ea10; ea08&ea26; ea09&ea13; ea09&ea22; ea10&ea14; ea12&ea18; ea12&ea20; ea12&ea21; ea12&ea27; ea14&ea16; ea16&ea21; ea18&ea20; ea19&ea21; ea19&ea22;  ea19&ea27;  ea21&ea27; ea24&ea28'
flagant7='ea01&ea02; ea01&ea04; ea01&ea05; ea01&ea13; ea01&ea18; ea01&ea20; ea02&ea13; ea02&ea14; ea04&ea05; ea04&ea09; ea04&ea11; ea04&ea12; ea04&ea17; ea04&ea18; ea04&ea20; ea04&ea27; ea02&ea18; ea05&ea19; ea05&ea20; ea07&ea10; ea08&ea26; ea09&ea12; ea09&ea20; ea09&ea22; ea09&ea27; ea10&ea14; ea10&ea16; ea11&ea25; ea11&ea27; ea12&ea17; ea12&ea18; ea12&ea20; ea12&ea22; ea14&ea16;ea16&ea21;  ea16&ea27; ea17&ea20; ea18&ea20; ea18&ea25; ea19&ea20;  ea19&ea25; ea19&ea27; ea20&ea22; ea20&ea25; ea20&ea27; ea21&ea27'

#flagdata(vis=calhs,field='0,1',mode='manual',spw='0',antenna=flagant0)
#flagdata(vis=calhs,field='0,1',mode='manual',spw='1',antenna=flagant1)
#flagdata(vis=calhs,field='0,1',mode='manual',spw='2',antenna=flagant2)
#flagdata(vis=calhs,field='0,1',mode='manual',spw='3',antenna=flagant3)
#flagdata(vis=calhs,field='0,1',mode='manual',spw='4',antenna=flagant4)
#flagdata(vis=calhs,field='0,1',mode='manual',spw='5',antenna=flagant5)
#flagdata(vis=calhs,field='0,1',mode='manual',spw='6',antenna=flagant6)
#flagdata(vis=calhs,field='0,1',mode='manual',spw='7',antenna=flagant7)


#antenna='';spwa='0:90~118,0:122~127,1:50~60,2:10~23,2:35~45,4:20~85,5:55~80'
#flagdata(vis=calhs,field='0,1',mode='manual',spw=spwa)

#scang='4,5,39,40,74,75,109,110,144,145,179,180,215,216,217,218'
sys.exit()
##calibration##
if docalib:
#### initial phase-only solution on the flux calibrator t-p#### 
#the ant ea23 RR in spw4~7 have unusual phase configuration, ant ea25 has very low amp in LL
    #flagdata(vis=calhs,field='3C48',spw='4~7',antenna='ea23',correlation='RR')
    #flagdata(vis=calhs,field='3C48',spw='0~7',antenna='ea25',correlation='LL')
    spw1g='0:20~24,1:20~24,2:55~60,3:28~32,4:10~14,5:28~32,6:28~32,7:28~32'
    gaincal(vis=calhs, caltable=calprefix+'.pha0', field='2',\
            spw=spw1g, gaintype='G',gaintable=[calprefix+'.rq'],refant=refant, calmode='p', solint='10s', minsnr=3)
    if plotipcal:
        plotcal(caltable=calprefix+'.pha0',xaxis='time',yaxis='phase',subplot=441,\
                spw='',iteration='antenna',plotrange=[-1,-1,-180,180])

#### delay calibration f-p####
#    spwd='0:2~16;20~28;32~34;40~45;54~60;65~70;87~94;110~125,1:2~15;20~39;65~75;85~95;110~125,2:2~16;55~80;85~125,\
#            3:2~6;12~125,4:2~28;44~46;68~77;85~95;102~106;115~124,5:15~55;80~125,6:2~14;16~35;42~78;83~115;121~125,7:2~45;75~82;90~125'
    spwd='0:0~16;20~28;32~46;54~60,1:0~5;10~49;55~59,2:0~6;17~62,\
            3:0~62,4:0~14;35~38;44~51;58~60,5:5~12;15~62,6:0~62,7:0~19;24~26;38~41;45~48;55~62'
#    gaincal(vis=calhs, caltable=calprefix+'.delay', field='0,1', spw=spwd, \
#            gaintype='K', gaintable=[calprefix+'.rq', calprefix+'.pha0'],\
#            gainfield=['','2'],refant=refant, \
#            combine='scan',solint='inf', minsnr=5)
    gaincal(vis=calhs, caltable=calprefix+'.delay', field='0,1', spw=spwd, \
            gaintype='K', gaintable=[calprefix+'.rq'],\
            gainfield=[''],refant=refant, \
            combine='scan',solint='inf', minsnr=5)
    plotcal(caltable=calprefix+'.delay',xaxis='freq',yaxis='delay',iteration='antenna',subplot=441,plotrange= [-1, -1, -10, 10])
    listcal(vis=calhs,caltable=calprefix+'.delay')
#### bandpass calibration ####
    #spwbp='0:0~15;21~46;51~63,1:0~39;45~62,2~3,4:0~13;22~26;35~38;43~63,5:4~63,6~7'
    spwbp='0:0~16;20~28;32~46;54~60,1:0~5;10~49;55~59,2:0~6;17~62,\
            3:0~62,4:0~14;35~38;44~51;58~60,5:5~12;15~62,6:0~62,7:0~19;24~26;38~41;45~48;55~62'
#    bandpass(vis=calhs,caltable=calprefix+'.bp',field='0,1',spw=spwbp,\
#             bandtype='B', refant=refant, fillgaps=10,minsnr=2.0,\
#             gaintable=[calprefix+'.rq',calprefix+'.pha0',calprefix+'.delay'],\
#             gainfield=['','2','0,1'],\
#             solnorm=False,combine='scan',solint='inf')
    bandpass(vis=calhs,caltable=calprefix+'.bp',field='0,1',spw=spwbp,\
             bandtype='B', refant=refant, fillgaps=10,minsnr=2.0,\
             gaintable=[calprefix+'.rq',calprefix+'.delay'],\
             gainfield=['','0,1'],\
             solnorm=False,combine='scan',solint='inf')
    if plotbp:
        plotbandpass(caltable=calprefix+'.bp',overlay='spw',spw='0~7',xaxis='freq',yaxis='phase',\
                     subplot=32,markersize=6,interactive=True,plotrange=[0,0,-80,80])
        plotbandpass(caltable=calprefix+'.bp',overlay='spw',spw='0~7',xaxis='freq',yaxis='amp',\
                     subplot=42,markersize=6,interactive=True,plotrange=[0,0,0,0])
        plotcal(caltable=calprefix+'.bp',spw='',
            xaxis='freq',yaxis='phase',field= '0,1',subplot=441, 
            iteration='antenna',plotrange=[-1,-1,-180,180])
        plotcal(caltable=calprefix+'.bp',
            xaxis='freq',yaxis='amp',field= '0,1',subplot=441, 
            iteration='antenna',plotrange=[-1,-1,0,15])
#### gain calibration ####
## phase only solution PER SCAN (to be used to calibrate solar scans) ##
    #restrict uvrange to avoid strong RFIs
    spw0='0:0~16;20~28;32~46;54~60,1:0~5;10~49;55~59,2:0~6;17~62,\
            3:0~62,4:0~14;35~38;44~51;58~60,5:5~12;15~62,6:0~62,7:0~19;24~26;38~41;45~48;55~62'
    gaincal(vis=calhs,caltable=calprefix+'.pha_inf',field='0,1',spw=spw0,gaintype='G',uvrange='>200m',\
            calmode='p',refant=refant,solint='inf',solnorm=False,\
            gaintable=[calprefix+'.rq',calprefix+'.delay',calprefix+'.bp'],\
            gainfield=['','0,1','0,1'])
#    plotcal(caltable=calprefix+'.pha_inf',xaxis='time',yaxis='phase',subplot=441,
#            spw='',iteration='antenna',plotrange=[-1,-1,-180,180])
## amp only solution PER SCAN ##
    gaincal(vis=calhs,caltable=calprefix+'.amp_inf',field='0,1',spw=spw0,gaintype='G',uvrange='>200m',\
            calmode='a',refant=refant,solint='inf',solnorm=False,\
            gaintable=[calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',calprefix+'.pha_inf'],\
            gainfield=['','0,1','0,1','0,1'])
#    plotcal(caltable=calprefix+'.amp_inf',xaxis='time',yaxis='amp',subplot=441,\
#            spw='0',iteration='antenna',plotrange=[-1,-1,-1,-1])

## amp only solution COMBING SCAN ##
    gaincal(vis=calhs,caltable=calprefix+'.amp_comb',field='0,1',spw=spw0,gaintype='G',uvrange='>200m',\
            calmode='a',refant=refant,solint='inf',combine='scan',solnorm=False,\
            gaintable=[calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',calprefix+'.pha_inf'],\
            gainfield=['','0,1','0,1','0,1'])

if apply:
#################################################
##### Apply calibrations to the calibrators #####
#################################################
#### 3C48 ####
    applycal(vis=calhs,field='0,1',
             gaintable=[calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',\
             calprefix+'.pha_inf',calprefix+'.amp_inf'],\
             gainfield=['','0,1','0,1','0,1','0,1'],\
             interp=['','','','',''],parang=False,flagbackup=True)

sys.exit()

## 20db
antennas='ea01,ea03,ea06,ea07,ea08,ea09,ea17,ea18,ea19,ea21,ea23,ea25,ea26,ea27'

par_mbd = [ -0.008,-0.005,-0.012,-0.006,0.012,0.164,0.009,0.001,0.121,0.120,0.000,0.000,0.019,0.020,0.004,0.012,0.099,0.104,0.000,0.000,-0.006,-0.006,-0.025,\
            -0.027,-0.023,-0.017,0.112,0.112]

par_pha = [ -1.0,0.0,-0.6,-0.5,-0.3,0.9,-1.5,-0.2,-1.6,-1.7,0.0,0.0,-3.5,-3.0,-3.3,-2.3,0.1,-0.1,0.0,0.0,5.4,5.5,-0.7,-0.5,-1.4,-1.6,-3.0,-3.0]

gencal(vis=calhs,caltable=calprefix+'.20db_mbd',caltype='mbd',spw='',\
       antenna=antennas,pol='R,L',parameter=par_mbd) 
gencal(vis=calhs,caltable=calprefix+'.20db_pha',caltype='ph',spw='',\
       antenna=antennas,pol='R,L',parameter=par_pha) 



## make image of sun
split(vis=msvis, outputvis='2050_L.1s.ms', timerange='20:35:00~20:50:00',datacolumn='data')
clearcal('2050_L.1s.ms')
applycal(vis='2050_L.1s.ms',field='SUN',\
         gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',\
         calprefix+'.pha_inf',calprefix+'.amp_inf',\
         calprefix+'.20db_mbd',calprefix+'.20db_pha'],\
         gainfield=['','','0','0','0','0','',''],\
         interp=['','','','','linear','linear','',''],\
         spwmap=[],applymode='',calwt=False, parang=False,flagbackup=True)

split(vis='2050_L.1s.ms', outputvis='2050_L.1s.cal.ms', datacolumn='corrected')

from suncasa.utils import dspec as ds
from suncasa.utils import svplot as sv
ds.plt_dspec(specdata="/srg/yjluo/xixi/20130423/sun_20130423_L_1s.med.spec.npz",dmin=0.005,dmax=0.3)
sv.svplot(vis='2050_L.1s.cal.ms',timerange='20:45:17~20:45:23',spw='*:1.03~1.93Ghz',specfile='sun_20130423_L_1s.med.spec.npz',dmin=0.005,dmax=0.3)


flagdata(vis='2050_L.1s.cal.ms',field='0',antenna='ea09&ea26')








#### SCRIPT for calibrating the 20130423 data at S Band ####
import os

mspath = '/srg/yjluo/xixi/20130423/'
msvis = '/srg/yjluo/xixi/20130423/20130423_S_1s.ms/'
vispre='20130423_S_1s.ms'
calhs='20130423_S_1s.ms.hs.ms'
calprefix='cal-S/calSUN_20130423-S'

# change directory to mspath
os.chdir(mspath)

# calibration controls
doinfo = 0
doflag = 0
precal=0
docalib = 0
plotipcal=0
plotbp=0
docal20db = 0
docalsun = 0
doimage = 0
doplot = 1


## information step##
if doinfo:
###### MS information  ######
#### listobs ####
    listobs(vis=msvis,listfile=mspath+vispre+'.listobs',overwrite=True)
#### plotants ####
    figfile=mspath+vispre+'.ants.png'
    plotants(vis=msvis,figfile=figfile)

## flag step ##
if doflag:
    ##### Data Flagging #####
    #### flagcmd: apply telescope flags ####
    flagcmd(vis=msvis, inpmode='list', inpfile="/srg/yjluo/xixi/20130423/20130423_S_50ms.ms.onlineflags")
    #### flagdata ####
    flagdata(vis=msvis, scan='1~6, 8, 18, 20, 30, 32, 42, 44, 54, 56, 66, 68')
    flagdata(vis=msvis, mode='quack', quackmode='beg', quackinterval=6.0)

if precal:
##### Calibrations ######
#### hanningsmooth ####
    hanningsmooth(vis=msvis,datacolumn='data',field='',outputvis=calhs)
#### check for antenna position ####
    gencal(vis=calhs,caltable=calprefix+'.antpos',caltype='antpos')
#### gain curve calibration ####
    gencal(vis=calhs,caltable=calprefix+'.gaincurve',caltype='gc')
    listcal(vis=calhs,caltable=calprefix+'.gaincurve')
# All solutions are unity, NOT going to use it
    gencal(vis=calhs,caltable=calprefix+'.rq',caltype='rq')

# setting flux density scale #
setjy(vis=calhs, field='3C48', model='3C48_S.im')

## try the CASA flag algorithm 
#flagmanager(vis=calhs, mode='save', versionname='before-auto-flag')
#flagmanager(vis=calhs, mode='restore', versionname='before-auto-flag')
flagdata(vis=calhs, mode='tfcrop', spw='',field='3C48',scan='7,19,31,43,55,67',\
         datacolumn='data', action='calculate',\
         display='both', flagbackup=False)


#flaginfo
flagInfo = flagdata(vis=calhs, mode='summary')

## Refant Antenna ##
refant='ea14'
antennas='ea02,ea04,ea05,ea10,ea11,ea12,ea13,ea14,ea15,ea16,ea22,ea24,ea28'

##calibration##
if docalib:
#### initial phase-only solution on the flux calibrator t-p#### 
#ant ea04 has very low amp in LL in spw0~8
    flagdata(vis=calhs,field='3C48',spw='0~8',antenna='ea04',correlation='LL')
    flagdata(vis=calhs,field='3C48',antenna='ea04')
    flagdata(vis=calhs,field='3C48',scan='43')
    flagdata(vis=calhs,field='3C48',scan='7')
    spw1='0:40~44,1:20~24,2:20~24,3:40~44,4:20~24,5:10~14,6:28~32,7:28~32,8:28~32,9:28~32,10:40~44,11:50~54,12~15:28~32'
    gaincal(vis=calhs, caltable=calprefix+'.pha0', field='3C48',\
            spw=spw1, gaintype='G',gaintable=[calprefix+'.antpos',calprefix+'.rq'],refant=refant, calmode='p', solint='10s', minsnr=3)
    if plotipcal:
        plotcal(caltable=calprefix+'.pha0',xaxis='time',yaxis='phase',subplot=441,\
                spw='',antenna=antennas,iteration='antenna',plotrange=[-1,-1,-180,180])

#### delay calibration f-p####
    spwd='0:4~60,1:2~30;45~60,2:2~33;55~60,3~4:2~60,5:2~38;42~60,6:2~37;41~48;52~60,7:2~60,\
            8:3~60,9~14:2~60,15:2~50'
    gaincal(vis=calhs, caltable=calprefix+'.delay', field='3C48', spw=spwd, \
            gaintype='K', gaintable=[calprefix+'.antpos',calprefix+'.rq', calprefix+'.pha0'],\
            gainfield=['','','3C48'],refant=refant, \
            combine='scan',solint='inf', minsnr=5)
    plotcal(caltable=calprefix+'.delay',xaxis='freq',yaxis='delay',iteration='antenna',subplot=441,plotrange= [-1, -1, -10, 10])
    listcal(vis=calhs,caltable=calprefix+'.delay')
#### bandpass calibration ####
    #spwbp='0:0~15;21~46;51~63,1:0~39;45~62,2~3,4:0~13;22~26;35~38;43~63,5:4~63,6~7'
    spwbp='0:4~62,1:0~30;45~62,2:0~33;55~62,3~4:0~62,5:0~38;42~62,6:0~37;41~48;52~62,7:0~62,\
            8:3~62,9~14:0~62,15:0~50'
    bandpass(vis=calhs,caltable=calprefix+'.bp',field='3C48',spw=spwbp,\
             bandtype='B', refant=refant, fillgaps=10,minsnr=2.0,\
             gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.pha0',calprefix+'.delay'],\
             gainfield=['','','3C48','3C48'],\
             solnorm=False,combine='scan',solint='inf')
    if plotbp:
        plotbandpass(caltable=calprefix+'.bp',overlay='spw',spw='0~7',xaxis='freq',yaxis='phase',\
                     subplot=32,markersize=6,interactive=True,plotrange=[0,0,-80,80])
        plotbandpass(caltable=calprefix+'.bp',overlay='spw',spw='0~7',xaxis='freq',yaxis='amp',\
                     subplot=42,markersize=6,interactive=True,plotrange=[0,0,0,0])
        plotcal(caltable=calprefix+'.bp',spw='',
            xaxis='freq',yaxis='phase',field= '3C48',subplot=441, 
            iteration='antenna',plotrange=[-1,-1,-180,180])
        plotcal(caltable=calprefix+'.bp',
            xaxis='freq',yaxis='amp',field= '3C48',subplot=441, 
            iteration='antenna',plotrange=[-1,-1,0,8])
#### gain calibration ####
## phase only solution PER SCAN (to be used to calibrate solar scans) ##
    #restrict uvrange to avoid strong RFIs
    spw0='0:4~60,1:2~30;45~60,2:2~33;55~60,3~4:2~60,5:2~38;42~60,6:2~37;41~48;52~60,7:2~60,\
            8:3~60,9~14:2~60,15:2~50'
    gaincal(vis=calhs,caltable=calprefix+'.pha_inf',field='3C48',spw=spw0,gaintype='G',\
            calmode='p',refant=refant,solint='inf',solnorm=False,\
            gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp'],\
            gainfield=['','','0','0'])
    plotcal(caltable=calprefix+'.pha_inf',xaxis='time',yaxis='phase',subplot=441,\
            spw='',antenna=antennas,iteration='antenna',plotrange=[-1,-1,-180,180])
## amp only solution PER SCAN ##
    gaincal(vis=calhs,caltable=calprefix+'.amp_inf',field='0',spw=spw0,gaintype='G',\
            calmode='a',refant=refant,solint='inf',solnorm=False,\
            gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',calprefix+'.pha_inf'],\
            gainfield=['','','0','0','0'])
    plotcal(caltable=calprefix+'.amp_inf',xaxis='time',yaxis='amp',subplot=441,\
            spw='0',antenna=antennas,iteration='antenna',plotrange=[-1,-1,-1,-1])

## amp only solution COMBING SCAN ##
    gaincal(vis=calhs,caltable=calprefix+'.amp_comb',field='0',spw=spw0,gaintype='G',\
            calmode='a',refant=refant,solint='inf',combine='scan',solnorm=False,\
            gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',calprefix+'.pha_inf'],\
            gainfield=['','','0','0','0'])

if apply:
#################################################
##### Apply calibrations to the calibrators #####
#################################################
#### 3C48 ####
    applycal(vis=calhs,field='0',
             gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',\
             calprefix+'.pha_inf',calprefix+'.amp_inf'],\
             gainfield=['','','0','0','0','0'],\
             interp=['','','','','',''],parang=False,flagbackup=True)



## 20db
antennas='ea02,ea04,ea05,ea10,ea11,ea12,ea13,ea14,ea15,ea16,ea22,ea24,ea28'

par_mbd = [ 0.119,0.233,0.012,0.010,0.095,0.077,0.031,-0.015,-0.014,-0.014,0.222,0.222,\
            0.012,0.003,0.003,-0.006,-0.008,-0.017,0.002,0.000,0.004,0.009,0.001,0.009,-0.023,-0.015]

par_pha = [ -0.7,-0.5,-1.9,-1.3,-2.6,61.5,-2.6,-2.0,-1.6,-1.5,-1.9,-1.6,-2.8,-2.2,-3.5,-1.7,-2.4,-1.8,0.3,0.2,-3.0,-2.6,-1.6,-2.0,-4.4,-4.3]

gencal(vis=calhs,caltable=calprefix+'.20db_mbd',caltype='mbd',spw='',\
       antenna=antennas,pol='R,L',parameter=par_mbd) 
gencal(vis=calhs,caltable=calprefix+'.20db_pha',caltype='ph',spw='',\
       antenna=antennas,pol='R,L',parameter=par_pha) 



## make image of sun
split(vis=msvis, outputvis='2050_S.1s.ms', timerange='20:35:00~20:50:00',datacolumn='data')
clearcal('2050_S.1s.ms')
applycal(vis='2050_S.1s.ms',field='SUN',\
         gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',\
         calprefix+'.pha_inf',calprefix+'.amp_inf',\
         calprefix+'.20db_mbd',calprefix+'.20db_pha'],\
         gainfield=['','','0','0','0','0','',''],\
         interp=['','','','','linear','linear','',''],\
         spwmap=[],applymode='',calwt=False, parang=False,flagbackup=True)

split(vis='2050_S.1s.ms', outputvis='2050_S.1s.cal.ms', datacolumn='corrected')

from suncasa.utils import dspec as ds
from suncasa.utils import svplot as sv
ds.plt_dspec(specdata="/srg/yjluo/xixi/20130423/sun_20130423_S_1s.med.spec.npz",dmin=0.005,dmax=0.3)
sv.svplot(vis='2050_S.1s.cal.ms',timerange='20:36:36~20:36:55',spw='*:2.7~3.95Ghz',specfile='sun_20130423_S_1s.med.spec.npz',dmin=0.005,dmax=0.3)




##split and test
#L-band
split(vis='/srg/yjluo/xixi/20130423/20130423_L_50ms.ms/', outputvis='20130423T2030-2050.L.50ms.ms' , timerange='20:30:00~20:50:00',datacolumn='data')
clearcal('20130423T2030-2050.L.50ms.ms')
applycal(vis='20130423T2030-2050.L.50ms.ms',field='SUN',\
         gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',\
         calprefix+'.pha_inf',calprefix+'.amp_inf',\
         calprefix+'.20db_mbd',calprefix+'.20db_pha'],\
         gainfield=['','','0','0','0','0','',''],\
         interp=['','','','','linear','linear','',''],\
         spwmap=[],applymode='',calwt=False, parang=False,flagbackup=True)

split(vis='20130423T2030-2050.L.50ms.ms', outputvis='20130423T2030-2050.L.50ms.cal.ms', datacolumn='corrected')
listobs(vis='20130423T2030-2050.L.50ms.ms',listfile='20130423T2030-2050.L.50ms.ms.listobs',overwrite=True)
from suncasa.utils import svplot as sv
sv.svplot(vis="/srg/yjluo/xixi/20130423/20130423T2030-2050.L.50ms.cal.ms/",spw=['0~1','2~3','4~5','6~7'],timerange='20:30:00~20:50:00',specfile="sun_20130423_L_1s.med.spec.npz",\
	overwrite=False,plotaia=True,mkmovie=True,twidth=240,plt_composite=False,xycen=[700, 300],fov=[512., 512.],dmin=0.000001,dmax=0.00001)

#S-band
split(vis='/srg/yjluo/xixi/20130423/20130423_S_50ms.ms/', outputvis='20130423T2030-2050.S.50ms.ms' , timerange='20:30:00~20:50:00',datacolumn='data')
clearcal('20130423T2030-2050.S.50ms.ms')
applycal(vis='20130423T2030-2050.S.50ms.ms',field='SUN',\
         gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',\
         calprefix+'.pha_inf',calprefix+'.amp_inf',\
         calprefix+'.20db_mbd',calprefix+'.20db_pha'],\
         gainfield=['','','0','0','0','0','',''],\
         interp=['','','','','linear','linear','',''],\
         spwmap=[],applymode='',calwt=False, parang=False,flagbackup=True)

split(vis='20130423T2030-2050.S.50ms.ms', outputvis='20130423T2030-2050.S.50ms.cal.ms', datacolumn='corrected')
listobs(vis='20130423T2030-2050.S.50ms.ms',listfile='20130423T2030-2050.S.50ms.ms.listobs',overwrite=True)
sv.svplot(vis="/srg/yjluo/xixi/20130423/20130423T2030-2050.S.50ms.cal.ms/",spw=['0~3','4~7','8~11','12~15'],timerange='20:30:00~20:50:00',specfile="sun_20130423_S_1s.med.spec.npz",\
	overwrite=False,plotaia=True,mkmovie=True,twidth=240,plt_composite=False,xycen=[700, 300],fov=[512., 512.],dmin=0.00001,dmax=0.0001)



## L-band image
from suncasa.utils import dspec as ds
from suncasa.utils import svplot as sv
ds.plt_dspec(specdata="/srg/yjluo/xixi/20130423/sun_20130423_L_1s.med.spec.npz",dmin=0.005,dmax=0.3)
sv.svplot(vis='/srg/yjluo/xixi/20130423/20130423T2030-2050.L.50ms.cal.ms/',timerange='20:45:17~20:45:23',spw='0~7:3~60',specfile='sun_20130423_L_1s.med.spec.npz',imagefile='L-band-image/burst2.image',dmin=0.005,dmax=0.3)
sv.svplot(vis='/srg/yjluo/xixi/20130423/20130423T2030-2050.L.50ms.cal.ms/',timerange='20:40:54~20:40:59',spw='0~7:3~60',specfile='sun_20130423_L_1s.med.spec.npz',imagefile='L-band-image/burst1.image',dmin=0.005,dmax=0.3)
from astropy.io import fits
from suncasa.utils import helioimage2fits as hf
vis='/srg/yjluo/xixi/20130423/20130423T2030-2050.L.50ms.cal.ms/'
mf=hf.read_msinfo(vis=vis)
hf.imreg(vis=vis,imagefile='L-band-image/burst2.image',fitsfile='L-band-image/burst2.fits',timerange='20:45:17~20:45:23',msinfo=mf,verbose=True,toTb=False)
hf.imreg(vis=vis,imagefile='L-band-image/burst1.image',fitsfile='L-band-image/burst1.fits',timerange='20:40:54~20:40:59',msinfo=mf,verbose=True,toTb=False)



##
ds.get_dspec(vis="/srg/yjluo/xixi/20130423/selfcal/20130423T2030-2050.L.50ms.selfcal.ms/",specfile='2030-2040L.med.spec.npz',domedian=True,spw='',timeran='')



#L-band
split(vis='/srg/yjluo/xixi/20130423/20130423_L_50ms.ms/', outputvis='20130423T2015-2030.L.50ms.ms' , timerange='20:15:00~20:30:00',datacolumn='data')
clearcal('20130423T2015-2030.L.50ms.ms')
applycal(vis='20130423T2015-2030.L.50ms.ms',field='SUN',\
         gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',\
         calprefix+'.pha_inf',calprefix+'.amp_inf',\
         calprefix+'.20db_mbd',calprefix+'.20db_pha'],\
         gainfield=['','','0','0','0','0','',''],\
         interp=['','','','','linear','linear','',''],\
         spwmap=[],applymode='',calwt=False, parang=False,flagbackup=True)

split(vis='20130423T2015-2030.L.50ms.ms', outputvis='20130423T2015-2030.L.50ms.cal.ms', datacolumn='corrected')
listobs(vis='20130423T2015-2030.L.50ms.cal.ms',listfile='20130423T2015-2030.L.50ms.cal.ms.listobs',overwrite=True)


#S-band
split(vis='/srg/yjluo/xixi/20130423/20130423_S_50ms.ms/', outputvis='20130423T2015-2030.S.50ms.ms' , timerange='20:15:00~20:30:00',datacolumn='data')
clearcal('20130423T2015-2030.S.50ms.ms')
applycal(vis='20130423T2015-2030.S.50ms.ms',field='SUN',\
         gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',\
         calprefix+'.pha_inf',calprefix+'.amp_inf',\
         calprefix+'.20db_mbd',calprefix+'.20db_pha'],\
         gainfield=['','','0','0','0','0','',''],\
         interp=['','','','','linear','linear','',''],\
         spwmap=[],applymode='',calwt=False, parang=False,flagbackup=True)

split(vis='20130423T2015-2030.S.50ms.ms', outputvis='20130423T2015-2030.S.50ms.cal.ms', datacolumn='corrected')
listobs(vis='20130423T2015-2030.S.50ms.cal.ms',listfile='20130423T2015-2030.S.50ms.cal.ms.listobs',overwrite=True)









