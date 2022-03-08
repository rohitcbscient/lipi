##Also set this as general script

#### SCRIPT for calibrating the 20130423 data at L Band ####
import os

mspath = '/srg/yjluo/xixi/20130423/'
msvis = '/srg/yjluo/xixi/20130423/20130423_L_1s.ms/'
vispre='20130423_L_1s.ms'
calhs='20130423_L_1s.ms.hs.ms'
calprefix='cal-L/calSUN_20130423-L'

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
    flagcmd(vis=msvis, inpmode='list', inpfile="/srg/yjluo/xixi/20130423/20130423_L_50ms.ms.onlineflags")
    #### flagdata ####
    flagdata(vis=msvis, scan='1~6, 8, 28, 30, 50, 52')
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
flagdata(vis=calhs, mode='tfcrop', spw='0',field='3C48',scan='7,29,51'
         datacolumn='data', action='calculate',
         display='both', flagbackup=False)


#flaginfo
flagInfo = flagdata(vis=calhs, mode='summary')

## Refant Antenna ##
refant='ea07'
antennas='ea01,ea03,ea06,ea07,ea08,ea09,ea17,ea18,ea19,ea21,ea23,ea25,ea26,ea27'

##calibration##
if docalib:
#### initial phase-only solution on the flux calibrator t-p#### 
#the ant ea23 RR in spw4~7 have unusual phase configuration, ant ea25 has very low amp in LL
    flagdata(vis=calhs,field='3C48',spw='4~7',antenna='ea23',correlation='RR')
    flagdata(vis=calhs,field='3C48',spw='0~7',antenna='ea25',correlation='LL')
    spw1='0:20~24,1:20~24,2:40~44,3:28~32,4:10~14,5:28~32,6:28~32,7:14~18'
    gaincal(vis=calhs, caltable=calprefix+'.pha0', field='3C48',\
            spw=spw1, gaintype='G',gaintable=[calprefix+'.antpos',calprefix+'.rq'],refant=refant, calmode='p', solint='10s', minsnr=3)
    if plotipcal:
        plotcal(caltable=calprefix+'.pha0',xaxis='time',yaxis='phase',subplot=441,\
                spw='',antenna=antennas,iteration='antenna',plotrange=[-1,-1,-180,180])

#### delay calibration f-p####
    spwd='0:2~16;20~28;32~46;54~60,1:2~5;10~49;55~59,2:2~6;17~60,\
            3:2~60,4:2~14;37~38;44~51;58~60,5:5~12;15~60,6:2~60,7:2~19;24~26;38~41;45~48;55~60'
    gaincal(vis=calhs, caltable=calprefix+'.delay', field='3C48', spw=spwd, \
            gaintype='K', gaintable=[calprefix+'.antpos',calprefix+'.rq', calprefix+'.pha0'],\
            gainfield=['','','3C48'],refant=refant, \
            combine='scan',solint='inf', minsnr=5)
    plotcal(caltable=calprefix+'.delay',xaxis='freq',yaxis='delay',iteration='antenna',subplot=441,plotrange= [-1, -1, -10, 10])
    listcal(vis=calhs,caltable=calprefix+'.delay')
#### bandpass calibration ####
    #spwbp='0:0~15;21~46;51~63,1:0~39;45~62,2~3,4:0~13;22~26;35~38;43~63,5:4~63,6~7'
    spwbp='0:0~16;20~28;32~46;54~60,1:0~5;10~49;55~59,2:0~6;17~62,\
            3:0~62,4:0~14;35~38;44~51;58~60,5:5~12;15~62,6:0~62,7:0~19;24~26;38~41;45~48;55~62'
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
            iteration='antenna',plotrange=[-1,-1,0,15])
#### gain calibration ####
## phase only solution PER SCAN (to be used to calibrate solar scans) ##
    #restrict uvrange to avoid strong RFIs
    spw0='0:0~16;20~28;32~46;54~60,1:0~5;10~49;55~59,2:0~6;17~62,\
            3:0~62,4:0~14;35~38;44~51;58~60,5:5~12;15~62,6:0~62,7:0~19;24~26;38~41;45~48;55~62'
    gaincal(vis=calhs,caltable=calprefix+'.pha_inf',field='3C48',spw=spw0,gaintype='G',uvrange='>200m',\
            calmode='p',refant=refant,solint='inf',solnorm=False,\
            gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp'],\
            gainfield=['','','0','0'])
        plotcal(caltable=calprefix+'.pha_inf',xaxis='time',yaxis='phase',subplot=441,\
                spw='',antenna=antennas,iteration='antenna',plotrange=[-1,-1,-180,180])
## amp only solution PER SCAN ##
    gaincal(vis=calhs,caltable=calprefix+'.amp_inf',field='0',spw=spw0,gaintype='G',uvrange='>200m',\
            calmode='a',refant=refant,solint='inf',solnorm=False,\
            gaintable=[calprefix+'.antpos',calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',calprefix+'.pha_inf'],\
            gainfield=['','','0','0','0'])
        plotcal(caltable=calprefix+'.amp_inf',xaxis='time',yaxis='amp',subplot=441,\
                spw='0',antenna=antennas,iteration='antenna',plotrange=[-1,-1,-1,-1])

## amp only solution COMBING SCAN ##
    gaincal(vis=calhs,caltable=calprefix+'.amp_comb',field='0',spw=spw0,gaintype='G',uvrange='>200m',\
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









