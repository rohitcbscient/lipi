# A script to run clean on a large number of time slices for a given spectral
# channel.
#
# Rohit 27Sept19


import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from suncasa.utils import helioimage2fits as hf

def define_TIMES(MSNAME,TIME_STEP):
#returns a list of time slices, taking steps defined in TIME_STEPS 
        ms.open(MSNAME)
        ms.selectinit(datadescid=0)
        tt = ms.getdata('time')
        ms.close()
        times = tt['time']#[::TIME_STEP]
        dt = float(tt['time'][1])-float(tt['time'][0])
        print "#### First time stamp %s" % (tt['time'][0])
        print "#### Integration time %f" % (dt)
        print "#### Number of time stamps %d" % (len(tt['time']))
        print "#### Last time stamp %s" % (tt['time'][len(tt['time'])-1])
        print "#### No. of time stamps to image %d" % (len(times))
        return times

###### im.open 
FBAND='L'
#MSNAME='sun_L_20160409T1844-1846UT.50ms.cal.ms'
MSNAME='sun_L_20160409T184000-184400UT.50ms.cal.ms'
###### im.defineimage 
IMSIZE=256
CELLSIZE='2.0arcsec'
STOKES='RR'
MODE='mfs'
NCHAN=1         # No. of channels in output image
START_CHAN=0    # Start ch. in the band
WIDTH=1        # How many channels to use from START_CHAN

UVRANGE_IN=''
###### im.weight
#TYPE='superuniform' 
WEIGHTING='natural'
###### im.clean 
SPW_='2'#'0:32~63'#'0:45' 
#ALGORITHM='cs' ## 'cs' was used when working with SanB
INTERPOLATION='linear'
NITER=500
GAIN=0.1
#MULTISCALE=[0,6,12]
#MULTISCALE=[0,5,10,20,30]
MULTISCALE=[]
INTERACTIVE=False
NPERCYCLE=1
NAME=MSNAME.split('.ms')[0]+'.pol.'+STOKES

#MASK=NAME[0]+'.mask.img'
#PSF=NAME[0]+'.psf.img'

TIME_STEP = 1
TIME_LIST=define_TIMES(MSNAME,TIME_STEP)# Get the time stamps from the MS
# It will be a better design to get all the time stamps and then later choose the 
# period and the sampling with which to do the imaging.


############ TIME
#TIME_LIST=TIME_LIST[0]+np.arange(120)

#i=0
#time_list=[0]*(TIME_LIST.shape[0])
#time_pair=[0]*(TIME_LIST.shape[0]-1)
#imagename=[0]*(TIME_LIST.shape[0]-1)
#for i in range(TIME_LIST.shape[0]):
#        timed=TIME_LIST[i]-TIME_LIST[0]
#        reftime=qa.time(qa.quantity(TIME_LIST[0],'s'),form="ymd")
#        reftime_sec=int(reftime[0].split('/')[3].split(':')[0])*3600+int(reftime[0].split('/')[3].split(':')[1])*60+int(reftime[0].split('/')[3].split(':')[2])
#        time_list[i]=time.strftime('%H:%M:%S', time.gmtime(np.round(timed+reftime_sec,2)))

tres=0.05
total_images=4800
TIME_LIST=TIME_LIST[0]+np.arange(total_images)*tres
first_time='18:40:00'
#first_time='18:44:00'
#first_time='22:04:00'
i=0
time_list=[0]*(TIME_LIST.shape[0])
time_pair=[0]*(TIME_LIST.shape[0]-1)
imagename=[0]*(TIME_LIST.shape[0]-1)
for i in range(TIME_LIST.shape[0]):
        timed=TIME_LIST[i]-TIME_LIST[0]
        #reftime=qa.time(qa.quantity(TIME_LIST[0],'s'),form="ymd")
        #reftime=qa.time(qa.quantity(1,'s'),form="ymd")
        #reftime_sec=int(reftime[0].split('/')[3].split(':')[0])*3600+int(reftime[0].split('/')[3].split(':')[1])*60+int(reftime[0].split('/')[3].split(':')[2])
        reftime_sec=int(first_time.split(':')[0])*3600+int(first_time.split(':')[1])*60+int(first_time.split(':')[2])
        #time_list[i]=time.strftime('%H:%M:%S', time.gmtime(np.round(timed+reftime_sec+1,2)))
        hour,min_,sec=int(first_time.split(':')[0]),int(first_time.split(':')[1]),int(first_time.split(':')[2])
        microsec=1.e6*(tres*i)
        if(microsec>=int(1.e6)):
            n=int(microsec/1.e6)
            microsec=int((microsec/1.e6-int(microsec/1.e6))*1.e6)
            sec=sec+n
            if(sec>=60):
                m=int(sec/60)
                sec=sec-m*60;min_=min_+m
        else:
            microsec=int(microsec)
        time_=datetime.datetime(2016, 4, 9, hour, min_, sec, microsec)
        time_list[i]=time_.strftime("%H:%M:%S.%f").rstrip('0')


#for i in range(TIME_LIST.shape[0]/2):
#        time_list[2*i+1]=time_list[2*i].split(':')[0]+':'+time_list[2*i].split(':')[1]+':'+str(float(time_list[2*i].split(':')[2])+0.5)
for i in range(TIME_LIST.shape[0]-1):
        time_pair[i]=time_list[i]+'~'+time_list[i+1]

############ FREQ
ave_chan=16
tot_chan=64
#tot_chan=128
fband_list=['0','1','2','3','4','5','6','7']
spw_list1=np.arange(tot_chan/ave_chan)*ave_chan
spw_list2=np.arange(tot_chan/ave_chan)*ave_chan+ave_chan-1
spw_list=[0]*len(spw_list1)*len(fband_list)
j=0
for fband in fband_list:
    for i in range(len(spw_list1)):
            spw_list[j]=fband+':'+str(spw_list1[i])+'~'+str(spw_list2[i])
            j=j+1

#TIME_LIST=TIME_LIST[62:]
#time_pair=time_pair[62:] # SELECT THE FLARE TIME FOR 1 SEC
#sys.exit()
#TIME_LIST=TIME_LIST[220:230]
#time_pair=time_pair[220:230]  # FOR AIA 12 sec average
rotate_crop=1
#phasecenter_='J2000 01h14m56.0  +8d00m47.084' # 12 sec
#phasecenter_='J2000 01h14m44.0  +7d56m00.0' # 1sec
#phasecenter_c='J2000 01h14m04.7s  +7d47m42.5s' # 1sec for the full map
#phasecenter_='J2000 01h14m43.863s  +07d56m01.107s'
#phasecenter_='J2000 01h14m43.820s  +7d56m02.709s' # 1sec for Tb submap (141,609) pixel: (-582",500")
#phasecenter_='J2000 01h14m43.446s  +7d55m59.701s' # 1sec for Tb submap (141,609) pixel: (-582",500")
phasecenter_='01h14m57.902s  +08d00m59.004s' # New MS
antennas=''#'ea02,ea03,ea04,ea05,ea06,ea07,ea08,ea09,\
date='2016/04/09'
mf=hf.read_msinfo(vis=MSNAME)
j=0#494*3
k=5
#sys.exit()
for i in range(TIME_LIST.shape[0]):
    ss=0
    for s in spw_list[4*k:4*(k+1)]:
    #for s in [spw_list[4*k+2]]:
            sname=s.split(':')[0]+'_'+s.split(':')[1].split('~')[0]+'-'+s.split(':')[1].split('~')[1]
            tname=time_pair[i].split('~')[0]+'-'+time_pair[i].split('~')[1]
            imagename[i]='images_RR_3/'+NAME+'.spw.'+sname+'.time.'+tname
            #imagename[i]='try/'+NAME+'.spw.'+s+'.time.'+time_pair[i]
            OUTNAME=imagename[i]+'.FITS'
            if (os.path.isfile(OUTNAME)==False):
                print 'Cleaning..'+imagename[i]
                tclean(vis=MSNAME,selectdata=True,field="",spw=s,timerange=time_pair[i], \
                        uvrange="",antenna=antennas,scan="",observation="",intent="",datacolumn="data", \
                        imagename=imagename[i],imsize=IMSIZE,cell=CELLSIZE,phasecenter=phasecenter_, \
                        stokes=STOKES,projection="SIN",startmodel="", \
                specmode="mfs",reffreq="",nchan=-1,start="",width="", \
                outframe="LSRK",veltype="radio",restfreq=[],interpolation="linear",gridder="standard",facets=1,chanchunks=1, \
                wprojplanes=1,vptable="",aterm=True,psterm=False,wbawp=True,conjbeams=True,cfcache="",computepastep=360.0, \
                rotatepastep=360.0,pblimit=0.2,normtype="flatnoise",deconvolver="hogbom",scales=[],nterms=2,smallscalebias=0.6, \
                restoration=True,restoringbeam='',pbcor=False,outlierfile="",weighting="natural",robust=0.5,npixels=0,uvtaper=[],niter=NITER,gain=0.1, \
                threshold="0.0mJy",nsigma=0.0,cycleniter=-1,cyclefactor=1.5,minpsffraction=0.05,maxpsffraction=0.8,interactive=False, \
                usemask="user",mask=[],pbmask=0.0,maskthreshold="",maskresolution="",nmask=0,sidelobethreshold=3.0,noisethreshold=5.0, \
                lownoisethreshold=1.5,negativethreshold=0.0,smoothfactor=1.0,minbeamfrac=0.3, \
                cutthreshold=0.01,growiterations=75,dogrowprune=True,minpercentchange=-1.0,verbose=False,restart=True,savemodel="none", \
                calcres=True,calcpsf=True,parallel=False)
                os.system('rm -rf '+imagename[i]+'.psf')   #remove files that are not wanted
                os.system('rm -rf '+imagename[i]+'.flux')
                os.system('rm -rf '+imagename[i]+'.model')
                os.system('rm -rf '+imagename[i]+'.mask')
                os.system('rm -rf '+imagename[i]+'.residual')
                os.system('rm -rf '+imagename[i]+'.pb')
                os.system('rm -rf '+imagename[i]+'.sumwt')
                if(rotate_crop==1):
                    time_fits=date+' '+time_pair[i].split('~')[0]+' ~ '+date+' '+time_pair[i].split('~')[1]
                    hf.imreg(vis=MSNAME,imagefile=imagename[i]+'.image',fitsfile=OUTNAME,timerange=time_fits ,msinfo=mf,verbose=True,toTb=True,scl100=True,p_ang=True,usephacenter=phasecenter_)
                    #os.system('rm -rf '+imagename[i]+'.image')
            ss=ss+1
