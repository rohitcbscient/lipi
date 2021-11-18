# A script to run clean on a large number of time slices for a given spectral
# channel.
#
# Divya 17Apr13

# TO DO
# Figure out how to skip over flagged time stamps in an automated manner

import matplotlib.pyplot as plt
import numpy as np
import os 
import sys

def define_TIMES(MSNAME,TIME_STEP):
#returns a list of time slices, taking steps defined in TIME_STEPS 
	ms.open(MSNAME)
	ms.select({'antenna1':[1],'antenna2':[1]})
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
MSNAME='1238472032_089-089.ms'
###### im.defineimage 
IMSIZE=4096
CELLSIZE='20arcsec'
STOKES='I'
MODE='mfs'
NCHAN=1 	# No. of channels in output image
START_CHAN=0	# Start ch. in the band
WIDTH=16	# How many channels to use from START_CHAN

#NCHAN_OUT=1
#START_OUT=0
###### im.selectvis 
#NCHAN_IN=16
#START_IN=0
#STEP_IN=1
####TIME_IN='04:11:02~04:16:58' ## The variable to iterate over
#UVRANGE_IN='10~1000lambda'
UVRANGE_IN=''#'10~800lambda'
###### im.weight
#TYPE='superuniform' 
WEIGHTING='superuniform' 
###### im.clean 
SPW=''#'0:32~63'#'0:45' 
#ALGORITHM='cs' ## 'cs' was used when working with SanB
INTERPOLATION='linear'
NITER=2000
GAIN=0.1
#MULTISCALE=[0,6,12]
#MULTISCALE=[0,5,10,20,30]
MULTISCALE=[]
INTERACTIVE=False
NPERCYCLE=1
NAME=MSNAME.split('.')

#MASK=NAME[0]+'.mask.img'
#PSF=NAME[0]+'.psf.img'

TIME_STEP = 1
TIME_LIST=define_TIMES(MSNAME,TIME_STEP)# Get the time stamps from the MS
# It will be a better design to get all the time stamps and then later choose the 
# period and the sampling with which to do the imaging.

i=0
time_list=[0]*(TIME_LIST.shape[0])
time_pair=[0]*(TIME_LIST.shape[0]-1)
imagename=[0]*(TIME_LIST.shape[0]-1)
for i in range(TIME_LIST.shape[0]):
	timed=TIME_LIST[i]-TIME_LIST[0]
	reftime=qa.time(qa.quantity(TIME_LIST[0],'s'),form="ymd")
	reftime_sec=int(reftime[0].split('/')[3].split(':')[0])*3600+int(reftime[0].split('/')[3].split(':')[1])*60+int(reftime[0].split('/')[3].split(':')[2])
	time_list[i]=time.strftime('%H:%M:%S', time.gmtime(np.round(timed+reftime_sec,2)))	
for i in range(TIME_LIST.shape[0]/2):
	time_list[2*i+1]=time_list[2*i].split(':')[0]+':'+time_list[2*i].split(':')[1]+':'+str(float(time_list[2*i].split(':')[2])+10)
for i in range(TIME_LIST.shape[0]-1):	
	time_pair[i]=time_list[i]+'~'+time_list[i+1]
#for i in range(TIME_LIST.shape[0]-2):
j=0#494*3
for i in range(TIME_LIST.shape[0]):
	index = '%04d' % (i+j)
	#imagename[i]='images/239MHz.'+index
	imagename[i]='images/'+NAME[0]+'.'+index
	print 'Cleaning..'+index+'..'+imagename[i]	
	clean(vis=MSNAME,imagename=imagename[i],outlierfile="",field="",spw=SPW,selectdata=True,timerange=time_pair[i],uvrange=UVRANGE_IN, antenna="",scan="",observation="",mode=MODE,gridmode="",wprojplanes=1,facets=1,cfcache="cfcache.dir",painc=360.0,epjtable="",interpolation=INTERPOLATION,niter=NITER,gain=GAIN,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean",ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=MULTISCALE,negcomponent=-1,smallscalebias=0.6,interactive=INTERACTIVE,mask=[],nchan=NCHAN,start=START_CHAN,width=WIDTH,outframe="",veltype="radio",imsize=[IMSIZE, IMSIZE],cell=[CELLSIZE, CELLSIZE],phasecenter="",restfreq="",stokes=STOKES,weighting=WEIGHTING,robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",npixels=0,npercycle=NPERCYCLE,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=True,flatnoise=True,allowchunk=False) 


