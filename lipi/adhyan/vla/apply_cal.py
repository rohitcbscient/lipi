import sys

MS='20120225_20min.ms'
split_=1
if(split_==1):
    OUTPUT_MS='2050.1s.cal.qs.ms'
    os.system('rm -rf 2050.1s.cal.qs.m*')
    timerange_='20:47:09~20:47:10'
    datacolumn_='corrected'
    print 'Splitting the MS....'+MS
    split(vis=MS,outputvis=OUTPUT_MS,keepmms=True, \
            field="",spw="",scan="",antenna="",correlation="",timerange=timerange_, \
            intent="",array="",uvrange="",observation="",feed="", \
            datacolumn=datacolumn_,keepflags=True,width=1, \
            timebin="0s",combine="")
    print 'Splitting Finished..'+MS

MS='2050.1s.cal.qs.ms'
origin_apply=0
if(origin_apply):
    print 'Applying clearcal....'
    clearcal(MS)
    print 'Applying original solutions.........'
    calprefix='yingjie_cal/calSUN_20120225'
    applycal(vis=MS,field='SUN',\
             gaintable=[calprefix+'.rq',calprefix+'.delay',calprefix+'.bp',\
             calprefix+'.pha_inf',calprefix+'.amp_inf',calprefix+'.amp_Finc',\
             calprefix+'.20db_mbd',calprefix+'.20db_pha'],\
             gainfield=['','3','3','0','0','0','',''],\
             interp=['','','','linear','linear','linear','',''],\
             spwmap=[],applymode='',calwt=False, parang=False,flagbackup=True)

apply_self=0
s=7
gaintable_=['selfsol/fsun.selfcal_soln.spw.7.cal.1.ap']
if(apply_self):
    print 'Applying selfcal solutions on spw:'+str(s)
    applycal(vis=MS,field="",spw=str(s),intent="",selectdata=True,timerange="",\
            uvrange="",antenna="",scan="",observation="",msselect="", \
            docallib=False,callib="", \
            gaintable=gaintable_,gainfield=[], \
            interp=[],spwmap=[],calwt=[True],parang=False, \
            applymode="calonly",flagbackup=True)

make_image=1
image_f='try'
if(make_image):
    print 'Making image....'
    clean(vis=MS,imagename=image_f, \
            outlierfile="",field="",spw=str(s),selectdata=True,timerange="", \
            uvrange="",antenna="",scan="",observation="",intent="",mode="channel",resmooth=False, \
            gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0, \
            aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="", \
            interpolation="linear",niter=200,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
            ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1, \
            smallscalebias=0.6,interactive=False,mask=[],nchan=32,start=0,width=4,outframe="",veltype="radio", \
            imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter="",restfreq="", \
            stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'], \
            modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", \
            npixels=0,npercycle=10,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=True, \
            flatnoise=True,allowchunk=False)

img_=imstat(imagename=image_f+".image/",axes=-1,region="",box="",chans="",stokes="",listit=True, \
        verbose=True,mask="",stretch=False,logfile="",append=True,algorithm="classic",fence=-1, \
        center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto",niter=3)
std=img_['rms']
max_=img_['max']
DR=max_/std
print DR,max_,std

sys.exit()

sys.exit()

s=5
gaintable_=['selfcal_soln/fsun.selfcal_soln.spw.5.cal.1.ap']
print 'Applying selfcal solutions on spw:'+str(s)
applycal(vis=MS,field="",spw=str(s),intent="",selectdata=True,timerange="",\
	uvrange="",antenna="",scan="",observation="",msselect="", \
	docallib=False,callib="", \
	gaintable=gaintable_,gainfield=[], \
	interp=[],spwmap=[],calwt=[True],parang=False, \
	applymode="calonly",flagbackup=True)

sys.exit()
s=6
gaintable_=['selfcal_soln/fsun.selfcal_soln.spw.6.cal.4.ap']
print 'Applying selfcal solutions on spw:'+str(s)
applycal(vis=MS,field="",spw=str(s),intent="",selectdata=True,timerange="",\
	uvrange="",antenna="",scan="",observation="",msselect="", \
	docallib=False,callib="", \
	gaintable=gaintable_,gainfield=[], \
	interp=[],spwmap=[],calwt=[True],parang=False, \
	applymode="calonly",flagbackup=True)
