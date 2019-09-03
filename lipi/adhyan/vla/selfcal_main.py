## GAINCAL: Divide the data column by model column
## APPLYCAL
import sys
import numpy as np
import os

MS='/nas08-data02/rohit/20120225_selfcal_yingie/2050.1s.cal.ms' 
split_=1
if(split_==1):
    OUTPUT_MS='2050.1s.cal.fl.ms'
    os.system('rm -rf 2050.1s.cal.fl.m*')
    timerange_='20:47:09~20:47:10'
    #timerange_='20:46:09~20:46:10'
    datacolumn_='data'
    print 'Splitting the MS....'+MS
    split(vis=MS,outputvis=OUTPUT_MS,keepmms=True, \
            field="",spw="",scan="",antenna="",correlation="",timerange=timerange_, \
            intent="",array="",uvrange="",observation="",feed="", \
            datacolumn=datacolumn_,keepflags=True,width=1, \
            timebin="0s",combine="")
    print 'Splitting Finished..'+MS
    
MS='2050.1s.cal.fl.ms'
prefix='selfsol/fsun.selfcal_soln'
spw_='7'

interactive=False
modes=['ap','ap','ap','ap','ap']
#modes=['ap']*1+['p']*2+['ap']
std=[0]*len(modes)
max_=[0]*len(modes)
DR=[0]*len(modes)
image_=[0]*len(modes)
image_f=[0]*len(modes)
for step in range(len(modes)):
    # Make model 
    image_[step]=prefix+'.spw.'+spw_+'.cal.'+str(step)+'.'+str(modes[step])
    niter_=5
    npercycle_=1
    make_model=1
    if(make_model):
        print 'Making a model from CLEAN..'
        if(interactive==True):
            clean(vis=MS,imagename=image_[step], \
            outlierfile="",field="",spw=spw_,selectdata=True,timerange="", \
            uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False, \
            gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0, \
            aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="", \
            interpolation="linear",niter=niter_,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
            ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1, \
            smallscalebias=0.6,interactive=True,mask=[],nchan=-1,start=0,width=1,outframe="",veltype="radio", \
            imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter="",restfreq="", \
            stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'], \
            modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", \
            npixels=0,npercycle=npercycle_,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
            flatnoise=True,allowchunk=False)
        if(interactive==False):
            clean(vis=MS,imagename=image_[step], \
            outlierfile="",field="",spw=spw_,selectdata=True,timerange="", \
            uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False, \
            gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0, \
            aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="", \
            interpolation="linear",niter=niter_,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
            ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1, \
            smallscalebias=0.6,interactive=False,mask=[],nchan=-1,start=0,width=1,outframe="",veltype="radio", \
            imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter="",restfreq="", \
            stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'], \
            modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", \
            npixels=0,npercycle=npercycle_,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
            flatnoise=True,allowchunk=False)
    # FT
    apply_model=1
    if(apply_model):
        print 'Applying model...'+image_[step]+'.model'
        #delmod(vis=MS)
        ft(vis=MS,field="",spw=spw_,model=image_[step]+'.model',nterms=1,reffreq="",complist="",incremental=True,usescratch=True)
    # gaincal
    cal_gaincal=1
    if(cal_gaincal):
        print 'Calculating solutions..'
        caltable_=image_[step]
        gaincal(vis=MS,caltable=caltable_,field="",spw=spw_,intent="",selectdata=True,timerange="",uvrange="", 
                antenna="",scan="",observation="",msselect="", \
                solint="inf",combine="",preavg=-1.0,refant="",refantmode="flex",minblperant=4,minsnr=3.0, \
                solnorm=False,gaintype="G",smodel=[],calmode=modes[step],append=True,splinetime=3600.0, \
                npointaver=3,phasewrap=180.0,docallib=False,callib="",gaintable=[],gainfield=[],interp=[],spwmap=[],parang=False)
    # applycal
    apply_cal=1
    if(apply_cal):
        print 'Applying solutions....',image_[0:step+1]
        #clearcal(vis=MS)
        applycal(vis=MS,field="",spw=spw_,intent="",selectdata=True,timerange="",uvrange="",antenna="", \
                scan="",observation="",msselect="",docallib=False,callib="", \
                #gaintable=[caltable_], \
                gaintable=image_[0:step+1], \
                gainfield=[],interp=[], \
                spwmap=[],calwt=[True],parang=False,applymode='calonly',flagbackup=True)
    make_image=1
    image_f[step]=image_[step]+'.final'
    if(make_image):
        print 'Making image....'
        clean(vis=MS,imagename=image_f[step], \
                outlierfile="",field="",spw=spw_,selectdata=True,timerange="", \
                uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False, \
                gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0, \
                aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="", \
                interpolation="linear",niter=200,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
                ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1, \
                smallscalebias=0.6,interactive=False,mask=[],nchan=-1,start=0,width=1,outframe="",veltype="radio", \
                imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter="",restfreq="", \
                stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'], \
                modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", \
                npixels=0,npercycle=10,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
                flatnoise=True,allowchunk=False)
    img_=imstat(imagename=image_f[step]+".image/",axes=-1,region="",box="",chans="",stokes="",listit=True, \
            verbose=True,mask="",stretch=False,logfile="",append=True,algorithm="classic",fence=-1, \
            center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto",niter=3)
    std[step]=img_['rms']
    max_[step]=img_['max']
    DR[step]=max_[step]/std[step]
    print DR[step],max_[step],std[step]




