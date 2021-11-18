## SCRIPT TO SELFCAL ON CALIBRATED VLA DATA

## GAINCAL: Divide the data column by model column
## APPLYCAL
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits

MS='/nas08-data02/vladata/20130423/L-Band/selfcal/20130423T2030-2050.L.50ms.selfcal.ms' 
split_=1
if(split_==1):
    OUTPUT_MS='20130423T203430-203431.L.50ms.ms'
    os.system('rm -rf 20130423T203430-203431.L.50ms.m*')
    #timerange_='20:45:00~20:45:30'
    timerange_='20:34:30~20:34:31' # MARINA timestamps
    #timerange_='20:46:09~20:46:10'
    datacolumn_='data'
    print 'Splitting the MS....'+MS
    split(vis=MS,outputvis=OUTPUT_MS,keepmms=True, \
            field="",spw="",scan="",antenna="",correlation="",timerange=timerange_, \
            intent="",array="",uvrange="",observation="",feed="", \
            datacolumn=datacolumn_,keepflags=True,width=1, \
            timebin="0s",combine="")
    print 'Splitting Finished..'+MS
   


# INPUT PARAMETERS
MS='20130423T203430-203431.L.50ms.ms'
DIR='selfsol'
prefix=DIR+'/fsun.selfcal_soln'
spw_list=['0','1','2','3','4','5','6','7']
interactive=False
modes=['p','p','p','a','p','p','p','p','p']#['ap','ap','ap','ap','ap','ap','ap']
niter_scale_factor=np.ones(len(spw_list))#np.linspace(2,1,len(spw_list))
#niter_scale_factor[1]=1.0#4.0
niter_base=np.array([100,200,400,400,400,400,400,400,400])#[,600,800,1000,1200,1300,1000,1000])
# RA:02 05 40.74 DEC:+12 44 01.6 20:34:00
#phase_centre='J2000 0.548505138896rad 0.222291421955rad'
phase_centre='J2000 02h05m40.74 +12d44m01.6'
do_organise=0
#modes=['ap']*1+['p']*2+['ap']

# SCRIPTS STARTS
ss=0
for spw_ in [spw_list[7]]:
    niter_list=np.array(niter_base/niter_scale_factor[ss],dtype=int)
    print '##### STARTING SELFCAL FOR SPW:'+str(spw_)+' #####'
    print 'Niter list: ', niter_list
    std=[0]*len(modes)
    max_=[0]*len(modes)
    DR=[0]*len(modes)
    image_=[0]*len(modes)
    model_=[0]*(len(modes)+1)
    image_f=[0]*len(modes)
    model_[0]=[]
    for step in range(len(modes)):
        # Make model 
        image_[step]=prefix+'.spw.'+spw_+'.cal.'+str(step)+'.'+str(modes[step])
        mod=image_[step]+'.model'
        model_[step+1]=mod#+model_[step]
        print model_
        niter_=niter_list[step]
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
                imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter=phase_centre,restfreq="", \
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
                smallscalebias=0.6,interactive=False,mask=['mask0.mask'],nchan=-1,start=0,width=1,outframe="",veltype="radio", \
                imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter=phase_centre,restfreq="", \
                stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'], \
                modelimage='',restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", \
                npixels=0,npercycle=npercycle_,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
                flatnoise=True,allowchunk=False)
            os.system('rm -rf '+image_[step]+'.residual')
            os.system('rm -rf '+image_[step]+'.mask')
            os.system('rm -rf '+image_[step]+'.flux')
            os.system('rm -rf '+image_[step]+'.psf')

        # FT
        apply_model=1
        if(apply_model):
            print 'Applying model...'+str(model_[step+1])+'.model'
            #delmod(vis=MS)
            ft(vis=MS,field="",spw=spw_,model=model_[step+1],nterms=1,reffreq="",complist="",incremental=False,usescratch=True)
        # gaincal
        cal_gaincal=1
        if(cal_gaincal):
            print 'Calculating solutions..'
            caltable_=image_[step]
            gaincal(vis=MS,caltable=caltable_,field="",spw=spw_,intent="",selectdata=True,timerange="",uvrange="", 
                    antenna="",scan="",observation="",msselect="", \
                    solint="inf",combine="",preavg=-1.0,refant="",refantmode="flex",minblperant=4,minsnr=3.0, \
                    solnorm=False,gaintype="G",smodel=[],calmode=modes[step],append=False,splinetime=3600.0, \
                    npointaver=3,phasewrap=180.0,docallib=False,callib="",gaintable=[],gainfield=[],interp=[],spwmap=[],parang=False)
        # plotcal
        plot_cal_amp=1
        if(plot_cal_amp):
            plotcal(caltable=caltable_,xaxis="antenna",yaxis="amp",poln="L",
                    field="",antenna="",spw="",timerange="",subplot=111,overplot=False,clearpanel="Auto",
                    iteration="",plotrange=[-1,-1,0,2],showflags=False,plotsymbol="o",plotcolor="blue",markersize=5.0,
                    fontsize=10.0,showgui=False,figfile=caltable_+'_amp_ant.png')
        plot_cal_phase=1
        if(plot_cal_phase):
            plotcal(caltable=caltable_,xaxis="antenna",yaxis="phase",poln="L",
                    field="",antenna="",spw="",timerange="",subplot=111,overplot=False,clearpanel="Auto",
                    iteration="",plotrange=[-1,-1,-180,180],showflags=False,plotsymbol="o",plotcolor="blue",markersize=5.0,
                    fontsize=10.0,showgui=False,figfile=caltable_+'_phase_ant.png')
        plot_cal_amp_chan=0
        if(plot_cal_amp_chan):
            plotcal(caltable=caltable_,xaxis="chan",yaxis="amp",poln="L",
                    field="",antenna="",spw="",timerange="",subplot=431,overplot=False,clearpanel="Auto",
                    iteration="antenna",plotrange=[-1,-1,0,1],showflags=False,plotsymbol="o",plotcolor="blue",markersize=5.0,
                    fontsize=10.0,showgui=False,figfile=caltable_+'_amp_chan.png')
        plot_cal_phase_chan=0
        if(plot_cal_phase_chan):
            plotcal(caltable=caltable_,xaxis="chan",yaxis="phase",poln="L",
                    field="",antenna="",spw="",timerange="",subplot=431,overplot=False,clearpanel="Auto",
                    iteration="antenna",plotrange=[-1,-1,0,1],showflags=False,plotsymbol="o",plotcolor="blue",markersize=5.0,
                    fontsize=10.0,showgui=False,figfile=caltable_+'_phase_chan.png')
        # Delmod
        del_mod=1
        if(del_mod):
            print 'Deleting old model...'
            delmod(vis=MS)
        # clearcal
        clear_cal=1
        if(clear_cal):
            clearcal(vis=MS)

        # Apply solutions
        apply_cal=1
        if(apply_cal):
            print 'Applying solutions....',[caltable_]#image_[0:step+1]
            #clearcal(vis=MS)
            applycal(vis=MS,field="",spw=spw_,intent="",selectdata=True,timerange="",uvrange="",antenna="", \
                    scan="",observation="",msselect="",docallib=False,callib="", \
                    gaintable=[caltable_], \
                    #gaintable=image_[0:step+1], \
                    gainfield=[],interp=[], \
                    spwmap=[],calwt=[False],parang=False,applymode='calonly',flagbackup=True)

        make_image=1
        image_f[step]=image_[step]+'.final'
        if(make_image):
            print 'Making image....'
            clean(vis=MS,imagename=image_f[step], \
                    outlierfile="",field="",spw=spw_,selectdata=True,timerange="", \
                    uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False, \
                    gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0, \
                    aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="", \
                    interpolation="linear",niter=800,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
                    ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1, \
                    smallscalebias=0.6,interactive=False,mask=['mask0.mask'],nchan=-1,start=0,width=1,outframe="",veltype="radio", \
                    imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter=phase_centre,restfreq="", \
                    stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'], \
                    modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", \
                    npixels=0,npercycle=10,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
                    flatnoise=True,allowchunk=False)
            os.system('rm -rf '+image_f[step]+'.model')
            os.system('rm -rf '+image_f[step]+'.residual')
            os.system('rm -rf '+image_f[step]+'.mask')
            os.system('rm -rf '+image_f[step]+'.flux')
            os.system('rm -rf '+image_f[step]+'.psf')
        img_=imstat(imagename=image_f[step]+".image/",axes=-1,region="",box="",chans="",stokes="",listit=True, \
                verbose=True,mask="",stretch=False,logfile="",append=True,algorithm="classic",fence=-1, \
                center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto",niter=3)
        img_std=imstat(imagename=image_f[step]+".image/",axes=-1,region="",box="10,10,200,200",chans="",stokes="",listit=True, \
                verbose=True,mask="",stretch=False,logfile="",append=True,algorithm="classic",fence=-1, \
                center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto",niter=3)
        std[step]=img_std['rms']
        max_[step]=img_['max']
        DR[step]=max_[step]/std[step]
        print DR[step],max_[step],std[step]
        # split the corrected data  
        split_corrdata=1
        if(split_corrdata):
            print 'Splitting the corrected data..'
            corrMS='corr_'+MS
            os.system('rm -rf '+corrMS)
            split(vis=MS,outputvis=corrMS,datacolumn='corrected')
        os.system('rm -rf '+MS)
        os.system('mv '+corrMS+' '+MS)
        ##############
        # Apply solutions to the final MS
        # Below steps are WRONG! As you have already applied solutions and previous solutions are already applied when splitting corrected data..
        do_clearcal=0
        if(do_clearcal):
            clearcal(MS)
        apply_cal=0
        if(apply_cal):
            print 'Applying solutions....',image_[0:step+1]
            #clearcal(vis=MS)
            applycal(vis=MS,field="",spw=spw_,intent="",selectdata=True,timerange="",uvrange="",antenna="", \
                    scan="",observation="",msselect="",docallib=False,callib="", \
                    #gaintable=image_, \
                    gaintable=image_[0:step+1], \
                    gainfield=[],interp=[], \
                    spwmap=[],calwt=[False],parang=False,applymode='calonly',flagbackup=True)
        make_images_edges=0
        image_lower=image_[step]+'.lower' # Lower in frequency
        if(make_image):
            print 'Making image....'
            clean(vis=MS,imagename=image_lower, \
                    outlierfile="",field="",spw=spw_+':0~3',selectdata=True,timerange="", \
                    uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False, \
                    gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0, \
                    aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="", \
                    interpolation="linear",niter=800,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
                    ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1, \
                    smallscalebias=0.6,interactive=False,mask=['mask0.mask'],nchan=-1,start=0,width=1,outframe="",veltype="radio", \
                    imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter=phase_centre,restfreq="", \
                    stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'], \
                    modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", \
                    npixels=0,npercycle=10,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
                    flatnoise=True,allowchunk=False)
            os.system('rm -rf '+image_lower+'.model')
            os.system('rm -rf '+image_lower+'.residual')
            os.system('rm -rf '+image_lower+'.mask')
            os.system('rm -rf '+image_lower+'.flux')
            os.system('rm -rf '+image_lower+'.psf')
        image_upper=image_[step]+'.upper' # Lower in frequency
        if(make_image):
            print 'Making image....'
            clean(vis=MS,imagename=image_upper, \
                    outlierfile="",field="",spw=spw_+':60~63',selectdata=True,timerange="", \
                    uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False, \
                    gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0, \
                    aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="", \
                    interpolation="linear",niter=800,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
                    ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1, \
                    smallscalebias=0.6,interactive=False,mask=['mask0.mask'],nchan=-1,start=0,width=1,outframe="",veltype="radio", \
                    imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter=phase_centre,restfreq="", \
                    stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'], \
                    modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", \
                    npixels=0,npercycle=10,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
                    flatnoise=True,allowchunk=False)
            os.system('rm -rf '+image_upper+'.model')
            os.system('rm -rf '+image_upper+'.residual')
            os.system('rm -rf '+image_upper+'.mask')
            os.system('rm -rf '+image_upper+'.flux')
            os.system('rm -rf '+image_upper+'.psf')

    plt.close()
    plt.plot(np.arange(len(DR))+1,DR,'o-')
    plt.xlabel('Iterations')
    plt.ylabel('Dynamic Range')
    plt.savefig(prefix+'_DR_spw'+spw_+'.png')
    plt.close()
    ss=ss+1

if(do_organise):
    os.system('rm -rf '+DIR+'/selfcal_plots')
    os.system('mkdir '+DIR+'/selfcal_plots')
    os.system('mv '+DIR+'/*.png '+DIR+'/selfcal_plots')
    os.system('rm -rf '+DIR+'/selfcal_images')
    os.system('mkdir '+DIR+'/selfcal_images')
    os.system('mv '+DIR+'/*upper.image '+DIR+'/selfcal_images')
    os.system('mv '+DIR+'/*lower.image '+DIR+'/selfcal_images')
    os.system('mv '+DIR+'/*final.image '+DIR+'/selfcal_images')
    os.system('rm -rf '+DIR+'/selfcal_models')
    os.system('mkdir '+DIR+'/selfcal_models')
    os.system('mv '+DIR+'/*.model '+DIR+'/selfcal_models')
    os.system('mv '+DIR+'/*.image '+DIR+'/selfcal_models')




