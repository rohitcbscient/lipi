## SCRIPT TO SELFCAL ON CALIBRATED VLA DATA

## GAINCAL: Divide the data column by model column
## APPLYCAL
import sys
import numpy as np
import os
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from suncasa.utils import helioimage2fits

#MOMMS='/nas08-data02/vladata/20130423/L-Band/selfcal/20130423T2030-2050.L.50ms.selfcal.ms' 
MOMMS='/media/rohit/VLA/20130423_sub_selfcal_new/flare1/20130423T202641.5-202642.5.L.50ms.cal.sub.ms' 
OUTPUT_MS=MOMMS
split_=0
if(split_==1):
    OUTPUT_MS='20130423T203430-203431.S.50ms.cal.ms'
    os.system('rm -rf '+OUTPUT_MS)
    #timerange_='20:40:00~20:46:00'
    timerange_='20:34:30~20:34:31' # MARINA timestamps
    #timerange_='20:46:09~20:46:10'
    datacolumn_='data'
    print 'Splitting the MS....'+MOMMS
    split(vis=MOMMS,outputvis=OUTPUT_MS,keepmms=True, \
            field="",spw="",scan="",antenna="",correlation="",timerange=timerange_, \
            intent="",array="",uvrange="",observation="",feed="", \
            datacolumn=datacolumn_,keepflags=True,width=1, \
            timebin="0s",combine="")
    print 'Splitting Finished..'+MOMMS
   


# INPUT PARAMETERS
MS_=OUTPUT_MS
DIR='selfsol'
prefix=DIR+'/fsun.selfcal_soln'
imagename_final='20130423T203430-203431_final'
spw_final='7:55~60'
timerange_final='20:34:30~20:34:31'
antennas_=''#'ea01,ea03,ea06,ea07,ea08,ea09,ea17,ea18,ea19,ea21,ea23,ea25,ea26,ea27'
refant_=''#'ea07'
# RA:02 05 i40.74 DEC:+12 44 01.6 20:34:00
#phase_centre='J2000 0.548505138896rad 0.222291421955rad'
solar_centre='J2000 02h05m40.74 +12d44m01.6' # For 203430-203431
phase_centre='J2000 02h04m51.211 12d43m12.47' # From Yingjie script
#phase_centre='J2000 02h05m41.97 +12d44m08.2' # For 204200
#phase_centre='J2000 02h05m42.44 +12d44m10.6' # For 20:45:00 solar centre
masks_=['mask0.mask']
#phase_centre='J2000 02h04m51.211 12d43m12.47'
initial_model_L_band='fsun.selfcal_soln.spw.7.cal.5.p.model'

# STANDARD PARAMETERS 
#spw_list=['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']
spw_list=['0','1','2','3','4','5','6','7']
chanlist=['3~60']*9
gainchan='28~32'
interactive=False
timerange_=''
modes=['p','p','a','p','p','p']#['p','p','p','p']#['ap','ap','ap','ap','ap','ap','ap']
snr_=[2,2,2,2,2,2,2]
caltable_all=[0]*len(spw_list)
MS=[0]*(len(modes)+1)
MS[0]=MS_
niter_scale_factor=np.ones(len(spw_list))#np.linspace(2,1,len(spw_list))
niter_scale_factor[1]=1.0#4.0
#niter_base=np.array([100,200,400,400,400,400,400,400,400])#[,600,800,1000,1200,1300,1000,1000])
niter_base=np.array([100,100,100,200,200,200,400,400,400])#[,600,800,1000,1200,1300,1000,1000])
#niter_base=np.ones(len(modes))*100
#uvrange_=['2000~7000lambda','500~7000lambda','','','','','','','','']
uvrange_=['','','','','','','','','','']
stokes_='I'
imsize_=[256,256]
cell_=['2.5arcsec','2.5arcsec']

# SCRIPT'S TAG 
first_clearcal=1
make_model=1
apply_model=1
use_initial_model=0
cal_gaincal=1 # Major Step to compute gaincal
plot_cal_amp=1 # plot amp soln
plot_cal_phase=1 # plot phase soln
make_image=1 # make image for each selfcal iteration
make_images_edges=0 # make images of edge channnels for spw
do_organise=0 # Put images and pngs and solns into separate folders
type1=1 # Use subsequent MS to apply computed soln
type2=0 # Use first MS to apply computed soln
apply_cal_type1=0 # Default initilisation of apply cal / apply soln to MS
apply_cal_type2=0
if(type1):
    apply_cal_type1=1
if(type2):
    apply_cal_type2=1
do_clearcal_momms=0 # Clear solns in Mom MS
apply_cal_momms=0 # Apply solns in Mom MS
plot_DR=1 # Plot Dynamic Range for each spw

#################### SCRIPTS STARTS ########################################
DR=[0]*len(spw_list)
max_=[0]*len(spw_list)
std=[0]*len(spw_list)
initial_model=[0]*len(spw_list)
initial_model[0]=initial_model_L_band
ss=0
for spw_ in spw_list[::-1]:
    niter_list=np.array(niter_base/niter_scale_factor[ss],dtype=int)
    print '##### STARTING SELFCAL FOR SPW:'+str(spw_)+' #####'
    print 'Niter list: ', niter_list
    std[ss]=[0]*len(modes)
    max_[ss]=[0]*len(modes)
    DR[ss]=[0]*len(modes)
    image_=[0]*len(modes)
    model_=[0]*(len(modes)+1)
    image_f=[0]*len(modes)
    model_[0]=[]
    if(first_clearcal):
        clearcal(vis=MS[0])
    ######## SELFCAL BEGINS #################################
    for step in range(len(modes)):
        # Make model 
        image_[step]=prefix+'.spw.'+spw_+'.cal.'+str(step)+'.'+str(modes[step])
        mod=image_[step]+'.model'
        model_[step+1]=mod#+model_[step]
        MS[step+1]=MS[0].split('.ms')[0]+'.slfcal.'+str(step)+'.ms'
        print model_
        niter_=niter_list[step]
        npercycle_=50
        if(make_model):
            print 'Making a model from CLEAN..'
            if(interactive==True):
                clean(vis=MS[step],imagename=image_[step], \
                        field="",spw=spw_+':'+chanlist[step],selectdata=True,timerange=timerange_, \
                uvrange=uvrange_[step],antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False, \
                interpolation="linear",niter=niter_,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
                interactive=True,mask=masks_, \
                imsize=imsize_,cell=cell_,phasecenter=phase_centre, \
                stokes=stokes_,weighting="natural", \
                modelimage='',restoringbeam=[''],usescratch=False, \
                npercycle=npercycle_) 
            if(interactive==False):
                clean(vis=MS[step],imagename=image_[step], \
                        field="",spw=spw_+':'+chanlist[step],selectdata=True,timerange=timerange_, \
                uvrange=uvrange_[step],antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False, \
                interpolation="linear",niter=niter_,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
                interactive=False,mask=masks_, \
                imsize=imsize_,cell=cell_,phasecenter=phase_centre, \
                stokes=stokes_,weighting="natural", \
                modelimage='',restoringbeam=[''],usescratch=False, \
                npercycle=npercycle_) 
            os.system('rm -rf '+image_[step]+'.residual')
            os.system('rm -rf '+image_[step]+'.mask')
            os.system('rm -rf '+image_[step]+'.flux')
            os.system('rm -rf '+image_[step]+'.psf')

        # FT
        if(apply_model):
            if(use_initial_model & (step==0)):
                ft(vis=MS[0],field="",spw=spw_,model=initial_model[ss],nterms=1,reffreq="",complist="",incremental=False,usescratch=True)
            else:
                print 'Applying model...'+str(model_[step+1])+'.model'
                #delmod(vis=MS)
                ft(vis=MS[step],field="",spw=spw_,model=model_[step+1],nterms=1,reffreq="",complist="",incremental=False,usescratch=True)
        # gaincal
        if(cal_gaincal):
            print 'Calculating solutions..'
            caltable_=image_[step]
            gaincal(vis=MS[step],caltable=caltable_,field="",spw=spw_+':'+gainchan,selectdata=True,timerange=timerange_,uvrange="", 
                    solint="inf",refant=refant_,minblperant=3,minsnr=snr_[step], \
                    gaintype="G",calmode=modes[step],append=True)
        # plotcal
        if(plot_cal_amp):
            plotcal(caltable=caltable_,xaxis="antenna",yaxis="amp",poln="RL",
                    field="",antenna=antennas_,spw="",timerange=timerange_,subplot=111,overplot=False,clearpanel="Current",
                    iteration="spw",plotrange=[-1,-1,0,2],showflags=False,plotsymbol="o",plotcolor="blue",markersize=5.0,
                    fontsize=10.0,showgui=False,figfile=caltable_+'_amp_ant.png')
        if(plot_cal_phase):
            plotcal(caltable=caltable_,xaxis="antenna",yaxis="phase",poln="RL",
                    field="",antenna=antennas_,spw="",timerange=timerange_,subplot=111,overplot=False,clearpanel="Current",
                    iteration="spw",plotrange=[-1,-1,-180,180],showflags=False,plotsymbol="o",plotcolor="blue",markersize=5.0,
                    fontsize=10.0,showgui=False,figfile=caltable_+'_phase_ant.png')
        plot_cal_amp_chan=0
        if(plot_cal_amp_chan):
            plotcal(caltable=caltable_,xaxis="chan",yaxis="amp",poln="L",
                    field="",antenna="",spw="",timerange="",subplot=431,overplot=False,clearpanel="Current",
                    iteration="antenna",plotrange=[-1,-1,0,1],showflags=False,plotsymbol="o",plotcolor="blue",markersize=5.0,
                    fontsize=10.0,showgui=False,figfile=caltable_+'_amp_chan.png')
        plot_cal_phase_chan=0
        if(plot_cal_phase_chan):
            plotcal(caltable=caltable_,xaxis="chan",yaxis="phase",poln="L",
                    field="",antenna="",spw="",timerange="",subplot=431,overplot=False,clearpanel="Current",
                    iteration="antenna",plotrange=[-1,-1,0,1],showflags=False,plotsymbol="o",plotcolor="blue",markersize=5.0,
                    fontsize=10.0,showgui=False,figfile=caltable_+'_phase_chan.png')
        # Delmod
        del_mod=0
        if(del_mod):
            print 'Deleting old model...'
            delmod(vis=MS)

        # Apply solutions
        if(apply_cal_type2):
            clearcal(vis=MS[0])
            print 'Applying solutions....',image_[0:step+1]
            applycal(vis=MS[0],field="",spw=spw_,intent="",selectdata=True,timerange="",uvrange="",antenna=antennas_, \
                    gaintable=image_[0:step+1], \
                    gainfield=[],interp='nearest', \
                    calwt=False,applymode='calonly',flagbackup=False)

        # Apply solutions
        if(apply_cal_type1):
            delmod(vis=MS[step])
            clearcal(vis=MS[step])
            print 'Applying solutions....',[caltable_]
            applycal(vis=MS[step],field="",spw=spw_,intent="",selectdata=True,timerange="",uvrange="",antenna=antennas_, \
                    gaintable=[caltable_], \
                    gainfield=[],interp='nearest', \
                    calwt=False,applymode='calonly',flagbackup=False)
        image_f[step]=image_[step]+'.final'
        if(make_image):
            if(type1):
                MSI=MS[step]
            if(type2):
                MSI=MS[0]
            print 'Making image with VIS='+MSI+' IMAGENAME= '+image_f[step]
            clean(vis=MSI,imagename=image_f[step], \
                    outlierfile="",field="",spw=spw_,selectdata=True,timerange="", \
                    uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False, \
                    gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0, \
                    aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="", \
                    interpolation="linear",niter=800,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
                    ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1, \
                    smallscalebias=0.6,interactive=False,mask=['mask0.mask'],nchan=-1,start=0,width=1,outframe="",veltype="radio", \
                    imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter=phase_centre,restfreq="", \
                    stokes='I',weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'], \
                    modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", \
                    npixels=0,npercycle=10,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
                    flatnoise=True,allowchunk=False)
            os.system('rm -rf '+image_f[step]+'.model')
            #os.system('rm -rf '+image_f[step]+'.residual')
            os.system('rm -rf '+image_f[step]+'.mask')
            os.system('rm -rf '+image_f[step]+'.flux')
            os.system('rm -rf '+image_f[step]+'.psf')
            img_=imstat(imagename=image_f[step]+".image/",axes=-1,region="",box="",chans="",stokes="",listit=True, \
                    verbose=True,mask="",stretch=False,logfile="",append=True,algorithm="classic",fence=-1, \
                    center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto",niter=3)
            img_std=imstat(imagename=image_f[step]+".residual/",axes=-1,region="",box="10,10,200,200",chans="",stokes="",listit=True, \
                    verbose=True,mask="",stretch=False,logfile="",append=True,algorithm="classic",fence=-1, \
                    center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto",niter=3)
            std[ss][step]=img_std['rms']
            max_[ss][step]=img_['max']
            DR[ss][step]=max_[ss][step]/std[ss][step]
            print DR[ss][step],max_[ss][step],std[ss][step]
        # split the corrected data  
        split_corrdata=1
        if(split_corrdata):
            if(type1):
                print 'Splitting the corrected data..'
                os.system('rm -rf '+MS[step+1])
                split(vis=MS[step],outputvis=MS[step+1],datacolumn='corrected')
            if(type2):
                print 'Splitting the corrected data..'
                os.system('rm -rf '+MS[step+1])
                split(vis=MS[0],outputvis=MS[step+1],datacolumn='corrected')
        caltable_all[ss]=image_[0:step+1]
    initial_model[ss+1]=mod
    print 'Initial model for next iteration: ',initial_model[ss+1]
    ss=ss+1
        ############## SELFCAL ENDS ####################

# Plotting Dynamic Range
if(plot_DR):
    i=0
    for spw_ in  spw_list:
        mpl.pyplot.plot(np.arange(len(DR[i]))+1,DR[i],'o-')
        mpl.pyplot.xlabel('Iterations')
        mpl.pyplot.ylabel('Dynamic Range')
        mpl.pyplot.savefig(prefix+'_DR_spw'+spw_+'.png')
        mpl.pyplot.close()
        i=i+1

# Apply solutions to the final MS
if(do_clearcal_momms):
    clearcal(MOMMS)
if(apply_cal_momms):
    caltable_apply=list(itertools.chain(*caltable_all))
    print 'Applying solutions....',caltable_apply
    #clearcal(vis=MS)
    applycal(vis=MOMMS,field="",spw=spw_,intent="",selectdata=True,timerange="",uvrange="",antenna="", \
            scan="",observation="",msselect="",docallib=False,callib="", \
            #gaintable=image_, \
            gaintable=caltable_apply, \
            gainfield=[],interp=[], \
            spwmap=[],calwt=[False],parang=False,applymode='calonly',flagbackup=True)

clean(vis=MOMMS,imagename=imagename_final, \
        outlierfile="",field="",spw=spw_final,selectdata=True,timerange=timerange_final, \
        uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False, \
        gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0, \
        aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="", \
        interpolation="linear",niter=800,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
        ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1, \
        smallscalebias=0.6,interactive=False,mask=['mask0.mask'],nchan=-1,start=0,width=1,outframe="",veltype="radio", \
        imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter=solar_centre,restfreq="", \
        stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'], \
        modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", \
        npixels=0,npercycle=10,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
        flatnoise=True,allowchunk=False)
#### MAKE EDGE IMAGES ####

for i in range(len(spw_list)):
    image_lower=image_[step]+'.spw.'+spw_list[i]+'.lower' # Lower in frequency
    if(make_images_edges):
        print 'Making Edge images....'
        clean(vis=MOMMS,imagename=image_lower, \
                outlierfile="",field="",spw=spw_+':4~9',selectdata=True,timerange="", \
                uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False, \
                gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0, \
                aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="", \
                interpolation="linear",niter=800,gain=0.1,threshold="0.0mJy",psfmode="hogbom",imagermode="csclean", \
                ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1, \
                smallscalebias=0.6,interactive=False,mask=['mask0.mask'],nchan=-1,start=0,width=1,outframe="",veltype="radio", \
                imsize=[1024, 1024],cell=['2.5arcsec', '2.5arcsec'],phasecenter=phase_centre,restfreq="", \
                stokes=stokes_,weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'], \
                modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", \
                npixels=0,npercycle=10,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
                flatnoise=True,allowchunk=False)
        os.system('rm -rf '+image_lower+'.model')
        os.system('rm -rf '+image_lower+'.residual')
        os.system('rm -rf '+image_lower+'.mask')
        os.system('rm -rf '+image_lower+'.flux')
        os.system('rm -rf '+image_lower+'.psf')
    image_upper=image_[step]+'.spw.'+spw_list[i]+'.upper' # Upper in frequency
    if(make_images_edges):
        print 'Making image....'
        clean(vis=MOMMS,imagename=image_upper, \
                outlierfile="",field="",spw=spw_+':55~60',selectdata=True,timerange="", \
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
    os.system('mkdir '+DIR+'/selfcal_ms')
    os.system('rm -rf '+MS[0].split('.ms')[0]+'.slfcal.*.ms')




