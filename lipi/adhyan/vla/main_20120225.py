from surya.utils import main as ut
from surya.plot import visual as vs
from surya.aia import main as aia
from surya.vla import main as vla
from surya.hmi import main as hmi
from surya.rhessi import main as rh
from surya.plot import movie as mov
from sunpy.map import Map
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
import glob
import matplotlib.pyplot as plt
import plot_20120225 as pl
import plot_20120225_old as plo
import os
import pickle

plt.style.use('/home/i4ds1807205/scripts/general/plt_style.py')
wavel=94
aiapath='/home/i4ds1807205/vla_data/20120225/aia/'
aiafile=aiapath+'aia_'+str(wavel)+'_sub_map.sav'
aiafile_all=sorted(glob.glob(aiapath+'aia*.sav'))
aiafile_all.sort(key=lambda ss: len(ss))
wavelength=wavel
wav=[94,131,171,193,211,304,335,1600,1700]
xl,xr,yl,yr=450,520,320,390
#xl,xr,yl,yr=420,540,280,420 # For Quiet Sun plotting
res=0.6
aiaxarray=np.linspace(xl,xr,100)
aiayarray=np.linspace(yl,yr,100)
#cmap=aia.get_submap(aiafile,res,wavelength)
cmap,aiatsec,aiatstring=aia.get_submap_crop(aiafile,res,wavelength,xl,xr,yl,yr)
cmap_all,aiatsec_all,aiatstring_all=aia.get_submap_all(aiafile_all,wav)
startt=' 20:46:00'
endt=' 20:50:00'
idxs,idxe=aia.get_time_idx(startt,endt,aiatsec)
ccmap=cmap[idxs:idxe]
aiactime=aiatstring[idxs:idxe]

aiamax=[0]*len(ccmap)
aiamax_xid=[0]*len(ccmap)
aiamax_yid=[0]*len(ccmap)
for j in range(len(ccmap)):
    aiamax[j]=np.max(ccmap[j])
    aiamax_xid[j]=aiaxarray[np.where(ccmap[j]==np.max(ccmap[j]))[0][0]]
    aiamax_yid[j]=aiayarray[np.where(ccmap[j]==np.max(ccmap[j]))[1][0]]

aiadata304=hmi.get_submap_('/media/rohit/VLA/20120225_aia/AIA20120225_205132_0304.fits',xl,xr,yl,yr)[0]
aiadata171=hmi.get_submap_('/media/rohit/VLA/20120225_aia/AIA20120225_204600_0171.fits',xl,xr,yl,yr)[0]
aiadata94=hmi.get_submap_('/media/rohit/VLA/20120225_aia/AIA20120225_204902_0094.fits',xl,xr,yl,yr)[0]

###################### RHESSI

rh_file='/home/i4ds1807205/vla_data/hsi_imagecube_15tx2e_20120225_204500_1min_aia_time.fits'
rhessi_time=['20:45:00','20:46:00','20:47:00','20:47:28','20:48:00','20:49:00','20:50:00','20:51:00','20:52:00','20:53:00','20:54:00','20:55:00','20:56:00','20:57:00','20:58:00']
rhres=1.0 # in arcsec
xrh,yrh,rhsize=480,363,129
rhsubmap_high,rhsubmap_low,rh_tsec=rh.get_submap(rh_file,rhessi_time,xrh,yrh,rhsize,rhres,xl,xr,yl,yr)
rhh,rhd=ut.read_fits(rh_file)
rhxl,rhxr,rhyl,rhyr=rhh['XCEN']-rhh['NAXIS1']*0.5,rhh['XCEN']+rhh['NAXIS1']*0.5,rhh['YCEN']-rhh['NAXIS2']*0.5,rhh['YCEN']+rhh['NAXIS2']*0.5
rhmap=Map(rhd[10][0],rhh)
rhw=WCS(rhh)
#rhsubmap=rh.get_submap_(rhmap,xl,xr,yl,yr,rhw)

#rhessi st: 20:08:46
###################### HMI
hmifile='/home/i4ds1807205/vla_data/20120225/hmi/hmi.m_45s.2012.02.25_20_51_45_TAI.magnetogram.fits'
#hmifile='/home/i4ds1807205/Downloads/try.fits'
#hmidata=hmi.get_submap(hmifile,xl,xr,yl,yr)
hmidata,xlh,ylh,xrh,yrh=hmi.get_submap_(hmifile,xl,xr,yl,yr)
#pl.plot_hmi_rhessi(hmifile,rhmap,rhxl,rhxr,rhyl,rhyr)

###################### VLA


tstart_vla='20:46:00'
tend_vla='20:50:00'

#t='20:47:09~20:47:10'
spw=['5','6','7']
pix=2.5
dec=-9.13056
x_offset,y_offset=vla.get_offset(dec)
vlapath='/media/rohit/VLA/20120225_sub/fits/'
#vlapath='/media/rohit/VLA/20120225_cube/'
filelist=sorted(glob.glob(vlapath+'hcs*.FITS'))
lev=[0.5,0.6,0.7,0.8,0.9]
lev=[0.7,0.8,0.9]
vlasubmap=[0]*len(filelist)
xc_mean=[0]*len(filelist)
yc_mean=[0]*len(filelist)
xc_std=[0]*len(filelist)
yc_std=[0]*len(filelist)
vlat=[0]*len(filelist)
vlatsec=[0]*len(filelist)
if(os.path.isfile(vlapath+'/centroid_20120225.p')==False):
    ss=0
    for f in filelist:
        s=f.split('.')[4]
        vlat[ss]=f.split('.')[-2].split('~')[0]
        vlatsec[ss]=ut.hms2sec_c(' '+vlat[ss])
        print f,vlatsec[ss],' '+vlat[ss]
        #vlafile_fl=vlapath+'hcs_2050.1s.cal.spw.'+s+'.time.'+t+'.FITS'
        #vlafile_qs=vlapath+'hcs_2050.1s.cal.spw.7.time.20:46:02~20:46:03.FITS'
        #vlafile_qs='/media/rohit/VLA/selfcal/selfcal_pf_qs/hcs_2050.1s.cal.spw.7.time.20:46:00~20:46:01.FITS'
        #vlafile_qs0='/media/rohit/VLA/yingie_data/20120225_images/time_studies/spw_7/flare_time/hcs_2050.1s.cal.spw.7:0~3.time.20:46:01~20:46:02.FITS'
        #vlafile_all=sorted(glob.glob(vlapath+'hcs*.FITS'))
        #vlafile_all.sort(key=lambda ss: len(ss))
        #vlamap,vlah=vla.get_map(vlafile_qs0)
        #vlamap,vlah=vla.get_map_cube(vlafile_fl)
        #vlasubmap,vlah=vla.get_submap(vlafile_qs0,xl,xr,yl,yr,pix,dec)
        vlasubmap[ss],vlah=vla.get_submap_cube(f,xl,xr,yl,yr,pix,dec)
        #xc,yc,w,l,ang=vla.fit_centroid_max(vlasubmap[0][0][::-1],lev[0])
        xc,yc,w,l,ang=vla.fit_centroid_map(vlasubmap[ss],lev,xl,yl,pix,dec)
        xc_mean[ss]=xc.mean(axis=0)
        xc_std[ss]=xc.std(axis=0)
        yc_mean[ss]=yc.mean(axis=0)
        yc_std[ss]=yc.std(axis=0)
        ss=ss+1
    vlasubmap=np.array(vlasubmap)
    tt=int(vlasubmap.shape[0]/3)
    vlasubmap=vlasubmap.reshape(3,tt,32,vlasubmap.shape[2],vlasubmap.shape[3]).swapaxes(0,1).reshape(tt,96,vlasubmap.shape[2],vlasubmap.shape[3])
    xc_mean=np.array(xc_mean)
    xc_mean=xc_mean.reshape(3,tt,32).swapaxes(0,1).reshape(tt,96)
    yc_mean=np.array(yc_mean)
    yc_mean=yc_mean.reshape(3,tt,32).swapaxes(0,1).reshape(tt,96)
    vlatsec=np.array(vlatsec).reshape(3,tt)
    freq=1700+4*np.arange(96) # Frequencies are always from low to high 0-->32
    pickle.dump([vlasubmap,xc_mean,yc_mean,freq,vlatsec,vlat],open(vlapath+'/centroid_20120225.p','wb'))


plot_check_cent=0
if(plot_check_cent):
    vlasubvs='/media/rohit/VLA/20120225_sub/fits/'
    fsub=open(vlasubvs+'/centroid_20120225.p','rb')
    subvs_cent=pickle.load(fsub)
    fsub.close()
    vlacube='/media/rohit/VLA/20120225_cube/'
    fcube=open(vlacube+'/centroid_20120225.p','rb')
    cube_cent=pickle.load(fcube)
    fcube.close()
    plt.plot(subvs_cent[1][68,:],subvs_cent[2][68,:],'o-')
    plt.plot(cube_cent[1][66,:],cube_cent[2][66,:],'o-')
    plt.show()

vlasubmap,xc_mean,yc_mean,freq,vlatsec,vlat=pickle.load(open(vlapath+'/centroid_20120225.p','r'))
noise,peak=pickle.load(open('/media/rohit/VLA/20120225_cube/20120225_noise.p','r'))
tt=int(vlasubmap.shape[0])
bmaj=np.concatenate((np.ones(96)*19.68,np.ones(96)*18.40,np.ones(96)*17.08),axis=0)
bmin=np.concatenate((np.ones(96)*12.30,np.ones(96)*11.46,np.ones(96)*10.69),axis=0)

contour_ref=0.9
region1=np.where(vlasubmap[68][-1]>contour_ref*np.max(vlasubmap[68][-1]))
Tb=[0]*vlasubmap.shape[0]
Tb1=[0]*vlasubmap.shape[0]
for i in range(vlasubmap.shape[0]):
    Tb[i]=[0]*vlasubmap.shape[1]
    Tb1[i]=[0]*vlasubmap.shape[1]
    Tb_noise=[0]*vlasubmap.shape[1]
    for j in range(vlasubmap.shape[1]):
        Tb[i][j]=ut.flux2Tb(vlasubmap[i,j]*1.e-4,bmaj[j],bmin[j],freq[j]/1000.)*100+1.e6 # Attenuation = 100
        Tb1[i][j]=np.mean(Tb[i][j][region1])
        Tb_noise[j]=ut.flux2Tb(noise[j]*1.e-4,bmaj[j],bmin[j],freq[j]/1000.)*100+1.e6 # Attenuation = 100
Tb=np.array(Tb)
Tb1=np.array(Tb1)
Tb_noise=np.array(Tb_noise)/3 # 1 sigma variation
Tb[85]=np.nan
Tb1[85]=np.nan
Tb[197]=np.nan
Tb1[197]=np.nan

cent_refx,cent_refy=484,350
rc_mean=np.sqrt((xc_mean-cent_refx*np.ones((tt,96)))**2 + (yc_mean-cent_refy*np.ones((tt,96)))**2)
rc_mean[85]=np.nan
rc_mean[197]=np.nan
xc_mean[85]=np.nan
xc_mean[197]=np.nan
yc_mean[85]=np.nan
yc_mean[197]=np.nan


pickle.dump([Tb1,Tb_noise],open('20120225_peaks.p','wb'))
np.savetxt('20120225_peaks_Tb.txt',Tb1)
np.savetxt('20120225_freq.txt',freq)
np.savetxt('20120225_error_Tb.txt',Tb_noise)

## Plotting
# rhtime (20:08:46)
#pl.euv_vla(cmap,vlasubmap[-1],xl,xr,yl,yr)
#pl.goes_vla_line_plot()
#pl.rhessi_vla_line_plot()
#pl.plot_ds()
#plo.Tb(Tb1.mean(axis=1),rc_mean.swapaxes(0,1),Tb1.swapaxes(0,1)/1.e6,freq)
#plo.plot_centroids(xc_mean[:,-3],yc_mean[:,-3])
lev_1=np.linspace(0.18,0.180001,2) #for 304
lev_2=np.linspace(0.2,0.2001,3) #for 171
#pl.composite_ma(hmidata,cmap_all,lev_all,430,540,300,410)
#pl.composite_map(hmidata,cmap_all,lev_all,xl,xr,yl,yr)
#pl.centroid_map(hmidata.data[::-1,::-1],cmap[30], xc_mean,yc_mean,lev_1,xl,xr,yl,yr)
#pl.hmi_map_inv_lines(hmidata.data[::-1,::-1],cmap[30], xc_mean,yc_mean,lev_all,xl,xr,yl,yr)
# Change to 131 \AA for the below plot
#pl.hmi_euv_map_paper(hmidata.data[::-1,::-1],aiadata304.data*1.0,aiadata171.data*1.0, xc_mean,yc_mean,lev_1,lev_2,xl,xr,yl,yr)
pl.euv_map_paper(aiadata94.data*1.0,aiadata171.data*1.0, xc_mean,yc_mean,lev_1,lev_2,xl,xr,yl,yr)
#plo.plot_spec(Tb1,Tb_noise,freq)
#plo.plot_spec_movie(Tb1,Tb_noise,freq)
sys.exit()
pngfiles=[0]*len(vlatsec[0])
#for i in range(1):
for i in range(len(vlatsec[0])):
    ii="%03d" %i 
    #pngfiles[i]=vlapath+'../png/'+'cube_rh'+ii+'.png'
    #pngfiles[i]='/media/rohit/VLA/20120225_sub/png/'+'subvis_rh'+ii+'_high.png'
    pngfiles[i]='/media/rohit/VLA/20120225_sub/aia/'+'aia'+str(wavel)+'_'+ii+'.png'
    aia_j=ut.find_predecessor(aiatsec,vlatsec[0][i])[0]
    rh_j=ut.find_predecessor(rh_tsec,vlatsec[0][i])[0]
    title='AIA: '+aiatstring[aia_j]+'   VLA: '+vlat[i]+'     RHESSI: '+rhessi_time[rh_j]
    rh_norm_low=rhsubmap_low[rh_j]/np.max(rhsubmap_low[rh_j])
    rh_norm_high=rhsubmap_high[rh_j]/np.max(rhsubmap_high[rh_j])
    rhlev=[0.65,0.75,0.85,0.95]
    #plo.euv_vla_rhessi_qs_centroids(cmap[aia_j],vlasubmap[i],rh_norm_low,rh_norm_high,rhlev,freq[::-1],xc_mean[i],yc_mean[i],xl,xr,yl,yr,21,13,-60,title)
    plo.euv(cmap[aia_j],xl,xr,yl,yr,'AIA '+str(wavel)+' $\AA$ : '+aiatstring[aia_j])
    #pl.euv_vla_qs_centroids(cmap[j],vlasubmap[i],freq,xc_mean[i],yc_mean[i],xl,xr,yl,yr,21,13,-60,title)
    plt.savefig(pngfiles[i])
    plt.close()

pngfiles.remove(pngfiles[85])
pngfiles.remove(pngfiles[196])
#mov.write_imagero(pngfiles,'/media/rohit/VLA/20120225_sub/cube_rh_20120225.gif',20)

#vs.xyz_line_plot(vlasubmap[70],freq,xc_mean[70],yc_mean[70],xl,xr,yl,yr,'')

