from surya.utils import main as ut
from surya.aia import main as aia
from surya.vla import main as vla
from surya.rhessi import main as rh
from surya.plot import movie as mov
from surya.plot import visual as vs
import numpy as np
import glob
import matplotlib.pyplot as plt
import plot_20120225 as pl
import os
import pickle

aiapath='/home/i4ds1807205/vla_data/20120225/aia/'
aiafile=aiapath+'aia_171_sub_map.sav'
aiafile_all=sorted(glob.glob(aiapath+'aia*.sav'))
aiafile_all.sort(key=lambda ss: len(ss))
wavelength=171
wav=[94,131,171,193,211,304,335,1600,1700]
xl,xr,yl,yr=450,510,320,380
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

###################### RHESSI

rh_file='/home/i4ds1807205/vla_data/hsi_imagecube_15tx2e_20120225_204500_1min_aia_time.fits'
rhessi_time=['20:45:00','20:46:00','20:47:00','20:47:28','20:48:00','20:49:00','20:50:00','20:51:00','20:52:00','20:53:00','20:54:00','20:55:00','20:56:00','20:57:00','20:58:00']
rhres=1.0 # in arcsec
xrh,yrh,rhsize=480,363,129
rhsubmap_high,rhsubmap_low,rh_tsec=rh.get_submap(rh_file,rhessi_time,xrh,yrh,rhsize,rhres,xl,xr,yl,yr)


###################### VLA

t='20:47:09~20:47:10'
spw=['5','6','7']
pix=2.5
dec=-9.13056
x_offset,y_offset=vla.get_offset(dec)
vlapath='/media/rohit/VLA/20120225_sub/fits/'
filelist=sorted(glob.glob(vlapath+'hcs*.FITS'))
lev=[0.5,0.6,0.7,0.8,0.9]
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
    vlasubmap=vlasubmap.reshape(3,223,32,vlasubmap.shape[2],vlasubmap.shape[3]).swapaxes(0,1).reshape(223,96,vlasubmap.shape[2],vlasubmap.shape[3])
    xc_mean=np.array(xc_mean)
    xc_mean=xc_mean.reshape(3,223,32).swapaxes(0,1).reshape(223,96)
    yc_mean=np.array(yc_mean)
    yc_mean=yc_mean.reshape(3,223,32).swapaxes(0,1).reshape(223,96)
    vlatsec=np.array(vlatsec).reshape(3,223)
    freq=1700+4*np.arange(96) # Frequencies are always from low to high 0-->32
    pickle.dump([vlasubmap,xc_mean,yc_mean,freq,vlatsec,vlat],open(vlapath+'/centroid_20120225.p','wb'))

vlasubmap,xc_mean,yc_mean,freq,vlatsec,vlat=pickle.load(open(vlapath+'/centroid_20120225.p','r'))
bmaj=np.concatenate((np.ones(96)*19.68,np.ones(96)*18.40,np.ones(96)*17.08),axis=0)
bmin=np.concatenate((np.ones(96)*12.30,np.ones(96)*11.46,np.ones(96)*10.69),axis=0)

contour_ref=0.8
region1=np.where(vlasubmap[70][-1]>contour_ref*np.max(vlasubmap[70][-1]))
Tb=[0]*vlasubmap.shape[0]
Tb1=[0]*vlasubmap.shape[0]
for i in range(vlasubmap.shape[0]):
    Tb[i]=[0]*vlasubmap.shape[1]
    Tb1[i]=[0]*vlasubmap.shape[1]
    for j in range(vlasubmap.shape[1]):
        Tb[i][j]=ut.flux2Tb(vlasubmap[i,j]*1.e-4,bmaj[j],bmin[j],freq[j]/1000.)*100 # Attenuation = 100
        Tb1[i][j]=np.mean(Tb[i][j][region1])
Tb=np.array(Tb)
Tb1=np.array(Tb1)
cent_refx,cent_refy=484,350
rc_mean=np.sqrt((xc_mean-cent_refx*np.ones((223,96)))**2 + (yc_mean-cent_refy*np.ones((223,96)))**2)

sys.exit()

## Plotting
#pl.euv_vla(cmap,vlasubmap[-1],xl,xr,yl,yr)

pngfiles=[0]*len(vlatsec[0])
for i in range(len(vlatsec[0])):
    i=68
    ii="%03d" %i 
    pngfiles[i]=vlapath+'../png/'+'subvis_rh'+ii+'.png'
    aia_j=ut.find_predecessor(aiatsec,vlatsec[0][i])[0]
    rh_j=ut.find_predecessor(rh_tsec,vlatsec[0][i])[0]
    title='AIA: '+aiatstring[aia_j]+' VLA:'+vlat[i]+' RH:'+rhessi_time[rh_j]
    rh_norm_low=rhsubmap_low[rh_j]/np.max(rhsubmap_low[rh_j])
    rh_norm_high=rhsubmap_high[rh_j]/np.max(rhsubmap_high[rh_j])
    rhlev=[0.65,0.75,0.85,0.95]
    pl.euv_vla_rhessi_qs_centroids(cmap[aia_j],vlasubmap[i],rh_norm_low,rh_norm_high,rhlev,freq[::-1],xc_mean[i],yc_mean[i],xl,xr,yl,yr,21,13,-60,title)
    #pl.euv_vla_qs_centroids(cmap[j],vlasubmap[i],freq,xc_mean[i],yc_mean[i],xl,xr,yl,yr,21,13,-60,title)
    #plt.savefig(pngfiles[i],dpi=40)
    plt.show()

pngfiles.remove(pngfiles[85])
pngfiles.remove(pngfiles[196])
#mov.write_imagero(pngfiles,vlapath+'../subvis_rh_20120225.gif',20)

#vs.xyz_line_plot(vlasubmap[70],freq,xc_mean[70],yc_mean[70],xl,xr,yl,yr,'')

