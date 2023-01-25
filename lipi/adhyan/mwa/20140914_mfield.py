import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import sys
import glob
from astropy.io import fits
from sunpy.map import Map
from astropy.coordinates import SkyCoord
import astropy.units as u
from surya.utils import main as ut
import pickle
from sunpy import sun
from dateutil import parser
from surya.utils import Bextrap
from reproject import reproject_interp
from sunpy.coordinates import frames
from astropy.wcs import WCS
from astropy.wcs import utils
import csv
from sunpy.coordinates import get_body_heliographic_stonyhurst
from surya.utils import model
from sunpy.map import make_fitswcs_header, Map
from scipy.io import readsav
import os
import sunpy
from tvtk.api import tvtk, write_data
import matplotlib.cm as cm
from surya.utils import model
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import pycwt as wavelet

rename=0
if(rename):
    lst=glob.glob('hmi.m_45s.*.fits')
    for l in lst:
        out=l.split('.')[0]+'_'+l.split('.')[2]+l.split('.')[3]+l.split('.')[4].split('_')[0]+'_'+l.split('.')[4].split('_')[1]+l.split('.')[4].split('_')[2]+l.split('.')[4].split('_')[3]+'_magnetogram.fits'
        os.system('mv '+l+' '+out)
        

hmifile='/sdata/20140914_hmi/jsoc/2014-09-14/hmi.M_720s.20140914_014800_TAI.magnetogram.fits'
#hmifile='/sdata/20140914_hmi/jsoc/2014-09-14/hmi.M_720s.20140914_022400_TAI.magnetogram.fits'
hmimap=Map(hmifile)
hmid=hmimap.data#[::-1,::-1]
hmid[np.where(hmid<-5000)]=1
hmimap=Map(hmid,hmimap.meta)
###########################
aiamap=Map('/media/rohit/MWA/20140914/EUV/fits/aia.lev1.193A_2014-09-14T02_24_33.45Z.image_lev1.fits')
aiamap=Map('/media/rohit/MWA/20140914/EUV/fits/aia.lev1.171A_2014-09-14T02_24_13.21Z.image_lev1.fits')
#aiamapb,aiamapd=np.load('aia.lev1.171A_2014-09-14T02_16_13.26Z.image_lev1.fits_bdiff.sav.sunpy.npy',allow_pickle=True)

blh = SkyCoord(625*u.arcsec, -425*u.arcsec, frame=hmimap.coordinate_frame)
trh = SkyCoord(775*u.arcsec, -275*u.arcsec, frame=hmimap.coordinate_frame)
bla = SkyCoord(625*u.arcsec, -350*u.arcsec, frame=aiamap.coordinate_frame)
tra = SkyCoord(775*u.arcsec, -200*u.arcsec, frame=aiamap.coordinate_frame)
hmisubmap=hmimap.submap(blh, top_right=trh)
aiasubmap=aiamap.submap(bla, top_right=tra)
###########################
out_hmi = hmimap.reproject_to(aiamap.wcs)

k=0 # Layer
expolB_=readsav('/sdata/20140914_hmi/2014-09-14/hmi.M_720s.20140914_014621.W98S17CR.CEA.NAS.sav');expolB=expolB_['box']
#expolB_=readsav('/data/Dropbox/20120225_VLA_work/gregory_gyrosynchroton/hmi.M_720s.20120225_203413.W133N16CR.CEA.NAS.sav');expolB=expolB_['box']
bx=expolB['bx'][0];by=expolB['by'][0];bz=expolB['bz'][0]
babs=np.sqrt(bx**2 + by**2 + bz**2)
ff=1
if(ff):
    index=expolB['index'][0];crval1=index['crval1'];crval2=index['crval2'];crpix1=index['crpix1'];crpix2=index['crpix2'];ctype1=index['ctype1']
    ctype2=index['ctype2'];cdelt1=index['cdelt1'];cdelt2=index['cdelt2'];cunit1=index['cunit1'];cunit2=index['cunit2']
    hdu = fits.PrimaryHDU(babs[k])
    list_all=list(index.dtype.names);list_all.remove('COMMENT');list_all.remove('HISTORY');list_all.remove('BITPIX');list_all.remove('NAXIS');list_all.remove('DATE_D$OBS')
    index['WCSNAME'],index['CTYPE1'],index['CUNIT1'],index['CTYPE2'],index['CUNIT2']=['Carrington-Heliographic'],['CRLN-CEA'],['deg'],['CRLT-CEA'],['deg']
    index['DATE_OBS']=['2014-09-14T01:55:30']
ii=0
for idx in list_all:
    #print idx
    #hdu.header.append((idx,index[list_all[ii]][0],[]))
    hdu.header.update({str(idx):index[list_all[ii]][0]})
    ii=ii+1

delta=1.400
hdu.data=babs[0]
hhdu=hdu.header
hdul = fits.HDUList([hdu])
mymap=Map(babs[k],hhdu)
hp_coord=mymap.reference_coordinate.transform_to(frames.Helioprojective(observer="earth"))
hp_hcc=mymap.reference_coordinate.transform_to(frames.Heliocentric(observer="earth"))
out_shape = (334, 334)
out_header = sunpy.map.make_fitswcs_header(mymap.data,hp_coord)
out_wcs = WCS(out_header)
#earth = get_body_heliographic_stonyhurst('earth', mymap.date)
#out_wcs.heliographic_observer = earth
output, footprint = reproject_interp(mymap, out_wcs, out_shape)
outmap = sunpy.map.Map((output, out_header))
outmap.plot_settings = mymap.plot_settings

bsmap,bcarrmap=Bextrap.get_gxs_sav2hpp('/sdata/20140914_hmi/2014-09-14/hmi.M_720s.20140914_014621.W98S17CR.CEA.NAS.sav','2014-09-14T01:46:21')
hmiaia_data, footprint = reproject_interp(mymap, hmimap.wcs, shape_out=(400,400));hmiaia_map=Map(hmiaia_data,hmimap.meta)
x0,y0,z0,bx0,by0,bz0=Bextrap.get_fieldlines('/sdata/20140914_hmi/2014-09-14/hmi.M_720s.20140914_014621.W98S17CR.CEA.NAS.sav_0000.vtk.x-ray.csv')
xs0,ys0,zs0,bsx0,bsy0,bsz0,bs_hp0,bs_hp_pix0,bs_carr0,bs_carr_pix0=Bextrap.transform_fieldlines(x0,y0,z0,bx0,by0,bz0,'2014/09/14T01:55:00',bsmap[0].wcs,bcarrmap[0])
x1,y1,z1,bx1,by1,bz1=Bextrap.get_fieldlines('/sdata/20140914_hmi/2014-09-14/hmi.M_720s.20140914_014621.W98S17CR.CEA.NAS.sav_0000.vtk.csv')
xs1,ys1,zs1,bsx1,bsy1,bsz1,bs_hp1,bs_hp_pix1,bs_carr1,bs_carr_pix1=Bextrap.transform_fieldlines(x1,y1,z1,bx1,by1,bz1,'2014/09/14T01:55:00',bsmap[0].wcs,bcarrmap[0])
x2,y2,z2,bx2,by2,bz2=Bextrap.get_fieldlines('/sdata/20140914_hmi/2014-09-14/hmi.M_720s.20140914_022221.W98S17CR.CEA.NAS.sav_0001.vtk.csv')
xs2,ys2,zs2,bsx2,bsy2,bsz2,bs_hp2,bs_hp_pix2,bs_carr2,bs_carr_pix2=Bextrap.transform_fieldlines(x2,y2,z2,bx2,by2,bz2,'2014/09/14T02:22:21',bsmap[0].wcs,bcarrmap[0])
x2,y2,z2,bx2,by2,bz2=Bextrap.get_fieldlines('/sdata/20140914_hmi/2014-09-14/hmi.M_720s.20140914_022221.W98S17CR.CEA.NAS.sav_0001.vtk.x-ray.csv')
xs2,ys2,zs2,bsx2,bsy2,bsz2,bs_hp2,bs_hp_pix2,bs_carr2,bs_carr_pix2=Bextrap.transform_fieldlines(x2,y2,z2,bx2,by2,bz2,'2014/09/14T02:22:21',bsmap[0].wcs,bcarrmap[0])
#corx=45;cory=50
corx=0;cory=0
bshp0x=bs_hp0.Tx.value-corx;bshp0y=bs_hp0.Ty.value-cory
bshp1x=bs_hp1.Tx.value-corx;bshp1y=bs_hp1.Ty.value-cory
bshp2x=bs_hp2.Tx.value-corx;bshp2y=bs_hp2.Ty.value-cory
bsabs1=np.sqrt(bsx1**2 +bsy1**2+bsz1**2)
rs1=np.sqrt(xs1**2 +ys1**2)


beam=np.array(([6.79143269857,7.19919993083],
[6.18038431804,6.55234578451],
[4.25709788005,7.0822555542],
[5.03817596436,5.4289469401],
[3.94947509766,4.26192092895],
[4.18227513631,4.35891418457],
[3.68575286865,3.90556131999],
[3.34071350098,3.54748077393],
[3.1191889445,3.2729019165]))


####### Radio Contours #############
freq=np.array([108, 120, 133, 145, 160, 179, 197, 217, 240])
flist=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
aa=pickle.load(open('/media/rohit/MWA/20140914/Tb_20140914_187-188_sub_test.p','rb'),encoding='bytes');Tb240=np.array(aa[0]);maxX=[0]*len(Tb240);maxY=[0]*len(Tb240);xc90_240=[0]*len(Tb240);yc90_240=[0]*len(Tb240)
aa=pickle.load(open('/media/rohit/MWA/20140914/Tb_20140914_125-126_sub_test.p','rb'),encoding='bytes');Tb160=np.array(aa[0]);xc90_160=[0]*len(Tb160);yc90_160=[0]*len(Tb160)
aa=pickle.load(open('/media/rohit/MWA/20140914/Tb_20140914_084-085_sub_test.p','rb'),encoding='bytes');Tb108=np.array(aa[0]);xc90_108=[0]*len(Tb108);yc90_108=[0]*len(Tb108)
xc90_240_,yc90_240_,xc90_160_,yc90_160_,xc90_108_,yc90_108_=0,0,0,0,0,0

Tball=[0]*9;xc90=[0]*9;yc90=[0]*9;maxX=[0]*9;maxY=[0]*9;Tbmax=[0]*9;xbsub50=[0]*9;ybsub50=[0]*9
Tbsuball=[0]*9;xcsub90=[0]*9;ycsub90=[0]*9;maxsubX=[0]*9;maxsubY=[0]*9
for j in range(9):
    print(freq[j],' MHz')
    bb=pickle.load(open('/media/rohit/MWA/20140914/pickle/Tb_20140914_'+flist[j]+'_sub_test.p','rb'),encoding='bytes')
    aa=pickle.load(open('/media/rohit/MWA/20140914/pickle/Tb_20140914_'+flist[j]+'.p','rb'),encoding='bytes')
    Tball[j]=np.array(aa[0]);xc90[j]=[0]*len(Tball[j]);yc90[j]=[0]*len(Tball[j]);maxX[j]=[0]*len(Tball[j]);maxY[j]=[0]*len(Tball[j])
    for i in range(len(Tball[j])):
        if(Tball[j][i].shape[0]==200):
            Tb_=Tball[j][i][50:-50,50:-50]
        else:
            Tb_=Tball[j][i];Tb_[np.isnan(Tb_)]=0
        y,x=np.where(Tb_==np.max(Tb_))
        maxX[j][i]=50.*(x[0]-50.);maxY[j][i]=50.*(y[0]-50.)
        if(np.nanmax(Tb_)!=0):
            bi=ut.get_bimage(Tb_,0.3);xc90_,yc90_,w,l,angle=ut.fitEllipse(bi)
        else:
            xc90_,yc90_=0,0
        xc90[j][i]=50.*(xc90_-50.);yc90[j][i]=50.*(yc90_-50.)
    #---- Sub-----------#########
    Tbsuball[j]=np.array(bb[0]);xcsub90[j]=[0]*len(Tbsuball[j]);ycsub90[j]=[0]*len(Tbsuball[j]);maxsubX[j]=[0]*len(Tbsuball[j]);maxsubY[j]=[0]*len(Tbsuball[j]);xbsub50[j]=[0]*len(Tbsuball[j]);ybsub50[j]=[0]*len(Tbsuball[j])
    for i in range(len(Tbsuball[j])):
        if(Tbsuball[j][i].shape[0]==200):
            Tb_=Tbsuball[j][i][50:-50,50:-50]
        else:
            Tb_=Tbsuball[j][i];Tb_[np.isnan(Tb_)]=0
        y,x=np.where(Tb_==np.max(Tb_))
        maxsubX[j][i]=50.*(x[0]-50.);maxsubY[j][i]=50.*(y[0]-50.)
        if(np.nanmax(Tb_)!=0):
            bi=ut.get_bimage(Tb_,0.9);xc90_,yc90_,w,l,angle=ut.fitEllipse(bi)
            bi50=ut.get_bimage(Tb_,0.6);xc50_,yc50_,w50,l50,angle=ut.fitEllipse(bi50)
        else:
            xc90_,yc90_=0,0;w50,l50=0,0
        xcsub90[j][i]=50.*(xc90_-50.);ycsub90[j][i]=50.*(yc90_-50.)
        xbsub50[j][i]=l50;ybsub50[j][i]=w50
    Tbmax[j]=Tball[j][0:1323].max(axis=(1,2))
    xc90[j]=xc90[j][0:1323];yc90[j]=yc90[j][0:1323]
#pa_angle=np.arctan((np.array(xcsub90[7])-np.array(xcsub90[1]))/(np.array(ycsub90[7])-np.array(ycsub90[1])))*180/np.pi
pa_angle=np.arctan((np.array(xcsub90)[5:].mean(axis=0)-np.array(xcsub90[:2]).mean(axis=0))/(np.array(ycsub90[5:]).mean(axis=0)-np.array(ycsub90[:2]).mean(axis=0)))*180/np.pi
Tball=np.array(Tball);xc90=np.array(xc90);yc90=np.array(yc90);xbsub50=np.array(xbsub50);ybsub50=np.array(ybsub50);maxX=np.array(maxX);maxY=np.array(maxY)
Tbsuball=np.array(Tbsuball);xcsub90=np.array(xcsub90);ycsub90=np.array(ycsub90);maxsubX=np.array(maxsubX);maxsubY=np.array(maxsubY)
Tbsuball[Tbsuball>5.e8]=np.nan
tmwa=np.array(aa[10]);tsubmwa=np.array(bb[5]);pa_angle=np.array(pa_angle);Tbmax=np.array(Tbmax);Tbsubmax=np.nanmax(Tbsuball,axis=(2,3))
xcsub90_xray=xcsub90-750;ycsub90_xray=ycsub90+300;rcsub90_xray=np.sqrt(xcsub90_xray**2 + ycsub90_xray**2);
pickle.dump([tmwa,tsubmwa,Tbmax,Tbsubmax,xcsub90,ycsub90,maxsubX,maxsubY,pa_angle],open('/media/rohit/MWA/20140914/Tb_centroid.p','wb'))

tmwa,tsubmwa,Tbmax,Tbsubmax,xcsub90,ycsub90,maxsubX,maxsubY,pa_angle = pickle.load(open('/media/rohit/MWA/20140914/Tb_centroid.p','rb'))
rcsub90=np.sqrt(xcsub90**2+ycsub90**2)



znewkirk=np.linspace(0,200,1000)
freq_arr=np.linspace(100,240,100)
r_arr=6.99e5*1.e-3*(model.nk_freq2r(freq_arr,1)[0]-1)
ne_arr=model.nk_freq2r(freq_arr,1)[1]
romegae_mwa=6.99e5*1.e-3*(model.nk_freq2r(freq,1)[0]-1)

bshp1x_mean=[0]*9;bshp1y_mean=[0]*9
for i in range(9):
    idx=np.where((romegae_mwa[i]<zs1[1:]) & (zs1[1:]<romegae_mwa[i]+1) & (bshp1x>700))
    bshp1x_mean[i] = bshp1x[idx].mean();bshp1y_mean[i] = bshp1y[idx].mean() 
rcsub90_mean_freq=np.sqrt((xcsub90.mean(axis=1)-bshp1x_mean)**2 +(ycsub90.mean(axis=1)-bshp1y_mean)**2)
rcsub90_std_freq=xcsub90[:,200:400].std(axis=1)+ycsub90[:,200:400].std(axis=1)
Tbsubstd=np.nanstd(Tbsuball[:,0:500,0:500,0:500],axis=(1,2,3))*5
#-----------------------------------------------

i=0;Tbmaxfit=[0]*len(Tbsubmax[0]);alpha=[0]*len(Tbsubmax[0]);ealpha=[0]*len(Tbsubmax[0])
for i in range(len(Tbsubmax[0])):
    freqfit=np.log10(freq)
    coef = np.polyfit(freqfit,np.log10(Tbsubmax[:,i]),1,cov=True);poly1d_fn = np.poly1d(coef[0]);alpha[i]=coef[0][0]
    ealpha[i]=3*np.sqrt(coef[1][0][0])
    #coef = np.polyfit(freqfit,np.log10(Tbsubmax[:,i]*2*1.38e-23*freqfit*freqfit*1.e12/(9.e16)*1.e-4),1);poly1d_fn = np.poly1d(coef);alpha[i]=coef[0]
    Tbmaxfit[i]=10**poly1d_fn(freqfit)

xbsub50_psf=xbsub50*0;ybsub50_psf=ybsub50*0;beam_psf=ybsub50*0
for i in range(899):
    xbsub50_psf[:,i]=xbsub50[:,i]/(beam[:,1]*60./50.)
    ybsub50_psf[:,i]=ybsub50[:,i]/(beam[:,0]*60./50.)
    beam_psf[:,i]=xbsub50_psf[:,i]*ybsub50_psf[:,i]



def get_pol(aa):
    XXamp=aa[5][0];YYamp=aa[5][3]
    XYamp=aa[5][1]
    XYphase=aa[6][1]
    YXamp=aa[5][2]
    YXphase=aa[6][2]
    V=(XYamp*np.sin(XYphase)-YXamp*np.sin(YXphase))/2
    I=(XXamp+YYamp)/2
    return I,V



ff1='/media/rohit/MWA/20140914/pickle/20140914_084-085.ms_T009-010.p'
aa=pickle.load(open(ff1,'rb'),encoding='bytes')
I108,V108=get_pol(aa)
ff2='/media/rohit/MWA/20140914/pickle/20140914_187-188.ms_T009-010.p'
bb=pickle.load(open(ff2,'rb'),encoding='bytes')
I240,V240=get_pol(bb)
filelist=sorted(glob.glob('/media/rohit/MWA/20140914/pickle/20140914_*.ms_T009-010.p'))
flist=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
freq=np.array([108,120,133,145,161,179,197,217,240]);I=[0]*9;V=[0]*9
for i in range(9):
    aa=pickle.load(open(filelist[i],'rb'),encoding='bytes')
    I_,V_=get_pol(aa);I[i]=I_[0];V[i]=V_[0]
I=np.array(I);V=np.array(V)

VbyI_sub=-1*V[:,270:270+899]/I[:,270:270+899]
plt.imshow(V/I*-1,aspect='auto',origin='lower',vmin=0,vmax=0.8,cmap='jet')
plt.xlabel('Time (HH:MM:SS UT)');plt.ylabel('Frequency (MHz)')
plt.yticks(np.arange(9),freq)
plt.xticks(np.arange(53)[::10]*25,aa[14][::10]);plt.colorbar(label="(V/I)")
plt.show()

#--- Plotting the centroids plot
hTb=np.nanmean(Tbsubmax[4:,:],axis=0)/1.e6;hTb[hTb<0] = np.nan
fTb=np.nanmean(Tbsubmax[0:4,:],axis=0)/1.e6;fTb[fTb<0] = np.nan
f,ax=plt.subplots(1,1)
ax.plot(hTb,color='red',label='161-240 MHz')
ax.plot(fTb,color='blue',label='108-145 MHz')
ax.axvline(x=379,color='k',linestyle='--',label='03:20:00 UT')
ax.legend();ax.set_ylim(0,200);ax.set_ylabel('$T_B$ (MK)');ax.set_xlabel('Time (HH:MM:SS)')
ax.set_xticks([0,200,400,600,800]);ax.set_xticklabels(['02:15:04','02:48:24','03:22:04','03:55:24','04:29:04'])
ax.text(10, 160, 'Fundamental Dominated', bbox=dict(facecolor='red', alpha=0.5))
ax.text(400, 160, 'Harmonic Dominated', bbox=dict(facecolor='green', alpha=0.5))
ax.axvspan(0, 379, color='red',alpha=0.2);ax.axvspan(379, 900, color='green', alpha=0.2)
plt.show()


hTb=np.nanmean(Tbsubmax[4:,:],axis=0)/1.e6;hTb[hTb<0] = np.nan
fTb=np.nanmean(Tbsubmax[0:4,:],axis=0)/1.e6;fTb[fTb<0] = np.nan
f,ax=plt.subplots(1,1)
ax.plot(hTb,color='red',label='161-240 MHz')
ax.plot(fTb,color='blue',label='108-145 MHz')
ax.axvline(x=404,color='brown',linestyle='--',label='03:22:24 UT') # 379 / 03:20:00 /Tb changed
ax.legend();ax.set_ylim(0,200);ax.set_ylabel('$T_B$ (MK)');ax.set_xlabel('Time (HH:MM:SS UT)')
ax.set_xticks([0,200,400,600,800]);ax.set_xticklabels(['02:15:04','02:48:24','03:22:04','03:55:24','04:29:04'])
ax.text(10, 180, 'Harmonic Dominated', bbox=dict(facecolor='red', alpha=0.5))
ax.text(428, 180, 'Fundamental Dominated', bbox=dict(facecolor='green', alpha=0.5))
ax.axvspan(0, 404, color='red',alpha=0.2);ax.axvspan(404, 900, color='green', alpha=0.2)
ax1=ax.twinx()
ax1.plot(np.arange(rcsub90.shape[1])[70:],np.mean(rcsub90[6:,70:],axis=0),'o-',color='k',markersize=2,linewidth=0.5);ax1.set_ylabel('Radial Coordinate (arcsec)')
ax1.set_ylim(750,1000)
plt.show()


f,ax=plt.subplots(1,1)
ax.imshow(Tbsubmax/1.e6,aspect='auto',origin='lower',vmin=0,vmax=200,cmap='coolwarm')
ax.contour(Tbsubmax/1.e6,[135],colors='k',linewidths=1)
ax.set_ylabel('$T_B$ (MK)');ax.set_xlabel('Time (HH:MM:SS)')
ax.set_xticks([0,200,300,400,500,600,800]);ax.set_xticklabels(['02:15:04','02:48:24','03:05:04','03:22:04','03:38:44','03:55:24','04:29:04'])
ax.set_yticks(np.arange(9));ax.set_yticklabels(freq);ax.set_xlim([190,610])
plt.show()

f,ax=plt.subplots(1,1)
ax.plot(freq,Tbsubmax[:,414]/1.e6,'o-',label='03:24:04 UT')
ax.legend();ax.set_ylabel('$T_B$ (MK)');ax.set_xlabel('Frequency (MHz)')
ax.axvline(x=161,linestyle='--',color='k')
ax.text(120, 161, 'Harmonic Emission', bbox=dict(facecolor='red', alpha=0.5))
ax.text(180, 161, 'Fundamental Emission', bbox=dict(facecolor='green', alpha=0.5))
ax.axvspan(108, 161, color='red',alpha=0.2);ax.axvspan(161, 240, color='green', alpha=0.2)
ax.set_xlim(100,250)
plt.show()

#############

rhessi0="/data/Dropbox/20140914/rhessi_analysis/hsi_20140914_010340_002.fits";r0=fits.open(rhessi0);r0_h,r0_d=r0[0].header,r0[0].data
rhessi1="/data/Dropbox/20140914/rhessi_analysis/hsi_20140914_021700_002.fits";r1=fits.open(rhessi1);r1_h,r1_d=r1[0].header,r1[0].data
rh_counts=readsav('/media/rohit/MWA/20140914/rhessi/counts_obs.sav')



############333
tmwa,tsubmwa,Tbmax,Tbsubmax,xcsub90,ycsub90,maxsubX,maxsubY,pa_angle=pickle.load(open('/media/rohit/MWA/20140914/Tb_centroid.p','rb'))


################ PA_angle Wavelet
def do_wavelet(Tbr, N):
    #datawave_=(Tbr-Tbr[0])/1.e6
    datawave_=Tbr
    #datawave_=np.convolve(datawave_, np.ones(N)/N, mode='valid')
    datawave=datawave_[int(N/2):int(-1*N/2+1)]-np.convolve(datawave_, np.ones(N)/N, mode='valid')
    std=np.std(datawave);datawave_std=datawave/std
    timewave=np.arange(len(datawave))*10
    #mother=wavelet.Morlet(6)
    mother=wavelet.MexicanHat()
    dt=10;s0 = 2 * dt;dj = 1. / 12 # Lowest scale s0
    J = 8. / dj # Number of scales -1; Largest scale: s0 * 2**(J * dj)
    var=std**2
    alpha, _, _ = wavelet.ar1(datawave_std)
    wave, scales, freqs, coi1, fft, fftfreqs = wavelet.cwt(datawave_std, dt, dj, s0, J,mother)
    iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std # Inverse CWT
    power = (np.abs(wave)) ** 2;fft_power = np.abs(fft) ** 2;period = 1 / freqs;power /= scales[:, None]
    signif, fft_theor = wavelet.significance(var, dt, scales, 0, alpha,
                                             significance_level=0.95,
                                             wavelet=mother)
    return datawave_std,power,iwave,coi1,period,signif,fft_theor,var

def plot_wavelet(Tb,t,dt,period,power,coi1,fil,tfreq,signif):
    plt.ioff();label='T$_B$';units='(deg)'
    figprops = dict(figsize=(20, 20), dpi=72)
    fig = plt.figure(**figprops)
    ax = plt.axes([0.1, 0.68, 0.67, 0.2])
    #ax.plot(t, iwave[i], '-', linewidth=1, color=[0.5, 0.5, 0.5])
    ax.plot(t,Tb,'o-', 'k', linewidth=1.5,markersize=1)
    #ax.axvline(x=t[195],linestyle='--',color='black')
    ax.grid(True)
    ax.set_ylabel(r'{}{}'.format(label, units))
    bx = plt.axes([0.1, 0.37, 0.65, 0.28], sharex=ax);levels = np.linspace(0.005,0.1,90)
    #bx.contourf(t, periodds, np.log2(powerds), np.log2(levels),extend='both', cmap=plt.cm.YlOrRd)
    bx.contourf(t, period, np.log2(power),np.log2(levels),extend='both', cmap=plt.cm.YlOrRd)
    extent = [t.min(), t.max(), 0, max(period)]
    bx.fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt, t[:1] - dt, t[:1] - dt]),np.concatenate([coi1, [1], period[-1:],period[-1:], [1]]),'k', alpha=0.3, hatch='x')
    #bx.set_title('b) Wavelet Power Spectrum: Freq:'+str(freq[i])+' GHz({})'.format(label, 'MORLET'))
    bx.set_ylabel('Period (sec)');bx.set_xlabel('Time (sec)');bx.set_yscale('log')
    bx.grid(True)
    cx = plt.axes([0.77, 0.37, 0.2, 0.28], sharey=bx)
    cx.plot(power.mean(axis=1), period, 'k-', linewidth=1.5)
    cx.set_title('Global Wavelet Spectrum')
    cx.set_xlabel('Power (deg$^2$)')
    #cx.set_ylim(([periodds.min(), periodds.max()]))
    cx.grid(True)
    cx.set_xscale('log');plt.title('Frequency: '+str(tfreq)+' GHz')
    cx.plot(signif, period, 'k--')
    plt.savefig(fil,dpi=100)
    plt.show()

#pa_angle[558]=pa_angle[555];pa_angle[559]=pa_angle[556]
pa_angle[628:632] = pa_angle[620];pa_angle[536:538]=pa_angle[532]
pa_angle1=pa_angle[0:897].reshape(3,299).mean(axis=0)
pa_angle2=pa_angle[120:315]
pa_angle3=pa_angle[340:860];pa_angle4=np.hstack((pa_angle2,pa_angle3))

N=40;dt=10;t2=np.arange(len(pa_angle2)-N+1)*dt
Tbrs_wave2,s_power2,s_iwave2,s_coi2,period2,signif2,fft_theor2,var2=do_wavelet(pa_angle2, N);ns_power2=s_power2/np.nanmax(s_power2)
plot_wavelet(Tbrs_wave2,t2,dt,period2,s_power2,s_coi2,'power0.png','',signif2)

N=40;dt=10;t3=np.arange(len(pa_angle3)-N+1)*dt
Tbrs_wave3,s_power3,s_iwave3,s_coi3,period3,signif3,fft_theor3,var3=do_wavelet(pa_angle3, N);ns_power3=s_power3/np.nanmax(s_power3)
plot_wavelet(Tbrs_wave3,t3,dt,period3,s_power3,s_coi3,'power1.png','',signif3)

f,ax=plt.subplots(1,1)
ax.plot(period,np.median(ns_power2,axis=1),'o-',label='Harmonic Phase')
ax.plot(period,np.median(ns_power3,axis=1),'o-',label='Fundamental Phase')
ax.set_ylabel('Normalised Power');ax.set_xlabel('Period (sec)');ax.legend()
plt.show()

############### AIA

ddiff_list=sorted(glob.glob('/sdata/fits/running_diff/*171*ddiff.fits'))
bdiff_list=sorted(glob.glob('/sdata/fits/running_diff/*171*bdiff.fits'))
linex=np.arange(3100,3500);liney=np.arange(1270,1670)[::-1]
line_check=1
if line_check:
    mapp=Map(ddiff_list[0]);d=mapp.data
    mapp.plot()
    plt.plot(linex,liney,color='k')
    plt.show()

intd171=[0]*len(ddiff_list)
time171=[0]*len(ddiff_list)
for i in range(len(ddiff_list)):
    mapp=Map(ddiff_list[i]);d=mapp.data;intd171[i]=[0]*len(linex)
    time171[i]=ddiff_list[i].split('T')[1].split('Z')[0]
    for j in range(len(linex)):
        intd171[i][j]=d[liney[j]-5:liney[j]+5,linex[j]-5:linex[j]+5].mean()
intd171=np.array(intd171).swapaxes(0,1)[:,::-1]

intb171=[0]*len(bdiff_list)
for i in range(len(bdiff_list)):
    mapp=Map(bdiff_list[i]);d=mapp.data
    intb171[i]=d[liney,linex]
intb171=np.array(intb171).swapaxes(0,1)[:,::-1]



ddiff211_list=sorted(glob.glob('/sdata/fits/running_diff/*211*ddiff.fits'))
bdiff211_list=sorted(glob.glob('/sdata/fits/running_diff/*211*bdiff.fits'))


sys.exit()

colors = iter(cm.jet(np.linspace(0, 1,9)))
f,(ax0,ax1)=plt.subplots(1,2)
for i in range(9):
    ax0.plot(xcsub90[i].mean()-bshp1x_mean[i],ycsub90[i].mean()-bshp1y_mean[i],'o',markersize=10,markeredgecolor='k',alpha=1.0,label=str(freq[i])+' MHz')
    ax0.plot(xcsub90[i]-bshp1x_mean[i],ycsub90[i]-bshp1y_mean[i],'o',markersize=2.0,color=next(colors),alpha=0.5)
ax0.legend();ax0.set_ylim(-200,300);ax0.set_xlim(-200,300);ax0.set_xlabel('X-Coordinate');ax0.set_ylabel('Y-Coordinate')
ax1.plot(freq,rcsub90_mean_freq,'o-');ax1.set_xlabel('Frequency (MHz)');ax1.set_ylabel('Separation (arcsec)')
ax1.errorbar(freq,rcsub90_mean_freq,yerr=rcsub90_std_freq)
plt.show()

f,ax=plt.subplots(2,1,sharex=True);ax0=ax[0];ax1=ax[1]
for i in [0,8]:
    ax0.plot(Tbsubmax[i]/1.e6,'o-',markersize=2,label=str(freq[i])+' MHz')
    ax1.plot(ycsub90[i],'o-',markersize=2,label=str(freq[i])+' MHz')
ax0.set_ylim(0,250);ax1.set_ylim(-600,-200)
ax0.legend();ax1.legend();ax1.set_ylabel('Y-Coordinate');ax0.set_ylabel('$T_B$ (MK)')
plt.show()


f,(ax0) = plt.subplots(1,1)
ax0.plot(np.array(xcsub90)[6:].mean(axis=0),'o',label=str(np.round(freq[6]))+'-240 MHz')
ax0.plot(np.array(xcsub90)[0:2].mean(axis=0),'o',label='108-'+str(np.round(freq[2]))+' MHz')
ax0.legend();ax0.set_ylim([700,1000]);ax0.set_ylabel('Y-Coordinate (arcsec)');ax0.set_xlabel('Time (HH:MM UT)')
ax0.set_xticks([0,100,200,300,400,500,600,700,800]);ax0.set_xticklabels(['02:15','02:35','02:55','03:15','03:35','03:55','04:15','04:35','04:55'])
plt.show()

f,ax=plt.subplots(2,1,sharex=True);ax0=ax[0];ax1=ax[1]
im0=ax0.imshow(Tbsubmax,origin='lower',aspect='auto',vmin=1.e6,vmax=2.e8)
ax0.contour(Tbsubmax,[5.e7],colors='white',linewidths=0.5)
ax0.contour(Tbsubmax,[1.e8],colors='k',linewidths=0.5)
ax0_divider = make_axes_locatable(ax0);cax0 = ax0_divider.append_axes("right", size="1%", pad="0%")
f.colorbar(im0,label='$T_B$ (K)',cax=cax0)
im1=ax1.imshow(rcsub90_xray,origin='lower',aspect='auto',vmin=0,vmax=200)
ax1.contour(rcsub90_xray,[100],colors='white',linewidths=0.5)
ax1.contour(rcsub90_xray,[150],colors='white',linewidths=0.5)
ax1.set_xlabel('Time');ax1.set_ylabel('Frequency (MHz)');ax0.set_ylabel('Frequency (MHz)')
ax1_divider = make_axes_locatable(ax1);cax1 = ax1_divider.append_axes("right", size="1%", pad="0%")
f.colorbar(im1,label='Radial Coordinate (arcsec)',cax=cax1)
plt.show()

f,ax=plt.subplots(2,1,sharex=True);ax0=ax[0];ax1=ax[1]
im0=ax0.imshow(Tbsubmax,origin='lower',aspect='auto',vmin=1.e6,vmax=2.e8)
ax0.contour(Tbsubmax,[5.e7],colors='white',linewidths=0.5)
ax0.contour(Tbsubmax,[1.e8],colors='k',linewidths=0.5)
ax0_divider = make_axes_locatable(ax0);cax0 = ax0_divider.append_axes("right", size="1%", pad="0%")
f.colorbar(im0,label='$T_B$ (K)',cax=cax0)
im1=ax1.imshow(ycsub90,origin='lower',aspect='auto',vmin=-500,vmax=-200)
ax1.contour(ycsub90,[-370],colors='white',linewidths=0.5)
ax1.contour(ycsub90,[-300],colors='white',linewidths=0.5)
ax1.set_xlabel('Time');ax1.set_ylabel('Frequency (MHz)');ax0.set_ylabel('Frequency (MHz)')
ax1_divider = make_axes_locatable(ax1);cax1 = ax1_divider.append_axes("right", size="1%", pad="0%")
f.colorbar(im1,label='Y-Coordinate (arcsec)',cax=cax1)
plt.show()

for i in range(len(Tbsubmax[0])):
    ii="%04d"%i
    f,ax=plt.subplots(1,1)
    ax.errorbar(freq,Tbsubmax[:,i]/1.e6,yerr=Tbsubstd/1.e6,color='k',label='MWA ($T_B$)')
    ax.plot(freq,Tbsubmax[:,i]/1.e6,'o-',color='k');ax.set_title(str(bb[4][i]))
    ax.plot(freq,Tbmaxfit[i]/1.e6,'--',color='k',label='$\\alpha=$'+str(np.round(alpha[i],2))+'$\\pm$'+str(np.round(ealpha[i],2)))
    ax.set_xlabel('Frequency (MHz)');ax.set_ylabel('$T_B$ (MK)');ax.set_ylim(1,250);ax.set_yscale('log');ax.set_xscale('log');ax.legend()
    f.savefig('/sdata/fits/pngs_spec/spec_'+str(ii)+'.png')
    plt.close()

f,ax=plt.subplots(1,1)
ax.plot(alpha,label='')
ax.set_ylabel('$\\alpha$');ax.set_xlabel('Time (12-sec)')
plt.show()

f,ax=plt.subplots(1,1)
ax1=ax.twinx()
ax.plot(tsubmwa,np.array(xcsub90[-1]),'o-',markersize=4.5,label='X-Coordinate')
ax.plot(tsubmwa,np.array(ycsub90[-1]),'o-',markersize=4.5,label='Y-Coordinate')
ax.plot(tmwa,np.array(maxX[-1]),'o-',markersize=0.5,label='X-Coordinate')
ax.plot(tmwa,np.array(maxY[-1]),'o-',markersize=0.5,label='Y-Coordinate')
ax1.plot(tmwa,Tball[-1].max(axis=(1,2)),'o-',markersize=4.5,color='k');ax.legend()
plt.show()


f,ax0=plt.subplots(1,1)
ax0.imshow(intd171[:,::-1],aspect='auto',origin='lower',vmin=-60,vmax=60,cmap='coolwarm')
ax0.set_xticks(np.arange(1311)[::200]);ax0.set_xticklabels(time171[::200]);ax0.set_xlabel('Time (HH:MM:SS)')
ax0.set_yticks(np.arange(400)[::50]);ax0.set_yticklabels(np.arange(400)[::50]*726/1.e3);ax0.set_ylabel('Radial Distance (Mm)')
#ax1=ax0.twinx()
#ax1.plot(Tbmax[1],color='magenta',label='120 MHz');ax1.plot(Tbmax[7],color='r',label='217 MHz');ax1.legend();ax1.set_ylim(0,80.e6)
ax2=ax0.twinx()
ax2.plot(-1*V[0]/I[0],'-',color='k')
ax2.axvline(x=460,linestyle='--',color='orange')
ax2.axvline(x=240,linestyle='--',color='orange')
ax2.axvline(x=180,linestyle='--',color='orange')
plt.show()

f,ax=plt.subplots(1,1)
ax.plot(freq,Tbsuball.max(axis=(2,3))[:,200:300].mean(axis=1)/1.e6,'o',color='blue')
ax.errorbar(freq,Tbsuball.max(axis=(2,3))[:,200:300].mean(axis=1)/1.e6,yerr=Tbsuball.max(axis=(2,3))[:,250:300].std(axis=1)/1.e6,label='Rising Phase',color='blue')
ax.plot(freq,Tbsuball.max(axis=(2,3))[:,600:650].mean(axis=1)/1.e6,'o',color='red')
ax.errorbar(freq,np.nanmean(Tbsuball.max(axis=(2,3))[:,600:650],axis=1)/1.e6,yerr=np.nanstd(Tbsuball.max(axis=(2,3))[:,600:650],axis=1)/1.e6,label='Decay Phase',color='red')
ax.legend();ax.set_xlabel('Frequency (MHz)');ax.set_ylabel('$T_B$ (MK)')
plt.show()

### Plot Spectrum
Tbsub_spec=Tbsubmax[:,250:550].mean(axis=1)
Tbsub_err = Tbsubmax[:,0:200].std(axis=1)
f,ax=plt.subplots(1,1)
ax.plot(freq,Tbsub_spec/1.e6,'o',color='k')
ax.errorbar(freq, Tbsub_spec/1.e6, yerr = Tbsub_err/1.e6,color='k')
ax.set_xlabel('Frequency (MHz)'); ax.set_ylabel('$T_B$ (MK)')
ax.axvline(x=freq[4],linestyle = ":")
ax.axvspan(100, freq[4], color='red',alpha=0.2);ax.axvspan(freq[4], 250, color='green', alpha=0.2)
ax.text(110, 130, 'Fundamental Emission', bbox=dict(facecolor='red', alpha=0.5))
ax.text(185, 90, 'Harmonic Emission', bbox=dict(facecolor='green', alpha=0.5))
ax.set_xlim(100,250)
plt.show()

f,(ax0,ax1,ax2)=plt.subplots(3,1,sharex=True)
im0=ax0.imshow(Tbmax[:,400:1299]/1.e6,aspect='auto',cmap='jet',origin='lower',vmin=1.e0,vmax=1.e2,interpolation='nearest')
ax0.set_xlabel('Time');ax0.set_ylabel('Frequency (MHz)');ax0_divider = make_axes_locatable(ax0);cax0 = ax0_divider.append_axes("right", size="1%", pad="0%")
ax0.contour(VbyI_sub,[0.4],colors='k',linewidths=2.0)
ax0.contour(VbyI_sub,[0.5],colors='k',linewidths=2.0)
ax0.contour(VbyI_sub,[0.6],colors='k',linewidths=2.0)
f.colorbar(im0,label='$T_B$ (MK)',cax=cax0)
im1=ax1.imshow(np.sqrt(xcsub90**2 + ycsub90**2),aspect='auto',origin='lower',vmin=800,vmax=1100,cmap='jet',interpolation='None')
ax1.set_xlabel('Time');ax1.set_ylabel('Frequency (MHz)');ax1_divider = make_axes_locatable(ax1);cax1 = ax1_divider.append_axes("right", size="1%", pad="0%")
f.colorbar(im1,label='Radial Coordinate (arcsec)',cax=cax1)
ax2.plot(pa_angle+45,'o');ax2.set_xlabel('Time (HH:MM)');ax2.set_ylabel('Position Angle (deg) w.r.t radial vector')
ax2.set_xticks([0,100,200,300,400,500,600,700,800]);ax2.set_xticklabels(['02:15','02:35','02:55','03:15','03:35','03:55','04:15','04:35','04:55'])
ax1.set_yticks([0,1,2,3,4,5,6,7,8]);ax1.set_yticklabels(freq)
ax0.set_yticks([0,1,2,3,4,5,6,7,8]);ax0.set_yticklabels(freq)
plt.show()

xcsub90[:,0:120]=np.nan;ycsub90[:,0:120]=np.nan

f,(ax0,ax1,ax2)=plt.subplots(3,1,sharex=True)
im0=ax0.imshow(Tbmax[:,270:270+899]/1.e6,aspect='auto',cmap='coolwarm',origin='lower',vmin=1.e0,vmax=1.e2,interpolation='nearest')
ax0.set_ylabel('Frequency (MHz)');ax0_divider = make_axes_locatable(ax0);cax0 = ax0_divider.append_axes("right", size="1%", pad="0%")
ax0.contour(VbyI_sub,[0.4],colors='k',linewidths=2.0)
ax0.contour(VbyI_sub,[0.5],colors='k',linewidths=2.0)
ax0.contour(VbyI_sub,[0.6],colors='k',linewidths=2.0)
ax0.axvline(x=372,linestyle='--',linewidth=3,color='green')
f.colorbar(im0,label='$T_B$ (MK)',cax=cax0)
im1=ax1.imshow(np.sqrt(xcsub90**2 + ycsub90**2),aspect='auto',origin='lower',vmin=750,vmax=1000,cmap='coolwarm',interpolation='None')
ax1.set_ylabel('Frequency (MHz)');ax1_divider = make_axes_locatable(ax1);cax1 = ax1_divider.append_axes("right", size="1%", pad="0%")
ax1.axvline(x=372,linestyle='--',linewidth=3,color='green')
f.colorbar(im1,label='Radial Coordinate (arcsec)',cax=cax1)
ax2.plot(pa_angle,'o',markersize=2);ax2.set_xlabel('Time (HH:MM)');ax2.set_ylabel('Orientation Angle (deg)')
ax2.set_xticks([0,120,240,360,480,600,720,840]);ax2.set_xticklabels(['02:15','02:35','02:55','03:15','03:35','03:55','04:15','04:35'])
ax2.axvline(x=372,linestyle='--',linewidth=3,color='green')
ax1.set_yticks([0,1,2,3,4,5,6,7,8]);ax1.set_yticklabels(freq)
ax0.set_yticks([0,1,2,3,4,5,6,7,8]);ax0.set_yticklabels(freq)
ax2.set_yticks([-50,-20,0,20,50]);ax2.set_ylim([-95,70])
ax0.text(130,7,'(A)');ax1.text(130,7,'(B)');ax2.text(130,40,'(C)')
plt.show()


bdiff_list=sorted(glob.glob('/sdata/fits/running_diff/*171*bdiff.fits'))[500:505]
for i in range(len(bdiff_list)):
    ii="%04d"%i
    mapp=Map(bdiff_list[i])
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=mapp);cc=['r','g','yellow','b']
    p0=mapp.plot(axes=ax,aspect='auto',vmin=-2000,vmax=2000,cmap='coolwarm')
    bl=SkyCoord(600*u.arcsec,-600*u.arcsec,frame=mapp.coordinate_frame)
    tr=SkyCoord(1000*u.arcsec,-200*u.arcsec,frame=mapp.coordinate_frame)
    blx=mapp.world_to_pixel(bl)[0].value;bly=mapp.world_to_pixel(bl)[1].value
    trx=mapp.world_to_pixel(tr)[0].value;rty=mapp.world_to_pixel(tr)[1].value
    ax.set_xlim([blx,trx]);ax.set_ylim([bly,rty])
    plt.savefig('/sdata/fits/pngs171_bdiff/aia171_bdiff_'+ii+'.png')

ddiff_list=sorted(glob.glob('/sdata/fits/running_diff/*94*bdiff.fits'))
for i in range(150,650):
    ii="%04d"%i
    mapp=Map(ddiff_list[i])
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=mapp);cc=['r','g','yellow','b']
    p0=mapp.plot(axes=ax,aspect='auto',vmin=-5,vmax=5,cmap='coolwarm')
    bl=SkyCoord(600*u.arcsec,-600*u.arcsec,frame=mapp.coordinate_frame)
    tr=SkyCoord(1000*u.arcsec,-200*u.arcsec,frame=mapp.coordinate_frame)
    blx=mapp.world_to_pixel(bl)[0].value;bly=mapp.world_to_pixel(bl)[1].value
    trx=mapp.world_to_pixel(tr)[0].value;rty=mapp.world_to_pixel(tr)[1].value
    ax.set_xlim([blx,trx]);ax.set_ylim([bly,rty])
    #------------------------------------------------------------------
    seeds2 = SkyCoord(710*u.arcsec, -314*u.arcsec,frame=mapp.coordinate_frame)
    ax.plot_coord(seeds2, color='k', marker='*', markersize=15.0,alpha=1.0,linestyle='None');l=0
    #-----------------
    tt=ddiff_list[i].split('T')[1].split('Z')[0].split('_');taiasec=int(tt[0])*3600+int(tt[1])*60+float(tt[2]);k=ut.find_nearest(bb[5],taiasec)[0]
    colors = iter(cm.jet(np.linspace(0, 1, 9)));pix=[0]*9
    for j in range(9):
        seeds0 = SkyCoord(np.array(xcsub90[j])[k]*u.arcsec, np.array(ycsub90[j])[k]*u.arcsec,frame=mapp.coordinate_frame)
        ax.plot_coord(seeds0, color=next(colors), marker='s', markersize=15,alpha=0.8,linestyle='None',markeredgecolor='white')
        pix_=mapp.wcs.world_to_pixel(seeds0);pix[j]=[float(pix_[0]),float(pix_[1])]
    pix=np.array(pix);qx=pix[7:9,0].mean();qy=pix[7:9,1].mean();qx1=pix[0:1,0].mean();qy1=pix[0:1,1].mean()
    ax.quiver(qx,qy,(-qx+qx1)/np.sqrt(qx**2 + qx1**2),(-qy+qy1)/np.sqrt(qy**2 + qy1**2), color='k', scale=0.2)
    plt.savefig('/sdata/fits/pngs94_ddiff/aia193_bdiff3_'+ii+'.png')
    plt.close()

ddiff_list=sorted(glob.glob('/sdata/fits/running_diff/*171*ddiff.fits'))
for i in range(150,450):
    ii="%04d"%i
    mapp=Map(ddiff_list[i])
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=mapp);cc=['r','g','yellow','b']
    p0=mapp.plot(axes=ax,aspect='auto',vmin=-200,vmax=200,cmap='coolwarm')
    bl=SkyCoord(600*u.arcsec,-600*u.arcsec,frame=mapp.coordinate_frame)
    tr=SkyCoord(1000*u.arcsec,-200*u.arcsec,frame=mapp.coordinate_frame)
    blx=mapp.world_to_pixel(bl)[0].value;bly=mapp.world_to_pixel(bl)[1].value
    trx=mapp.world_to_pixel(tr)[0].value;rty=mapp.world_to_pixel(tr)[1].value
    ax.set_xlim([blx,trx]);ax.set_ylim([bly,rty])
    #------------------------------------------------------------------
    hp_lon0 = bshp1x * u.arcsec;hp_lat0 = bshp1y * u.arcsec
    seeds0 = SkyCoord(hp_lon0.ravel(), hp_lat0.ravel(),frame=mapp.coordinate_frame)
    hp_lon4 = bshp2x * u.arcsec;hp_lat4 = bshp2y * u.arcsec
    seeds4 = SkyCoord(hp_lon4.ravel(), hp_lat4.ravel(),frame=mapp.coordinate_frame)
    seeds2 = SkyCoord(710*u.arcsec, -314*u.arcsec,frame=mapp.coordinate_frame)
    ax.plot_coord(seeds0, color='white', marker='s', markersize=0.1,alpha=0.5,linestyle='None')
    #ax.plot_coord(seeds4, color='tab:brown', marker='s', markersize=0.5,alpha=0.8,linestyle='None')
    ax.plot_coord(seeds2, color='tab:cyan', marker='*', markersize=10.0,alpha=1.0,linestyle='None');l=0
    #-----------------
    for k in [1,7]:
        colors = iter(cm.jet(np.linspace(0, 1, 9)));pix=[0]*9
        for i in range(9):
            seeds0 = SkyCoord(np.array(xcsub90[i])[100*k:100*(k+1)].mean()*u.arcsec, np.array(ycsub90[i])[100*k:100*(k+1)].mean()*u.arcsec,frame=mapp.coordinate_frame)
            ax.plot_coord(seeds0, color=next(colors), marker='s', markersize=10,alpha=0.8,linestyle='None',markeredgecolor='white')
            pix[i]=mapp.wcs.world_to_pixel(seeds0)
        #ax.quiver(pix[0][0],pix[0][1],pix[8][0]-pix[0][0],pix[8][1]-pix[0][1], color='r', scale=21)
        ax.quiver(pix[8][0],pix[8][1],4*(-pix[8][0]+pix[0][0])/np.sqrt(pix[8][0]**2 + pix[0][0]**2),4*(-pix[8][1]+pix[0][1])/np.sqrt(pix[8][0]**2 + pix[0][0]**2), color=cc[l], scale=1)
        l=l+1
    colors = iter(cm.jet(np.linspace(0, 1, 9)));pix=[0]*9
    for i in range(9):
        idx=np.where((romegae_mwa[i]<zs1[1:]) & (zs1[1:]<romegae_mwa[i]+1) & (bshp1x>700))
        hp_lon0 = bshp1x[idx].mean() * u.arcsec;hp_lat0 = bshp1y[idx].mean() * u.arcsec
        seeds0 = SkyCoord(hp_lon0.ravel(), hp_lat0.ravel(),frame=mapp.coordinate_frame)
        ax.plot_coord(seeds0, color=next(colors), marker='o', markersize=10,alpha=0.5,linestyle='None',markeredgecolor='black')
    plt.savefig('/sdata/fits/pngs171_ddiff/aia171_ddiff2_'+ii+'.png')
    plt.show()
    plt.close()


f,((ax0,ax1),(ax2,ax3))=plt.subplots(2,2)
#im=ax.scatter(x1,y1,c=bsz1,cmap='jet',vmin=0,vmax=40)
im0=ax0.scatter(rs1*delta,z1*delta,c=bsabs1,cmap='jet',vmin=0,vmax=40,alpha=0.5)
ax0.set_xlabel('Radial Coordinate (arcsec)');ax0.set_ylabel('Z-Coordinate (arcsec)');ax0_divider = make_axes_locatable(ax0)
cax0 = ax0_divider.append_axes("right", size="7%", pad="2%")
f.colorbar(im0,label='|B| (G)',cax=cax0)
ax1.plot(r_arr,freq_arr,'o-',label='$\omega_e$ (Newkirk)')
ax1.plot(np.arange(300)*1.4,babs[:,200,200]/2.8*s,'o-',label='$\Omega_e$ (s=2)');ax1.legend();ax1.set_xlim(10,200);ax1.set_ylabel('Frequency (MHz)'),ax1.set_xlabel('Coronal Height (Mm)')
im2=ax2.scatter(rs1*delta,z1*delta,c=bsabs1/2.8*s,cmap='jet',vmin=100,vmax=240,alpha=0.5);ax2_divider = make_axes_locatable(ax2)
cax2 = ax2_divider.append_axes("right", size="7%", pad="2%")
f.colorbar(im2,label='$\Omega_e$ (MHz) (s=2)',cax=cax2)
ax3.plot(r_arr,ne_arr,'o-');ax3.set_xlabel('Coronal Height (Mm)');ax3.set_ylabel('$n_e$ cm$^{-3}$')
plt.show()

f,ax0=plt.subplots()
im0=ax0.scatter(rs1*delta,z1*delta,c=bsabs1,cmap='cool',vmin=0,vmax=40,alpha=0.5)
ax0.set_xlabel('Radial Coordinate (arcsec)');ax0.set_ylabel('Z-Coordinate (arcsec)');ax0_divider = make_axes_locatable(ax0)
cax0 = ax0_divider.append_axes("right", size="7%", pad="2%")
f.colorbar(im0,label='|B| (G)',cax=cax0);c = iter(cm.jet(np.linspace(0, 1, 9)));pix=[0]*9
for i in range(9):
    ax0.axhline(y=romegae_mwa[i],color=next(c),label=str(freq[i])+' MHz')
ax0.legend()
plt.show()

for i in range(len(Tb240)):
    if(Tb240[i].shape[0]==200):
        Tb240_=Tb240[i][50:-50,50:-50]
    else:
        Tb240_=Tb240[i]
    y,x=np.where(Tb240_==np.max(Tb240_))
    maxX[i]=50.*(x[0]-50.);maxY[i]=50.*(y[0]-50.)
    if(np.nanmax(Tb240_)!=0):
        bi=ut.get_bimage(Tb240_,0.9);xc90_240_,yc90_240_,w,l,angle=ut.fitEllipse(bi)
    else:
        xc90_240_,yc90_240_=0,0
    xc90_240[i]=50.*(xc90_240_-50.);yc90_240[i]=50.*(yc90_240_-50.)
for i in range(len(Tb160)):
    if(Tb160[i].shape[0]==200):
        Tb160_=Tb160[i][50:-50,50:-50]
    else:
        Tb160_=Tb160[i]
    if(np.nanmax(Tb160_)!=0):
        bi=ut.get_bimage(Tb160_,0.9);xc90_160_,yc90_160_,w,l,angle=ut.fitEllipse(bi)
    else:
        xc90_160_,yc90_160_=0,0
    xc90_160[i]=50.*(xc90_160_-50.);yc90_160[i]=50.*(yc90_160_-50.)
for i in range(len(Tb108)):
    if(Tb108[i].shape[0]==200):
        Tb108_=Tb108[i][50:-50,50:-50]
    else:
        Tb108_=Tb108[i]
    if(np.nanmax(Tb108_)!=0):
        bi=ut.get_bimage(Tb108_,0.9);xc90_108_,yc90_108_,w,l,angle=ut.fitEllipse(bi)
    else:
        xc90_108_,yc90_108_=0,0
    xc90_108[i]=50.*(xc90_108_-50.);yc90_108[i]=50.*(yc90_108_-50.)
maxX=np.array(maxX);maxY=np.array(maxY)
tl=300;tr=1200
xc90_240=np.array(xc90_240)[tl:tr];yc90_240=np.array(yc90_240)[tl:tr];xc90_160=np.array(xc90_160)[tl:tr]
yc90_160=np.array(yc90_160)[tl:tr];xc90_108=np.array(xc90_108)[tl:tr];yc90_108=np.array(yc90_108)[tl:tr]

import matplotlib.cm as cm

plot_one=1
if plot_one:
    k=4;l=0
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=aiamap);cc=['r','g','yellow','b']
    #------- 
    hp_lon0 = bshp1x * u.arcsec;hp_lat0 = bshp1y * u.arcsec
    seeds0 = SkyCoord(hp_lon0.ravel(), hp_lat0.ravel(),frame=aiamap.coordinate_frame)
    hp_lon4 = bshp2x * u.arcsec;hp_lat4 = bshp2y * u.arcsec
    seeds4 = SkyCoord(hp_lon4.ravel(), hp_lat4.ravel(),frame=aiamap.coordinate_frame)
    seeds2 = SkyCoord(710*u.arcsec, -314*u.arcsec,frame=aiamap.coordinate_frame)
    ax.plot_coord(seeds0, color='white', marker='s', markersize=0.1,alpha=0.5,linestyle='None')
    #ax.plot_coord(seeds4, color='tab:brown', marker='s', markersize=0.5,alpha=0.8,linestyle='None')
    ax.plot_coord(seeds2, color='tab:cyan', marker='*', markersize=10.0,alpha=1.0,linestyle='None')
    #-----------------
    for k in [1,7]:
        p0=aiamap.plot(axes=ax,aspect='auto')
        colors = iter(cm.jet(np.linspace(0, 1, 9)));pix=[0]*9
        for i in range(9):
            seeds0 = SkyCoord(np.array(xcsub90[i])[100*k:100*(k+1)].mean()*u.arcsec, np.array(ycsub90[i])[100*k:100*(k+1)].mean()*u.arcsec,frame=aiamap.coordinate_frame)
            ax.plot_coord(seeds0, color=next(colors), marker='s', markersize=10,alpha=0.8,linestyle='None',markeredgecolor='white')
            pix[i]=aiamap.wcs.world_to_pixel(seeds0)
        #ax.quiver(pix[0][0],pix[0][1],pix[8][0]-pix[0][0],pix[8][1]-pix[0][1], color='r', scale=21)
        ax.quiver(pix[8][0],pix[8][1],4*(-pix[8][0]+pix[0][0])/np.sqrt(pix[8][0]**2 + pix[0][0]**2),4*(-pix[8][1]+pix[0][1])/np.sqrt(pix[8][0]**2 + pix[0][0]**2), color=cc[l], scale=1)
        l=l+1
    colors = iter(cm.jet(np.linspace(0, 1, 9)));pix=[0]*9
    for i in range(9):
        idx=np.where((romegae_mwa[i]<zs1[1:]) & (zs1[1:]<romegae_mwa[i]+1) & (bshp1x>700))
        hp_lon0 = bshp1x[idx].mean() * u.arcsec;hp_lat0 = bshp1y[idx].mean() * u.arcsec
        seeds0 = SkyCoord(hp_lon0.ravel(), hp_lat0.ravel(),frame=aiamap.coordinate_frame)
        ax.plot_coord(seeds0, color=next(colors), marker='o', markersize=10,alpha=0.5,linestyle='None',markeredgecolor='black')
    #seeds1 = SkyCoord(xc90[0][350:]*u.arcsec, yc90[0][350:]*u.arcsec,frame=aiamap.coordinate_frame)
    #ax.plot_coord(seeds1, color='tab:green', marker='s', markersize=5,alpha=0.8,linestyle='None')
    #bl=SkyCoord(-850*u.arcsec,300*u.arcsec,frame=aiamap.coordinate_frame)
    #tr=SkyCoord(-550*u.arcsec,600*u.arcsec,frame=aiamap.coordinate_frame)
    #blx=aiamap.world_to_pixel(bla)[0].value;bly=aiamap.world_to_pixel(bla)[1].value
    #trx=aiamap.world_to_pixel(tra)[0].value;rty=aiamap.world_to_pixel(tra)[1].value
    #ax.set_xlim([blx,trx]);ax.set_ylim([bly,rty])
    plt.show()

plot_one=1
if plot_one:
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=aiamap)
    p0=aiamap.plot(axes=ax,aspect='auto')
    colors = iter(cm.cool(np.linspace(0, 1, len(xc90[0][350:]))))
    for i in range(len(xc90[0][350:])):
        seeds0 = SkyCoord(xc90[0][350:][i]*u.arcsec, yc90[0][350:][i]*u.arcsec,frame=aiamap.coordinate_frame)
        ax.plot_coord(seeds0, color=next(colors), marker='s', markersize=5,alpha=0.8,linestyle='None')
    colors = iter(cm.hot(np.linspace(0, 1, len(xc90[8][350:]))))
    for i in range(len(xc90[3][350:])):
        seeds0 = SkyCoord(xc90[3][350:][i]*u.arcsec, yc90[3][350:][i]*u.arcsec,frame=aiamap.coordinate_frame)
        ax.plot_coord(seeds0, color=next(colors), marker='s', markersize=5,alpha=0.8,linestyle='None')
    #seeds1 = SkyCoord(xc90[0][350:]*u.arcsec, yc90[0][350:]*u.arcsec,frame=aiamap.coordinate_frame)
    #ax.plot_coord(seeds1, color='tab:green', marker='s', markersize=5,alpha=0.8,linestyle='None')
    #bl=SkyCoord(-850*u.arcsec,300*u.arcsec,frame=aiamap.coordinate_frame)
    #tr=SkyCoord(-550*u.arcsec,600*u.arcsec,frame=aiamap.coordinate_frame)
    #blx=aiamap.world_to_pixel(bla)[0].value;bly=aiamap.world_to_pixel(bla)[1].value
    #trx=aiamap.world_to_pixel(tra)[0].value;rty=aiamap.world_to_pixel(tra)[1].value
    #ax.set_xlim([blx,trx]);ax.set_ylim([bly,rty])
    plt.show()

plot_one=1
if plot_one:
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=bcarrmap[0])
    p0=bcarrmap[0].plot(axes=ax,aspect='auto')
    hp_lon = bs_carr.lon
    hp_lat = bs_carr.lat
    seeds0 = SkyCoord(hp_lon.ravel(), hp_lat.ravel(),frame=bcarrmap[0].coordinate_frame)
    ax.plot_coord(seeds0, color='tab:red', marker='s', markersize=0.5,alpha=0.8,linestyle='None')
    #bl=SkyCoord(-850*u.arcsec,300*u.arcsec,frame=aiamap.coordinate_frame)
    #tr=SkyCoord(-550*u.arcsec,600*u.arcsec,frame=aiamap.coordinate_frame)
    #blx=aiamap.world_to_pixel(bla)[0].value;bly=aiamap.world_to_pixel(bla)[1].value
    #trx=aiamap.world_to_pixel(tra)[0].value;rty=aiamap.world_to_pixel(tra)[1].value
    #ax.set_xlim([blx,trx]);ax.set_ylim([bly,rty])
    plt.show()

#####- PAPER

plot_one=1
if plot_one:
    aiamap = Map('/media/rohit/MWA/20140914/EUV/fits/aia.lev1.171A_2014-09-14T02_00_35.34Z.image_lev1.fits')
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=aiamap)
    p0=aiamap.plot(axes=ax,aspect='auto',vmin=100,vmax=8000)
    hp_lon0 =  bshp0x* u.arcsec;hp_lat0 = bshp0y* u.arcsec
    hp_lon1 = bshp1x * u.arcsec;hp_lat1 = bshp1y* u.arcsec
    #seeds00 = SkyCoord(hp_lon0.ravel(), hp_lat0.ravel(),frame=aiamap.coordinate_frame)
    #ax.plot_coord(seeds00, color='brown', marker='s', markersize=0.5,alpha=0.8,linestyle='None')
    seeds01 = SkyCoord(hp_lon1.ravel(), hp_lat1.ravel(),frame=aiamap.coordinate_frame)
    ax.plot_coord(seeds01, color='tab:blue', marker='s', markersize=0.5,alpha=0.5,linestyle='None')
    seeds1 = SkyCoord(xc90[8][400:600]*u.arcsec, yc90[8][400:600]*u.arcsec,frame=aiamap.coordinate_frame)
    ax.plot_coord(seeds1, color='tab:red', marker='s', markersize=4.0,alpha=0.6,linestyle='None')
    #-------------------
    rhessi_image_list = sorted(glob.glob('/media/rohit/MWA/20140914/rhessi/fitsimage_15-25keV_*.fits'))
    img = fits.open(rhessi_image_list[85]);imgdata = img[0].data;imghead = img[0].header
    aiahead1=aiamap.meta
    aiahead1['naxis1']=imghead['NAXIS1'];aiahead1['naxis2']=imghead['NAXIS2']
    aiahead1['CRPIX1']=imghead['CRPIX1'];aiahead1['CRPIX2']=imghead['CRPIX2']
    aiahead1['CRVAL1']=imghead['CRVAL1'];aiahead1['CRVAL2']=imghead['CRVAL2']
    aiahead1['CDELT1']=imghead['CDELT1'];aiahead1['CDELT2']=imghead['CDELT2']
    img = fits.open(rhessi_image_list[trhessi_idx]);imgdata = img[0].data;imghead = img[0].header
    rhessi_map=Map(imgdata,aiahead1)
    frac_r1=0.9;frac_r2=0.8;frac_r3=0.85
    lev_r1=np.nanmax(rhessi_map.data)*frac_r1
    lev_r2 = np.nanmax(rhessi_map.data) * frac_r2
    lev_r3 = np.nanmax(rhessi_map.data) * frac_r3
    c1=rhessi_map.contour(level=lev_r1* u.ct);ax.plot_coord(c1[0],color='green')
    c2 = rhessi_map.contour(level=lev_r2 * u.ct);ax.plot_coord(c2[0], color='green')
    c3 = rhessi_map.contour(level=lev_r3 * u.ct);    ax.plot_coord(c3[0], color='green')
    #seeds2 = SkyCoord(xc90[4][400:1000]*u.arcsec, yc90[4][400:1000]*u.arcsec,frame=aiamap.coordinate_frame)
    #ax.plot_coord(seeds2, color='tab:green', marker='s', markersize=1.0,alpha=0.6,linestyle='None')
    #seeds3 = SkyCoord(xc90[0][400:1000]*u.arcsec, yc90[0][400:1000]*u.arcsec,frame=aiamap.coordinate_frame)
    #ax.plot_coord(seeds3, color='tab:blue', marker='s', markersize=1.0,alpha=1.0,linestyle='None')
    #seeds4 = SkyCoord(710*u.arcsec, -314*u.arcsec,frame=hmimap.coordinate_frame)
    #ax.plot_coord(seeds4, color='tab:cyan', marker='o', markersize=10.0,alpha=1.0,linestyle='None')
    #bl=SkyCoord(-850*u.arcsec,300*u.arcsec,frame=aiamap.coordinate_frame)
    #tr=SkyCoord(-550*u.arcsec,600*u.arcsec,frame=aiamap.coordinate_frame)
    #blx=aiamap.world_to_pixel(bla)[0].value;bly=aiamap.world_to_pixel(bla)[1].value
    #trx=aiamap.world_to_pixel(tra)[0].value;rty=aiamap.world_to_pixel(tra)[1].value
    #ax.set_xlim([blx,trx]);ax.set_ylim([bly,rty])
    ax.set_xlim([2700,3700]);ax.set_ylim([800,1800])
    plt.show()

plot_hmi=1
if plot_hmi:
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=hmimap)
    p0=hmimap.plot(axes=ax,aspect='auto',vmin=-1250,vmax=1250)
    hp_lon0 = bshp0x * u.arcsec;hp_lat0 = bshp0y * u.arcsec
    seeds0 = SkyCoord(hp_lon0.ravel(), hp_lat0.ravel(),frame=hmimap.coordinate_frame)
    hp_lon4 = bshp2x * u.arcsec;hp_lat4 = bshp2y * u.arcsec
    seeds4 = SkyCoord(hp_lon4.ravel(), hp_lat4.ravel(),frame=hmimap.coordinate_frame)
    seeds1 = SkyCoord(xc90_240.mean()*u.arcsec, yc90_240.mean()*u.arcsec,frame=hmimap.coordinate_frame)
    seeds2 = SkyCoord(710*u.arcsec, -314*u.arcsec,frame=hmimap.coordinate_frame)
    ax.plot_coord(seeds0, color='black', marker='s', markersize=0.5,alpha=0.8,linestyle='None')
    ax.plot_coord(seeds4, color='tab:brown', marker='s', markersize=0.5,alpha=0.8,linestyle='None')
    ax.plot_coord(seeds1, color='tab:red', marker='s', markersize=10.0,alpha=1.0,linestyle='None')
    ax.plot_coord(seeds2, color='tab:cyan', marker='o', markersize=10.0,alpha=1.0,linestyle='None')
    plt.show()


tmin=np.arange(len(xc90_240))*12
plot_xc=1
if(plot_xc):
    f,ax=plt.subplots(2,1,sharex=True);ax0=ax[0];ax1=ax[1]
    ax0.plot(tmin,xc90_240[1:],'o-',label='240 MHz');ax0.plot(tmin,xc90_160,'o-',label='160 MHz');ax0.plot(tmin,xc90_108,'o-',label='108 MHz')
    ax1.plot(tmin,yc90_240[1:],'o-',label='240 MHz');ax1.plot(tmin,yc90_160,'o-',label='160 MHz');ax1.plot(tmin,yc90_108,'o-',label='108 MHz')
    ax0.set_ylabel('X-Coordinate');ax0.legend();ax0.set_ylim(700,1050)
    ax1.set_ylabel('Y-Coordinate');ax1.legend();ax1.set_ylim(-600,-200);ax1.set_xlabel('Time (sec)')
    plt.show()

#fig = plt.figure(figsize=(12, 5))
#ax1 = fig.add_subplot(1, 2, 1, projection=map_aia)


#expolB_=readsav('/sdata/20140914_hmi/20140914_magnetic_extrapolation_1arcsec.sav');expolB=expolB_['box']
#expolB_1=readsav('/sdata/20140914_hmi/20140914_magnetic_extrapolation_024450.sav');expolB1=expolB_1['box']
listB=sorted(glob.glob('/sdata/20140914_hmi/2014-09-14/hmi*NAS.sav'))
for i in range(len(listB)):
    ii="%04d" % i
    #expolB_=readsav(listB[i]);expolB=expolB_['bout']
    expolB_=readsav(listB[i]);expolB=expolB_['box']
    bx,by,bz=expolB['bx'][0],expolB['by'][0],expolB['bz'][0]
    babs=np.sqrt(bx*bx+by*by+bz*bz)
    # Generate the grid
    dim=bx.shape
    zz,yy,xx=np.mgrid[0:dim[0],0:dim[1],0:dim[2]]
    pts = np.empty(bx.shape + (3,), dtype=np.int32)
    pts[..., 0] = xx
    pts[..., 1] = yy
    pts[..., 2] = zz
    vectors = np.empty(bx.shape + (3,), dtype=np.float64)
    vectors[..., 0] = bx
    vectors[..., 1] = by
    vectors[..., 2] = bz
    # We reorder the points and vectors so this is as per VTK's
    # requirement of x first, y next and z last.
    pts = pts.transpose(2, 1, 0, 3).copy()
    pts.shape = pts.size // 3, 3
    vectors = vectors.transpose(2, 1, 0, 3).copy()
    vectors.shape = vectors.size // 3, 3
    sg = tvtk.StructuredGrid(dimensions=xx.shape, points=pts)
    sg.point_data.vectors = vectors
    sg.point_data.vectors.name = 'Magnetic Field'
    absv=np.sqrt(vectors[:,0]**2+vectors[:,1]**2+vectors[:,2]**2)
    sg.point_data.scalars = absv
    sg.point_data.scalars.name = 'Magnitude'
    write_data(sg, listB[i]+'_'+str(ii)+'.vtk')


#source_height=1400;hdu.header.update({'NAXIS':3});hdu.header.append(('CRPIX3',0));hdu.header.append(('CRVAL3',695700.0));hdu.header.append(('CTYPE3','HECH'));hdu.header.append(('CUNIT3','km'));hdu.header.append(('CDELT1',1400))
#hdu.header.append(('RSUN_REF',source_height))

############################





omegap=2.8*babs[:,160,172]/1000 # in GHz
h=np.arange(200)*1400/1.e3 # in Mm
hscale=5000;k=10
ne=1.16e17*np.exp(-1000*h/hscale)+1.e9
fp=9000*np.sqrt(ne)/1.e9
ne5=1.16e17*(1+(h*1000/(10*500)))**(-10+1)
fp5=9000*np.sqrt(ne5)/1.e9
ne2=1.16e17*(1+(h*1000/(10*200)))**(-10+1)
fp2=9000*np.sqrt(ne2)/1.e9
ne3=1.16e17*(1+(h*1000/(10*300)))**(-10+1)
fp3=9000*np.sqrt(ne3)/1.e9
plt.plot(h,omegap,'o-',label='Gyroresonce Freqeuncy')
plt.plot(h,fp,'-',label='Plasma Freqeuncy ($1.e9+1.16e17*exp(-h/h0)$)')
plt.plot(h,fp5,'-',label='Plasma Freqeuncy ($1.16e17*(1+h/(k*h0))^{-k+1}$) h0=500 km')
plt.plot(h,fp2,'-',label='h0=200 km')
plt.plot(h,fp3,'-',label='h0=300 km')
plt.axhline(y=1,color='k',linewidth=2)
plt.xlabel('Coronal Height (Mm)');plt.ylabel('Freqeuncy (GHz)'),plt.legend()
plt.yscale('log');plt.xlim([0,80]);plt.ylim([0.01,100])
plt.show()

plt.plot(h,omegap,'o-',label='Gyroresonce Freqeuncy')
plt.axhline(y=1.4);plt.axhline(y=0.99)
plt.show()
### CANS
d=readsav('/media/rohit/VLA/20160409_sim/code/output/test_param.sav')
x=d['x']*2.e7/1.e8;ro=d['ro']*1.e17;temp=d['te']*1.e4
fpro=9000*np.sqrt(ro)/1.e9
f,ax=plt.subplots(1,1)
ax.plot(x,fpro[0],'o-',label='Plasma Frequency')
ax.plot(h,omegap,'o-',label='Gyro-resonance Freqeuncy')
ax.set_xlabel('Coronal Height (Mm)');ax.set_ylabel('Freqeuncy (GHz)'),ax.legend()
plt.show()

### FORWARD

ff=readsav('/media/rohit/VLA/20160409_EUV/forward/RT_params.sav')
dens=ff['densall'];r=ff['r3dall'];temp=ff['tempall']
plt.plot((r[81,154]-1)*700,dens[81,154],'o-')
plt.xlabel('Coronal Height (Mm)');plt.ylabel('Density ($cm^{-3}$)')
plt.show()
plt.plot((r[81,154]-1)*700,9000*np.sqrt(dens[81,154])/1.e6,'o-')
plt.xlabel('Coronal Height (Mm)');plt.ylabel('Plasma Freqeuncy (MHz)')
plt.show()

f,ax=plt.subplots(1,1)
ax.plot(x,ro[0],'o-',label='Density')
ax1=ax.twinx()
ax1.plot(x,temp[0]/1.e6,'o-',label='Temperature',color='k')
ax.set_ylabel('Density (cm^-3)');ax1.set_ylabel('Temperature (MK)');ax1.legend()
ax.set_xlabel('Coronal Height (Mm)');ax.set_yscale('log')
plt.show()


dem=readsav('/media/rohit/VLA/20160409/demmap.sav')
xc=dem['demmap']['xc'][0];yc=dem['demmap']['yc'][0];densdata=dem['demmap']['data']
demptemp=[0.5,1.0,1.5,2.0,3.0,4.0,6.0,8.0,11.0,14.0]

dem_region=[0]*10
dem_region1=[0]*10
for i in range(10):
    dem_region[i]=densdata[i][80:130,20:35].mean()
    dem_region1[i]=densdata[i][120:160,40:55].mean()

#nk=model.nk_freq2r(np.arange(100)+1000,1)
R=np.linspace(2000,280.e3,1000)
nkne=4.2e4 *10**(4.32/((695700.+R)/695700.))
th=10*np.pi/180.
saitone=(3.09e8*(1-0.5*np.sin(th))/((695700.+R)/695700.)**16)+(1.58e8*(1-0.95*np.sin(th))/((695700.+R)/695700.)**6)+(0.0251e8*(1-np.sqrt(np.sin(th)))/((695700.+R)/695700.)**2.5)
nkfp=9000*np.sqrt(nkne)
saitofp=9000*np.sqrt(saitone)

#plt.plot(R,nkne)
plt.plot(R,nkfp/1.e6,label='Newkirk')
plt.plot(R,saitofp/1.e6,label='SAITO')
plt.legend()
#plt.xlabel('Coronal Height (km)');plt.ylabel('Electron Density (ccm)')
plt.xlabel('Coronal Height (km)');plt.ylabel('Plasma Frequency (MHz)')
plt.show()

ne25=1.e9+1.16e17*np.exp(-1*(h*1000/500.))
fp25=9000*np.sqrt(ne25)/1.e9
plt.plot(h,omegap,'o-',label='Gyroresonce Freqeuncy')
#plt.plot(h,fpu,'-',label='Plasma Freqeuncy ($1.16e17*(1+h/(k*h0))^{-k+1}$) h0=400 km')
plt.plot(h,fp25,'-',label='Plasma Freqeuncy ($1.e6+1.16e17*e^{(-h/h0)}$) h0=110 km')
#plt.axvline(x=25,color='k',linewidth=4,label='25 Mm')
#plt.axhline(y=1.0,color='r',linewidth=4, label='1 GHz')
#plt.axhline(y=0.4,color='g',linewidth=4, label='400 MHz')
plt.xlabel('Coronal Height (Mm)');plt.ylabel('Freqeuncy (GHz)'),plt.legend()
plt.yscale('log');plt.xlim([0,50]);plt.ylim([0.01,100])
plt.show()

plt.imshow(densdata[3],aspect='auto',origin=0,extent=[xc-0.6*83,xc+83,yc-0.6*83,yc+0.6*83])
#plt.plot(xcrmax[0],ycrmax[0],'o',color='k')
plt.title('2.0-3.0 MK')
plt.colorbar(label='EM (cm$^{-5}$)')
plt.xlabel('Solar X (arcsec)');plt.ylabel('Solar Y (arcsec)')
plt.show()


h=np.arange(200)*1400/1.e3 # in Mm
h1=np.arange(200)*500/1.e3 # in Mm
sh=np.linspace(10,500,100)
ux=np.linspace(0,50,50)
dem_arr=[0]*100;ne=[0]*100;ne_dem=[0]*100;wpe_dem=[0]*100;omegap_dem=[0]*100
tem=2.6e27

for i in range(100):
    dem_arr[i]=[0]*200;ne[i]=[0]*200
    #ne[i]=1.16e17*(1+(h*1000/(10*sh[i])))**(-10+1)
    ne[i]=1.16e17*(np.exp(-1*h1*1000/sh[i]))+9.5e8*np.exp(-1*h1*1000/31.e3)#+4.2e4 *10**(4.32/((695700.+h*1000)/695700.))
    #ne[i]=4.2e4 *10**(4.32/((695700.+h*1000)/695700.))+1.16e17*(np.exp(-1*h*1000/sh[i]))
    #ne[i]=1.16e17*(np.exp(-1*h*1000/sh[i]))+1.e6
    #ne[i]=1.16e17*(np.exp(-1*h*1000/sh[i]))+(3.09e8*(1-0.5*np.sin(th))/((695700.+h*1000)/695700.)**16)+(1.58e8*(1-0.95*np.sin(th))/((695700.+h*1000)/695700.)**6)+(0.0251e8*(1-np.sqrt(np.sin(th)))/((695700.+h*1000)/695700.)**2.5)
    for j in range(200):
        dem_arr[i][j]=np.sum(ne[i][j:]*ne[i][j:])*(h[1]-h[0])*1.e8
    idx=ut.find_nearest(dem_arr[i],tem)[0]
    ne_dem[i]=ne[i][idx];wpe_dem[i]=9000*np.sqrt(ne_dem[i])/1.e9
    omegap_dem[i]=omegap[idx]

dem_arr=np.array(dem_arr);ne=np.array(ne);ne_dem=np.array(ne_dem);wpe_dem=np.array(wpe_dem);omegap_dem=np.array(omegap_dem)
#dem_arr[np.where((dem_arr<2.e27) & (dem_arr>1.e27))]=np.nan

s=2;sf=np.math.factorial(s);e=1.e-19;me=9.1e-31;c=3.e8;nu=1.e9;mu=1;th=45
lb=0.10e6;ratio=wpe_dem/omegap_dem
def optical_depth(s,ratio):
    Ftho=(np.sin(th*np.pi/180)*np.sin(th*np.pi/180)*(np.sin(th*np.pi/180)**2+4*np.cos(th*np.pi/180)**2 +np.sqrt(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2))**2)/(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2+np.sin(th*np.pi/180)*np.sin(th*np.pi/180)*np.sqrt(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2))
    Fthx=(np.sin(th*np.pi/180)*np.sin(th*np.pi/180)*(np.sin(th*np.pi/180)**2+4*np.cos(th*np.pi/180)**2 -np.sqrt(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2))**2)/(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2-np.sin(th*np.pi/180)*np.sin(th*np.pi/180)*np.sqrt(np.sin(th*np.pi/180)**4 + 16*np.cos(th*np.pi/180)**2))
    tauo=np.pi*e*e*ne[28]*1.e9*lb*s*s*s**(2*(s-1))*np.sin(th*np.pi/180)*Ftho/(2*me*c*nu*sf*(2*mu)**(s-1))
    taux=np.pi*e*e*ne[28]*1.e9*lb*s*s*s**(2*(s-1))*np.sin(th*np.pi/180)*Fthx/(2*me*c*nu*sf*(2*mu)**(s-1))
    return tauo,taux

tauo2,taux2=optical_depth(2,ratio)
tauo3,taux3=optical_depth(3,ratio)
tauo4,taux4=optical_depth(5,ratio)


plt.imshow(np.log10(dem_arr),origin=0,aspect='auto',extent=[h[0],h[-1],sh[0],sh[-1]])
plt.colorbar(label='Log(EM(cm$^{-5}$))')
plt.ylabel('Scale Height (km)');plt.xlabel('Source Height (Mm)')
plt.show()

plt.plot(sh,wpe_dem/omegap_dem,'o-')
plt.ylabel('$\omega_p/\Omega_e$');plt.xlabel('Scale Height (km)')
plt.show()

plt.plot(sh,dem_arr[:,10],'o-')
plt.axhline(y=tem,color='k')
plt.xscale('log');plt.yscale('log')
plt.xlabel('Scale Height (km)');plt.ylabel('TEM ($cm^{-5}$)')
plt.show()


f,ax=plt.subplots(2,2)
ax[0,0].plot(h,babs[:,160,172],'o-',label='Magnetic Field')
ax[0,0].set_ylabel('|B| (Gauss)')
ax[0,1].plot(h,ne[28],'o-',label='Electron density')
ax[0,1].set_ylabel(' n$_e$ ($cm^{-3}$)')
ax[1,0].plot(h,omegap,'o-',label='Gyroresonce Frequency')
ax[1,0].set_ylabel('Frequency (GHz)')
ax[1,1].plot(h,9000*np.sqrt(ne[28])/1.e9,'-',label='Plasma Freqeuncy h0=150 km')
ax[1,1].plot(h,omegap,'o-',label='Gyroresonce Frequency')
ax[1,1].set_ylabel('Frequency (GHz)')
ax[1,0].set_xlabel('Coronal Height (Mm)');ax[1,1].set_xlabel('Coronal Height (Mm)')
ax[0,0].legend()
ax[0,1].legend()
ax[1,0].legend()
ax[1,1].legend()
plt.show()

f,ax=plt.subplots(1,1)
ax.plot(h,9000*np.sqrt(ne[28])/1.e9,'-',label='Plasma Freqeuncy h0=150 km')
ax.plot(h,omegap,'o-',label='Gyroresonce Frequency')
ax1=ax.twinx();ax.legend()
ax1.plot(h,tauo3,'o--',color='k',label='O-mode, s=3')
ax1.plot(h,taux3,'s--',color='k',label='X-mode, s=3')
ax1.axhline(y=1.0,color='k')
ax1.legend(loc=2);ax.set_xscale('log');ax.set_yscale('log');ax1.set_xscale('log');ax1.set_yscale('log')
ax.set_xlabel('Coronal Heights (Mm)');ax.set_ylabel('Freqeuncy (GHz)');ax1.set_ylabel('$\\tau$')
plt.show()

vA=babs[:,160,172]/np.sqrt(4*np.pi*ne[28]*9.1e-28)/1.e5
f,ax=plt.subplots(1,1)
ax.plot(h,9000*np.sqrt(ne[28])/1.e9,'-',label='Plasma Freqeuncy h0=150 km')
ax.plot(h,omegap,'o-',label='Gyroresonce Frequency')
ax1=ax.twinx();ax.legend()
ax1.plot(h,vA,'-',color='k')
ax.set_xlabel('Coronal Heights (Mm)');ax.set_ylabel('Freqeuncy (GHz)');ax1.set_ylabel('V${_A}$ km/s')
ax1.legend(loc=2);ax.set_xscale('log');ax.set_yscale('log');ax1.set_xscale('log');ax1.set_yscale('log')
plt.show()

plt.plot(h,omegap,'o-',label='Gyroresonce Frequency')
plt.plot(h1,9000*np.sqrt(ne[28])/1.e9,'-',label='Plasma Frequency h0=150 km')
#plt.axhline(y=1.4,color='k');plt.axhline(y=0.99,color='k')
#plt.axvline(x=20.0,color='k');plt.axvline(x=26.0,color='k')
plt.axvline(x=hh2_0[0],color='k');plt.axvline(x=hh2_0[16],color='k');plt.axvline(x=hh2_0[30],color='k')
plt.xlabel('Coronal Height (Mm)');plt.ylabel('Frequency (GHz)'),plt.legend()
plt.yscale('log');plt.xscale('log');plt.xlim([0,100]);plt.ylim([0.01,100])
plt.show()


f,ax=plt.subplots(1,1)
ax.plot(h,ne[28],'o-',label='Inside Density')
ax.plot(h,ne[28]*(2/2.6),'o-',label='Outside Density')
ax.legend();ax.set_ylim(1.e8,1.e9);ax.set_xscale('log')
ax.set_xlim(1.e1,1.e2)
ax.set_ylabel('Density ($cm^{-3}$)');ax.set_xlabel('Coronal Height (Mm)')
plt.show()
