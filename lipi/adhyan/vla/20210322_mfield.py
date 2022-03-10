import numpy as np
from astropy.io import fits
from surya.utils import Bextrap
from astropy import units as u
import matplotlib.pyplot as plt
import pickle
from astropy.coordinates import SkyCoord
import sys

euifile='/data/Dropbox/20210321/EUI/solo_L2_eui-fsi174-image_20210322T203913290_V02.fits'
euimap=Map(euifile)
euiwcs=euimap.wcs



sys.exit()

file_name_fr="/media/rohit/VLA/paraview/radio_loop_new.csv"
file_name_ls="/media/rohit/VLA/paraview/large_scale.csv"
file_name_euv="/media/rohit/VLA/paraview/EUV_loop2.csv"
file_name_radio='/media/rohit/VLA/paraview/EUV_loop_radio.csv'
x,y,z,bx,by,bz=Bextrap.get_fieldlines(file_name_fr)
ff='/media/rohit/VLA/20160409_EUV/hmi.M_720s.20160409_183417.E18N10CR.CEA.NAS.sav'
bsmap,bcarrmap=Bextrap.get_gxs_sav2hpp(ff,'2016-04-09T18:34:17.00')
x,y,z,bx,by,bz,b_hp,b_hp_pix,b_carr,b_carr_pix=Bextrap.transform_fieldlines(x,y,z,bx,by,bz,'2016/04/09T18:45:00',bsmap[0].wcs,bcarrmap[0])

#####
xcrmax,ycrmax,xcr90,ycr90,maxTbr,Tbr_r1,Tbr_r2,eTbr=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r.p','rb'),encoding = 'latin1')
qsx,qsy,qsxcr90,qsycr90,qsmaxTbr,qsTbr_r1,qsTbr_r2,qsarear50,qstimevla=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r_qs.p','rb'),encoding = 'latin1')
qsmaxTbr=np.array(qsmaxTbr)
x0,y0,sigx0,sigy0,rot0,tb0,x1,y1,sigx1,sigy1,rot1,tb1=np.load('/media/rohit/VLA/20160409/blob/all_params.npz.npy')
corr_t=np.load('/media/rohit/VLA/20160409/correlation_t.npz.npy')
xlaia=-948;xraia=-648;ylaia=70;yraia=370

plot_fields=1
if(plot_fields):
    plt.close()
    fig = plt.figure(figsize=(16, 8))
    ax0 = fig.add_subplot(121,projection=bsmap[0]);ax1 = fig.add_subplot(122,projection=bcarrmap[0])
    p0=bsmap[0].plot(axes=ax0,aspect='auto')
    p1=bcarrmap[0].plot(axes=ax1,aspect='auto')
    hp_lon = b_hp.Tx.value * u.arcsec
    hp_lat = b_hp.Ty.value* u.arcsec
    seeds0 = SkyCoord(hp_lon.ravel(), hp_lat.ravel(),frame=bsmap[0].coordinate_frame)
    ax0.plot_coord(seeds0, color='tab:olive', marker='o', markersize=10)
    #ca_lon = b_proj[0] * u.deg
    #ca_lat = b_proj[1] * u.deg
    seeds1 = SkyCoord(b_carr.lon.ravel(), b_carr.lat.ravel(),frame=bcarrmap[0].coordinate_frame)
    ax1.plot_coord(seeds1, color='tab:olive', marker='o', markersize=10)
    #plt.plot(xcr90[15][1000:1400],ycr90[15][1000:1400],'o',alpha=0.2)
    #plt.plot(xcr90[10][1000:1400],ycr90[10][1000:1400],'o',alpha=0.2)
    #plt.plot(xcr90[5][1000:1400],ycr90[5][1000:1400],'o',alpha=0.2)
    #plt.plot(xcr90[1][1000:1400],ycr90[1][1000:1400],'o',alpha=0.2)
    #plt.plot(xcr90[12][500:600],ycr90[12][500:600],'o',alpha=0.2)
    #plt.plot(xcr90[10][500:600],ycr90[10][500:600],'o',alpha=0.2)
    #plt.plot(xcr90[5][500:600],ycr90[5][500:600],'o',alpha=0.2)
    #plt.plot(xcr90[1][500:600],ycr90[1][500:600],'o',alpha=0.2)
    #ax0.set_xlim(-860,-680);ax0.set_ylim(180,360)
    plt.show()

hp_lon = b_hp.Tx.value * u.arcsec
hp_lat = b_hp.Ty.value* u.arcsec
fig = plt.figure()
ax = plt.subplot(projection=bsmap[0])
bsmap[0].plot(axes=ax)
seeds = SkyCoord(hp_lon.ravel(), hp_lat.ravel(),frame=bsmap[0].coordinate_frame)
ax.plot_coord(seeds, color='white', marker='o', linewidth=0)

