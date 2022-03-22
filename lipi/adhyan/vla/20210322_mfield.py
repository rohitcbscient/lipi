import numpy as np
from astropy.io import fits
from surya.utils import Bextrap
from astropy import units as u
import matplotlib.pyplot as plt
import pickle
from astropy.coordinates import SkyCoord
import sys
from sunpy.map import Map
from reproject import reproject_interp

euifile='/data/Dropbox/20210322/EUI/solo_L2_eui-fsi174-image_20210322T203913290_V02.fits'
hmifile='/data/Dropbox/20210322/EUI/hmi.m_45s.2021.03.22_20_39_00_TAI.magnetogram.fits'
aiafile='/media/rohit/VLA/20210322_EUV/aia.lev1.171A_2021-03-22T20_30_33.35Z.image_lev1.fits'
euimap=Map(euifile);hmimap=Map(hmifile);aiamap=Map(aiafile)
euiwcs=euimap.wcs;hmiwcs=hmimap.wcs;aiawcs=aiamap.wcs
hmieui_map, footprint = reproject_interp(hmimap, euiwcs, shape_out=(2048,2048))
aiaeui_map, footprint = reproject_interp(aiamap, euiwcs, shape_out=(2048,2048))
aiahmi_map, footprint = reproject_interp(aiamap, hmiwcs, shape_out=(2048,2048))
euihmi_map, footprint = reproject_interp(euimap, hmiwcs, shape_out=(2048,2048))
euiaia_map, footprint = reproject_interp(euimap, aiawcs, shape_out=(2048,2048))

sun_to_stereo = aiamap.observer_coordinate.transform_to('hcrs')
stereo_to_sun = SkyCoord(-sun_to_stereo.spherical, obstime=sun_to_stereo.obstime,frame='hcrs')
tbl_crds = SkyCoord(ra=result[0]['RA_ICRS'],dec=result[0]['DE_ICRS'],distance=Distance(parallax=u.Quantity(result[0]['Plx'])),pm_ra_cosdec=result[0]['pmRA'],pm_dec=result[0]['pmDE'],radial_velocity=result[0]['RV'],frame='icrs',obstime=Time(result[0]['Epoch'], format='jyear'))
tbl_crds = tbl_crds.apply_space_motion(new_obstime=cor2.date)
mars = get_body_heliographic_stonyhurst('mars', cor2.date, observer=cor2.observer_coordinate)

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

