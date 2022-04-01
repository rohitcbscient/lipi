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
hmieui_data, footprint = reproject_interp(hmimap, euiwcs, shape_out=(2048,2048));hmieui_map=Map(hmieui_data,euimap.meta)
aiaeui_data, footprint = reproject_interp(aiamap, euiwcs, shape_out=(2048,2048));aiaeui_map=Map(aiaeui_data,euimap.meta)
aiahmi_data, footprint = reproject_interp(aiamap, hmiwcs, shape_out=(4096,4096));aiahmi_map=Map(aiahmi_data,hmimap.meta)
euihmi_data, footprint = reproject_interp(euimap, hmiwcs, shape_out=(4096,4096));euihmi_map=Map(euihmi_data,hmimap.meta)
euiaia_data, footprint = reproject_interp(euimap, aiawcs, shape_out=(4096,4096));euiaia_map=Map(euiaia_data,aiamap.meta)

stereo=0
if stereo:
    sun_to_stereo = aiamap.observer_coordinate.transform_to('hcrs')
    stereo_to_sun = SkyCoord(-sun_to_stereo.spherical, obstime=sun_to_stereo.obstime,frame='hcrs')
    tbl_crds = SkyCoord(ra=result[0]['RA_ICRS'],dec=result[0]['DE_ICRS'],distance=Distance(parallax=u.Quantity(result[0]['Plx'])),pm_ra_cosdec=result[0]['pmRA'],pm_dec=result[0]['pmDE'],radial_velocity=result[0]['RV'],frame='icrs',obstime=Time(result[0]['Epoch'], format='jyear'))
    tbl_crds = tbl_crds.apply_space_motion(new_obstime=cor2.date)
    mars = get_body_heliographic_stonyhurst('mars', cor2.date, observer=cor2.observer_coordinate)

#file_name_fr='/data/Dropbox/proposals/vla_stix_psp_2022/fieldmodel1.csv'
file_name_fr='/data/Dropbox/proposals/vla_stix_psp_2022/20210322_fieldlines.csv'
x,y,z,bx,by,bz=Bextrap.get_fieldlines(file_name_fr)
ff='/sdata/20210322/hmi.M_720s.20210322_202236.E115N20CR.CEA.NAS.sav'
bsmap,bcarrmap=Bextrap.get_gxs_sav2hpp(ff,'2021-03-22T20:22:36.00')
xs,ys,zs,bsx,bsy,bsz,bs_hp,bs_hp_pix,bs_carr,bs_carr_pix=Bextrap.transform_fieldlines(x,y,z,bx,by,bz,'2021/03/22T20:22:36',bsmap[0].wcs,bcarrmap[0])
xe,ye,ze,bxe,bye,bze,be_hp,be_hp_pix,be_carr,be_carr_pix=Bextrap.transform_fieldlines_2SO(x,y,z,bx,by,bz,'2021/03/22T20:22:36',euiwcs,bcarrmap[0],euimap)
xa,ya,za,bxa,bya,bza,ba_hp,ba_hp_pix,ba_carr,ba_carr_pix=Bextrap.transform_fieldlines_2SO(x,y,z,bx,by,bz,'2021/03/22T20:22:36',aiawcs,bcarrmap[0],aiamap)


#####
xlaia=800;xraia=1200;ylaia=200;yraia=600

plot_fields_euvi=1
if(plot_fields_euvi):
    plt.close()
    fig = plt.figure(figsize=(16, 8))
    ax0 = fig.add_subplot(231,projection=euimap);ax1 = fig.add_subplot(232,projection=bcarrmap[0])
    ax2 = fig.add_subplot(233,projection=euimap);ax3 = fig.add_subplot(234,projection=bcarrmap[0])
    ax4 = fig.add_subplot(235,projection=aiamap);ax5 = fig.add_subplot(236,projection=aiamap)
    p0=euimap.plot(axes=ax0,aspect='auto')
    p2=euimap.plot(axes=ax2,aspect='auto')
    p1=bcarrmap[0].plot(axes=ax1,aspect='auto')
    p3=bcarrmap[0].plot(axes=ax3,aspect='auto')
    p4=aiamap.plot(axes=ax4,aspect='auto')
    p5=aiamap.plot(axes=ax5,aspect='auto')
    ####
    hp_lon = be_hp.Tx.value * u.arcsec
    hp_lat = be_hp.Ty.value* u.arcsec
    hpa_lon = ba_hp.Tx.value * u.arcsec
    hpa_lat = ba_hp.Ty.value* u.arcsec
    seeds0 = SkyCoord(hp_lon.ravel(), hp_lat.ravel(),frame=euimap.coordinate_frame)
    ax0.plot_coord(seeds0, color='tab:red', marker='s', markersize=0.1,alpha=0.5,linestyle='None')
    seeds1 = SkyCoord(be_carr.lon.ravel(), be_carr.lat.ravel(),frame=bcarrmap[0].coordinate_frame)
    ax1.plot_coord(seeds1, color='tab:red', marker='s', markersize=0.1,alpha=0.5,linestyle='None')
    seeds2 = SkyCoord(hpa_lon.ravel(), hpa_lat.ravel(),frame=aiamap.coordinate_frame)
    ax4.plot_coord(seeds2, color='tab:red', marker='s', markersize=0.1,alpha=0.5,linestyle='None')
    ####
    bl=SkyCoord(900*u.arcsec,200*u.arcsec,frame=euimap.coordinate_frame)
    tr=SkyCoord(1300*u.arcsec,600*u.arcsec,frame=euimap.coordinate_frame)
    pixbl=euimap.world_to_pixel(bl);pixtr=euimap.world_to_pixel(tr) 
    blx=euimap.world_to_pixel(bl)[0].value;bly=euimap.world_to_pixel(bl)[1].value
    trx=euimap.world_to_pixel(tr)[0].value;rty=euimap.world_to_pixel(tr)[1].value
    ax0.set_xlim([blx,trx]);ax0.set_ylim([bly,rty])
    ax2.set_xlim([blx,trx]);ax2.set_ylim([bly,rty])
    ####
    bl=SkyCoord(-850*u.arcsec,300*u.arcsec,frame=aiamap.coordinate_frame)
    tr=SkyCoord(-550*u.arcsec,600*u.arcsec,frame=aiamap.coordinate_frame)
    pixbl=aiamap.world_to_pixel(bl);pixtr=aiamap.world_to_pixel(tr) 
    blx=aiamap.world_to_pixel(bl)[0].value;bly=aiamap.world_to_pixel(bl)[1].value
    trx=aiamap.world_to_pixel(tr)[0].value;rty=aiamap.world_to_pixel(tr)[1].value
    ax4.set_xlim([blx,trx]);ax4.set_ylim([bly,rty])
    ax5.set_xlim([blx,trx]);ax5.set_ylim([bly,rty])
    ####
    plt.show()

plot_one=1
if plot_one:
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=aiamap)
    p0=aiamap.plot(axes=ax,aspect='auto')
    hp_lon = ba_hp.Tx.value * u.arcsec
    hp_lat = ba_hp.Ty.value* u.arcsec
    seeds0 = SkyCoord(hp_lon.ravel(), hp_lat.ravel(),frame=aiamap.coordinate_frame)
    ax.plot_coord(seeds0, color='tab:red', marker='s', markersize=0.1,alpha=0.5,linestyle='None')
    bl=SkyCoord(-850*u.arcsec,300*u.arcsec,frame=aiamap.coordinate_frame)
    tr=SkyCoord(-550*u.arcsec,600*u.arcsec,frame=aiamap.coordinate_frame)
    pixbl=aiamap.world_to_pixel(bl);pixtr=aiamap.world_to_pixel(tr) 
    blx=aiamap.world_to_pixel(bl)[0].value;bly=aiamap.world_to_pixel(bl)[1].value
    trx=aiamap.world_to_pixel(tr)[0].value;rty=aiamap.world_to_pixel(tr)[1].value
    ax.set_xlim([blx,trx]);ax.set_ylim([bly,rty])
    plt.show()

plot_one=1
if plot_one:
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=euimap)
    p0=euimap.plot(axes=ax,aspect='auto')
    hp_lon = be_hp.Tx.value * u.arcsec
    hp_lat = be_hp.Ty.value* u.arcsec
    seeds0 = SkyCoord(hp_lon.ravel(), hp_lat.ravel(),frame=euimap.coordinate_frame)
    ax.plot_coord(seeds0, color='tab:red', marker='s', markersize=0.1,alpha=0.5,linestyle='None')
    bl=SkyCoord(900*u.arcsec,200*u.arcsec,frame=euimap.coordinate_frame)
    tr=SkyCoord(1300*u.arcsec,600*u.arcsec,frame=euimap.coordinate_frame)
    pixbl=euimap.world_to_pixel(bl);pixtr=euimap.world_to_pixel(tr) 
    blx=euimap.world_to_pixel(bl)[0].value;bly=euimap.world_to_pixel(bl)[1].value
    trx=euimap.world_to_pixel(tr)[0].value;rty=euimap.world_to_pixel(tr)[1].value
    ax.set_xlim([blx,trx]);ax.set_ylim([bly,rty])
    plt.show()

plot_fields_hp=1
if(plot_fields_hp):
    plt.close()
    fig = plt.figure(figsize=(16, 8))
    ax0 = fig.add_subplot(121,projection=bsmap[0]);ax1 = fig.add_subplot(122,projection=bcarrmap[0])
    p0=bsmap[0].plot(axes=ax0,aspect='auto')
    p1=bcarrmap[0].plot(axes=ax1,aspect='auto')
    hp_lon = b_hp.Tx.value * u.arcsec
    hp_lat = b_hp.Ty.value* u.arcsec
    seeds0 = SkyCoord(hp_lon.ravel(), hp_lat.ravel(),frame=bsmap[0].coordinate_frame)
    ax0.plot_coord(seeds0, color='tab:olive', marker='s', markersize=5)
    #ca_lon = b_proj[0] * u.deg
    #ca_lat = b_proj[1] * u.deg
    seeds1 = SkyCoord(b_carr.lon.ravel(), b_carr.lat.ravel(),frame=bcarrmap[0].coordinate_frame)
    ax1.plot_coord(seeds1, color='tab:olive', marker='s', markersize=5)
    plt.show()

hp_lon = b_hp.Tx.value * u.arcsec
hp_lat = b_hp.Ty.value* u.arcsec
fig = plt.figure()
ax = plt.subplot(projection=bsmap[0])
bsmap[0].plot(axes=ax)
seeds = SkyCoord(hp_lon.ravel(), hp_lat.ravel(),frame=bsmap[0].coordinate_frame)
ax.plot_coord(seeds, color='white', marker='o', linewidth=0)

