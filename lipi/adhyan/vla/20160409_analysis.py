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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import os
from scipy.io import readsav
from astropy.time import Time
import sunpy
import sunpy.map as smap
from sunpy import sun
from dateutil import parser
from scipy import stats
import numpy as np

def idl2sunpy_hmi(mapsav):
    mapstruc = readsav(mapsav)
    keys = list(mapstruc.keys())
    mp = mapstruc[keys[0]][0]
    dat = np.nan_to_num(mp[0]) + 1e-6 # just to get rid of zeros
    header = {}
    header['EXPTIME'] = mp['dur']
    header['CRVAL1'] = mp[1]
    header['CRVAL2'] = mp[2]
    header['CDELT1'] = mp[3]
    header['CDELT2'] = mp[4]
    header['DATE-OBS'] = Time(parser.parse(mp[5])).isot
    header['CUNIT1'] = 'arcsec'
    header['CUNIT2'] = 'arcsec'
    header['CTYPE1'] = 'HPLN-TAN'
    header['CTYPE2'] = 'HPLT-TAN'
    header['CRPIX2'], header['CRPIX1'] = (np.array(dat.shape) + 1.0) / 2 + 0.5
    header['TELESCOP'] = 'SDO/HMI'
    header['P_ANGLE'] = 0.0
    header['NAXIS'] = 2
    header['NAXIS2'], header['NAXIS1'] = dat.shape
    header['HGLN_OBS'] = 0.
    #header['HGLT_OBS'] = sunpy.coordinates.get_sun_B0(header['DATE-OBS']).value
    header['HGLT_OBS'] = sunpy.coordinates.sun.B0(header['DATE-OBS']).value
    header['RSUN_REF'] = sunpy.sun.constants.radius.value
    header['RSUN_OBS'] = 958.11#sunpy.coordinates.sun.angular_radius(header['DATE-OBS']).value
    header['DSUN_OBS'] = sunpy.coordinates.sun.earth_distance(header['DATE-OBS']).to(u.meter).value
    smp = smap.Map(dat, header)
    return smp


def idl2sunpy_sdo(mapsav,wave,inst):
    mapstruc = readsav(mapsav)
    keys = list(mapstruc.keys())
    mp = mapstruc[keys[0]][0]
    dat = np.nan_to_num(mp[0]) + 1e-6 # just to get rid of zeros
    header = {}
    header['EXPTIME'] = mp['dur']
    header['CRVAL1'] = mp[1]
    header['CRVAL2'] = mp[2]
    header['CDELT1'] = mp[3]
    header['CDELT2'] = mp[4]
    header['DATE-OBS'] = Time(parser.parse(mp[5])).isot
    header['CUNIT1'] = 'arcsec'
    header['CUNIT2'] = 'arcsec'
    header['CTYPE1'] = 'HPLN-TAN'
    header['CTYPE2'] = 'HPLT-TAN'
    header['CRPIX2'], header['CRPIX1'] = (np.array(dat.shape) + 1.0) / 2 + 0.5
    header['TELESCOP'] = 'SDO'
    header['INSTRUME'] = inst
    header['WAVELNTH'] = wave
    header['waveunit'] = 'angstrom'
    header['P_ANGLE'] = 0.0
    header['NAXIS'] = 2
    header['NAXIS2'], header['NAXIS1'] = dat.shape
    header['HGLN_OBS'] = 0.
    #header['HGLT_OBS'] = sunpy.coordinates.get_sun_B0(header['DATE-OBS']).value
    header['HGLT_OBS'] = sunpy.coordinates.sun.B0(header['DATE-OBS']).value
    header['RSUN_REF'] = sunpy.sun.constants.radius.value
    header['RSUN_OBS'] = 958.11#sunpy.coordinates.sun.angular_radius(header['DATE-OBS']).value
    header['DSUN_OBS'] = sunpy.coordinates.sun.earth_distance(header['DATE-OBS']).to(u.meter).value
    smp = smap.Map(dat, header)
    return smp


def produce_tstring(f):
    mapp=Map(f);d=mapp.date
    date1=mapp.date.ymdhms#.date
    #hhmmss=' '+str(date.hour)+':'+str(date.minute)+':'+str(date.second)+'.'+str(date.microsecond/1.e6).split('.')[1]
    hhmmss=' '+str(date1[3])+':'+str(date1[4])+':'+str(date1[5])
    sec=ut.hms2sec_c(hhmmss)
    print(d, date1,sec,hhmmss)
    return sec

def get_sunpy_maps_rot(f,wave,inst):
    print('Reading...'+f[0])
    n=len(f);maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(n):
        if(inst=='HMI'):
            maplist[i]=idl2sunpy_hmi(f[i])
        if(inst=='AIA'):
            maplist[i]=idl2sunpy_sdo(f[i],wave,inst)
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time

def get_sunpy_maps(f):
    print('Reading...'+f[0])
    n=len(f);maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(n):
        maplist[i]=Map(f[i])
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time

def get_sunpy_basediff_maps(f,wave,inst):
    print('Reading...'+f[0])
    n=len(f);maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(n):
        if(inst=='AIA'):
            aa=idl2sunpy_sdo(f[i],wave,inst)
            ab=idl2sunpy_sdo(f[0],wave,inst)
        if(inst=='HMI'):
            aa=idl2sunpy_hmi(f[i])
            ab=idl2sunpy_hmi(f[0])
        h,d=aa.meta,aa.data
        hb,db=ab.meta,ab.data
        if(d.shape[0]==db.shape[0]):
            d=d-db
        else:
            if(d.shape[0]<db.shape[0]):
                d=d-db[0:d.shape[0],0:d.shape[1]]
            else:
                d=d[0:db.shape[0],0:db.shape[1]]-db
        maplist[i]=Map(d,h)
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time


def get_sunpy_rundiff_maps(f,wave,inst):
    print('Reading...'+f[0])
    n=len(f)-1;maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(1,n):
        if(inst=='AIA'):
            ab=idl2sunpy_sdo(f[i-1],wave,inst);hb,db=ab.meta,ab.data
            aa=idl2sunpy_sdo(f[i],wave,inst);h,d=aa.meta,aa.data
        if(inst=='HMI'):
            ab=idl2sunpy_hmi(f[i-1]);hb,db=ab.meta,ab.data
            aa=idl2sunpy_hmi(f[i]);h,d=aa.meta,aa.data
        if(d.shape[0]==db.shape[0]):
            d=d-db
        else:
            if(d.shape[0]<db.shape[0]):
                d=d-db[0:d.shape[0],0:d.shape[1]]
            else:
                d=d[0:db.shape[0],0:db.shape[1]]-db
        maplist[i]=Map(d,h)
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time

#############################################################
def get_hmi_basediff_maps(f):
    print('Reading...'+f[0])
    n=len(f);maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(n):
        aa=Map(f[i])
        ab=Map(f[0])
        h,d=aa.meta,aa.data
        hb,db=ab.meta,ab.data
        if(d.shape[0]==db.shape[0]):
            d=d-db
        else:
            if(d.shape[0]<db.shape[0]):
                d=d-db[0:d.shape[0],0:d.shape[1]]
            else:
                d=d[0:db.shape[0],0:db.shape[1]]-db
        maplist[i]=Map(d[::-1,::-1],aa.meta)
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time
def get_hmi_rundiff_maps(f):
    print('Reading...'+f[0])
    n=len(f)-1;maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(1,n):
        ab=Map(f[i-1])
        hb,db=ab.meta,ab.data
        aa=Map(f[i])
        h,d=aa.meta,aa.data
        if(d.shape[0]==db.shape[0]):
            d=d-db
        else:
            if(d.shape[0]<db.shape[0]):
                d=d-db[0:d.shape[0],0:d.shape[1]]
            else:
                d=d[0:db.shape[0],0:db.shape[1]]-db
        maplist[i]=Map(d[::-1,::-1],h)
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time


def get_hmi_submap(f,xbl,ybl,xtr,ytr):
    print('Reading...'+f[0])
    #xcen=(xbl+xtr)*0.5;ycen=(ybl+ytr)*0.5
    n=len(f);maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(n):
        g=fits.open(f[i],ignore_missing_simple=True)
        map_1=Map(f[i])
        map_=Map(map_1.data[::-1,::-1],map_1.meta)
        #cent_pix=map_.world_to_pixel(SkyCoord(xcen*u.arcsec, ycen*u.arcsec, frame=map_.coordinate_frame)) 
        bl = SkyCoord(xbl*u.arcsec, ybl*u.arcsec, frame=map_.coordinate_frame)
        tr = SkyCoord(xtr*u.arcsec, ytr*u.arcsec, frame=map_.coordinate_frame)
        maplist[i]=map_.submap(bl,tr)
        maplist[i]=map_
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time

def get_evla_submap(f,xbl,ybl,xtr,ytr):
    print('Reading...'+f[0])
    #xcen=(xbl+xtr)*0.5;ycen=(ybl+ytr)*0.5
    n=len(f);maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(n):
        g=fits.open(f[i])
        #g[0].header['CRVAL1']=0;g[0].header['CRVAL2']=0
        #g[0].header['CRVAL1']=-770;g[0].header['CRVAL2']=220
        #g[0].header['CRVAL1']=-731.12;g[0].header['CRVAL2']=243.5
        #g[0].data[np.isnan(g[0].data)]=0
        map_=Map(g[0].data[0][0],g[0].header)
        #cent_pix=map_.world_to_pixel(SkyCoord(xcen*u.arcsec, ycen*u.arcsec, frame=map_.coordinate_frame)) 
        bl = SkyCoord(xbl*u.arcsec, ybl*u.arcsec, frame=map_.coordinate_frame)
        tr = SkyCoord(xtr*u.arcsec, ytr*u.arcsec, frame=map_.coordinate_frame)
        #maplist[i]=map_.submap(bl,tr)
        maplist[i]=map_
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i].date)
    time=np.array(time)
    return maplist,datalist,time

def change_pointing(f,xc,yc):
    g=fits.open(f)
    g[0].header['CRVAL1']=xc;g[0].header['CRVAL2']=yc
    g.writeto('/media/rohit/VLA/20160409/sun_L_20160409T1844-1846UT.50ms.cal.pol.LL.spw.0_16-31.time.18:45:10.-18:45:10.05.fits')
    g.close()

sys.exit()

dump_submaps=1
if(dump_submaps):
    # NOT USED
    list1600=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/1600/*rot.sav'))
    list1700=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/1700/*rot.sav'))
    list335=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/335/*rot.sav'))
    list304=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/304/*rot.sav'))
    list131=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/131/*rot.sav'))
    list193=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/193/*rot.sav'))
    list171=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/171/*rot.sav'))
    list211=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/211/*rot.sav'))
    list94=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/94/*rot.sav'))
    map94,data94,time94=get_sunpy_maps_rot(list94,'94','AIA')
    map131,data131,time131=get_sunpy_maps_rot(list131,'131','AIA')
    map171,data171,time171=get_sunpy_maps_rot(list171,'171','AIA')
    map193,data193,time193=get_sunpy_maps_rot(list193,'193','AIA')
    map211,data211,time211=get_sunpy_maps_rot(list211,'211','AIA')
    map304,data304,time304=get_sunpy_maps_rot(list304,'304','AIA')
    map335,data335,time335=get_sunpy_maps_rot(list335,'335','AIA')
    map1600,data1600,time1600=get_sunpy_maps_rot(list1600,'1600','AIA')
    map1700,data1700,time1700=get_sunpy_maps_rot(list1700,'1700','AIA')
    ################
    #map94,data94,time94=get_sunpy_maps(list94)
    #map131,data131,time131=get_sunpy_maps(list131)
    #map335,data335,time335=get_sunpy_maps(list335)
    #map171,data171,time171=get_sunpy_maps(list171)
    #map211,data211,time211=get_sunpy_maps(list211)
    #map193,data193,time193=get_sunpy_maps(list193)
    #map1600,data1600,time1600=get_sunpy_maps(list1600)
    #map1700,data1700,time1700=get_sunpy_maps(list1700)
    allaiamaps={};allaiamaps['aia94']={'map94':map94,'data94':data94,'time94':time94}
    allaiamaps['aia131']={'map131':map131,'data131':data131,'time131':time131}
    allaiamaps['aia171']={'map171':map171,'data171':data171,'time171':time171}
    allaiamaps['aia193']={'map193':map193,'data193':data193,'time193':time193}
    allaiamaps['aia211']={'map211':map211,'data211':data211,'time211':time211}
    allaiamaps['aia304']={'map304':map304,'data304':data304,'time304':time304}
    allaiamaps['aia335']={'map335':map335,'data335':data335,'time335':time335}
    allaiamaps['aia1600']={'map1600':map1600,'data1600':data1600,'time1600':time1600}
    allaiamaps['aia1700']={'map1700':map1700,'data1700':data1700,'time1700':time1700}
    print("Writing..")
    pickle.dump(allaiamaps,open('/media/rohit/VLA/20160409/20160409_submap_aia_50ms.p','wb'))
    sys.exit()
    ##########################
    mapd94,datad94,timed94=get_sunpy_rundiff_maps(list94,'94','AIA')
    mapd131,datad131,timed131=get_sunpy_rundiff_maps(list131,'131','AIA')
    mapd304,datad304,timed304=get_sunpy_rundiff_maps(list304,'304','AIA')
    mapd335,datad335,timed335=get_sunpy_rundiff_maps(list335,'335','AIA')
    mapd171,datad171,timed171=get_sunpy_rundiff_maps(list171,'171','AIA')
    mapd211,datad211,timed211=get_sunpy_rundiff_maps(list211,'211','AIA')
    mapd193,datad193,timed193=get_sunpy_rundiff_maps(list193,'193','AIA')
    mapd1600,datad1600,timed1600=get_sunpy_rundiff_maps(list1600,'1600','AIA')
    mapd1700,datad1700,timed1700=get_sunpy_rundiff_maps(list1700,'1700','AIA')
    allmapsd={};allmapsd['aiad94']={'mapd94':mapd94,'datad94':datad94,'timed94':timed94}
    allmapsd['aiad131']={'mapd131':mapd131,'datad131':datad131,'timed131':timed131}
    allmapsd['aiad171']={'mapd171':mapd171,'datad171':datad171,'timed171':timed171}
    allmapsd['aiad193']={'mapd193':mapd193,'datad193':datad193,'timed193':timed193}
    allmapsd['aiad211']={'mapd211':mapd211,'datad211':datad211,'timed211':timed211}
    allmapsd['aiad304']={'mapd304':mapd304,'datad304':datad304,'timed304':timed304}
    allmapsd['aiad335']={'mapd335':mapd335,'datad335':datad335,'timed335':timed335}
    allmapsd['aiad1600']={'mapd1600':mapd1600,'datad1600':datad1600,'timed1600':timed1600}
    allmapsd['aiad1700']={'mapd1700':mapd1700,'datad1700':datad1700,'timed1700':timed1700}
    print("Writing..")
    pickle.dump(allmapsd,open('/media/rohit/VLA/20160409/20160409_submap_aia_diff_50ms.p','wb'))
    ##########################
    mapb94,datab94,timeb94=get_sunpy_basediff_maps(list94,'94','AIA')
    mapb131,datab131,timeb131=get_sunpy_basediff_maps(list131,'131','AIA')
    mapb335,datab335,timeb335=get_sunpy_basediff_maps(list335,'335','AIA')
    mapb304,datab304,timeb304=get_sunpy_basediff_maps(list304,'304','AIA')
    mapb171,datab171,timeb171=get_sunpy_basediff_maps(list171,'171','AIA')
    mapb193,datab193,timeb193=get_sunpy_basediff_maps(list193,'193','AIA')
    mapb211,datab211,timeb211=get_sunpy_basediff_maps(list211,'211','AIA')
    mapb1600,datab1600,timeb1600=get_sunpy_basediff_maps(list1600,'1600','AIA')
    mapb1700,datab1700,timeb1700=get_sunpy_basediff_maps(list1700,'1700','AIA')
    allmapsb={};allmapsb['aiab94']={'mapb94':mapb94,'datab94':datab94,'timeb94':timeb94}
    allmapsb['aiab131']={'mapb131':mapb131,'datab131':datab131,'timeb131':timeb131}
    allmapsb['aiab171']={'mapb171':mapb171,'datab171':datab171,'timeb171':timeb171}
    allmapsb['aiab193']={'mapb193':mapb193,'datab193':datab193,'timeb193':timeb193}
    allmapsb['aiab211']={'mapb211':mapb211,'datab211':datab211,'timeb211':timeb211}
    allmapsb['aiab304']={'mapb304':mapb304,'datab304':datab304,'timeb304':timeb304}
    allmapsb['aiab335']={'mapb335':mapb335,'datab335':datab335,'timeb335':timeb335}
    allmapsb['aiab1600']={'mapb1600':mapb1600,'datab1600':datab1600,'timeb1600':timeb1600}
    allmapsb['aiab1700']={'mapb1700':mapb1700,'datab1700':datab1700,'timeb1700':timeb1700}
    print("Writing..")
    pickle.dump(allmapsb,open('/media/rohit/VLA/20160409/20160409_submap_aia_base_50ms.p','wb'))
    ###########################


dump_submaps=0
if(dump_submaps):
    #mapvla,datavla,timevla=get_evla_submap(listvla,-1000,100,-600,400)
    #listvla=sorted(glob.glob('/media/rohit/VLA/20160409/images_1s/sun*FITS'))
    listvla0=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_0/*.spw.0_*FITS'))
    listvla1=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_1/*.spw.1_*FITS'))
    listvla2=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_2/*.spw.2_*FITS'))
    listvla3=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_3/*.spw.3_*FITS'))
    listvla4=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_4/*.spw.4_*FITS'))
    listvla5=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_5/*.spw.5_*FITS'))
    listvla6=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_6/*.spw.6_*FITS'))
    listvla7=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_7/*.spw.7_*FITS'))
    mapvla=[0]*2399;datavla=[0]*2399
    mapvla0=[0]*2399;mapvla1=[0]*2399;mapvla2=[0]*2399;mapvla3=[0]*2399;mapvla4=[0]*2399;mapvla5=[0]*2399;mapvla6=[0]*2399;mapvla7=[0]*2399
    datavla0=[0]*2399;datavla1=[0]*2399;datavla2=[0]*2399;datavla3=[0]*2399;datavla4=[0]*2399;datavla5=[0]*2399;datavla6=[0]*2399;datavla7=[0]*2399
    for i in range(0,2399):
        print(str(listvla0[i].split('.')[8]),i)
        if(os.path.isfile('/media/rohit/VLA/20160409/vlamaps_RR/vlamap_'+str(listvla0[i].split('.')[8])+'_'+"%04d"%i+'.p')==False):
            mapvla0[i]=[0]*4;mapvla1[i]=[0]*4;mapvla2[i]=[0]*4;mapvla3[i]=[0]*4;mapvla4[i]=[0]*4;mapvla5[i]=[0]*4;mapvla6[i]=[0]*4;mapvla7[i]=[0]*4
            datavla0[i]=[0]*4;datavla1[i]=[0]*4;datavla2[i]=[0]*4;datavla3[i]=[0]*4;datavla4[i]=[0]*4;datavla5[i]=[0]*4;datavla6[i]=[0]*4;datavla7[i]=[0]*4
            for k in range(4):
                mapvla0[i][k],datavla0[i][k],timevla=get_evla_submap([listvla0[i+2399*k]],0,-1,0,-1)
                mapvla1[i][k],datavla1[i][k],timevla=get_evla_submap([listvla1[i+2399*k]],0,-1,0,-1)
                mapvla2[i][k],datavla2[i][k],timevla=get_evla_submap([listvla2[i+2399*k]],0,-1,0,-1)
                mapvla3[i][k],datavla3[i][k],timevla=get_evla_submap([listvla3[i+2399*k]],0,-1,0,-1)
                mapvla4[i][k],datavla4[i][k],timevla=get_evla_submap([listvla4[i+2399*k]],0,-1,0,-1)
                mapvla5[i][k],datavla5[i][k],timevla=get_evla_submap([listvla5[i+2399*k]],0,-1,0,-1)
                mapvla6[i][k],datavla6[i][k],timevla=get_evla_submap([listvla6[i+2399*k]],0,-1,0,-1)
                mapvla7[i][k],datavla7[i][k],timevla=get_evla_submap([listvla7[i+2399*k]],0,-1,0,-1)
            mapvla[i]=mapvla0[i]+mapvla1[i]+mapvla2[i]+mapvla3[i]+mapvla4[i]+mapvla5[i]+mapvla6[i]+mapvla7[i]
            print(len(mapvla[i]))
            datavla[i]=datavla0[i]+datavla1[i]+datavla2[i]+datavla3[i]+datavla4[i]+datavla5[i]+datavla6[i]+datavla7[i]
            #allmaps={};allmaps['vla']={'mapvla':mapvla,'datavla':datavla,'timevla':timevla}
            allmaps={};allmaps['vla']={'mapvla':mapvla[i],'timevla':timevla}
            pickle.dump(allmaps,open('/media/rohit/VLA/20160409/vlamaps_RR/vlamap_'+str(listvla0[i].split('.')[8])+'_'+"%04d"%i+'.p','wb'))
    sys.exit()
    mapvla0,datavla0,timevla0=get_evla_submap(listvla0,0,-1,0,-1)
    mapvla1,datavla1,timevla1=get_evla_submap(listvla1,0,-1,0,-1)
    mapvla2,datavla2,timevla2=get_evla_submap(listvla2,0,-1,0,-1)
    mapvla3,datavla3,timevla3=get_evla_submap(listvla3,0,-1,0,-1)
    mapvla4,datavla4,timevla4=get_evla_submap(listvla4,0,-1,0,-1)
    mapvla5,datavla5,timevla5=get_evla_submap(listvla5,0,-1,0,-1)
    mapvla6,datavla6,timevla6=get_evla_submap(listvla6,0,-1,0,-1)
    mapvla7,datavla7,timevla7=get_evla_submap(listvla7,0,-1,0,-1)
    allmaps={};allvlamaps0={};allvlamaps1={};allvlamaps2={};allvlamaps3={};allvlamaps4={};allvlamaps5={};allvlamaps6={};allvlamaps7={}
    allvlamaps0['vla0']={'mapvla':mapvla0,'datavla':datavla0,'timevla':timevla0}
    pickle.dump(allvlamaps0,open('/media/rohit/VLA/20160409/20160409_vla_spw_0_50ms.p','wb'))
    allvlamaps1['vla1']={'mapvla':mapvla1,'datavla':datavla1,'timevla':timevla1}
    pickle.dump(allvlamaps1,open('/media/rohit/VLA/20160409/20160409_vla_spw_1_50ms.p','wb'))
    allvlamaps2['vla2']={'mapvla':mapvla2,'datavla':datavla2,'timevla':timevla2}
    pickle.dump(allvlamaps2,open('/media/rohit/VLA/20160409/20160409_vla_spw_2_50ms.p','wb'))
    allvlamaps3['vla3']={'mapvla':mapvla3,'datavla':datavla3,'timevla':timevla3}
    pickle.dump(allvlamaps3,open('/media/rohit/VLA/20160409/20160409_vla_spw_3_50ms.p','wb'))
    allvlamaps4['vla4']={'mapvla':mapvla4,'datavla':datavla4,'timevla':timevla4}
    pickle.dump(allvlamaps4,open('/media/rohit/VLA/20160409/20160409_vla_spw_4_50ms.p','wb'))
    allvlamaps5['vla5']={'mapvla':mapvla5,'datavla':datavla5,'timevla':timevla5}
    pickle.dump(allvlamaps5,open('/media/rohit/VLA/20160409/20160409_vla_spw_5_50ms.p','wb'))
    allvlamaps6['vla6']={'mapvla':mapvla6,'datavla':datavla6,'timevla':timevla6}
    pickle.dump(allvlamaps6,open('/media/rohit/VLA/20160409/20160409_vla_spw_6_50ms.p','wb'))
    allvlamaps7['vla7']={'mapvla':mapvla7,'datavla':datavla7,'timevla':timevla7}
    pickle.dump(allvlamaps7,open('/media/rohit/VLA/20160409/20160409_vla_spw_7_50ms.p','wb'))

def get_mapI(ml,mr):
    meta=ml.meta;nd=(ml.data+mr.data)*0.5
    mi=Map(nd,meta)
    return mi

sys.exit()

################ GOES

goes=readsav('/media/rohit/VLA/20160409/20160409_idlsave_goes.sav')
gflux=goes['yclean']
gtime=goes['tarray']


################ Only radio analysis

get_submap_pol=1
if(get_submap_pol):
    listvla_v=sorted(glob.glob('/media/rohit/VLA/20160409/images_V/*FITS'))
    listvla_i=sorted(glob.glob('/media/rohit/VLA/20160409/images_I/*FITS'))
    mapvla_v=[0]*32;datavla_v=[0]*32;mapvla_i=[0]*32;datavla_i=[0]*32
    for k in range(32):
        mapvla_v[k],datavla_v[k],timevla=get_evla_submap([listvla_v[k]],0,-1,0,-1)
        mapvla_i[k],datavla_i[k],timevla=get_evla_submap([listvla_i[k]],0,-1,0,-1)

get_max=1
if(get_max):
    spc=['0-15','16-31','32-47','48-63']
    maxTbv=[0]*8;xcvmax=[0]*8;ycvmax=[0]*8;maxTbi=[0]*8;xcimax=[0]*8;ycimax=[0]*8
    maxTbl=[0]*8;xclmax=[0]*8;yclmax=[0]*8;maxTbr=[0]*8;xcrmax=[0]*8;ycrmax=[0]*8
    xci90=[0]*8;yci90=[0]*8;xcv90=[0]*8;ycv90=[0]*8;xcl90=[0]*8;ycl90=[0]*8;xcr90=[0]*8;ycr90=[0]*8
    Tbi_r1=[0]*8;Tbi_r2=[0]*8;Tbv_r1=[0]*8;Tbv_r2=[0]*8;Tbl_r1=[0]*8;Tbl_r2=[0]*8;Tbr_r1=[0]*8;Tbr_r2=[0]*8
    eTbi=[0]*8;eTbl=[0]*8;eTbr=[0]*8;eTbv=[0]*8;areai50=[0]*8
    for i in range(8):
        maxTbv[i]=[0]*4;xcvmax[i]=[0]*4;ycvmax[i]=[0]*4;maxTbi[i]=[0]*4;xcimax[i]=[0]*4;ycimax[i]=[0]*4;maxTbl[i]=[0]*4;xclmax[i]=[0]*4;yclmax[i]=[0]*4;maxTbr[i]=[0]*4;xcrmax[i]=[0]*4;ycrmax[i]=[0]*4
        Tbi_r1[i]=[0]*4;Tbi_r2[i]=[0]*4;Tbv_r1[i]=[0]*4;Tbv_r2[i]=[0]*4;Tbl_r1[i]=[0]*4;Tbl_r2[i]=[0]*4;Tbr_r1[i]=[0]*4;Tbr_r2[i]=[0]*4;areai50[i]=[0]*4;eTbi[i]=[0]*4;eTbl[i]=[0]*4;eTbr[i]=[0]*4;eTbv[i]=[0]*4
        xci90[i]=[0]*4;yci90[i]=[0]*4;xcv90[i]=[0]*4;ycv90[i]=[0]*4;xcl90[i]=[0]*4;ycl90[i]=[0]*4;xcr90[i]=[0]*4;ycr90[i]=[0]*4
        for j in range(4):
            listvla_l=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_LL/spw_'+str(i)+'/*spw.*'+spc[j]+'*.FITS'))[0:2000]
            listvla_r=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_'+str(i)+'/*spw.*'+spc[j]+'*.FITS'))[0:2000]
            maxTbi[i][j]=[0]*2000;xcimax[i][j]=[0]*2000;ycimax[i][j]=[0]*2000;Tbi_r1[i][j]=[0]*2000;Tbi_r2[i][j]=[0]*2000;eTbi[i][j]=[0]*2000;eTbl[i][j]=[0]*2000;eTbr[i][j]=[0]*2000;eTbv[i][j]=[0]*2000
            listvla_v=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_V/spw_'+str(i)+'/*spw.*'+spc[j]+'*.FITS'))[0:2000]
            maxTbv[i][j]=[0]*2000;xcvmax[i][j]=[0]*2000;ycvmax[i][j]=[0]*2000;Tbv_r1[i][j]=[0]*2000;Tbv_r2[i][j]=[0]*2000
            maxTbl[i][j]=[0]*2000;xclmax[i][j]=[0]*2000;yclmax[i][j]=[0]*2000;Tbl_r1[i][j]=[0]*2000;Tbl_r2[i][j]=[0]*2000
            maxTbr[i][j]=[0]*2000;xcrmax[i][j]=[0]*2000;ycrmax[i][j]=[0]*2000;Tbr_r1[i][j]=[0]*2000;Tbr_r2[i][j]=[0]*2000;areai50[i][j]=[0]*2000
            xci90[i][j]=[0]*2000;yci90[i][j]=[0]*2000;xcv90[i][j]=[0]*2000;ycv90[i][j]=[0]*2000;xcl90[i][j]=[0]*2000;ycl90[i][j]=[0]*2000;xcr90[i][j]=[0]*2000;ycr90[i][j]=[0]*2000
            for k in range(2000):
                ml,dl,tl=get_evla_submap([listvla_l[k]],0,-1,0,-1)
                mr,dr,tr=get_evla_submap([listvla_r[k]],0,-1,0,-1)
                hl=ml[0];hr=mr[0]
                Tbi=(hl.data+hr.data)*0.5;maxTbi[i][j][k]=np.nanmax(Tbi)
                yc,xc=np.where(Tbi==np.nanmax(Tbi))
                xcimax[i][j][k]=hl.reference_coordinate.Tx.value+(xc-(hl.reference_pixel.x.value-1))*hl.scale.axis1.value
                ycimax[i][j][k]=hl.reference_coordinate.Ty.value+(yc-(hl.reference_pixel.y.value-1))*hl.scale.axis2.value
                Tbi_r1[i][j][k]=Tbi[138:148,120:130].mean();Tbi_r2[i][j][k]=Tbi[124:136,116:126].mean();eTbi[i][j][k]=Tbi[100:200,15:85].std()
                areai50[i][j][k]=len(np.where(Tbi>np.nanmax(Tbi)*0.5)[0])*2.0
                Tbi_f=Tbi*1.0;Tbi_f[np.isnan(Tbi_f)]=0
                bi=ut.get_bimage(Tbi_f,0.95)
                xcf,ycf=ut.fitEllipse(bi)[0:2]
                xci90[i][j][k]=hl.reference_coordinate.Tx.value+(xcf-(hl.reference_pixel.x.value-1))*hl.scale.axis1.value
                yci90[i][j][k]=hl.reference_coordinate.Ty.value+(ycf-(hl.reference_pixel.y.value-1))*hl.scale.axis2.value
                ######
                mv,dv,tv=get_evla_submap([listvla_v[k]],0,-1,0,-1)
                h=mv[0]
                Tbv=h.data;maxTbv[i][j][k]=np.nanmax(Tbv)
                yc,xc=np.where(Tbv==np.nanmax(Tbv))
                xcvmax[i][j][k]=h.reference_coordinate.Tx.value+(xc-(h.reference_pixel.x.value-1))*h.scale.axis1.value
                ycvmax[i][j][k]=h.reference_coordinate.Ty.value+(yc-(h.reference_pixel.y.value-1))*h.scale.axis2.value
                Tbv_r1[i][j][k]=Tbv[138:148,120:130].mean();Tbv_r2[i][j][k]=Tbv[124:136,116:126].mean();eTbv[i][j][k]=Tbv[100:200,15:85].std()
                Tbv_f=Tbv*1.0;Tbv_f[np.isnan(Tbv_f)]=0
                bi=ut.get_bimage(Tbv_f,0.95)
                xcf,ycf=ut.fitEllipse(bi)[0:2]
                xcv90[i][j][k]=h.reference_coordinate.Tx.value+(xcf-(h.reference_pixel.x.value-1))*h.scale.axis1.value
                ycv90[i][j][k]=h.reference_coordinate.Ty.value+(ycf-(h.reference_pixel.y.value-1))*h.scale.axis2.value
                #######
                h=ml[0]
                Tbl=h.data;maxTbl[i][j][k]=np.nanmax(Tbl)
                yc,xc=np.where(Tbl==np.nanmax(Tbl))
                xclmax[i][j][k]=h.reference_coordinate.Tx.value+(xc-(h.reference_pixel.x.value-1))*h.scale.axis1.value
                yclmax[i][j][k]=h.reference_coordinate.Ty.value+(yc-(h.reference_pixel.y.value-1))*h.scale.axis2.value
                Tbl_r1[i][j][k]=Tbl[138:148,120:130].mean();Tbl_r2[i][j][k]=Tbl[124:136,116:126].mean();eTbl[i][j][k]=Tbl[100:200,15:85].std()
                Tbl_f=Tbl*1.0;Tbl_f[np.isnan(Tbl_f)]=0
                bi=ut.get_bimage(Tbl_f,0.95)
                xcf,ycf=ut.fitEllipse(bi)[0:2]
                xcl90[i][j][k]=h.reference_coordinate.Tx.value+(xcf-(h.reference_pixel.x.value-1))*h.scale.axis1.value
                ycl90[i][j][k]=h.reference_coordinate.Ty.value+(ycf-(h.reference_pixel.y.value-1))*h.scale.axis2.value
                #######
                h=mr[0]
                Tbr=h.data;maxTbr[i][j][k]=np.nanmax(Tbr)
                yc,xc=np.where(Tbr==np.nanmax(Tbr))
                xcrmax[i][j][k]=h.reference_coordinate.Tx.value+(xc-(h.reference_pixel.x.value-1))*h.scale.axis1.value
                ycrmax[i][j][k]=h.reference_coordinate.Ty.value+(yc-(h.reference_pixel.y.value-1))*h.scale.axis2.value
                Tbr_r1[i][j][k]=Tbr[138:148,120:130].mean();Tbr_r2[i][j][k]=Tbr[124:136,116:126].mean();eTbr[i][j][k]=Tbr[100:200,15:85].std()
                Tbr_f=Tbr*1.0;Tbr_f[np.isnan(Tbr_f)]=0
                bi=ut.get_bimage(Tbr_f,0.95)
                xcf,ycf=ut.fitEllipse(bi)[0:2]
                xcr90[i][j][k]=h.reference_coordinate.Tx.value+(xcf-(h.reference_pixel.x.value-1))*h.scale.axis1.value
                ycr90[i][j][k]=h.reference_coordinate.Ty.value+(ycf-(h.reference_pixel.y.value-1))*h.scale.axis2.value
    areai50=np.array(areai50);areai50=areai50.reshape(32,2000)
    maxTbi=np.array(maxTbi);maxTbi=maxTbi.reshape(32,2000);eTbi=np.array(eTbi);eTbi=eTbi.reshape(32,2000);eTbv=np.array(eTbv);eTbv=eTbv.reshape(32,2000);eTbl=np.array(eTbl);eTbl=eTbl.reshape(32,2000);eTbr=np.array(eTbr);eTbr=eTbr.reshape(32,2000)
    xcimax=np.array(xcimax);xcimax=xcimax.reshape(32,2000);ycimax=np.array(ycimax);ycimax=ycimax.reshape(32,2000)
    xci90=np.array(xci90);xci90=xci90.reshape(32,2000);yci90=np.array(yci90);yci90=yci90.reshape(32,2000)
    xcv90=np.array(xcv90);xcv90=xcv90.reshape(32,2000);ycv90=np.array(ycv90);ycv90=ycv90.reshape(32,2000)
    xcl90=np.array(xcl90);xcl90=xcl90.reshape(32,2000);ycl90=np.array(ycl90);ycl90=ycl90.reshape(32,2000)
    xcr90=np.array(xcr90);xcr90=xcr90.reshape(32,2000);ycr90=np.array(ycr90);ycr90=ycr90.reshape(32,2000)
    Tbi_r1=np.array(Tbi_r1);Tbi_r1=Tbi_r1.reshape(32,2000)
    Tbi_r2=np.array(Tbi_r2);Tbi_r2=Tbi_r2.reshape(32,2000)
    pickle.dump([xcimax,ycimax,xci90,yci90,maxTbi,Tbi_r1,Tbi_r2,areai50,eTbi],open('/media/rohit/VLA/20160409/vlamax_loc_i.p','wb'))
    maxTbv=np.array(maxTbv);maxTbv=maxTbv.reshape(32,2000)
    xcvmax=np.array(xcvmax);xcvmax=xcvmax.reshape(32,2000)
    ycvmax=np.array(ycvmax);ycvmax=ycvmax.reshape(32,2000)
    Tbv_r1=np.array(Tbv_r1);Tbv_r1=Tbv_r1.reshape(32,2000)
    Tbv_r2=np.array(Tbv_r2);Tbv_r2=Tbv_r2.reshape(32,2000)
    pickle.dump([xcvmax,ycvmax,xcv90,ycv90,maxTbv,Tbv_r1,Tbv_r2,eTbv],open('/media/rohit/VLA/20160409/vlamax_loc_v.p','wb'))
    maxTbl=np.array(maxTbl);maxTbl=maxTbl.reshape(32,2000)
    xclmax=np.array(xclmax);xclmax=xclmax.reshape(32,2000)
    yclmax=np.array(yclmax);yclmax=yclmax.reshape(32,2000)
    Tbl_r1=np.array(Tbl_r1);Tbl_r1=Tbl_r1.reshape(32,2000)
    Tbl_r2=np.array(Tbl_r2);Tbl_r2=Tbl_r2.reshape(32,2000)
    pickle.dump([xclmax,yclmax,xcl90,ycl90,maxTbl,Tbl_r1,Tbl_r2,eTbl],open('/media/rohit/VLA/20160409/vlamax_loc_l.p','wb'))
    maxTbr=np.array(maxTbr);maxTbr=maxTbr.reshape(32,2000)
    xcrmax=np.array(xcrmax);xcrmax=xcrmax.reshape(32,2000)
    ycrmax=np.array(ycrmax);ycrmax=ycrmax.reshape(32,2000)
    Tbr_r1=np.array(Tbr_r1);Tbr_r1=Tbr_r1.reshape(32,2000)
    Tbr_r2=np.array(Tbr_r2);Tbr_r2=Tbr_r2.reshape(32,2000)
    pickle.dump([xcrmax,ycrmax,xcr90,ycr90,maxTbr,Tbr_r1,Tbr_r2,eTbr],open('/media/rohit/VLA/20160409/vlamax_loc_r.p','wb'))


xcvmax,ycvmax,xcv90,ycv90,maxTbv,Tbv_r1,Tbv_r2,eTbv=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_v.p','rb'),encoding='latin1')
#Tbmayn[np.where(ycmax<255)]=np.nan;Tbmays[np.where(ycmax>255)]=np.nan
xcimax,ycimax,xci90,yci90,maxTbi,Tbi_r1,Tbi_r2,areai50,eTbi=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_i.p','rb'),encoding='latin1')
xclmax,yclmax,xcl90,ycl90,maxTbl,Tbl_r1,Tbl_r2,eTbl=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_l.p','rb'),encoding='latin1')
xcrmax,ycrmax,xcr90,ycr90,maxTbr,Tbr_r1,Tbr_r2,eTbr=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r.p','rb'),encoding='latin1')
Tbmayn=maxTbr*1.0;Tbmays=maxTbr*1.0
#Tbmayn=maxTbi*1.0;Tbmays=maxTbi*1.0
#Tbmayn[np.where(ycmax<255)]=np.nan;Tbmays[np.where(ycmax>255)]=np.nan

qtime=146
plot_polarisation=1
if(plot_polarisation):
    plt.plot(freq,maxTbv[:,qtime]/maxTbi[:,qtime]*100,'o',color='blue')
    plt.errorbar(freq,maxTbv[:,qtime]/maxTbi[:,qtime]*100,yerr=eTbv[:,qtime]/maxTbv[:,qtime]*100,color='blue',label='Continuum (18:44:07.5-18:44:07.55 UT)')
    plt.plot(freq,maxTbv[:,858]/maxTbi[:,858]*100,'o',color='red')
    plt.errorbar(freq,maxTbv[:,858]/maxTbi[:,858]*100,yerr=eTbv[:,858]/maxTbv[:,858]*100,color='red',label='Bursts (18:44:42.9-18:44:42.95 UT)')
    plt.legend();plt.xlabel('Frequency (GHz)');plt.ylabel("Degree of Polarisation (V/I %)");plt.ylim(0,100)
    plt.show()

xx=np.log10(freq*1.e9);ycn=np.log10(maxTbr[:,qtime]);ybs=np.log10(maxTbr[:,858])
zc=np.polyfit(xx,ycn, 1, cov=True);zb=np.polyfit(xx,ybs, 1, cov=True);pc=np.poly1d(zc[0]);pb=np.poly1d(zb[0]);ycn_fit=10**pc(xx);ybs_fit=10**pb(xx)
ezc=np.sqrt(zc[1][0][0]);ezb=np.sqrt(zb[1][0][0])
plot_Tb=1
if(plot_polarisation):
    plt.plot(freq,maxTbr[:,qtime],'o',color='blue')
    plt.errorbar(freq,maxTbr[:,qtime],yerr=eTbr[:,10],color='blue',label='Continuum (18:44:07.5-18:44:07.55 UT)')
    plt.plot(freq,maxTbr[:,858],'o',color='red')
    plt.errorbar(freq,maxTbr[:,858],yerr=eTbr[:,858],color='red',label='Bursts (18:44:42.9-18:44:42.95 UT)')
    plt.plot(freq,ycn_fit,'--',color='blue',label='$\\beta$='+str(np.round(zc[0][0],2))+'$\\pm$'+str(np.round(ezc,2)))
    plt.plot(freq,ybs_fit,'--',color='red',label='$\\beta$='+str(np.round(zb[0][0],2))+'$\\pm$'+str(np.round(ezb,2)))
    plt.legend();plt.xlabel('Frequency (GHz)');plt.ylabel("$T_B$ (K)");plt.xscale('log');plt.yscale('log')#;plt.ylim(0,100)
    plt.show()

do_pearson=1
if(do_pearson):
    pr=[0]*32;pv=[0]*32
    for i in range(32):
        pr[i],pv[i]=stats.pearsonr(maxTbi[0],maxTbi[i])


get_submap_pol=1
if(get_submap_pol):
    listvla_v=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_V/spw_0/*spw.0_16*.FITS'))[0:2000]
    listvla_ll=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_LL/spw_0/*FITS'))[0:2000]
    listvla_rr=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_0/*FITS'))[0:2000]
    mapvla_v=[0]*2000;datavla_v=[0]*2000;mapvla_i=[0]*2000;datavla_i=[0]*2000;timevla=[0]*2000;mapvla_l=[0]*2000;mapvla_r=[0]*2000
    for k in range(2000):
        mapvla_v[k],datavla_v[k],timevla[k]=get_evla_submap([listvla_v[k]],0,-1,0,-1)
        ml,d,t=get_evla_submap([listvla_ll[k]],0,-1,0,-1)
        mr,d,t=get_evla_submap([listvla_rr[k]],0,-1,0,-1)
        mapvla_i[k]=get_mapI(ml[0],mr[0])
        mapvla_l[k]=ml;mapvla_r[k]=mr
vlamax_v=np.nanmax(np.array(datavla_v),axis=(1,2,3))

get_qs=1
if(get_qs):
    listqs=sorted(glob.glob('/media/rohit/VLA/20160409/1844/fits_5_0-15/*spw.*FITS'))
    xcrmax=[0]*len(listqs);ycrmax=[0]*len(listqs);maxTbr=[0]*len(listqs);Tbr_r1=[0]*len(listqs);Tbr_r2=[0]*len(listqs);arear50=[0]*len(listqs);qstimevla=[0]*len(listqs)
    qsxcr90=[0]*len(listqs);qsycr90=[0]*len(listqs)
    for k in range(len(listqs)):
        m,d,qstimevla[k]=get_evla_submap([listqs[k]],0,-1,0,-1)
        hl=m[0]
        Tbr=hl.data;maxTbr[k]=np.nanmax(Tbr)
        yc,xc=np.where(Tbr==np.nanmax(Tbr))
        xcrmax[k]=hl.reference_coordinate.Tx.value+(xc-(hl.reference_pixel.x.value-1))*hl.scale.axis1.value
        ycrmax[k]=hl.reference_coordinate.Ty.value+(yc-(hl.reference_pixel.y.value-1))*hl.scale.axis2.value
        Tbr_r1[k]=Tbr[138:148,120:130].mean();Tbr_r2[k]=Tbr[124:136,116:126].mean()
        arear50[k]=len(np.where(Tbr>np.nanmax(Tbr)*0.5)[0])*2.0
        Tbr_f=Tbr*1.0;Tbr_f[np.isnan(Tbr_f)]=0
        bi=ut.get_bimage(Tbr_f,0.95)
        xcf,ycf=ut.fitEllipse(bi)[0:2]
        qsxcr90[k]=hl.reference_coordinate.Tx.value+(xcf-(hl.reference_pixel.x.value-1))*hl.scale.axis1.value
        qsycr90[k]=hl.reference_coordinate.Ty.value+(ycf-(hl.reference_pixel.y.value-1))*hl.scale.axis2.value
    pickle.dump([xcrmax,ycrmax,qsxcr90,qsycr90,maxTbr,Tbr_r1,Tbr_r2,arear50,qstimevla],open('/media/rohit/VLA/20160409/vlamax_loc_r_qs_5.p','wb'))
qsx,qsy,qsxcr90,qsycr90,qsmaxTbr,qsTbr_r1,qsTbr_r2,qsarear50,qstimevla=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r_qs_3.p','rb'),encoding='latin1')
qsx,qsy,qsxcr90,qsycr90,qsmaxTbr,qsTbr_r1,qsTbr_r2,qsarear50,qstimevla=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r_qs_5.p','rb'),encoding='latin1')
qsx,qsy,qsxcr90,qsycr90,qsmaxTbr,qsTbr_r1,qsTbr_r2,qsarear50,qstimevla=pickle.load(open('/media/rohit/VLA/20160409/vlamax_loc_r_qs.p','rb'),encoding='latin1')
TbLLr1,TbLLr2,TbLLr3,TbLLr4,TbLLr5,xcLLmax,ycLLmax,xcLL90,ycLL90,timevla_all=pickle.load(open('/media/rohit/VLA/20160409/MaxLL.p','rb'),encoding='latin1')
timevla_all1=np.hstack((np.array(qstimevla)[:,0],timevla_all)) 
timevla_all1=np.hstack((np.array(qstimevla)[:,0],timevla_all[0:2000])) 
Tbr1_rall=np.hstack((np.array(qsTbr_r1),Tbr_r1[1]))
xcrmax_all=np.hstack((np.array(qsx)[:,0],xcrmax[1]));ycrmax_all=np.hstack((np.array(qsy)[:,0],ycrmax[1]))

############### Get Spatial Correlation

def do_wavelet(Tbr, N):
    datawave_=Tbr/1.e6
    datawave=datawave_[int(N/2):int(-1*N/2+1)]-np.convolve(datawave_, np.ones(N)/N, mode='valid')
    std=np.std(datawave);datawave_std=datawave/std
    return datawave/std

get_qs=1
if(get_qs):
    ss=['0','3','5'];cc0_t=[0]*3;cc0_s=[0]*3
    for i in range(3):
        listqs=sorted(glob.glob('/media/rohit/VLA/20160409/1844/fits_'+ss[i]+'_0-15/*spw.*FITS'))[0:4799]
        listqs0=sorted(glob.glob('/media/rohit/VLA/20160409/1844/fits_3_0-15/*spw.*FITS'))[0:4799]
        ts1=[0]*4799;ts0=[0]*4799
        for k in range(4799):
            mr0,dr0,tr0=get_evla_submap([listqs0[k]],0,-1,0,-1);ts0[k]=dr0[0][78:178,78:178]
            mr1,dr1,tr1=get_evla_submap([listqs[k]],0,-1,0,-1);ts1[k]=dr1[0][78:178,78:178]
        ts1=np.array(ts1);ts0=np.array(ts0);cc0_s[i]=[0]*100;cc0_t[i]=[0]*100
        for k1 in range(100):
            cc0_s[i][k1]=[0]*100;cc0_t[i][k1]=[0]*100
            for k2 in range(100):
                tt1=ts0[:,135-78,124-78]/np.sum(ts0[:,135-78,124-78]);tt2=ts1[:,k1,k2]/np.sum(ts1[:,k1,k2])
                tt1=do_wavelet(tt1,400);tt2=do_wavelet(tt2,400);cc0_s[i][k1][k2]=np.correlate(tt2,tt1, "same");cc0_s[i][k1][k2]=cc0_s[i][k1][k2]/np.max(cc0_s[i][k1][k2])
                t1=ut.find_nearest(cc0_s[i][k1][k2][2100:2200],np.max(cc0_s[i][k1][k2][2100:2300])*0.4)[0];t2=ut.find_nearest(cc0_s[i][k1][k2][2200:2300][::-1],np.max(cc0_s[i][k1][k2][2100:2300])*0.4)[0]
                cc0_t[i][k1][k2]=t2-t1
    cc0_s=np.array(cc0_s).reshape(3,100,100,4400);cc0_t=np.array(cc0_t).reshape(3,100,100)
    np.save('/media/rohit/VLA/20160409/correlation_t_q.npz',cc0_t)
    np.save('/media/rohit/VLA/20160409/correlation_s_q.npz',cc0_s)


get_max=1
if(get_max):
    spc=['0-15','16-31','32-47','48-63'];cc0_s=[0]*8;cc0_t=[0]*8;pr0_s=[0]*8;pv0_s=[0]*8
    for i in range(1):
        i=2
        cc0_s[i]=[0]*4;cc0_t[i]=[0]*4;pr0_s[i]=[0]*100;pv0_s[i]=[0]*100
        for j in range(1):
            j=2
            listvla_r=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_'+str(i)+'/*spw.*'+spc[j]+'*.FITS'))[0:2000]
            mr0,dr0,tr0=get_evla_submap([listvla_r[0]],0,-1,0,-1);ts1=[0]*2000
            for k in range(2000):
                ts1[k]=[0]*100
                mr1,dr1,tr1=get_evla_submap([listvla_r[k]],0,-1,0,-1);ts1[k]=dr1[0][78:178,78:178]
            ts1=np.array(ts1)
            cc0_s[i][j]=[0]*100;cc0_t[i][j]=[0]*100;pr0_s[i][j]=[0]*100;pv0_s[i][j]=[0]*100
            for k1 in range(100):
                cc0_s[i][j][k1]=[0]*100;cc0_t[i][j][k1]=[0]*100;pr0_s[i][j][k1]=[0]*100;pv0_s[i][j][k1]=[0]*100
                for k2 in range(100):
                    #tt1=ts1[:,124-78,135-78]/np.sum(ts1[:,124-78,135-78]);tt2=ts1[:,k1,k2]/np.sum(ts1[:,k1,k2])
                    tt1=ts1[:,135-78,124-78]/np.sum(ts1[:,135-78,124-78]);tt2=ts1[:,k1,k2]/np.sum(ts1[:,k1,k2])
                    tt1=do_wavelet(tt1,400);tt2=do_wavelet(tt2,400)
                    cc0_s[i][j][k1][k2]=np.correlate(tt2,tt1, "same");cc0_s[i][j][k1][k2]=cc0_s[i][j][k1][k2]/np.max(cc0_s[i][j][k1][k2])
                    t1=ut.find_nearest(cc0_s[i][j][k1][k2][750:800],np.max(cc0_s[i][j][k1][k2][750:850])*0.4)[0];t2=ut.find_nearest(cc0_s[i][j][k1][k2][800:850][::-1],np.max(cc0_s[i][j][k1][k2][750:850])*0.4)[0]
                    cc0_t[i][j][k1][k2]=t2-t1
                    pr0_s[i][j][k1][k2],pv0_s[i][j][k1][k2]=stats.pearsonr(tt2,tt1)
                    #cc0_t[i][j][k1][k2]=ut.find_nearest(cc0_s[i][j][k1][k2],np.max(cc0_s[i][j][k1][k2]))[0]-1000
            sys.exit()
    cc0_s=np.array(cc0_s).reshape(32,100,100,1601);cc0_t=np.array(cc0_t).reshape(32,100,100)
    np.save('/media/rohit/VLA/20160409/correlation_t_0.npz',cc0_t)
    np.save('/media/rohit/VLA/20160409/correlation_s_0.npz',cc0_s)
                
get_max=1
if(get_max):
    spc=['0-15','16-31','32-47','48-63'];cc0_s=[0]*8;cc0_t=[0]*8
    for i in range(1):
        i=2
        cc0_s[i]=[0]*4;cc0_t[i]=[0]*4
        for j in range(1):
            j=2
            listvla_r=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_'+str(i)+'/*spw.*'+spc[j]+'*.FITS'))[0:2000]
            #listvla_r0=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_0/*spw.*0-15*.FITS'))[0:2000]
            listvla_r0=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_2/*spw.*32-47*.FITS'))[0000:2000]
            ts1=[0]*2000;ts0=[0]*2000
            for k in range(2000):
                mr1,dr1,tr1=get_evla_submap([listvla_r[k]],0,-1,0,-1);ts1[k]=dr1[0][78:178,78:178]
                mr0,dr0,tr0=get_evla_submap([listvla_r0[k]],0,-1,0,-1);ts0[k]=dr0[0][78:178,78:178]
            ts1=np.array(ts1);ts0=np.array(ts0)
            cc0_s[i][j]=[0]*100;cc0_t[i][j]=[0]*100
            for k1 in range(100):
                cc0_s[i][j][k1]=[0]*100;cc0_t[i][j][k1]=[0]*100
                for k2 in range(100):
                    #tt1=ts1[:,124-78,135-78]/np.sum(ts1[:,124-78,135-78]);tt2=ts1[:,k1,k2]/np.sum(ts1[:,k1,k2])
                    tt1=ts0[:,135-78,124-78]/np.sum(ts0[:,135-78,124-78]);tt2=ts1[:,k1,k2]/np.sum(ts1[:,k1,k2])
                    tt1=do_wavelet(tt1,400);tt2=do_wavelet(tt2,400)
                    cc0_s[i][j][k1][k2]=np.correlate(tt2,tt1, "same");cc0_s[i][j][k1][k2]=cc0_s[i][j][k1][k2]/np.max(cc0_s[i][j][k1][k2])
                    t1=ut.find_nearest(cc0_s[i][j][k1][k2][750:800],np.max(cc0_s[i][j][k1][k2][750:850])*0.4)[0];t2=ut.find_nearest(cc0_s[i][j][k1][k2][800:850][::-1],np.max(cc0_s[i][j][k1][k2][750:850])*0.4)[0]
                    cc0_t[i][j][k1][k2]=t2-t1
                    #cc0_t[i][j][k1][k2]=ut.find_nearest(cc0_s[i][j][k1][k2],np.max(cc0_s[i][j][k1][k2]))[0]-1000
            sys.exit()
    cc0_s=np.array(cc0_s).reshape(32,100,100,1601);cc0_t=np.array(cc0_t).reshape(32,100,100)
    np.save('/media/rohit/VLA/20160409/correlation_t_spw2.npz',cc0_t)
    np.save('/media/rohit/VLA/20160409/correlation_s_spw2.npz',cc0_s)

############### HMI
hmilist=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/hmi/*rot.sav'))
maphmi=[0]*len(hmilist);datahmi=[0]*len(hmilist)
#for k in range(len(hmilist)):
#    maphmi[k],datahmi[k],timevla=get_hmi_submap([hmilist[k]],0,-1,0,-1)
maphmi,datahmi,timehmi=get_sunpy_maps_rot(hmilist,'','HMI')
maphmib,datalistb,time=get_sunpy_basediff_maps(hmilist,'','HMI')
maphmir,datalistr,time=get_sunpy_rundiff_maps(hmilist,'','HMI')

hmifile='/media/rohit/VLA/20160409_EUV/hmi.m_45s.2016.04.09_18_45_00_TAI.magnetogram.fits'
hmimap=Map(hmifile)
hmid=hmimap.data[::-1,::-1]
hmid[np.where(hmid<-5000)]=0
hmimap=Map(hmid,hmimap.meta)

# Extrapolation
from scipy.io import readsav
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
mapex_babs=pickle.load(open('/media/rohit/VLA/20160409_EUV/20160409_babs.p','rb'))
mapex_bz=pickle.load(open('/media/rohit/VLA/20160409_EUV/20160409_bz.p','rb'))
mapex_bx=pickle.load(open('/media/rohit/VLA/20160409_EUV/20160409_bx.p','rb'))
mapex_by=pickle.load(open('/media/rohit/VLA/20160409_EUV/20160409_by.p','rb'))
dx_km=1400
freq=np.round(np.linspace(0.994,2.006,32),3)

######################## Extrapolated Magnetic Fields
pmap=allmaps['aia171']['map171'][0]
b_hp_fr,b_hp_ls,b_hp_euv,hmiwcs=pickle.load(open('/media/rohit/VLA/paraview/bproj.p','rb'))
from astropy.wcs import utils
b_proj_fr = utils.skycoord_to_pixel(b_hp_fr, pmap.wcs)
b_proj_ls = utils.skycoord_to_pixel(b_hp_ls, pmap.wcs)
b_proj_euv = utils.skycoord_to_pixel(b_hp_euv, pmap.wcs)



plot_ds=1
if(plot_ds):
    f,(ax1,ax2,ax3)=plt.subplots(3,1,figsize=(15,15))
    ax1.plot(gtime,gflux[0])
    ax1.plot(gtime,gflux[1])
    ax2.imshow(TbLLr1.swapaxes(0,1),aspect='auto')
    plt.show()



plot_position=1
i=0
for i in range(32):
    f,ax=plt.subplots(1,1,figsize=(15,6))
    #ax.plot(np.arange(2000)*0.05,Tbmayn[i]/1.e6,'o',color='blue',label='North Source')
    ax.plot(np.arange(2000)*0.05,Tbmays[i]/1.e6,'o',color='red',label='South Source')
    ax.set_ylabel('$T_B$ (MK)');ax.set_ylabel('$T_B$ (MK)');ax.set_xlabel('Time (sec)')
    ax.legend();ax.set_title(str(freq[i])+' GHz');ax.set_ylim(0,40);ax.set_xlim(0,100)
    plt.savefig('/media/rohit/VLA/20160409/spec_max/max_poss_'+"%02d"%i+'.png')
    plt.close()

plot_dcp_hmi=1
if(plot_dcp_hmi):
    for i in range(32):
        f=plt.figure(figsize=(10,10))
        ax0 = f.add_subplot(111)
        #xlaia=cc.center.Tx.value-0.6*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.6*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.6*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.6*int(cc.data.shape[0]/2)
        xlaia=-1037.6;xraia=1031.6;ylaia=-1031.2;yraia=1031.9
        p=hmimap.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-2000,vmax=2000)
        #cc.draw_grid()
        #lev1=(1.5e7/dd0.data.max())*np.array([20,30,40,50,60,70,80,90])*u.percent
        lev0=np.array([50,55,60,65,70,75,80,85,90])*u.percent
        lev1=np.array([10,20,30,40,50,60,70,80,90])*u.percent
        lev2=np.array([1,5,10,20,50,60,70,80,90])*u.percent
        dcp_data=mapvla_v[i][0].data*100/mapvla_i[i][0].data;dcp_data[np.isnan(dcp_data)]=0;dcp_data[np.where(dcp_data<0)]=0;dcp_data[np.where(dcp_data>100)]=0
        dd0=Map(dcp_data,mapvla_i[i][0].meta)
        xlvla=dd0.center.Tx.value-2.0*int(dd0.data.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(dd0.data.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(dd0.data.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(dd0.data.shape[0]/2)
        dd0.draw_contours(levels=lev0,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        #mapvla_i[i][0].draw_contours(levels=lev1,colors='white',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        #mapvla_v[i][0].draw_contours(levels=lev1,colors='yellow',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        #ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
        #ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
        ax0.set_xlim([-850,-700]);ax0.set_ylim([200,350])
        ax0.text(-840,210,str(np.round(dd0.meta['crval3']/1.e9,4))+' GHz',color='blue')
        ax0.text(-840,330,'Contours (%): 50, 55, 60, 65, 70, 75, 80, 85, 90',color='yellow')
        plt.savefig('pngs_dcp_hmi/dcp_hmi_'+"%02d"%i+'.png')
        plt.close()

plot_vi_hmi=1
if(plot_vi_hmi):
    for i in range(32):
        f=plt.figure(figsize=(10,10))
        ax0 = f.add_subplot(111)
        #xlaia=cc.center.Tx.value-0.6*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.6*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.6*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.6*int(cc.data.shape[0]/2)
        xlaia=-1037.6;xraia=1031.6;ylaia=-1031.2;yraia=1031.9
        p=hmimap.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-2000,vmax=2000)
        #cc.draw_grid()
        #lev1=(1.5e7/dd0.data.max())*np.array([20,30,40,50,60,70,80,90])*u.percent
        lev0=np.array([50,55,60,65,70,75,80,85,90])*u.percent
        lev1=np.array([10,20,30,40,50,60,70,80,90])*u.percent
        lev2=np.array([1,5,10,20,50,60,70,80,90])*u.percent
        dcp_data=mapvla_v[i][0].data*100/mapvla_i[i][0].data;dcp_data[np.isnan(dcp_data)]=0;dcp_data[np.where(dcp_data<0)]=0;dcp_data[np.where(dcp_data>100)]=0
        dd0=Map(dcp_data,mapvla_i[i][0].meta)
        xlvla=dd0.center.Tx.value-2.0*int(dd0.data.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(dd0.data.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(dd0.data.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(dd0.data.shape[0]/2)
        #dd0.draw_contours(levels=lev0,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mapvla_i[i][0].data[np.isnan(mapvla_i[i][0].data)]=0;mapvla_v[i][0].data[np.isnan(mapvla_v[i][0].data)]=0
        mapvla_i[i][0].draw_contours(levels=lev1,colors='white',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mapvla_v[i][0].draw_contours(levels=lev1,colors='yellow',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        #ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
        #ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
        ax0.set_xlim([-850,-700]);ax0.set_ylim([200,350])
        ax0.text(-840,210,str(np.round(dd0.meta['crval3']/1.e9,4))+' GHz',color='blue')
        ax0.text(-840,330,'Contours (%): 10,20,30,40,50,60,70,80,90',color='yellow')
        plt.savefig('pngs_vi_hmi/vi_hmi_'+"%02d"%i+'.png')
        plt.close()

plot_dcp_2d=1
if(plot_dcp_2d):
    f=plt.figure(figsize=(10,10))
    ax0 = f.add_subplot(111)
    ax0.plot(freq,maxTbv[:,770:870].mean(axis=1)/maxTbi[:,770:870].mean(axis=1)*100,'o-')
    ax0.errorbar(freq,maxTbv[:,770:870].mean(axis=1)/maxTbi[:,770:870].mean(axis=1)*100,yerr=eTbv[:,770]/maxTbi[:,770]*100)
    ax0.set_ylabel('DCP (%)');ax0.set_xlabel('Frequency (GHz)')
    plt.show()


plot_dcp=1
if(plot_dcp):
    for i in range(32):
        f=plt.figure(figsize=(10,10))
        ax0 = f.add_subplot(111)
        tidx_171=ut.find_predecessor(allmaps['aia171']['time171'],allmaps['aia171']['time171'][0])[0]
        cc=allmaps['aia171']['map171'][tidx_171]
        #xlaia=cc.center.Tx.value-0.6*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.6*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.6*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.6*int(cc.data.shape[0]/2)
        xlaia=-1230;xraia=-569;ylaia=-47.9;yraia=572
        p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
        #cc.draw_grid()
        #lev1=(1.5e7/dd0.data.max())*np.array([20,30,40,50,60,70,80,90])*u.percent
        lev0=np.array([50,55,60,65,70,75,80,85,90])*u.percent
        lev1=np.array([10,20,30,40,50,60,70,80,90])*u.percent
        lev2=np.array([1,5,10,20,50,60,70,80,90])*u.percent
        dcp_data=mapvla_v[i][0].data*100/mapvla_i[i][0].data;dcp_data[np.isnan(dcp_data)]=0;dcp_data[np.where(dcp_data<0)]=0;dcp_data[np.where(dcp_data>100)]=0
        dd0=Map(dcp_data,mapvla_i[i][0].meta)
        xlvla=dd0.center.Tx.value-2.0*int(dd0.data.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(dd0.data.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(dd0.data.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(dd0.data.shape[0]/2)
        #dd0.draw_contours(levels=lev0,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mapvla_i[i][0].data[np.isnan(mapvla_i[i][0].data)]=0;mapvla_v[i][0].data[np.isnan(mapvla_v[i][0].data)]=0
        mapvla_i[i][0].draw_contours(levels=lev1,colors='white',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mapvla_v[i][0].draw_contours(levels=lev1,colors='yellow',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        #ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
        #ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
        ax0.set_xlim([-850,-700]);ax0.set_ylim([200,350])
        ax0.text(-840,210,str(np.round(dd0.meta['crval3']/1.e9,4))+' GHz',color='blue')
        ax0.text(-840,330,'Contours (%): 50, 55, 60, 65, 70, 75, 80, 85, 90',color='yellow')
        plt.savefig('pngs_vi_aia171/vi_aia171_'+"%02d"%i+'.png')
        plt.close()



plot_dcp=1
if(plot_dcp):
    for i in range(2000):
        f,ax=plt.subplots(2,1,figsize=(8,15));ax0=ax[0];ax1=ax[1]
        tidx_171=ut.find_predecessor(allmaps['aia171']['time171'],timevla[i])[0]
        cc=allmaps['aia171']['map171'][tidx_171]
        #xlaia=cc.center.Tx.value-0.6*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.6*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.6*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.6*int(cc.data.shape[0]/2)
        xlaia=-1230;xraia=-569;ylaia=-47.9;yraia=572
        p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
        #cc.draw_grid()
        #lev1=(1.5e7/dd0.data.max())*np.array([20,30,40,50,60,70,80,90])*u.percent
        lev0=np.array([50,55,60,65,70,75,80,85,90])*u.percent
        lev1=np.array([30,40,50,60,70,80,90])*u.percent
        lev2=np.array([1,5,10,20,50,60,70,80,90])*u.percent
        dcp_data=mapvla_v[i][0].data*100/mapvla_i[i].data*(mapvla_v[i][0].data/np.nanmax(mapvla_v[i][0].data));dcp_data[np.isnan(dcp_data)]=0;dcp_data[np.where(dcp_data<0)]=0;dcp_data[np.where(dcp_data>100)]=0
        dd0=Map(dcp_data,mapvla_i[i].meta)
        xlvla=dd0.center.Tx.value-2.0*int(dd0.data.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(dd0.data.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(dd0.data.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(dd0.data.shape[0]/2)
        dd0.draw_contours(levels=lev0,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mapvla_i[i].data[np.isnan(mapvla_i[i].data)]=0;mapvla_v[i][0].data[np.isnan(mapvla_v[i][0].data)]=0
        #mapvla_i[i].draw_contours(levels=lev1,colors='white',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        #mapvla_v[i][0].draw_contours(levels=lev1,colors='yellow',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        #ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
        #ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
        ax0.set_xlim([-850,-700]);ax0.set_ylim([200,350])
        ax0.text(-840,210,str(np.round(dd0.meta['crval3']/1.e9,4))+' GHz',color='blue')
        ax0.text(-840,330,'Contours (%): 50, 55, 60, 65, 70, 75, 80, 85, 90',color='yellow')
        ax1.plot(np.arange(2000)*0.05,vlamax_v,'-',label='Stokes V');ax1.axvline(x=i*0.05,color='k');ax1.set_xlabel('Time (sec)');ax1.set_ylabel('T$_{B}$');ax1.legend()
        plt.savefig('pngs_v_aia171_time/dcp_aia171_'+"%02d"%i+'.png')
        plt.close()


plot_dcp=1
if(plot_dcp):
    #for i in range(2000):
    for i in range(800,1200):
        f,ax=plt.subplots(2,1,figsize=(8,15));ax0=ax[0];ax1=ax[1]
        tidx_171=ut.find_predecessor(allmaps['aia171']['time171'],timevla[i])[0]
        cc=allmaps['aia171']['map171'][tidx_171]
        #cc0=allmaps['aia171']['map171'][tidx_171]
        #cc1=allmaps['aia171']['map171'][tidx_171-1]
        #cc=Map(cc0.data-cc1.data,cc0.meta)
        #xlaia=cc.center.Tx.value-0.6*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.6*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.6*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.6*int(cc.data.shape[0]/2)
        xlaia=-1230;xraia=-569;ylaia=-47.9;yraia=572
        p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
        #p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-80,vmax=80)
        #cc.draw_grid()
        #lev1=(1.5e7/dd0.data.max())*np.array([20,30,40,50,60,70,80,90])*u.percent
        lev1=np.array([30,40,50,60,70,80,90])/np.nanmax(mapvla_v[i][0].data)*(4.e7)*u.percent
        lev2=np.array([1,5,10,20,50,60,70,80,90])*u.percent
        dcp_data=mapvla_v[i][0].data*100/mapvla_i[i].data*(mapvla_v[i][0].data/np.nanmax(mapvla_v[i][0].data));dcp_data[np.isnan(dcp_data)]=0;dcp_data[np.where(dcp_data<0)]=0;dcp_data[np.where(dcp_data>100)]=0
        dd0=Map(dcp_data,mapvla_i[i].meta)
        lev0=np.array([50,55,60,65,70,75,80,85,90])/100*np.nanmax(dd0.data)*u.percent
        xlvla=dd0.center.Tx.value-2.0*int(dd0.data.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(dd0.data.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(dd0.data.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(dd0.data.shape[0]/2)
        #dd0.draw_contours(levels=lev0,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mapvla_i[i].data[np.isnan(mapvla_i[i].data)]=0;mapvla_v[i][0].data[np.isnan(mapvla_v[i][0].data)]=0
        #mapvla_i[i].draw_contours(levels=lev1,colors='white',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mapvla_v[i][0].draw_contours(levels=lev1,colors='yellow',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        #ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
        #ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
        ax0.set_xlim([-850,-700]);ax0.set_ylim([200,350])
        ax0.text(-840,210,str(np.round(dd0.meta['crval3']/1.e9,4))+' GHz',color='blue')
        #ax0.text(-840,330,'Contours (%): 50, 55, 60, 65, 70, 75, 80, 85, 90',color='blue')
        ax0.text(-840,330,'Contours (%): 20, 22, 24, 26, 28, 30, 32, 34, 36 MK',color='green')
        ax1.plot(np.arange(2000)*0.05,vlamax_v,'o-',label='Stokes V');ax1.axvline(x=i*0.05,color='k');ax1.set_xlabel('Time (sec)');ax1.set_ylabel('T$_{B}$');ax1.legend()
        ax1.set_xlim(40,60)
        #plt.savefig('pngs_v_aia171_time/dcp_run_aia171_'+"%02d"%i+'.png')
        plt.savefig('pngs_v_aia171_time/v_abs_aia171_'+"%02d"%i+'.png')
        plt.close()


xcvmax=[0]*2000;ycvmax=[0]*2000
xclmax=[0]*2000;yclmax=[0]*2000
xcrmax=[0]*2000;ycrmax=[0]*2000
xcvmed=[0]*2000;ycvmed=[0]*2000
xcimax=[0]*2000;ycimax=[0]*2000
xcimed=[0]*2000;ycimed=[0]*2000
Tbi_r1=[0]*2000;Tbi_r2=[0]*2000
Tbv_r1=[0]*2000;Tbv_r2=[0]*2000
Tbl_r1=[0]*2000;Tbl_r2=[0]*2000
Tbr_r1=[0]*2000;Tbr_r2=[0]*2000
for k in range(2000):
    h=mapvla_v[k][0]
    Tbv=h.data
    yc,xc=np.where(Tbv==np.nanmax(Tbv))
    xcvmax[k]=h.reference_coordinate.Tx.value+(xc-(h.reference_pixel.x.value-1))*h.scale.axis1.value
    ycvmax[k]=h.reference_coordinate.Ty.value+(yc-(h.reference_pixel.y.value-1))*h.scale.axis2.value
    Tbv_r1[k]=h.data[124:137,121:124].mean();Tbv_r2[k]=h.data[138:141,123:126].mean()
    Tbv[np.where(Tbv<Tbv.max()*0.5)]=np.nan
    idx=np.where(np.abs(Tbv-np.nanmedian(Tbv))==np.nanmin(np.abs(Tbv-np.nanmedian(Tbv))))
    if(len(idx[0])==1):
        yc,xc=idx[1],idx[0]
    else:
        yc,xc=idx[0][1],idx[0][0]
    xcvmed[k]=h.reference_coordinate.Tx.value+(xc-(h.reference_pixel.x.value-1))*h.scale.axis1.value
    ycvmed[k]=h.reference_coordinate.Ty.value+(yc-(h.reference_pixel.y.value-1))*h.scale.axis2.value
    hi=mapvla_i[k][0]
    Tbi=hi.data
    yc,xc=np.where(Tbi==np.nanmax(Tbi))
    xcimax[k]=hi.reference_coordinate.Tx.value+(xc-(hi.reference_pixel.x.value-1))*hi.scale.axis1.value
    ycimax[k]=hi.reference_coordinate.Ty.value+(yc-(hi.reference_pixel.y.value-1))*hi.scale.axis2.value
    Tbi_r1[k]=hi.data[124:137,121:124].mean();Tbi_r2[k]=hi.data[138:141,123:126].mean()
    Tbi[np.where(Tbi<Tbi.max()*0.5)]=np.nan
    idx=np.where(np.abs(Tbi-np.nanmedian(Tbi))==np.nanmin(np.abs(Tbi-np.nanmedian(Tbi))))
    if(len(idx[0])==1):
        yc,xc=idx[1],idx[0]
    else:
        yc,xc=idx[0][1],idx[0][0]
    xcimed[k]=hi.reference_coordinate.Tx.value+(xc-(hi.reference_pixel.x.value-1))*hi.scale.axis1.value
    ycimed[k]=hi.reference_coordinate.Ty.value+(yc-(hi.reference_pixel.y.value-1))*hi.scale.axis2.value
    hl=mapvla_l[k][0]
    Tbl=hl.data
    yc,xc=np.where(Tbl==np.nanmax(Tbl))
    xclmax[k]=hl.reference_coordinate.Tx.value+(xc-(hl.reference_pixel.x.value-1))*hl.scale.axis1.value
    yclmax[k]=hl.reference_coordinate.Ty.value+(yc-(hl.reference_pixel.y.value-1))*hl.scale.axis2.value
    Tbl_r1[k]=hl.data[124:137,121:124].mean();Tbl_r2[k]=hl.data[138:141,123:126].mean()
    hr=mapvla_r[k][0]
    Tbr=hr.data
    yc,xc=np.where(Tbl==np.nanmax(Tbl))
    xclmax[k]=hl.reference_coordinate.Tx.value+(xc-(hl.reference_pixel.x.value-1))*hl.scale.axis1.value
    yclmax[k]=hl.reference_coordinate.Ty.value+(yc-(hl.reference_pixel.y.value-1))*hl.scale.axis2.value
    Tbr_r1[k]=hi.data[124:137,121:124].mean();Tbr_r2[k]=hi.data[138:141,123:126].mean()
    


i=0
for i in range(32):
    f,ax=plt.subplots(2,1,figsize=(15,8),sharex=True);ax0=ax[0];ax1=ax[1]
    ax0.plot(np.arange(2000)*0.05,np.array(Tbi_r1[i])/1.e6,'-',label='R1 (North)') #
    ax0.plot(np.arange(2000)*0.05,np.array(Tbi_r2[i])/1.e6,'-',label='R2 (South)');ax0.legend()
    ax1.plot(np.arange(2000)*0.05,(np.array(Tbi_r2[i])-np.array(Tbi_r1[i]))/1.e6,'-',label='R2-R1');ax1.legend()
    ax1.set_ylim(-8,12);ax0.set_ylim(0,50);ax0.legend();ax0.set_title(str(freq[i])+' GHz')
    ax1.set_ylabel('$T_B$ (MK)');ax0.set_ylabel('T$_B$ (MK)');ax1.set_xlabel('time (sec)')
    plt.savefig('/media/rohit/VLA/20160409/pngs_spec/Tbi_regions_'+"%02d"%i+'.png')
    plt.close()

i=0
for i in range(32):
    f,ax=plt.subplots(2,1,figsize=(15,8),sharex=True);ax0=ax[0];ax1=ax[1]
    ax0.plot(np.arange(2000)*0.05,np.array(Tbi_r1[i])/1.e6,'-',label='R1 (North)') #
    ax0.plot(np.arange(2000)*0.05,np.array(Tbi_r2[i])/1.e6,'-',label='R2 (South)');ax0.legend()
    ax1.plot(np.arange(2000)*0.05,(np.array(Tbv_r1[i])/np.array(Tbi_r1[i]))*100,'-',label='R1 (North)');ax1.legend()
    ax1.plot(np.arange(2000)*0.05,(np.array(Tbv_r2[i])/np.array(Tbi_r2[i]))*100,'-',label='R2 (South)');ax1.legend()
    ax1.set_ylim(0,100);ax0.set_ylim(0,50);ax0.legend();ax0.set_title(str(freq[i])+' GHz')
    ax1.set_ylabel('DCP (V/I) (\%)');ax0.set_ylabel('T$_b$ (MK)');ax1.set_xlabel('time (sec)')
    plt.savefig('/media/rohit/VLA/20160409/pngs_spec/Tbi_dcp_regions_'+"%02d"%i+'.png')
    plt.close()

i=0
for i in range(2000):
    f,ax=plt.subplots(2,1,figsize=(15,8));ax0=ax[0];ax1=ax[1]
    ax0.plot(np.arange(2000)*0.05,np.array(Tbi_r1[0])/1.e6,'-',label='R1 (North)') #
    ax0.plot(np.arange(2000)*0.05,np.array(Tbi_r2[0])/1.e6,'-',label='R2 (South)');ax0.legend()
    ax0.axvline(x=np.arange(2000)[i]*0.05,color='k')
    ax1.plot(freq,(np.array(Tbv_r1[:,i])/np.array(Tbi_r1[:,i]))*100,'o-',label='R1 (North)');ax1.legend()
    ax1.plot(freq,(np.array(Tbv_r2[:,i])/np.array(Tbi_r2[:,i]))*100,'o-',label='R2 (South)');ax1.legend()
    ax1.set_ylim(0,100);ax0.set_ylim(0,50);ax0.legend();ax0.set_title(str(freq[0])+' GHz')
    ax1.set_ylabel('DCP (V/I) (\%)');ax0.set_ylabel('T$_b$ (MK)');ax0.set_xlabel('time (sec)');ax1.set_xlabel('Frequency (MHz)')
    plt.savefig('/media/rohit/VLA/20160409/pngs_spec/dcp_spec_regions_'+"%02d"%i+'.png')
    plt.close()

f,ax=plt.subplots(2,1,figsize=(8,15),sharex=True);ax0=ax[0];ax1=ax[1]
ax0.plot(np.arange(2000)*0.05,np.array(Tb_r1)/1.e6,'o-',label='R1')
ax0.plot(np.arange(2000)*0.05,np.array(Tb_r2)/1.e6,'o-',label='R2');ax0.legend()
ax1.plot(np.arange(2000)*0.05,(np.array(Tb_r2)-np.array(Tb_r1))/1.e6,'o-',label='R2-R1');ax1.legend()
ax1.set_ylabel('$T_B$ (MK)');ax0.set_ylabel('T$_B$ (MK)');ax1.set_xlabel('Time (sec)')
plt.show()

f,ax=plt.subplots(2,1,figsize=(8,15),sharex=True);ax0=ax[0];ax1=ax[1]
y=np.hstack((np.array(qsy)[:,0],ycrmax[1]))
ax0.plot(np.arange(len(y))*0.05,y,'o-')
ax1.plot(np.arange(len(y))*0.05,np.hstack((qsTbr_r1,Tbr_r1[1]))/1.e6,'o-')
ax1.set_ylabel('$T_B$ (MK)');ax0.set_ylabel('Y-Coordinate (arcsec)');ax1.set_xlabel('Time (sec)')
plt.show()

f,ax=plt.subplots(2,1,figsize=(8,15),sharex=True);ax0=ax[0];ax1=ax[1]
ax0.plot(np.arange(2000)*0.05,ycvmax,'o-')
ax1.plot(np.arange(2000)*0.05,vlamax_v/1.e6,'o-')
ax1.set_ylabel('$T_B$ (MK)');ax0.set_ylabel('Y-Coordinate (arcsec)');ax1.set_xlabel('Time (sec)')
plt.show()

idx1=np.where(ycrmax_all<255)[0];idx2=np.where(ycrmax_all>255)[0]
t=np.arange(len(ycrmax_all))*0.05
t1=t[idx1];ycrmax1=ycrmax_all[idx1];t2=t[idx2];ycrmax2=ycrmax_all[idx2];xcrmax1=xcrmax_all[idx1];t2=t[idx2];xcrmax2=xcrmax_all[idx2];vlamax_r1=Tbr1_rall[idx1];vlamax_r2=Tbr1_rall[idx2]
f,ax=plt.subplots(3,1,figsize=(8,15),sharex=True);ax0=ax[0];ax1=ax[1];ax2=ax[2]
ax0.plot(t1,ycrmax1,'o',color='r');ax0.plot(t2,ycrmax2,'o',color='g')
ax1.plot(t1,xcrmax1,'o',color='r');ax1.plot(t2,xcrmax2,'o',color='g')
ax2.plot(t,Tbr1_rall/1.e6,'-',color='k')
ax2.plot(t1,vlamax_r1/1.e6,'o',markersize=8,color='r');ax2.plot(t2,vlamax_r2/1.e6,'o',markersize=8,color='g');ax1.set_ylim(-750,-800)
ax2.set_ylabel('$T_B$ (MK)');ax0.set_ylabel('Y-Coordinate (arcsec)');ax2.set_xlabel('Time (sec)')
plt.show()

t=np.arange(2000)*0.05
f,ax=plt.subplots(2,1,figsize=(8,15));ax0=ax[0];ax1=ax[1]
#cc=allmaps['aia171']['map171'][20]
#p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
xlaia=-1037.6;xraia=1031.6;ylaia=-1031.2;yraia=1031.9;p=hmimap.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-2000,vmax=2000)
f.colorbar(p,label='B (G)')
ax0.plot(xcvmax[781],ycvmax[781],'o',color='g')
ax0.plot(xcvmax[858],ycvmax[858],'o',color='g')
ax0.plot(xcvmax[883],ycvmax[883],'o',color='g')
ax0.plot(xcvmax[954],ycvmax[954],'o',color='r')
ax0.plot(xcvmax[966],ycvmax[966],'o',color='r')
ax0.plot(xcvmax[1058],ycvmax[1058],'o',color='g')
ax0.plot(xcvmax[1153],ycvmax[1153],'o',color='r')
#ax0.plot(xcvmax[880:885],ycvmax[880:885],'o',color='g')
#ax0.plot(xcvmax[944:950],ycvmax[944:950],'o',color='r')
#ax0.plot(xcvmax[1057:1062],ycvmax[1057:1062],'o',color='g')
ax1.plot(t,vlamax_v,'o-',label='Stokes V');ax1.axvline(x=i*0.05,color='k');ax1.set_xlabel('Time (sec)');ax1.set_ylabel('T$_{B}$');ax1.legend(loc=4)
ax1.plot(t[781],vlamax_v[781],'o',color='g')
ax1.plot(t[858],vlamax_v[858],'o',color='g')
ax1.plot(t[883],vlamax_v[883],'o',color='g')
ax1.plot(t[954],vlamax_v[954],'o',color='r')
ax1.plot(t[966],vlamax_v[966],'o',color='r')
ax1.plot(t[1058],vlamax_v[1058],'o',color='g')
ax1.plot(t[1153],vlamax_v[1153],'o',color='r')
#ax1.plot(t[855:860],vlamax_v[855:860],'o',color='g')
#ax1.plot(t[964:970],vlamax_v[964:970],'o',color='r')
ax1.set_xlim(35,60)
plt.show()

t=np.arange(2000)*0.05
f,ax=plt.subplots(2,1,figsize=(12,15));ax0=ax[0];ax1=ax[1]
#cc=allmaps['aia171']['map171'][20]
#p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
xlaia=-949;xraia=-649.1;ylaia=70;yraia=369.4;p=mapex_bz[0].plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-750,vmax=750)
f.colorbar(p,label='B (G)')
ax0.plot(xcvmax[781],ycvmax[781],'o',color='g')
ax0.plot(xcvmax[858],ycvmax[858],'o',color='g')
ax0.plot(xcvmax[883],ycvmax[883],'o',color='g')
ax0.plot(xcvmax[954],ycvmax[954],'o',color='r')
ax0.plot(xcvmax[966],ycvmax[966],'o',color='r')
ax0.plot(xcvmax[1058],ycvmax[1058],'o',color='g')
ax0.plot(xcvmax[1153],ycvmax[1153],'o',color='r')
#ax0.plot(xcvmax[880:885],ycvmax[880:885],'o',color='g')
#ax0.plot(xcvmax[944:950],ycvmax[944:950],'o',color='r')
#ax0.plot(xcvmax[1057:1062],ycvmax[1057:1062],'o',color='g')
ax1.plot(t,vlamax_v,'o-',label='Stokes V');ax1.axvline(x=i*0.05,color='k');ax1.set_xlabel('Time (sec)');ax1.set_ylabel('T$_{B}$');ax1.legend(loc=4)
ax1.plot(t[781],vlamax_v[781],'o',color='g')
ax1.plot(t[858],vlamax_v[858],'o',color='g')
ax1.plot(t[883],vlamax_v[883],'o',color='g')
ax1.plot(t[954],vlamax_v[954],'o',color='r')
ax1.plot(t[966],vlamax_v[966],'o',color='r')
ax1.plot(t[1058],vlamax_v[1058],'o',color='g')
ax1.plot(t[1153],vlamax_v[1153],'o',color='r')
#ax1.plot(t[855:860],vlamax_v[855:860],'o',color='g')
#ax1.plot(t[964:970],vlamax_v[964:970],'o',color='r')
ax1.set_xlim(35,60)
plt.show()

vlalistLL=sorted(glob.glob('/media/rohit/VLA/20160409/vlamaps_LL/vla*.p'))
vlalistRR=sorted(glob.glob('/media/rohit/VLA/20160409/vlamaps_RR/vla*.p'))
re=[[124,134,120,130],[140,150,90,100],[95,115,70,90],[105,125,30,50],[95,115,95,115]]
TbLLr1=[0]*len(vlalistLL);TbLLr2=[0]*len(vlalistLL);TbLLr3=[0]*len(vlalistLL);TbLLr4=[0]*len(vlalistLL);TbLLr5=[0]*len(vlalistLL)
xcLLr1=[0]*len(vlalistLL);xcLLr2=[0]*len(vlalistLL);xcLLr3=[0]*len(vlalistLL);xcLLr4=[0]*len(vlalistLL);xcLLr5=[0]*len(vlalistLL)
ycLLr1=[0]*len(vlalistLL);ycLLr2=[0]*len(vlalistLL);ycLLr3=[0]*len(vlalistLL);ycLLr4=[0]*len(vlalistLL);ycLLr5=[0]*len(vlalistLL)
TbRRr1=[0]*len(vlalistRR);TbRRr2=[0]*len(vlalistRR);TbRRr3=[0]*len(vlalistRR);TbRRr4=[0]*len(vlalistRR);TbRRr5=[0]*len(vlalistRR)
xcRRr1=[0]*len(vlalistRR);xcRRr2=[0]*len(vlalistRR);xcRRr3=[0]*len(vlalistRR);xcRRr4=[0]*len(vlalistRR);xcRRr5=[0]*len(vlalistRR)
ycRRr1=[0]*len(vlalistRR);ycRRr2=[0]*len(vlalistRR);ycRRr3=[0]*len(vlalistRR);ycRRr4=[0]*len(vlalistRR);ycRRr5=[0]*len(vlalistRR)
xcRRmax=[0]*len(vlalistRR);ycRRmax=[0]*len(vlalistRR);xcLLmax=[0]*len(vlalistLL);ycLLmax=[0]*len(vlalistLL)
xcRR90=[0]*len(vlalistRR);ycRR90=[0]*len(vlalistRR);xcLL90=[0]*len(vlalistLL);ycLL90=[0]*len(vlalistLL)
timevla_all=[0]*len(vlalistLL)
freq=np.linspace(0.994,2.006,32)
v=0
for vfile in vlalistLL:
    d=pickle.load(open(vfile,'rb'))
    timevla_all[v]=d['vla']['timevla'][0]
    TbLLr1[v]=[0]*32;TbLLr2[v]=[0]*32;TbLLr3[v]=[0]*32;TbLLr4[v]=[0]*32;TbLLr5[v]=[0]*32;xcLLmax[v]=[0]*32;ycLLmax[v]=[0]*32;xcLL90[v]=[0]*32;ycLL90[v]=[0]*32
    print(vfile)
    for k in range(32):
        h=d['vla']['mapvla'][k][0]
        TbLLv=h.data
        TbLLv_f=TbLLv*1.0
        TbLLv_f[np.isnan(TbLLv_f)]=0
        bi=ut.get_bimage(TbLLv_f,0.95)
        xcf,ycf=ut.fitEllipse(bi)[0:2]
        xcLL90[v][k]=h.reference_coordinate.Tx.value+(xcf-(h.reference_pixel.x.value-1))*h.scale.axis1.value
        ycLL90[v][k]=h.reference_coordinate.Ty.value+(ycf-(h.reference_pixel.y.value-1))*h.scale.axis2.value
        yc,xc=np.where(TbLLv==np.nanmax(TbLLv))
        xcLLmax[v][k]=h.reference_coordinate.Tx.value+(xc-(h.reference_pixel.x.value-1))*h.scale.axis1.value
        ycLLmax[v][k]=h.reference_coordinate.Ty.value+(yc-(h.reference_pixel.y.value-1))*h.scale.axis2.value
        TbLLr1[v][k]=np.nanmax(TbLLv[re[0][0]:re[0][1],re[0][2]:re[0][3]])
        TbLLr2[v][k]=np.nanmax(TbLLv[re[1][0]:re[1][1],re[1][2]:re[1][3]])
        TbLLr3[v][k]=np.nanmax(TbLLv[re[2][0]:re[2][1],re[2][2]:re[2][3]])
        TbLLr4[v][k]=np.nanmax(TbLLv[re[3][0]:re[3][1],re[3][2]:re[3][3]])
        TbLLr5[v][k]=np.nanmax(TbLLv[re[4][0]:re[4][1],re[4][2]:re[4][3]])
    v=v+1


v=0
for vfile in vlalistRR:
    d=pickle.load(open(vfile,'rb'))
    TbRRr1[v]=[0]*32;TbRRr2[v]=[0]*32;TbRRr3[v]=[0]*32;TbRRr4[v]=[0]*32;TbRRr5[v]=[0]*32;xcRRmax[v]=[0]*32;ycRRmax[v]=[0]*32;xcRR90[v]=[0]*32;ycRR90[v]=[0]*32
    print(vfile)
    for k in range(32):
        h=d['vla']['mapvla'][k][0]
        TbRRv=h.data
        TbRRv_f=TbRRv*1.0
        TbRRv_f[np.isnan(TbRRv_f)]=0
        bi=ut.get_bimage(TbRRv_f,0.95)
        xcf,ycf=ut.fitEllipse(bi)[0:2]
        xcRR90[v][k]=h.reference_coordinate.Tx.value+(xcf-(h.reference_pixel.x.value-1))*h.scale.axis1.value
        ycRR90[v][k]=h.reference_coordinate.Ty.value+(ycf-(h.reference_pixel.y.value-1))*h.scale.axis2.value
        yc,xc=np.where(TbRRv==np.nanmax(TbRRv))
        xcRRmax[v][k]=h.reference_coordinate.Tx.value+(xc-(h.reference_pixel.x.value-1))*h.scale.axis1.value
        ycRRmax[v][k]=h.reference_coordinate.Ty.value+(yc-(h.reference_pixel.y.value-1))*h.scale.axis2.value
        TbRRr1[v][k]=np.nanmax(TbRRv[re[0][0]:re[0][1],re[0][2]:re[0][3]])
        TbRRr2[v][k]=np.nanmax(TbRRv[re[1][0]:re[1][1],re[1][2]:re[1][3]])
        TbRRr3[v][k]=np.nanmax(TbRRv[re[2][0]:re[2][1],re[2][2]:re[2][3]])
        TbRRr4[v][k]=np.nanmax(TbRRv[re[3][0]:re[3][1],re[3][2]:re[3][3]])
        TbRRr5[v][k]=np.nanmax(TbRRv[re[4][0]:re[4][1],re[4][2]:re[4][3]])
    v=v+1
TbLLr1=np.array(TbLLr1);TbLLr2=np.array(TbLLr2);TbLLr3=np.array(TbLLr3);TbLLr4=np.array(TbLLr4);TbLLr5=np.array(TbLLr5)
TbRRr1=np.array(TbRRr1);TbRRr2=np.array(TbRRr2);TbRRr3=np.array(TbRRr3);TbRRr4=np.array(TbRRr4);TbRRr5=np.array(TbRRr5)
xcRRmax=np.array(xcRRmax);ycRRmax=np.array(ycRRmax);xcLLmax=np.array(xcLLmax);ycLLmax=np.array(ycLLmax);xcLL90=np.array(xcLL90);ycLL90=np.array(ycLL90)
xcRR90=np.array(xcRR90);ycRR90=np.array(ycRR90)

pickle.dump([TbLLr1,TbLLr2,TbLLr3,TbLLr4,TbLLr5,xcLLmax,ycLLmax,xcLL90,ycLL90,timevla_all],open('MaxLL.p','wb'))
pickle.dump([TbRRr1,TbRRr2,TbRRr3,TbRRr4,TbRRr5,xcRRmax,ycRRmax,xcRR90,ycRR90],open('MaxRR.p','wb'))

#### READ

TbRRr1,TbRRr2,TbRRr3,TbRRr4,TbRRr5,xcRRmax,ycRRmax,xcRR90,ycRR90=pickle.load(open('/media/rohit/VLA/20160409/MaxRR.p','rb'),encoding='latin1')
TbLLr1,TbLLr2,TbLLr3,TbLLr4,TbLLr5,xcLLmax,ycLLmax,xcLL90,ycLL90,timevla_all=pickle.load(open('/media/rohit/VLA/20160409/MaxLL.p','rb'),encoding='latin1')
dcp_r1=(TbRRr1-TbLLr1)/(TbLLr1+TbRRr1)
ff=11

fig,ax=plt.subplots(1,1,sharex=True,figsize=(10,10))
ax1=ax.twinx()
ax1.plot(np.linspace(0,120,2399),TbRRr1[:,0]/1.e6,'o-',color='black',label='T$_B$')
ax.plot(np.linspace(0,120,2399),ycRR90.std(axis=1),'o-',color='red',label='Y-Centroid')
ax.set_xlabel('Time (sec)');ax.set_ylabel('Y-Centroid (arcsec)');ax1.set_ylabel('T$_B$ (MK)')
ax.legend();ax1.legend(loc=2)
plt.show()

fig,ax=plt.subplots(1,1,sharex=True,figsize=(10,10))
ax1=ax.twinx()
ax1.plot(np.linspace(0,120,2399),TbRRr1[:,0]/1.e6,'o-',color='black',label='T$_B$')
ax.plot(np.linspace(0,120,2399),xcRR90.mean(axis=1)-xcLL90.mean(axis=1),'o-',color='red',label='X-Centroid (RR-LL)')
ax.set_xlabel('Time (sec)');ax.set_ylabel('X-Centroid (arcsec)');ax1.set_ylabel('T$_B$ (MK)')
ax.legend();ax1.legend(loc=2)
plt.show()

x=np.linspace(0,120,2399)
for ff in range(32):
    loc_diff_LR=np.sqrt((np.concatenate(xcLLmax[:,ff])[:2399]-xcRRmax[:2399,ff,0])**2 + (np.concatenate(ycLLmax[:,ff])[:2399]-ycRRmax[:2399,ff,0])**2)
    fig,ax=plt.subplots(3,1,sharex=True,figsize=(10,13))
    ax[0].plot(x,(TbLLr1[:,ff]+TbRRr1[:,ff])/1.e6,'o-',label=str(np.round(freq[ff],3))+' GHz')
    ax[0].set_ylabel('Stokes I $\ T_B$ (MK)')
    ax[0].legend();ax[0].set_ylim(0,150)
    ax[1].plot(x,dcp_r1[:,ff],'o-')
    ax[1].set_ylabel('Stokes (V/I)')
    ax[1].set_ylim(0,1.0)
    ax[2].plot(x,loc_diff_LR,'o-')
    ax[2].set_ylabel('Position (arcsec)')
    ax[2].set_xlabel('Time (sec)')
    ax[2].set_ylim(0,50);ax[2].set_xlim(40,80)
    plt.title('Frequency: '+str(np.round(freq[ff],3))+' GHz')
    plt.savefig('pngs_dcp/dcp_'+"%02d"%ff+'_zoom.png')
    plt.close()
plt.show()


from mpl_toolkits.axes_grid1 import make_axes_locatable
fig, ax = plt.subplots(2, 2)
im1=ax[0,0].imshow(TbLLr1.swapaxes(0,1),aspect='auto',origin=0,vmin=1.e6,vmax=10.e6)
ax[1,0].imshow(TbLLr5.swapaxes(0,1),aspect='auto',origin=0,vmin=1.e6,vmax=10.e6)
ax[0,1].imshow(TbLLr3.swapaxes(0,1),aspect='auto',origin=0,vmin=1.e6,vmax=10.e6)
ax[1,1].imshow(TbLLr4.swapaxes(0,1),aspect='auto',origin=0,vmin=1.e6,vmax=10.e6)
ax[0,0].set_ylabel('Frequency (GHz)');ax[1,0].set_ylabel('Frequency (GHz)');ax[0,1].set_ylabel('Frequency (GHz)');ax[1,1].set_ylabel('Frequency (GHz)')
ax[0,0].set_xlabel('Time (sec)');ax[1,0].set_xlabel('Time (sec)');ax[0,1].set_xlabel('Time (sec)');ax[1,1].set_xlabel('Time (sec)')
ax[0,0].text(25,25,' Region 1',color='white')
ax[0,1].text(25,25,' Region 5',color='white')
ax[1,0].text(25,25,' Region 3',color='white')
ax[1,1].text(25,25,' Region 4',color='white')
ax[0,0].set_xticks([0,480,960,1440,1920,2400]);ax[0,0].set_xticklabels(['0','24','48','72','98','120'])
ax[0,1].set_xticks([0,480,960,1440,1920,2400]);ax[0,1].set_xticklabels(['0','24','48','72','98','120'])
ax[1,0].set_xticks([0,480,960,1440,1920,2400]);ax[1,0].set_xticklabels(['0','24','48','72','98','120'])
ax[1,1].set_xticks([0,480,960,1440,1920,2400]);ax[1,1].set_xticklabels(['0','24','48','72','98','120'])
ax[0,0].set_yticks(np.arange(32)[::4]);ax[0,0].set_yticklabels(np.round(np.linspace(0.994,2.006,32)[::4],2))
ax[1,0].set_yticks(np.arange(32)[::4]);ax[1,0].set_yticklabels(np.round(np.linspace(0.994,2.006,32)[::4],2))
ax[0,1].set_yticks(np.arange(32)[::4]);ax[0,1].set_yticklabels(np.round(np.linspace(0.994,2.006,32)[::4],2))
ax[1,1].set_yticks(np.arange(32)[::4]);ax[1,1].set_yticklabels(np.round(np.linspace(0.994,2.006,32)[::4],2))
divider = make_axes_locatable(ax[0,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im1, cax=cax, orientation='vertical')
plt.show()

plt.plot(TbLLr1[860],'o-')
plt.plot(TbRRr1[860],'o-')

Tb_r=[0]*5
for i in range(5):
    Tb_r[i]=[0]*2399
    for j in range(2399):
        Tb_r[i][j]=np.mean(Tb1[0][j][re[i][0]:re[i][1],re[i][2]:re[i][3]])
Tb_r=np.array(Tb_r)


for i in range(2399):
    dLL=pickle.load(open(vlalistLL[i],'rb'))
    vmLL=dLL['vla']['mapvla'][0][0]
    dRR=pickle.load(open(vlalistRR[i],'rb'))
    vmRR=dRR['vla']['mapvla'][0][0]
    f=plt.figure(figsize=(6,10))
    ax0 = f.add_subplot(211)
    tidx_171=ut.find_predecessor(allmaps['aia171']['time171'],timevla_all[i])[0]
    cc=allmaps['aia171']['map171'][tidx_171]
    xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
    p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    #lev1=(1.5e7/dd0.data.max())*np.array([20,30,40,50,60,70,80,90])*u.percent
    lev0=np.array([40,50,60,70,90])*u.percent
    dd0=vmLL;dd0.data[np.isnan(dd0.data)]=0
    xlvla=dd0.center.Tx.value-2.0*int(dd0.data.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(dd0.data.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(dd0.data.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(dd0.data.shape[0]/2)
    dd0.draw_contours(levels=lev0,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    lev1=np.array([40,50,60,70,90])*u.percent
    dd1=vmRR;dd1.data[np.isnan(dd1.data)]=0
    xlvla=dd1.center.Tx.value-2.0*int(dd1.data.shape[0]/2);xrvla=dd1.center.Tx.value+2.0*int(dd1.data.shape[0]/2);ylvla=dd1.center.Ty.value-2.0*int(dd1.data.shape[1]/2);yrvla=dd1.center.Ty.value+2.0*int(dd1.data.shape[0]/2)
    dd1.draw_contours(levels=lev1,colors='b',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    #ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
    #ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
    ax0.set_xlim([-1000,-700]);ax0.set_ylim([50,400])
    ax0.text(-1000,60,str(np.round(dd0.meta['crval3']/1.e9,4))+' GHz',color='white')
    ax0.add_patch(patches.Rectangle((-781,233),40,40,linewidth=5,edgecolor='cyan',facecolor='none'))
    ax0.add_patch(patches.Rectangle((-841,265),20,20,linewidth=5,edgecolor='g',facecolor='none'))
    ax0.add_patch(patches.Rectangle((-881,175),40,40,linewidth=5,edgecolor='red',facecolor='none'))
    ax0.add_patch(patches.Rectangle((-961,195),40,40,linewidth=5,edgecolor='b',facecolor='none'))
    ax0.add_patch(patches.Rectangle((-831,175),40,40,linewidth=5,edgecolor='magenta',facecolor='none'))
    ax1 = f.add_subplot(212)
    ax1.plot(freq,TbLLr1[i]/1.e6,'o-',markersize=2,label='Region 1 (LL)',color='cyan')
    ax1.plot(freq,TbLLr2[i]/1.e6,'o-',markersize=2,label='Region 2 (LL)',color='green')
    ax1.plot(freq,TbLLr3[i]/1.e6,'o-',markersize=2,label='Region 3 (LL)',color='red')
    ax1.plot(freq,TbLLr4[i]/1.e6,'o-',markersize=2,label='Region 4 (LL)',color='blue')
    ax1.plot(freq,TbLLr5[i]/1.e6,'o-',markersize=2,label='Region 5 (LL)',color='magenta')
    ax1.plot(freq,TbRRr1[i]/1.e6,'o--',markersize=2,label='Region 1 (RR)',color='cyan')
    ax1.plot(freq,TbRRr2[i]/1.e6,'o--',markersize=2,label='Region 2 (RR)',color='green')
    ax1.plot(freq,TbRRr3[i]/1.e6,'o--',markersize=2,label='Region 3 (RR)',color='red')
    ax1.plot(freq,TbRRr4[i]/1.e6,'o--',markersize=2,label='Region 4 (RR)',color='blue')
    ax1.plot(freq,TbRRr5[i]/1.e6,'o--',markersize=2,label='Region 5 (RR)',color='magenta')
    ax1.legend(fontsize=10)#;ax1.set_ylim(0,100)
    ax1.set_yscale('log');ax1.set_ylim(0,300);ax1.grid(True)
    ax1.set_xlabel('Frequency (GHz)');ax1.set_ylabel('$T_B$ (MK)')
    plt.savefig('pngs_spec_movie/aia171_'+"%02d"%i+'.png')
    plt.close()


#for i in range(2399):
for i in range(1):
    i=860
    dLL=pickle.load(open(vlalistLL[i],'rb'))
    vmLL=dLL['vla']['mapvla'][0][0]
    dRR=pickle.load(open(vlalistRR[i],'rb'))
    vmRR=dRR['vla']['mapvla'][0][0]
    f=plt.figure(figsize=(6,10))
    ax0 = f.add_subplot(211)
    tidx_171=ut.find_predecessor(allmaps['aia171']['time171'],timevla_all[i])[0]
    cc=allmaps['aia171']['map171'][tidx_171]
    xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
    p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    #lev1=(1.5e7/dd0.data.max())*np.array([20,30,40,50,60,70,80,90])*u.percent
    lev0=np.array([40,50,60,70,90])*u.percent
    dd0=vmLL;dd0.data[np.isnan(dd0.data)]=0
    xlvla=dd0.center.Tx.value-2.0*int(dd0.data.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(dd0.data.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(dd0.data.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(dd0.data.shape[0]/2)
    #dd0.draw_contours(levels=lev0,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    lev1=np.array([40,50,60,70,90])*u.percent
    dd1=vmRR;dd1.data[np.isnan(dd1.data)]=0
    xlvla=dd1.center.Tx.value-2.0*int(dd1.data.shape[0]/2);xrvla=dd1.center.Tx.value+2.0*int(dd1.data.shape[0]/2);ylvla=dd1.center.Ty.value-2.0*int(dd1.data.shape[1]/2);yrvla=dd1.center.Ty.value+2.0*int(dd1.data.shape[0]/2)
    dd1.draw_contours(levels=lev1,colors='white',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    #ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
    #ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
    ax0.set_xlim([-900,-700]);ax0.set_ylim([150,350])
    ax0.text(-1000,60,str(np.round(dd0.meta['crval3']/1.e9,4))+' GHz',color='white')
    #ax0.add_patch(patches.Rectangle((-781,233),40,40,linewidth=5,edgecolor='cyan',facecolor='none'))
    #ax0.add_patch(patches.Rectangle((-841,265),20,20,linewidth=5,edgecolor='g',facecolor='none'))
    #ax0.add_patch(patches.Rectangle((-881,175),40,40,linewidth=5,edgecolor='red',facecolor='none'))
    #ax0.add_patch(patches.Rectangle((-961,195),40,40,linewidth=5,edgecolor='b',facecolor='none'))
    #ax0.add_patch(patches.Rectangle((-831,175),40,40,linewidth=5,edgecolor='magenta',facecolor='none'))
    ax1 = f.add_subplot(212)
    ax1.plot(Tbi_r1[0:8].mean(axis=0)/1.e6,'o-',markersize=2,label='Region 1 (Stokes-I)',color='k')
    ax1.legend(fontsize=10)#;ax1.set_ylim(0,100)
    ax1.set_yscale('linear');ax1.set_ylim(0,50);ax1.grid(True)
    ax1.set_xticks([0,600,1200,1800]);ax1.set_xticklabels(['18:44:00','18:44:30','18:45:00','18:45:30'])
    ax1.set_xlabel('Time (HH:MM:SS UT)');ax1.set_ylabel('$T_B$ (MK)')
    #plt.savefig('pngs_spec_movie/aia171_'+"%02d"%i+'.png')
    plt.show()
    plt.figure(figsize=(5,5))
    plt.imshow(Tbi_r1[0:16]/1.e6,aspect='auto',origin=0,cmap='YlOrRd')
    plt.colorbar(label='$T_B$ (MK)')
    plt.xticks([0,600,1200,1800],['18:44:00','18:44:30','18:45:00','18:45:30'])
    plt.yticks([0,4,8,12,15],[freq[0],freq[4],freq[8],freq[12],freq[15]])
    plt.xlabel('Time (HH:MM:SS UT)');plt.ylabel('Frequency (GHz)')
    plt.show()


sys.exit()

#################
vlafileLL='/media/rohit/VLA/20160409/vlamaps_LL/vlamap_18:44:43_0860.p'
vlafileRR='/media/rohit/VLA/20160409/vlamaps_RR/vlamap_18:44:43_0860.p'
vlaLL=pickle.load(open(vlafileLL,'rb'),encoding='latin1');vlaRR=pickle.load(open(vlafileRR,'rb'),encoding='latin1')
tidx_94=ut.find_predecessor(allmaps['aia94']['time94'],vlaLL['vla']['timevla'][0])[0]
tidx_131=ut.find_predecessor(allmaps['aia131']['time131'],vlaLL['vla']['timevla'][0])[0]
tidx_171=ut.find_predecessor(allmaps['aia171']['time171'],vlaLL['vla']['timevla'][0])[0]
tidx_193=ut.find_predecessor(allmaps['aia193']['time193'],vlaLL['vla']['timevla'][0])[0]
tidx_1600=ut.find_predecessor(allmaps['aia1600']['time1600'],vlaLL['vla']['timevla'][0])[0]


for i in range(32):
    f=plt.figure(figsize=(6,6))
    ax0 = f.add_subplot(111)
    cc=allmaps['aia171']['map171'][tidx_171]
    xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
    p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    #lev1=(1.5e7/dd0.data.max())*np.array([20,30,40,50,60,70,80,90])*u.percent
    lev0=np.array([40,50,60,70,90])*u.percent
    dd0=vlaLL['vla']['mapvla'][i][0];dd0.data[np.isnan(dd0.data)]=0
    xlvla=dd0.center.Tx.value-2.0*int(dd0.data.shape[0]/2);xrvla=dd0.center.Tx.value+2.0*int(dd0.data.shape[0]/2);ylvla=dd0.center.Ty.value-2.0*int(dd0.data.shape[1]/2);yrvla=dd0.center.Ty.value+2.0*int(dd0.data.shape[0]/2)
    dd0.draw_contours(levels=lev0,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    lev1=np.array([40,50,60,70,90])*u.percent
    dd1=vlaRR['vla']['mapvla'][i][0];dd1.data[np.isnan(dd1.data)]=0
    xlvla=dd1.center.Tx.value-2.0*int(dd1.data.shape[0]/2);xrvla=dd1.center.Tx.value+2.0*int(dd1.data.shape[0]/2);ylvla=dd1.center.Ty.value-2.0*int(dd1.data.shape[1]/2);yrvla=dd1.center.Ty.value+2.0*int(dd1.data.shape[0]/2)
    dd1.draw_contours(levels=lev1,colors='b',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    #ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
    ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
    ax0.text(-1200,0,str(np.round(dd0.meta['crval3']/1.e9,4))+' GHz',color='green')
    plt.savefig('pngs_spec/aia171_'+"%02d"%i+'.png')
    plt.close()
    #ax0.text(-1200,50,'Contours: 3, 4.5, 6, 7.5, 9, 10, 12, 13 MK',color='yellow')
    #ax0.text(-1200,50,'Contours: 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%',color='yellow')

################ Make Time String ##########
listvla=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_LL/spw_0/sun*.spw.0_16-31*FITS'))
list1600=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/1600/*fits'))
list171=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/full_sun/171/*fits'))
timstr_vla=[0]*len(listvla);timstr_1600=[0]*len(list1600);timstr_171=[0]*len(list171)
for i in range(len(listvla)):
    timstr_vla[i]=listvla[i].split('time.')[1].split('.FITS')[0].split('-')[0]
for i in range(len(list1600)):
    timstr_1600[i]= list1600[i].split('T')[-1].split('Z')[0]
for i in range(len(list171)):
    timstr_171[i]= list171[i].split('T')[-1].split('Z')[0]
    


################ Read submaps ##############
allmaps=pickle.load(open('/media/rohit/VLA/20160409/20160409_submap_aia_50ms.p','rb'),encoding='latin-1')
vla0=pickle.load(open('/media/rohit/VLA/20160409/20160409_vla_spw_0_50ms.p','rb'),encoding='latin-1')
#vla1=pickle.load(open('/media/rohit/VLA/20160409/20160409_vla_spw_1_50ms.p','rb'))
#vla2=pickle.load(open('/media/rohit/VLA/20160409/20160409_vla_spw_2_50ms.p','rb'))
#vla3=pickle.load(open('/media/rohit/VLA/20160409/20160409_vla_spw_3_50ms.p','rb'))
print("Reading 5 and 7..")
#vla4=pickle.load(open('/media/rohit/VLA/20160409/20160409_vla_spw_4_50ms.p','rb'))
#vla5=pickle.load(open('/medi8a/rohit/VLA/20160409/20160409_vla_spw_5_50ms.p','rb'))
#vla6=pickle.load(open('/media/rohit/VLA/20160409/20160409_vla_spw_6_50ms.p','rb'))
#vla7=pickle.load(open('/media/rohit/VLA/20160409/20160409_vla_spw_7_50ms.p','rb'))

timevla=vla0['vla0']['timevla']

#########################
time94=allmaps['aia94']['time94'];time131=allmaps['aia131']['time131'];time335=allmaps['aia335']['time335']
time1600=allmaps['aia1600']['time1600'];time1700=allmaps['aia1700']['time1700'];time171=allmaps['aia171']['time171']
timed94=allmaps['aiad94']['timed94'];timed131=allmaps['aiad131']['timed131'];time335=allmaps['aiad335']['timed335']
timed1600=allmaps['aiad1600']['timed1600'];timed1700=allmaps['aiad1700']['timed1700'];time171=allmaps['aiad171']['timed171']
freq=np.linspace(0.997,1.245,32)
################# Analysis #################
#Tb=np.array(allmaps['vla']['datavla']).reshape(32,119,150,200)
Tb0=np.array(vla0['vla0']['datavla']).reshape(1,2399,256,256)
#Tb1=np.array(vla1['vla1']['datavla']).reshape(1,2399,256,256)
#Tb2=np.array(vla2['vla2']['datavla']).reshape(1,2399,256,256)
#Tb3=np.array(vla3['vla3']['datavla']).reshape(1,2399,256,256)
#Tb4=np.array(vla4['vla4']['datavla']).reshape(1,2399,256,256)
#Tb5=np.array(vla5['vla5']['datavla']).reshape(1,2399,256,256)
#Tb6=np.array(vla6['vla6']['datavla']).reshape(1,2399,256,256)
Tb7=np.array(vla7['vla7']['datavla']).reshape(1,2399,256,256)
Tbmax_ds=Tb0.max(axis=(2,3))
#Tb_mean_r1=Tb[:,:,30:80,115:145].mean(axis=(2,3))
Tb_mean_r1=Tb0[:,:,115:145,115:145].mean(axis=(2,3))
m=len(timevla_all1)
tidx94=[0]*m;tidx131=[0]*m;tidx171=[0]*m
tidx335=[0]*m;tidx1600=[0]*m;tidx1700=[0]*m
for i in range(m):
    tidx94[i]=ut.find_predecessor(allmaps['aia94']['time94'],timevla_all1[i])[0]
    tidx131[i]=ut.find_predecessor(allmaps['aia131']['time131'],timevla_all1[i])[0]
    tidx171[i]=ut.find_predecessor(allmaps['aia171']['time171'],timevla_all1[i])[0]
    tidx335[i]=ut.find_predecessor(allmaps['aia335']['time335'],timevla_all1[i])[0]
    tidx1600[i]=ut.find_predecessor(allmaps['aia1600']['time1600'],timevla_all1[i])[0]
    tidx1700[i]=ut.find_predecessor(allmaps['aia1700']['time1700'],timevla_all1[i])[0]
    #tidx94[i]=ut.find_predecessor(allmaps['aia94']['time94'],vla0['vla0']['timevla'][i])[0]
    #tidx131[i]=ut.find_predecessor(allmaps['aia131']['time131'],vla0['vla0']['timevla'][i])[0]
    #tidx171[i]=ut.find_predecessor(allmaps['aia171']['time171'],vla0['vla0']['timevla'][i])[0]
    #tidx335[i]=ut.find_predecessor(allmaps['aia335']['time335'],vla0['vla0']['timevla'][i])[0]
    #tidx1600[i]=ut.find_predecessor(allmaps['aia1600']['time1600'],vla0['vla0']['timevla'][i])[0]
    #tidx1700[i]=ut.find_predecessor(allmaps['aia1700']['time1700'],vla0['vla0']['timevla'][i])[0]

#vlafreqidx=np.array(list(np.arange(32))*119).reshape(119,32).swapaxes(0,1).flatten()
vlafreqidx=np.array(list(np.arange(1))*2399).reshape(2399,1).swapaxes(0,1).flatten()

lev=[0.5,0.6,0.7,0.9]
vlasize=[0]*len(lev)
vlasize_fwhm=[0]*len(lev)
for k in range(len(lev)):
    vlasize[k]=[0]*Tb0.shape[0]
    vlasize_fwhm[k]=[0]*Tb0.shape[0]
    for i in range(Tb0.shape[0]):
        vlasize[k][i]=[0]*Tb0.shape[1]
        vlasize_fwhm[k][i]=[0]*Tb0.shape[1]
        for j in range(Tb0.shape[1]):
            bimage=ut.get_bimage(Tb0[i][j],lev[k]*3.e7/np.nanmax(Tb0[i][j]))
            bimage_fwhm=ut.get_bimage(Tb0[i][j],lev[k])
            vlasize[k][i][j]=bimage[np.isfinite(bimage)].sum()*2.0
            vlasize_fwhm[k][i][j]=bimage_fwhm[np.isfinite(bimage_fwhm)].sum()*2.0
vlasize=np.array(vlasize)
vlasize_fwhm=np.array(vlasize_fwhm)
vlasize_fwhm[np.isnan(vlasize_fwhm)]=0


re=[[124,134,120,130],[140,150,90,100],[95,115,70,90],[105,125,30,50],[95,115,95,115]]

Tb_r=[0]*5
for i in range(5):
    Tb_r[i]=[0]*2000
    for j in range(2000):
        #Tb_r[i][j]=np.mean(Tb1[0][j][re[i][0]:re[i][1],re[i][2]:re[i][3]])
        Tb_r[i][j]=np.mean(datavla_v[j][0][re[i][0]:re[i][1],re[i][2]:re[i][3]])
Tb_r=np.array(Tb_r)

def get_ts(d,xl,xr,yl,yr):
    ts=[0]*len(d)
    for i in range(len(d)):
        #ts[i]=np.nanmax(d[i][xl:xr,yl:yr])
        ts[i]=np.nanmean(d[i][xl:xr,yl:yr])
    ts=np.array(ts)/np.nanmax(np.array(ts))
    return ts


data94=allmaps['aia94']['data94']
data131=allmaps['aia131']['data131']
data171=allmaps['aia171']['data171']
data335=allmaps['aia335']['data335']
data1600=allmaps['aia1600']['data1600']
data1700=allmaps['aia1700']['data1700']
#Timeseries
ts94=get_ts(data94,160,190,300,330) # 300-330,160-190
ts131=get_ts(data131,160,190,300,330)
ts335=get_ts(data335,160,190,300,330)
ts1600=get_ts(data1600,160,190,300,330)
ts1700=get_ts(data1700,160,190,300,330)
time94=allmaps['aia94']['time94'];time131=allmaps['aia131']['time131'];time335=allmaps['aia335']['time335']
time1600=allmaps['aia1600']['time1600'];time1700=allmaps['aia1700']['time1700'];time171=allmaps['aia171']['time171']

#############################################

tidx94=[0]*2000;tidx131=[0]*2000;tidx171=[0]*2000;tidx193=[0]*2000;tidx335=[0]*2000;tidx1600=[0]*2000;tidx1700=[0]*2000
for i in range(2000):
    tidx94[i]=ut.find_predecessor(allmaps['aia94']['time94'],timevla[i][0])[0]
    tidx131[i]=ut.find_predecessor(allmaps['aia131']['time131'],timevla[i][0])[0]
    tidx171[i]=ut.find_predecessor(allmaps['aia171']['time171'],timevla[i][0])[0]
    tidx335[i]=ut.find_predecessor(allmaps['aia335']['time335'],timevla[i][0])[0]
    tidx1600[i]=ut.find_predecessor(allmaps['aia1600']['time1600'],timevla[i][0])[0]
    tidx1700[i]=ut.find_predecessor(allmaps['aia1700']['time1700'],timevla[i][0])[0]


################ FERMI & GOES & VLA DS #######################
goes=readsav('/data/Dropbox/20160409/fermi/idlsave_goes.sav')
gf0540=goes['lx'][1];gf1080=goes['lx'][0];gtime=goes['tarray']


fm=fits.open('/data/Dropbox/20160409/fermi/glg_cspec_n5_160409_v00_data.fits')
fm0=fm[0];fm1=fm[1]
fmrate=fm1.data['RATE'];fmtime=fm1.data['TIME']
fmtime=fmtime-fmtime[0]
plt.plot(fmtime,fmrate,'o-')
plt.show()



############### MAKE DS #####################
ds1=np.load('/media/rohit/VLA/20160409/sun_L_20160409T184000-184400UT.50ms.cal.ms.dspec.npz')
specfile='/media/rohit/VLA/20160409/sun_L_20160409T1844-1846UT.50ms.cal.ms.dspec.median.npz'
data=np.load(specfile)
tim_=data['tim'];freq=data['freq'];ds_=data['spec']
tim=[0]*8
for i in range(8):
    tim[i]=np.hstack((tim_[2089*i:2089*(i+1)],tim_[16712+i*311:16712+(i+1)*311]))
tim=np.array(tim)
ds_=ds_[:,0,:,:]
ds=[0]*2
for i in range(2):
    ds[i]=[0]*64
    for j in range(64):
        ds[i][j]=[0]*8
        for k in range(8):
            ds[i][j][k]=np.hstack((ds_[i,j,2089*k:2089*(k+1)],ds_[i,j,16712+k*311:16712+(k+1)*311]))

ds=np.array(ds)
ds_LL=ds[0].swapaxes(0,1).reshape(512,2400)
ds_RR=ds[1].swapaxes(0,1).reshape(512,2400)
freq=freq.swapaxes(0,1).flatten()/1.e9
ds1s=np.subtract(ds1['spec'][0][0].swapaxes(0,1),ds1['spec'][0][0][:,0]).swapaxes(0,1)/3 
ds_LL1=np.hstack((ds1s,np.subtract(ds1s[:,-1],-1*ds_LL[6:512].swapaxes(0,1)+28).swapaxes(0,1)))
ds_RR1=np.hstack((ds1s,np.subtract(ds1s[:,-1],-1*ds_RR[6:512].swapaxes(0,1)+28).swapaxes(0,1)))


sys.exit()
################ PLOTTING ###################

plot_euv_time=1
if(plot_euv_time):
    dsplt=ds_LL1[0:128].mean(axis=0)/120; dsplt[2130:2170]=np.nan; dsplt[4530:4570]=np.nan; dsplt[6930:6970]=np.nan
    f,ax=plt.subplots(3,1,sharex=True,gridspec_kw={'height_ratios': [1, 3, 1]})
    im0=ax[0].imshow(ds_LL1/300*7,aspect='auto',origin='lower',cmap='coolwarm',vmin=0.1,vmax=2,extent=[67200.0,67560.,freq[0],freq[-1]])
    ax[0].set_ylabel('Frequency (GHz)')
    divider = make_axes_locatable(ax[0]);cax = divider.append_axes('right', size='0.5%', pad=0.01)
    f.colorbar(im0, cax=cax, orientation='vertical',label='(SFU)')
    ax[1].axhline(y=0.5,linestyle='--',color='gray')
    ax[1].axvline(x=67315, linestyle='--', color='gray')
    ax[1].axvline(x=67430, linestyle='--', color='gray')
    ax[1].axvline(x=67443, linestyle='--', color='gray')
    ax[1].axvline(x=67484, linestyle='--', color='gray')
    ax[1].text(67315,2.4,"18:41:55",color='blue');ax[1].text(67415,2.5,"18:43:50",color='blue');ax[1].text(67443,2.5,"18:44:03",color='blue');ax[1].text(67484,2.5,"18:44:44",color='blue')
    ax[1].plot(time94,ts94-ts94[1]+0.5,'o-',markersize=2,label='AIA 94 $\AA$')
    ax[1].plot(time131,ts131-ts131[1]+0.5,'o-',markersize=2,label='AIA 131 $\AA$')
    ax[1].plot(time335,ts335-ts335[1]+0.5,'o-',markersize=2,label='AIA 335 $\AA$')
    ax[1].plot(time1600,ts1600-ts1600[1]+0.5,'o-',markersize=2,label='AIA 1600 $\AA$')
    ax[1].plot(qstimevla[0][0]+np.arange(7200)*0.05,dsplt,'-',color='k',label='0.99-1.25 GHz')
    ax[1].plot(fmtime[14244:14332],fmrate[14244:14332:,0:20].mean(axis=1)/200,'o-',markersize=2,color='orange',label='FERMI (10-25 keV)')
    #ax1.plot(np.hstack((np.array(qstimevla)[:,0],np.array(timevla_all)[0:2000])),np.hstack((np.array(qsTbr_r1),Tbr_r1[1]))/1.e8,'-',color='k',label='1.077 GHz')
    ax[1].set_xticks(np.arange(11)*60+67212-11)
    ax[1].set_xticklabels(['18:40','18:41','18:42','18:43','18:44','18:45','18:46','18:47','18:48','18:49','18:50'])
    #plt.ylabel('T$_{B}$ ($\\times10^{5}$ (K))');plt.xlabel('Time (HH:MM UT)')
    ax[1].legend(loc=2,prop={'size':15})
    ax[2].plot(gtime[292:585]-600+67200,gf0540[292:585]*1.e-23,'o-',color='g',label='0.5-4.0$\AA$')
    ax2=ax[2].twinx()
    ax2.plot(gtime[292:585]-600+67200,gf1080[292:585]*1.e-23,'o-',label='1.0-8.0$\AA$');ax[2].legend(loc=4);ax2.legend(loc=2)
    ax[1].set_ylabel('Amplitude');ax[2].set_xlabel('Time (HH:MM UT)')#;ax1.set_ylabel('')
    ax[2].set_xlim([67200.,67560.]);ax[2].set_ylabel('Flux ($\\times10^{-23}$W/m$^2$)');ax2.set_ylabel('Flux ($\\times10^{-23}$W/m$^2$)')
    plt.show()

f,ax=plt.subplots(1,1)
ax.plot(gtime-600+67200,gf0540*1.e-23,'o-',color='g',label='0.5-4.0$\AA$')
ax1=ax.twinx()
ax1.plot(gtime-600+67200,gf1080*1.e-23,'o-',color='b',label='1.0-8.0$\AA$')
#ax.plot(gtime[292:585]-600+67200,gf0540[292:585],'o-',color='g',label='0.5-4.0$\AA$')
ax.axvline(x=67200,linestyle='--',color='k')
ax.axvline(x=67798,linestyle='--',color='k')
ax.set_xticks(np.arange(11)*180+66601)
ax.set_xticklabels(['18:30','18:33','18:36','18:39','18:42','18:45','18:48','18:51','18:54','18:57','19:00'])
ax.set_ylabel('Flux ($\\times10^{-23}$W/m$^2$)');ax1.set_ylabel('Flux ($\\times10^{-23}$W/m$^2$)');ax.set_xlabel('Time (HH:MM UT)')
ax.legend(loc=4);ax1.legend(loc=2)
plt.show()

listvla_rr=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_0/*spw.0_16-31*FITS'));v=Map(listvla_rr[860]);v.data[np.isnan(v.data)]=0
listvla_rr=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_2/*spw.2_16-31*FITS'));v2=Map(listvla_rr[860]);v2.data[np.isnan(v2.data)]=0
listvla_rr=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_4/*spw.4_16-31*FITS'));v4=Map(listvla_rr[860]);v4.data[np.isnan(v4.data)]=0
listvla_rr=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_6/*spw.6_16-31*FITS'));v6=Map(listvla_rr[860]);v6.data[np.isnan(v6.data)]=0
f=plt.figure(figsize=(10,10))
ax0 = f.add_subplot(111)
cc=allmaps['aia171']['map171'][10]
dd=v;dd.data[np.isnan(dd.data)]=0
xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
#lev1=(1.5e7/dd.data.max())*np.array([60,70,80,90])*u.percent
lev1=np.array([60,70,80,90])*u.percent
xlvla=dd.center.Tx.value-2.0*int(dd.data.shape[0]/2);xrvla=dd.center.Tx.value+2.0*int(dd.data.shape[0]/2);ylvla=dd.center.Ty.value-2.0*int(dd.data.shape[1]/2);yrvla=dd.center.Ty.value+2.0*int(dd.data.shape[0]/2)
dd.draw_contours(levels=lev1,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
ax0.set_title('AIA 171 $\AA$:18:42:10 UT VLA: 18:44:43.00-18:44:43.05 UT')
ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
#ax0.text(-1200,0,str(np.round(dd.meta['crval3']/1.e9,4))+' GHz',color='r')
#ax0.text(-1200,50,'Contours: 3, 4.5, 6, 7.5, 9, 10, 12, 13 MK',color='yellow')
ax0.set_xlim([-900,-700]),ax0.set_ylim([100,300])
plt.show()


f=plt.figure(figsize=(10,10))
ax0 = f.add_subplot(111)
cc=allmaps['aia171']['map171'][10]
dd=v;dd.data[np.isnan(dd.data)]=0
xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
ax0.plot(b_hp_euv[0],b_hp_euv[1],'o',markersize=2,color='cyan')
ax0.plot(b_hp_fr[0],b_hp_fr[1],'o',markersize=2,color='magenta')
ax0.plot(b_hp_ls[0],b_hp_ls[1],'o',markersize=2,color='orange')
#lev1=(1.5e7/dd.data.max())*np.array([60,70,80,90])*u.percent
lev1=np.array([60,70,80,90])*u.percent
xlvla=dd.center.Tx.value-2.0*int(dd.data.shape[0]/2);xrvla=dd.center.Tx.value+2.0*int(dd.data.shape[0]/2);ylvla=dd.center.Ty.value-2.0*int(dd.data.shape[1]/2);yrvla=dd.center.Ty.value+2.0*int(dd.data.shape[0]/2)
dd.draw_contours(levels=lev1,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
ax0.set_title('AIA 171 $\AA$:18:42:10 UT VLA: 18:44:43.00-18:44:43.05 UT')
#ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
#ax0.plot(b_hp_fr.Tx,b_hp_fr.Ty,'o',markersize=3)
#ax0.text(-1200,0,str(np.round(dd.meta['crval3']/1.e9,4))+' GHz',color='r')
#ax0.text(-1200,50,'Contours: 3, 4.5, 6, 7.5, 9, 10, 12, 13 MK',color='yellow')
ax0.set_xlim([-900,-700]),ax0.set_ylim([100,300])
plt.show()


from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Ellipse
ncolors = 256
color_array = plt.get_cmap('gist_rainbow')(range(ncolors))
color_array[:,-1] = np.linspace(0.0,1.0,ncolors)
map_object = LinearSegmentedColormap.from_list(name='rainbow_alpha',colors=color_array)
plt.register_cmap(cmap=map_object)

cc=maphmi[-1]
dd=v;dd.data[np.isnan(dd.data)]=0
ee=allmaps['aia171']['map171'][10]
f=plt.figure(figsize=(10,10))
ax0 = f.add_subplot(111,projection=cc)
xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
p=cc.plot(axes=ax0,aspect='auto')#,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
#lev1=(1.5e7/dd.data.max())*np.array([60,70,80,90])*u.percent
lev1=np.array([70,75,85,95])*u.percent
xlvla=dd.center.Tx.value-2.0*int(dd.data.shape[0]/2);xrvla=dd.center.Tx.value+2.0*int(dd.data.shape[0]/2);ylvla=dd.center.Ty.value-2.0*int(dd.data.shape[1]/2);yrvla=dd.center.Ty.value+2.0*int(dd.data.shape[0]/2)
dd.draw_contours(levels=lev1,colors=color_array[0][0:3],linewidths=2)#,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
v2.draw_contours(levels=lev1,colors=color_array[30][0:3],linewidths=2)#,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
v4.draw_contours(levels=lev1,colors=color_array[40][0:3],linewidths=2)#,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
v6.draw_contours(levels=lev1,colors=color_array[50][0:3],linewidths=2)#,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
p=ee.plot(axes=ax0,aspect='auto',cmap='rainbow_alpha',vmin=100,vmax=8000,alpha=0.7)#,extent=[xlaia,xraia,ylaia,yraia])
ax0.set_title('HMI: 18:50:47 UT VLA: 18:44:43.00-18:44:43.05 UT')
ax0.text(400,100,str(np.round(dd.meta['restfrq']/1.e9,3))+' GHz',fontsize=18,color=color_array[0][0:3])
ax0.text(400,120,str(np.round(v2.meta['restfrq']/1.e9,3))+' GHz',fontsize=18,color=color_array[30][0:3])
ax0.text(400,140,str(np.round(v4.meta['restfrq']/1.e9,3))+' GHz',fontsize=18,color=color_array[40][0:3])
ax0.text(400,160,str(np.round(v6.meta['restfrq']/1.e9,3))+' GHz',fontsize=18,color=color_array[50][0:3])
beam = Ellipse(xy=(340, 80), width=24, height=40, edgecolor='blue', fc='None', lw=2)
ax0.add_patch(beam)
#ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
#ax0.text(-1200,0,str(np.round(dd.meta['crval3']/1.e9,4))+' GHz',color='r')
#ax0.text(-1200,50,'Contours: 3, 4.5, 6, 7.5, 9, 10, 12, 13 MK',color='yellow')
#ax0.set_xlim([-900,-700]),ax0.set_ylim([100,300])
plt.show()


plot_vlasize=1
if(plot_vlasize):
    fig,ax=plt.subplots(3,1,sharex=True)
    im0=ax[0].imshow(Tb_mean_r1/1.e6,origin=0,cmap='jet',interpolation='None')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im0,cax=cax,label='T$_{B}$ (MK)')
    ax[0].set_yticks([0,10,20,30]);ax[0].set_yticklabels([freq[0],freq[10],freq[20],freq[30]])
    ax[1].plot(Tb_mean_r1[10]/1.e6,'o-',label='1.077 GHz');ax[1].legend()
    ax[2].plot(vlasize[0,10],'o-',label='15 MK (1.077 GHz)')
    #ax[2].plot(vlasize[1,10],'o-',label='18 MK')
    ax[2].plot(vlasize[2,10],'o-',label='21 MK')
    ax[2].plot(vlasize[3,10],'o-',label='27 MK')
    ax[2].legend()
    ax[2].set_xticks([0,30,60,90,120]);ax[2].set_xticklabels(['18:44:00','18:44:30','18:45:00','18:45:30','18:46:00'])
    ax[2].set_xlabel('Time (HH:MM:SS UT)')
    ax[0].set_ylabel('Frequency (GHz)');ax[1].set_ylabel('T$_{B}$ (MK)');ax[2].set_ylabel('Source size (arcsec$^{2}$)')
    plt.show()

plot_vlasize_fwhm=1
if(plot_vlasize_fwhm):
    fig,ax=plt.subplots(3,1,sharex=True)
    im0=ax[0].imshow(Tb_mean_r1/1.e6,origin=0,cmap='jet',interpolation='None')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im0,cax=cax,label='T$_{B}$ (MK)')
    ax[0].set_yticks([0,10,20,30]);ax[0].set_yticklabels([freq[0],freq[10],freq[20],freq[30]])
    ax[1].plot(Tb_mean_r1[10]/1.e6,'o-',label='1.077 GHz');ax[1].legend()
    ax[2].plot(vlasize_fwhm[0,10],'o',label='50% Contours')
    ax[2].plot(vlasize_fwhm[2,10],'o',label='70% Contours')
    ax[2].plot(vlasize_fwhm[3,10],'o',label='90% Contours')
    ax[2].legend()
    ax[2].set_xticks([0,30,60,90,120]);ax[2].set_xticklabels(['18:44:00','18:44:30','18:45:00','18:45:30','18:46:00'])
    ax[2].set_xlabel('Time (HH:MM:SS UT)')
    ax[0].set_ylabel('Frequency (GHz)');ax[1].set_ylabel('T$_{B}$ (MK)');ax[2].set_ylabel('Source size (arcsec$^{2}$)')
    plt.show()
plot_vlasize_freq=1
if(plot_vlasize_freq):
    fig,ax=plt.subplots(2,1,sharex=True)
    ax[0].plot(freq,vlasize_fwhm[0,:,46],'o-',label='18:44:46')
    ax[0].plot(freq,vlasize_fwhm[0,:,51],'o-',label='18:44:51')
    ax[0].legend();ax[0].set_ylabel('Source size (arcsec$^{2}$)');ax[0].set_xlabel('Frequency (GHz)')
    ax[1].plot(freq,Tb_mean_r1[:,46],'o-',label='18:44:46')
    ax[1].plot(freq,Tb_mean_r1[:,51],'o-',label='18:44:51')
    ax[1].legend();ax[1].set_ylabel('Source size (arcsec$^{2}$)');ax[1].set_xlabel('Frequency (GHz)')
    #ax[1].set_xticks([0,10,20,30]);ax[1].set_xticklabels([freq[0],freq[10],freq[20],freq[30]])
    plt.show()

plot_aia_movie=1
if(plot_aia_movie):
    for i in range(len(allmaps['aia131']['data131'])):
        f=plt.figure(figsize=(10,10))
        ax0 = f.add_subplot(221)
        cc=allmaps['aia131']['map131'][i]
        xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
        p=cc.plot(extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
        ax0.set_xlim([-1000,-700]);ax0.set_ylim([50,350])
        ax0 = f.add_subplot(222)
        cc=allmaps['aia171']['map171'][i]
        xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
        p=cc.plot(extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
        ax0.set_xlim([-1000,-700]);ax0.set_ylim([50,350])
        ax0 = f.add_subplot(223)
        cc=allmaps['aia94']['map94'][i]
        xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
        p=cc.plot(extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
        ax0.set_xlim([-1000,-700]);ax0.set_ylim([50,350])
        ax0 = f.add_subplot(224)
        cc=allmaps['aia335']['map335'][i]
        xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
        p=cc.plot(extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
        ax0.set_xlim([-1000,-700]);ax0.set_ylim([50,350])
        plt.savefig('/media/rohit/VLA/20160409/pngs/aia_'+str("%03d"%i)+'.png')
        plt.close()

plot_run_hmi_movie=1
if(plot_base_hmi_movie):
    for i in range(1,len(maphmir)):
        f=plt.figure(figsize=(10,10))
        ax0 = f.add_subplot(111)
        cc=maphmir[i]
        xlaia=cc.center.Tx.value-0.5*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.5*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.5*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.5*int(cc.data.shape[0]/2)
        p=cc.plot(extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-80,vmax=80,cmap='binary')
        ax0.set_xlim([-800,-730]);ax0.set_ylim([230,300]);ax0.plot([-778,-778],[245,268],'o',color='cyan')
        f.colorbar(p,orientation='horizontal',label='B (G)')
        plt.savefig('/media/rohit/VLA/20160409/pngs/hmir2_'+str("%03d"%i)+'.png')
        plt.close()

plot_base_hmi_movie=1
if(plot_base_hmi_movie):
    for i in range(len(maphmib)):
        f=plt.figure(figsize=(10,10))
        ax0 = f.add_subplot(111)
        cc=maphmib[i]
        xlaia=cc.center.Tx.value-0.5*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.5*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.5*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.5*int(cc.data.shape[0]/2)
        p=cc.plot(extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-150,vmax=150,cmap='binary')
        ax0.set_xlim([-800,-730]);ax0.set_ylim([200,270]);ax0.plot([-778,-778],[245,268],'o',color='cyan')
        f.colorbar(p,orientation='horizontal',label='B (G)')
        plt.savefig('/media/rohit/VLA/20160409/pngs/hmib1_'+str("%03d"%i)+'.png')
        plt.close()

plot_base_movie=1
if(plot_base_movie):
    for i in range(len(allmapsb['aiab171']['datab171'])):
        f=plt.figure(figsize=(10,10))
        ax0 = f.add_subplot(221)
        cc=allmapsb['aiab131']['mapb131'][i]
        xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
        p=cc.plot(extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-20,vmax=20,cmap='coolwarm')
        ax0.set_xlim([-800,-730]);ax0.set_ylim([200,270]);ax0.plot([-778,-778],[245,268],'o',color='k')
        ax0 = f.add_subplot(222)
        cc=allmapsb['aiab171']['mapb171'][i]
        xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
        p=cc.plot(extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-150,vmax=150,cmap='coolwarm')
        ax0.set_xlim([-800,-730]);ax0.set_ylim([200,270]);ax0.plot([-778,-778],[245,268],'o',color='k')
        ax0 = f.add_subplot(223)
        cc=allmapsb['aiab94']['mapb94'][i]
        xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
        p=cc.plot(extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-10,vmax=10,cmap='coolwarm')
        ax0.set_xlim([-800,-730]);ax0.set_ylim([200,270]);ax0.plot([-778,-778],[245,268],'o',color='k')
        ax0 = f.add_subplot(224)
        cc=allmapsb['aiab335']['mapb335'][i]
        xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
        p=cc.plot(extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-10,vmax=10,cmap='coolwarm')
        ax0.set_xlim([-800,-730]);ax0.set_ylim([200,270]);ax0.plot([-778,-778],[245,268],'o',color='k')
        plt.savefig('/media/rohit/VLA/20160409/pngs/base_zoom_aia_'+str("%03d"%i)+'.png')
        plt.close()

plot_euv1600_movie=1
if(plot_euv1600_movie):
    for i in range(len(allmaps['aia1600']['data1600'])):
        p=allmaps['aia1600']['map1600'][i].plot()
        plt.xlim([400,900]);plt.ylim([50,550])
        plt.savefig('/media/rohit/VLA/20160409/pngs/aia1600_'+str("%03d"%i)+'.png')
        plt.close()


freq_id=[0,15,30]
for i in range(1):
    i=1400
    p=allmaps['aia1600']['map1600'][tidx131[i]].plot()
    v1=allmaps['vla']['mapvla'][2399*freq_id[0]+i]
    v1.data[np.isnan(v1.data)]=0
    #v2=allmaps['vla']['mapvla'][2399*freq_id[1]+i]
    #v3=allmaps['vla']['mapvla'][2399*freq_id[2]+i]
    #lev1=(1.e7/v1.data.max())*np.array([50,60,70,90])*u.percent
    lev1=np.array([50,60,70,90])*u.percent
    #lev2=(3.e7/v2.data.max())*np.array([50,60,70,90])*u.percent
    #lev3=(3.e7/v3.data.max())*np.array([50,60,70,90])*u.percent
    v1.draw_contours(levels=lev1,colors='r',linewidths=3,extent=p.get_extent())
    #v2.draw_contours(levels=lev2,colors='g',linewidths=3,extent=p.get_extent())
    #v3.draw_contours(levels=lev3,colors='b',linewidths=3,extent=p.get_extent())
    #plt.xlim([400,900]);plt.ylim([50,550])
    plt.text(401,250,str(np.round(v1.meta['crval3']/1.e9,4))+' GHz',color='r')
    #plt.text(401,200,str(np.round(v2.meta['crval3']/1.e9,4))+' GHz',color='g')
    #plt.text(401,150,str(np.round(v3.meta['crval3']/1.e9,4))+' GHz',color='b')
    plt.text(401,50,'Contours: 1.5,1.8,2.1,2.7 MK',color='yellow')
    plt.title('VLA: '+v1.meta['date-obs'])
    plt.show()
    plt.savefig('/media/rohit/VLA/20160409/pngs_50ms/aia131_contour_'+str("%03d"%i)+'.png')
    plt.close()

#xl = SkyCoord(dd.center.Tx.value-0.5*(2.0/0.6)*dd.data.shape[0],dd.center.Ty.value-0.5*(2.0/0.6)*dd.data.shape[1], frame=cc.coordinate_frame, unit=(u.arcsec, u.arcsec))
#xlpix=cc.world_to_pixel(xl).x.value;ylpix=cc.world_to_pixel(xl).y.value
#yr = SkyCoord(dd.center.Tx.value+0.5*(2.0/0.6)*dd.data.shape[0],dd.center.Ty.value+0.5*(2.0/0.6)*dd.data.shape[1], frame=cc.coordinate_frame, unit=(u.arcsec, u.arcsec))
#xrpix=cc.world_to_pixel(yr).x.value;yrpix=cc.world_to_pixel(yr).y.value

i=0
for i in range(2399):
    f=plt.figure(figsize=(6, 10))
    ax0 = f.add_subplot(211)
    cc=allmaps['aia131']['map131'][tidx131[i]]
    dd=allmaps['vla']['mapvla'][i];dd.data[np.isnan(dd.data)]=0
    xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
    p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    lev1=(1.5e7/dd.data.max())*np.array([20,30,40,50,60,70,80,90])*u.percent
    #lev1=np.array([20,30,40,50,60,70,80,90])*u.percent
    xlvla=dd.center.Tx.value-2.0*int(dd.data.shape[0]/2);xrvla=dd.center.Tx.value+2.0*int(dd.data.shape[0]/2);ylvla=dd.center.Ty.value-2.0*int(dd.data.shape[1]/2);yrvla=dd.center.Ty.value+2.0*int(dd.data.shape[0]/2)
    dd.draw_contours(levels=lev1,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
    ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
    ax0.text(-1200,0,str(np.round(dd.meta['crval3']/1.e9,4))+' GHz',color='r')
    ax0.text(-1200,50,'Contours: 3, 4.5, 6, 7.5, 9, 10, 12, 13 MK',color='yellow')
    #ax0.text(-1200,50,'Contours: 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%',color='yellow')
    ax1 = f.add_subplot(212)
    ax1.plot(np.arange(2400)*0.05,ds[0].mean(axis=0)[0],'-',label='1.077 GHz')
    ax1.axvline(x=i*0.05,color='k')
    ax1.set_xlabel('Time (sec)');ax1.set_ylabel('Median Amplitude')
    plt.savefig('/media/rohit/VLA/20160409/pngs_50ms/aia131_contour_Tb_'+str("%03d"%i)+'.png')
    plt.close()


i=0
for i in range(32):
    x=np.arange(2000)*0.05
    f=plt.figure(figsize=(20, 10))
    ax0 = f.add_subplot(211)
    ax0.plot(x,ycimax[i],'o',color='red',label='Stokes I')
    ax0.plot(x,ycvmax[i],'o',color='blue',label='Stokes V')
    ax1 = f.add_subplot(212)
    ax1.plot(x,ycimax[i]-ycvmax[i],'o')
    ax1.set_xlabel('Time (sec)');ax1.set_ylabel('Solar-Y (arcsec)');ax0.set_ylabel('Solar-Y (arcsec)')
    ax1.set_ylim(-30,30);ax0.set_ylim(235,275);ax0.legend();ax0.set_title(str(freq[i])+' GHz')
    ax0.axhline(y=255,linestyle='--',color='black')
    ax2=ax1.twinx();ax2.plot(x,maxTbi[i]/1.e6,'-',color='black');ax2.set_ylabel('T$_B$ (MK)')
    plt.savefig('/media/rohit/VLA/20160409/spec_max/stokes_I_'+str("%04d"%i)+'.png')
    f.clear()
    plt.close()

################ Plot centroids on AIA

i=0
for i in range(1900):
    x=np.arange(2000)[i:i+100]*0.05
    f=plt.figure(figsize=(10, 15))
    ax0 = f.add_subplot(211)
    #cc=allmapsb['aiab171']['mapb171'][tidx171[i]]
    cc=allmaps['aia171']['map171'][tidx171[i]]
    #cc=allmaps['aia171']['map171'][tidx171[i]]
    #dd=mapvla_v[i+50][0];dd.data[np.isnan(dd.data)]=0
    xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
    p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    #im0=ax0.scatter(xcl90[i+50].reshape(1,32).mean(axis=0)[0:16],ycl90[i+50].reshape(1,32).mean(axis=0)[0:16],c=freq.reshape(1,32).mean(axis=0)[0:16],s=40,cmap=plt.cm.get_cmap('winter'),edgecolors='none')
    im0=ax0.scatter(xcr90[i+50].reshape(1,32).mean(axis=0)[0:16],ycr90[i+50].reshape(1,32).mean(axis=0)[0:16],c=freq.reshape(1,32).mean(axis=0)[0:16],s=40,cmap=plt.cm.get_cmap('spring'),edgecolors='none')
    ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
    #ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
    #ax0.set_xlim([-1000,-700]);ax0.set_ylim([50,350])
    ax0.set_xlim([-810,-750]);ax0.set_ylim([220,280])
    f.colorbar(im0,label='(RR) Frequency (GHz)')
    ax1 = f.add_subplot(212)
    ax1.plot(x,TbRRr1[i:i+100,0]/1.e6,'o-',markersize=5,label='Region 1 (0.994 GHz)')
    ax1.plot(x[50],TbRRr1[i+50][0]/1.e6,'o',color='red',markersize=10)
    #ax1.plot(x,Tb_r[1]/1.e6,'-',label='Region 2')
    #ax1.plot(x,Tb_r[2]/1.e6,'-',label='Region 3')
    #ax1.plot(x,Tb_r[3]/1.e6,'-',label='Region 4')
    #ax1.plot(x,Tb_r[4]/1.e6,'-',label='Region 5')
    ax1.legend(loc=2)
    ax1.axvline(x=x[50],color='k');ax1.set_ylim(0,100)
    ax1.set_xlabel('Time (sec)');ax1.set_ylabel('$T_B$ (MK)')
    plt.savefig('/media/rohit/VLA/20160409/pngs_50ms_spw0/centroid_llrr_'+str("%04d"%i)+'.png')
    f.clear()
    plt.close()

i=0
for i in range(2399,6799):
    x=np.arange(6799)[i:i+100]*0.05
    y=np.hstack((np.array(qsTbr_r1),Tbr_r1[1]))
    centx=np.hstack((np.array(qsxcr90),xcr90[1]));centy=np.hstack((np.array(qsycr90),ycr90[1]))
    f=plt.figure(figsize=(10, 15))
    cc=allmaps['aia171']['map171'][tidx171[i]]
    xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
    ax0 = f.add_subplot(211)
    p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    ax0.plot(centx[i+40:i+50],centy[i+40:i+50],'o',color='red')
    #ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
    ax0.set_xlim([-810,-750]);ax0.set_ylim([220,280])
    ax1 = f.add_subplot(212)
    ax1.plot(x,y[i:i+100]/1.e6,'o-',markersize=5,label='Region 1 (0.994 GHz)')
    ax1.plot(x[50],y[i+50]/1.e6,'o',color='red',markersize=10)
    ax1.legend(loc=2)
    ax1.axvline(x=x[50],color='k');ax1.set_ylim(0,100)
    ax1.set_xlabel('Time (sec)');ax1.set_ylabel('$T_B$ (MK)')
    plt.savefig('/media/rohit/VLA/20160409/pngs_50ms_spw0_qs/centroid_llrr_'+str("%04d"%i)+'.png')
    f.clear()
    plt.close()


f=plt.figure(figsize=(15, 15))
cc=allmaps['aia171']['map171'][10]
xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
ax0 = f.add_subplot(111)
p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
ax0.plot(centx[0:2000],centy[0:2000],'o',color='b',alpha=0.2,label='18:44:00 to 18:45:40 UT')
ax0.plot(centx[2000:5000],centy[2000:5000],'o',color='g',alpha=0.2,label='18:45:40 to 18:48:10 UT')
ax0.plot(centx[5000:],centy[5000:],'o',color='r',alpha=0.2,label='18:48:10 to 18:49:40 UT')
ax0.legend()
ax0.set_xlim([-810,-750]);ax0.set_ylim([220,280])
plt.show()


i=0
for i in range(1900):
    x=np.arange(2000)[i:i+100]*0.05
    f=plt.figure(figsize=(10, 15))
    ax0 = f.add_subplot(211)
    #cc=allmapsb['aiab171']['mapb171'][tidx171[i]]
    cc=allmaps['aia171']['map171'][tidx171[i]]
    #cc=allmaps['aia171']['map171'][tidx171[i]]
    #dd=mapvla_v[i+50][0];dd.data[np.isnan(dd.data)]=0
    xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
    p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    im0=ax0.errorbar(xcLL90[i+50].reshape(1,32).mean(axis=0)[0:6].mean(),ycLL90[i+50].reshape(1,32).mean(axis=0)[0:6].mean(),xerr=xcLL90[i+50].reshape(1,32).mean(axis=0)[0:6].std(),yerr=ycLL90[i+50].reshape(1,32).mean(axis=0)[0:6].std(),color='blue')
    im0=ax0.errorbar(xcRR90[i+50].reshape(1,32).mean(axis=0)[0:16].mean(),ycRR90[i+50].reshape(1,32).mean(axis=0)[0:16].mean(),xerr=xcRR90[i+50].reshape(1,32).mean(axis=0)[0:16].std(),yerr=ycRR90[i+50].reshape(1,32).mean(axis=0)[0:16].std(),color='magenta')
    ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
    #ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
    #ax0.set_xlim([-1000,-700]);ax0.set_ylim([50,350])
    ax0.set_xlim([-810,-750]);ax0.set_ylim([220,280])
    ax1 = f.add_subplot(212)
    ax1.plot(x,TbRRr1[i:i+100,0]/1.e6,'o-',markersize=5,label='Region 1 (0.994 GHz)')
    ax1.plot(x[50],TbRRr1[i+50][0]/1.e6,'o',color='red',markersize=10)
    #ax1.plot(x,Tb_r[1]/1.e6,'-',label='Region 2')
    #ax1.plot(x,Tb_r[2]/1.e6,'-',label='Region 3')
    #ax1.plot(x,Tb_r[3]/1.e6,'-',label='Region 4')
    #ax1.plot(x,Tb_r[4]/1.e6,'-',label='Region 5')
    ax1.legend(loc=2)
    ax1.axvline(x=x[50],color='k');ax1.set_ylim(0,100)
    ax1.set_xlabel('Time (sec)');ax1.set_ylabel('$T_B$ (MK)')
    plt.savefig('/media/rohit/VLA/20160409/pngs_50ms_spw0/centroid_llrr1_'+str("%04d"%i)+'.png')
    f.clear()
    plt.close()

i=0
for i in range(1900):
    x=np.arange(2000)[i:i+100]*0.05
    f=plt.figure(figsize=(10, 15))
    ax0 = f.add_subplot(211)
    #cc=allmapsb['aiab171']['mapb171'][tidx171[i]]
    cc=allmaps['aia171']['map171'][tidx171[i]]
    #cc=allmaps['aia171']['map171'][tidx171[i]]
    dd=mapvla_v[i+50][0];dd.data[np.isnan(dd.data)]=0
    xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
    p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    im0=ax0.scatter(xcvmax[:25,i],ycvmax[:25,i],c=freq[:25],s=40,cmap=plt.cm.get_cmap('RdYlBu'),edgecolors='none')
    ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
    #ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
    #ax0.set_xlim([-1000,-700]);ax0.set_ylim([50,350])
    ax0.set_xlim([-830,-730]);ax0.set_ylim([200,300])
    f.colorbar(im0,label='Frequency (GHz)')
    ax1 = f.add_subplot(212)
    ax1.plot(x,Tb_r[0][i:i+100]/1.e6,'o-',markersize=5,label='Region 1 (0.994 GHz)')
    ax1.plot(x[50],Tb_r[0][i+50]/1.e6,'o',color='red',markersize=10)
    #ax1.plot(x,Tb_r[1]/1.e6,'-',label='Region 2')
    #ax1.plot(x,Tb_r[2]/1.e6,'-',label='Region 3')
    #ax1.plot(x,Tb_r[3]/1.e6,'-',label='Region 4')
    #ax1.plot(x,Tb_r[4]/1.e6,'-',label='Region 5')
    ax1.legend(loc=2)
    ax1.axvline(x=x[50],color='k');ax1.set_ylim(0,50)
    ax1.set_xlabel('Time (sec)');ax1.set_ylabel('$T_B$ (MK)')
    plt.savefig('/media/rohit/VLA/20160409/pngs_50ms_spw0/centroid_'+str("%04d"%i)+'.png')
    f.clear()
    plt.close()

i=0
for i in range(1900):
    x=np.arange(2000)[i:i+100]*0.05
    f=plt.figure(figsize=(6, 10))
    ax0 = f.add_subplot(211)
    #cc=allmapsb['aiab171']['mapb171'][tidx171[i]]
    cc=allmaps['aia171']['map171'][tidx171[i]]
    #cc=allmaps['aia171']['map171'][tidx171[i]]
    dd=mapvla_v[i+50][0];dd.data[np.isnan(dd.data)]=0
    xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
    p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    #p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto',vmin=-150,vmax=150,cmap='coolwarm')
    #lev1=(1.5e7/dd.data.max())*np.array([20,30,40,50,60,70,80,90])*u.percent
    lev1=np.array([50,60,70,80,90,99])*u.percent
    xlvla=dd.center.Tx.value-2.0*int(dd.data.shape[0]/2);xrvla=dd.center.Tx.value+2.0*int(dd.data.shape[0]/2);ylvla=dd.center.Ty.value-2.0*int(dd.data.shape[1]/2);yrvla=dd.center.Ty.value+2.0*int(dd.data.shape[0]/2)
    dd.draw_contours(levels=lev1,colors='black',linewidths=1,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
    #ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
    #ax0.set_xlim([-1000,-700]);ax0.set_ylim([50,350])
    ax0.set_xlim([-830,-730]);ax0.set_ylim([200,300])
    ax0.text(-830,210,str(np.round(dd.meta['crval3']/1.e9,4))+' GHz',color='white')
    #ax0.text(-1200,50,'Contours: 0.6,0.75,0.9,1.0,1.2,1.3 MK',color='yellow')
    #ax0.text(-990,64,'Contours: 20%, 30%, 40%,50%, 60%, 70%, 80%, 90%',color='yellow')
    #ax0.text(-990,64,'Contours: 3, 4.5, 6, 7.5, 9, 10, 12, 13 MK',color='yellow')
    #ax0.add_patch(patches.Rectangle((-781,233),40,40,linewidth=5,edgecolor='b',facecolor='none'))
    #ax0.add_patch(patches.Rectangle((-841,265),20,20,linewidth=5,edgecolor='g',facecolor='none'))
    #ax0.add_patch(patches.Rectangle((-881,175),40,40,linewidth=5,edgecolor='r',facecolor='none'))
    #ax0.add_patch(patches.Rectangle((-961,195),40,40,linewidth=5,edgecolor='cyan',facecolor='none'))
    #ax0.add_patch(patches.Rectangle((-831,175),40,40,linewidth=5,edgecolor='magenta',facecolor='none'))
    ax1 = f.add_subplot(212)
    ax1.plot(x,Tb_r[0][i:i+100]/1.e6,'o-',markersize=5,label='Region 1')
    ax1.plot(x[50],Tb_r[0][i+50]/1.e6,'o',color='red',markersize=10)
    #ax1.plot(x,Tb_r[1]/1.e6,'-',label='Region 2')
    #ax1.plot(x,Tb_r[2]/1.e6,'-',label='Region 3')
    #ax1.plot(x,Tb_r[3]/1.e6,'-',label='Region 4')
    #ax1.plot(x,Tb_r[4]/1.e6,'-',label='Region 5')
    ax1.legend(loc=2)
    ax1.axvline(x=x[50],color='k');ax1.set_ylim(0,50)
    ax1.set_xlabel('Time (sec)');ax1.set_ylabel('$T_B$ (MK)')
    plt.savefig('/media/rohit/VLA/20160409/pngs_50ms_spw0/regions_contour_base_'+str("%04d"%i)+'.png')
    f.clear()
    plt.close()


i=0
for i in range(2399):
    f=plt.figure(figsize=(6, 10))
    ax0 = f.add_subplot(211)
    cc=allmaps['aia171']['map171'][tidx131[i]]
    dd=allmaps['vla']['mapvla'][i];dd.data[np.isnan(dd.data)]=0
    xlaia=cc.center.Tx.value-0.61*int(cc.data.shape[0]/2);xraia=cc.center.Tx.value+0.61*int(cc.data.shape[0]/2);ylaia=cc.center.Ty.value-0.61*int(cc.data.shape[1]/2);yraia=cc.center.Ty.value+0.61*int(cc.data.shape[0]/2)
    p=cc.plot(axes=ax0,extent=[xlaia,xraia,ylaia,yraia],aspect='auto')
    #lev1=(1.5e7/dd.data.max())*np.array([50,60,70,80,90])*u.percent
    lev1=np.array([50,60,70,80,90])*u.percent
    xlvla=dd.center.Tx.value-2.0*int(dd.data.shape[0]/2);xrvla=dd.center.Tx.value+2.0*int(dd.data.shape[0]/2);ylvla=dd.center.Ty.value-2.0*int(dd.data.shape[1]/2);yrvla=dd.center.Ty.value+2.0*int(dd.data.shape[0]/2)
    dd.draw_contours(levels=lev1,colors='r',linewidths=2,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    ax0.set_title('AIA 171 $\AA$:'+timstr_171[tidx171[i]]+' VLA: '+timstr_vla[i]+' UT')
    ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
    ax0.text(-1200,0,str(np.round(dd.meta['crval3']/1.e9,4))+' GHz',color='r')
    #ax0.text(-1200,50,'Contours: 0.6,0.75,0.9,1.0,1.2,1.3 MK',color='yellow')
    ax0.text(-1200,50,'Contours: 50%, 60%, 70%, 80%, 90%',color='yellow')
    ax1 = f.add_subplot(212)
    ax1.plot(np.arange(2400)*0.05,ds[0].mean(axis=0)[0],'-',label='1.077 GHz')
    ax1.axvline(x=i*0.05,color='k')
    ax1.set_xlabel('Time (sec)');ax1.set_ylabel('Median Amplitude')
    plt.savefig('/media/rohit/VLA/20160409/pngs_50ms/aia171_contour_'+str("%03d"%i)+'.png')
    plt.close()

for i in range(2399):
    cc=allmaps['aia1600']['map1600'][tidx131[i]]
    dd=allmaps['vla']['mapvla'][i];dd.data[np.isnan(dd.data)]=0
    p=cc.plot()
    lev1=(1.5e7/dd.data.max())*np.array([40,50,60,70,80,90])*u.percent
    xl = SkyCoord(dd.center.Tx.value-0.5*(2.0/0.6)*dd.data.shape[0],dd.center.Ty.value-0.5*(2.0/0.6)*dd.data.shape[1], frame=cc.coordinate_frame, unit=(u.arcsec, u.arcsec))
    xlpix=cc.world_to_pixel(xl).x.value;ylpix=cc.world_to_pixel(xl).y.value
    yr = SkyCoord(dd.center.Tx.value+0.5*(2.0/0.6)*dd.data.shape[0],dd.center.Ty.value+0.5*(2.0/0.6)*dd.data.shape[1], frame=cc.coordinate_frame, unit=(u.arcsec, u.arcsec))
    xrpix=cc.world_to_pixel(yr).x.value;yrpix=cc.world_to_pixel(yr).y.value
    dd.draw_contours(levels=lev1,colors='r',linewidths=3,extent=[xlpix,xrpix,ylpix,yrpix])
    plt.text(401,250,str(np.round(dd.meta['crval3']/1.e9,4))+' GHz',color='r')
    plt.text(401,50,'Contours: 0.6,0.75,0.9,1.0,1.2,1.3 MK',color='yellow')
    plt.title('AIA 1600 $\AA$:'+timstr_1600[tidx131[i]]+' VLA: '+timstr_vla[i]+' UT')
    plt.xlim([400,900]);plt.ylim([50,550])
    plt.savefig('/media/rohit/VLA/20160409/pngs_50ms/aia1600_contour_'+str("%03d"%i)+'.png')
    plt.close()

for i in range(1):
    cc=allmaps['aia1600']['map1600'][tidx131[i]]
    dd=allmaps['vla']['mapvla'][i];dd.data[np.isnan(dd.data)]=0
    p=cc.plot()
    lev1=(1.5e7/dd.data.max())*np.array([40,50,60,70,80,90])*u.percent
    xl = SkyCoord(dd.center.Tx.value-0.5*(2.0/0.6)*dd.data.shape[0],dd.center.Ty.value-0.5*(2.0/0.6)*dd.data.shape[1], frame=cc.coordinate_frame, unit=(u.arcsec, u.arcsec))
    xlpix=cc.world_to_pixel(xl).x.value;ylpix=cc.world_to_pixel(xl).y.value
    yr = SkyCoord(dd.center.Tx.value+0.5*(2.0/0.6)*dd.data.shape[0],dd.center.Ty.value+0.5*(2.0/0.6)*dd.data.shape[1], frame=cc.coordinate_frame, unit=(u.arcsec, u.arcsec))
    xrpix=cc.world_to_pixel(yr).x.value;yrpix=cc.world_to_pixel(yr).y.value
    dd.draw_contours(levels=lev1,colors='r',linewidths=3,extent=[xlpix,xrpix,ylpix,yrpix])
    plt.text(401,250,str(np.round(dd.meta['crval3']/1.e9,4))+' GHz',color='r')
    plt.text(401,50,'Contours: 0.6,0.75,0.9,1.0,1.2,1.3 MK',color='yellow')
    plt.xlim([400,900]);plt.ylim([50,550])
    plt.show()
    plt.savefig('/media/rohit/VLA/20160409/pngs_50ms/aia1600_contour_'+str("%03d"%i)+'.png')
    plt.close()

for i in range(1):
    i=1400
    compmap=Map(allmaps['aia171']['map171'][tidx171[i]],allmaps['vla']['mapvla'][i],composite=True)
    compmap.set_levels(1,[30,40,50,60,70,80,90],percent=True)
    compmap.plot()
    plt.text(401,250,str(np.round(v1.meta['crval3']/1.e9,4))+' GHz',color='r')
    plt.text(401,50,'Contours: 1.5,1.8,2.1,2.7 MK',color='yellow')
    #plt.xlim([400,900]);plt.ylim([50,550])
    plt.show()
    plt.savefig('/media/rohit/VLA/20160409/pngs_50ms/aia131_contour_'+str("%03d"%i)+'.png')
    plt.close()

for i in range(119):
    p=allmaps['aia1600']['map1600'][tidx1600[i]].plot()
    v1=allmaps['vla']['mapvla'][119*freq_id[0]+i]
    v2=allmaps['vla']['mapvla'][119*freq_id[1]+i]
    v3=allmaps['vla']['mapvla'][119*freq_id[2]+i]
    lev1=(3.e7/v1.data.max())*np.array([50,60,70,90])*u.percent
    lev2=(3.e7/v2.data.max())*np.array([50,60,70,90])*u.percent
    lev3=(3.e7/v3.data.max())*np.array([50,60,70,90])*u.percent
    v1.draw_contours(levels=lev1,colors='r',linewidths=3,extent=p.get_extent())
    v2.draw_contours(levels=lev2,colors='g',linewidths=3,extent=p.get_extent())
    v3.draw_contours(levels=lev3,colors='b',linewidths=3,extent=p.get_extent())
    plt.xlim([400,900]);plt.ylim([50,550])
    plt.text(401,250,str(np.round(v1.meta['crval3']/1.e9,4))+' GHz',color='r')
    plt.text(401,200,str(np.round(v2.meta['crval3']/1.e9,4))+' GHz',color='g')
    plt.text(401,150,str(np.round(v3.meta['crval3']/1.e9,4))+' GHz',color='b')
    plt.text(401,50,'Contours: 15,18,21,27 MK',color='yellow')
    plt.title('VLA: '+v1.meta['date-obs'])
    plt.savefig('/media/rohit/VLA/20160409/pngs/aia1600_contour_'+str("%03d"%i)+'.png')
    plt.close()
    


compmap=Map(map131[10],mapvla[0],composite=True)
compmap.set_levels(1,[10,20,30,40,50,60,70,80,90],percent=True)
compmap.set_plot_settings(1,{'cmap':'gray','norm':mpl.colors.Normalize(vmin=-1.,vmax=1.),'origin':0})
compmap.plot()
plt.show()


sys.exit()
# Time from 17:59:56-22:57:52
plt.style.use('/home/rohit/my_git/lipi/plt_style.py')
data=np.load('/media/rohit/VLA/20160409/sun_L_20160409.1s.ms.dspec.npz')
ds=data['spec']
freq=data['freq']/1.e9
time=data['tim']-data['tim'][0]
ds_LL=ds[0][0]
ds_RR=ds[1][0]
ds_I=0.5*(ds_LL+ds_RR)
#ds_I=ds_I[:,771:]
#for i in range(120):
#   ds_I[:,i*120]=np.nan

# For feature 1
#ds_I1=ds_I[:,:2680]
#ds_I2=ds_I[:,2050:2200]
# Feature 2
# 18:40:00 to 19:00:00
time_feature2=time[2585:2585+10*60]
ds_feature2=ds_I[:,2585:2585+10*60]

#plt.imshow(np.log(ds_I1),aspect='auto',extent=[time[0],time[-1],freq[0],freq[-1]],cmap='jet',interpolation=None,origin=0,vmin=-3,vmax=1)

f=plt.figure()
ax1=f.add_subplot(211)
ax2=f.add_subplot(212,sharex=ax1)
ax1.imshow(ds_feature2,origin=0,vmin=0.01,vmax=0.8,aspect='auto')
ax2.plot(ds_feature2[10],'o-')
ax2.set_xticks([0,60,120,180,240,300,360,420,480,540,600])
ax2.set_xticklabels(['18:40','18:41','18:42','18:43','18:44','18:45','18:46','18:47','18:48','18:49','18:50'])
ax1.set_yticks([0,100,200,300,400,500])
ax1.set_yticklabels([1.0,1.2,1.4,1.6,1.8,2.0])
ax2.set_xlim([0,600])
ax1.set_ylabel('Frequency (GHz)')
ax2.set_xlabel('Time (HH:MM) UT')
plt.show()

plot_full_ds=0
if(plot_full_ds):
    f=plt.figure()
    ax1=f.add_subplot(211)
    ax2=f.add_subplot(212)
    ax1.imshow(np.log(ds_I2),aspect='auto',cmap='jet',extent=[time[0],7.5,freq[0],1.122],interpolation=None,origin=0,vmin=-3,vmax=1)
    ax2.set_xlabel('Time (s)')
    ax1.set_ylabel('Frequency (GHz)')
    ax2.plot(np.arange(150)/20.0,np.nanmean(ds_I2[10:20,:],axis=0),'o-')
    ax2.set_ylabel('Amp')
    plt.show()

