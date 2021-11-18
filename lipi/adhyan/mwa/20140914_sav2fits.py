from astropy.io import fits
import glob
import pickle
from sunpy.map import Map
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from surya.utils import main as ut
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import os
from scipy.io import readsav
from astropy.time import Time
import sunpy
import sunpy.map as smap
from sunpy import sun
from dateutil import parser
import numpy as np


def idl2sunpy_hmi(mapsav,n):
    mapstruc = readsav(mapsav)
    keys = list(mapstruc.keys())
    mp = mapstruc[keys[0]][0]
    dat = np.nan_to_num(mp[0]) + 1e-6 # just to get rid of zeros
    header = {}
    header['EXPTIME'] = mp['dur']
    header['CRVAL1'] = mp[1]
    header['CRVAL2'] = mp[2]
    header['CDELT1'] = mp[3]*n
    header['CDELT2'] = mp[4]*n
    header['DATE-OBS'] = Time(parser.parse(mp[5])).isot
    header['CUNIT1'] = 'arcsec'
    header['CUNIT2'] = 'arcsec'
    header['CTYPE1'] = 'HPLN-TAN'
    header['CTYPE2'] = 'HPLT-TAN'
    header['CRPIX2'], header['CRPIX1'] = (np.array(dat.shape)/n + 1.0) / 2 + 0.5
    header['TELESCOP'] = 'SDO/HMI'
    header['P_ANGLE'] = 0.0
    header['NAXIS'] = 2
    header['NAXIS2'], header['NAXIS1'] = np.array(dat.shape)/n
    header['HGLN_OBS'] = 0.
    header['HGLT_OBS'] = sunpy.coordinates.get_sun_B0(header['DATE-OBS']).value
    header['RSUN_REF'] = sunpy.sun.constants.radius.value
    header['RSUN_OBS'] = 958.11#sunpy.coordinates.sun.angular_radius(header['DATE-OBS']).value
    header['DSUN_OBS'] = sunpy.coordinates.get_sunearth_distance(header['DATE-OBS']).to(u.meter).value
    d=reduce_size(dat,n)
    smp = smap.Map(d, header)
    return smp


def idl2sunpy_sdo(mapsav,n,wave,inst):
    mapstruc = readsav(mapsav)
    keys = list(mapstruc.keys())
    mp = mapstruc[keys[0]][0]
    dat = np.nan_to_num(mp[0]) + 1e-6 # just to get rid of zeros
    header = {}
    header['EXPTIME'] = mp['dur']
    header['CRVAL1'] = mp[1]
    header['CRVAL2'] = mp[2]
    header['CDELT1'] = mp[3]*n
    header['CDELT2'] = mp[4]*n
    header['DATE-OBS'] = Time(parser.parse(mp[5])).isot
    header['CUNIT1'] = 'arcsec'
    header['CUNIT2'] = 'arcsec'
    header['CTYPE1'] = 'HPLN-TAN'
    header['CTYPE2'] = 'HPLT-TAN'
    header['CRPIX2'], header['CRPIX1'] = (np.array(dat.shape)/n + 1.0) / 2 + 0.5
    header['TELESCOP'] = 'SDO'
    header['INSTRUME'] = inst
    header['WAVELNTH'] = wave
    header['waveunit'] = 'angstrom'
    header['P_ANGLE'] = 0.0
    header['NAXIS'] = 2
    header['NAXIS2'], header['NAXIS1'] = np.array(dat.shape)/n
    header['HGLN_OBS'] = 0.
    header['HGLT_OBS'] = sunpy.coordinates.get_sun_B0(header['DATE-OBS']).value
    header['RSUN_REF'] = sunpy.sun.constants.radius.value
    header['RSUN_OBS'] = 958.11#sunpy.coordinates.sun.angular_radius(header['DATE-OBS']).value
    header['DSUN_OBS'] = sunpy.coordinates.get_sunearth_distance(header['DATE-OBS']).to(u.meter).value
    d=reduce_size(dat,n)
    smp = smap.Map(d, header)
    return smp


def produce_tstring(mapp):
    date=mapp.date
    hhmmss=' '+str(date.hour)+':'+str(date.minute)+':'+str(date.second)+'.'+str(date.microsecond/1.e6).split('.')[1]
    sec=ut.hms2sec_c(hhmmss)
    return sec

def get_sunpy_maps_rot(f,m,wave,inst,filename):
    print 'Reading...'+f[0]
    n=len(f)
    for i in range(n):
        ii="%04d"%i
        if(inst=='HMI'):
            maplist=idl2sunpy_hmi(f[i],m)
        if(inst=='AIA'):
            maplist=idl2sunpy_sdo(f[i],m,wave,inst)
        tt=maplist.meta['date-obs'].split('T')[1];ti=produce_tstring(maplist)
        pickle.dump(maplist,open(filename+'_'+tt+'_'+str(int(ti))+'.sunpy','wb'))

def get_sunpy_maps(f):
    print 'Reading...'+f[0]
    n=len(f);maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(n):
        maplist[i]=Map(f[i])
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time

def get_sunpy_basediff_maps(f,m,wave,inst,filename):
    print 'Reading...'+f[0]
    n=len(f)
    for i in range(n):
        ii="%04d"%i
        if(inst=='AIA'):
            aa=idl2sunpy_sdo(f[i],m,wave,inst)
            ab=idl2sunpy_sdo(f[0],m,wave,inst)
        if(inst=='HMI'):
            aa=idl2sunpy_hmi(f[i],m)
            ab=idl2sunpy_hmi(f[0],m)
        h,d=aa.meta,aa.data
        hb,db=ab.meta,ab.data
        if(d.shape[0]==db.shape[0]):
            d=d-db
        else:
            if(d.shape[0]<db.shape[0]):
                d=d-db[0:d.shape[0],0:d.shape[1]]
            else:
                d=d[0:db.shape[0],0:db.shape[1]]-db
        maplist=Map(d,h)
        tt=maplist.meta['date-obs'].split('T')[1];ti=produce_tstring(maplist)
        pickle.dump(maplist,open(filename+'_'+tt+'_'+str(int(ti))+'.sunpy','wb'))


def get_sunpy_rundiff_maps(f,m,wave,inst,filename):
    print 'Reading...'+f[0]
    n=len(f)-1
    for i in range(1,n):
        ii="%04d"%i
        if(inst=='AIA'):
            ab=idl2sunpy_sdo(f[i-1],m,wave,inst);hb,db=ab.meta,ab.data
            aa=idl2sunpy_sdo(f[i],m,wave,inst);h,d=aa.meta,aa.data
        if(inst=='HMI'):
            ab=idl2sunpy_hmi(f[i-1],m);hb,db=ab.meta,ab.data
            aa=idl2sunpy_hmi(f[i],m);h,d=aa.meta,aa.data
        if(d.shape[0]==db.shape[0]):
            d=d-db
        else:
            if(d.shape[0]<db.shape[0]):
                d=d-db[0:d.shape[0],0:d.shape[1]]
            else:
                d=d[0:db.shape[0],0:db.shape[1]]-db
        maplist=Map(d,h)
        tt=maplist.meta['date-obs'].split('T')[1];ti=produce_tstring(maplist)
        pickle.dump(maplist,open(filename+'_'+tt+'_'+str(int(ti))+'.sunpy','wb'))

def get_evla_submap(f,xbl,ybl,xtr,ytr):
    print 'Reading...'+f[0]
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
        maplist[i]=map_.submap(bl,tr)
        maplist[i]=map_
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time

def change_pointing(f,xc,yc):
    g=fits.open(f)
    g[0].header['CRVAL1']=xc;g[0].header['CRVAL2']=yc
    g.writeto('/media/rohit/VLA/20160409/sun_L_20160409T1844-1846UT.50ms.cal.pol.LL.spw.0_16-31.time.18:45:10.-18:45:10.05.fits')
    g.close()

def reduce_size(d,n):
    return d.reshape(int(d.shape[0]/n),n,int(d.shape[0]/n),n).mean(axis=(1,3))
sys.exit()
dump_submaps=1
if(dump_submaps):
    list193=sorted(glob.glob('/media/rohit/MWA/20140914/EUV/fits/*193*rot.sav'))[0:5]
    list171=sorted(glob.glob('/media/rohit/MWA/20140914/EUV/fits/*171*.sav'))[0:5]
    hmilist=sorted(glob.glob('/media/rohit/MWA/20140914/EUV/fits/hmi*rot.sav'))[0:5]
    get_sunpy_maps_rot(list171,1,'171','AIA','/media/rohit/MWA/20140914/sunpy_maps/map171')
    get_sunpy_maps_rot(list193,1,'193','AIA','/media/rohit/MWA/20140914/sunpy_maps/map193')
    get_sunpy_basediff_maps(list171,1,'171','AIA','/media/rohit/MWA/20140914/sunpy_maps/map171b')
    get_sunpy_basediff_maps(list193,1,'193','AIA','/media/rohit/MWA/20140914/sunpy_maps/map193b')
    get_sunpy_rundiff_maps(list171,1,'171','AIA','/media/rohit/MWA/20140914/sunpy_maps/map171r')
    get_sunpy_rundiff_maps(list193,1,'193','AIA','/media/rohit/MWA/20140914/sunpy_maps/map193r')
    get_sunpy_maps_rot(hmilist,1,'','HMI','/media/rohit/MWA/20140914/sunpy_maps/hmi')
    get_sunpy_basediff_maps(hmilist,1,'','HMI','/media/rohit/MWA/20140914/sunpy_maps/hmib')
    get_sunpy_rundiff_maps(hmilist,1,'','HMI','/media/rohit/MWA/20140914/sunpy_maps/hmir')
    sys.exit()
    ###########################

list171fits=sorted(glob.glob('/media/rohit/MWA/20140914/sunpy_maps/map171_*.sunpy'))
list171bfits=sorted(glob.glob('/media/rohit/MWA/20140914/sunpy_maps/map171b_*.sunpy'))
list171rfits=sorted(glob.glob('/media/rohit/MWA/20140914/sunpy_maps/map171r_*.sunpy'))
list193fits=sorted(glob.glob('/media/rohit/MWA/20140914/sunpy_maps/map193_*.sunpy'))
list193bfits=sorted(glob.glob('/media/rohit/MWA/20140914/sunpy_maps/map193b_*.sunpy'))
list193rfits=sorted(glob.glob('/media/rohit/MWA/20140914/sunpy_maps/map193r_*.sunpy'))
listhmifits=sorted(glob.glob('/media/rohit/MWA/20140914/sunpy_maps/hmi_*.sunpy'))
listhmibfits=sorted(glob.glob('/media/rohit/MWA/20140914/sunpy_maps/hmib_*.sunpy'))
listhmirfits=sorted(glob.glob('/media/rohit/MWA/20140914/sunpy_maps/hmir_*.sunpy'))

tidx171=[0]*len(list171fits);tidx171_str=[0]*len(list171fits)
tidx193=[0]*len(list193fits);tidx193_str=[0]*len(list193fits)
tidxhmi=[0]*len(listhmifits);tidxhmi_str=[0]*len(listhmifits)
for i in range(len(list171fits)):
    tidx171[i]=int(list171fits[i].split('_')[3].split('.')[0]);tidx171_str[i]=list171fits[0].split('_')[2]
for i in range(len(list193fits)):
    tidx193[i]=int(list193fits[i].split('_')[3].split('.')[0]);tidx193_str[i]=list193fits[0].split('_')[2]
for i in range(len(listhmifits)):
    tidxhmi[i]=int(listhmifits[i].split('_')[3].split('.')[0]);tidxhmi_str[i]=listhmifits[0].split('_')[2]
tidx171=np.array(tidx171);tidx193=np.array(tidx193);tidxhmi=np.array(tidxhmi)

######################## MWA #################################3


listTb=sorted(glob.glob('/media/rohit/MWA/20140914/Tb*.p'))
aia=fits.open('/media/rohit/MWA/20140914/EUV/aia.lev1.171A_2014-09-14T03_06_11.34Z.image_lev1.fits')
freq=[108,161,240]
for j in range(3):
    aa=pickle.load(open(listTb[j],'rb'));Tb=aa[0];time=aa[9]
    for i in range(len(aa[0])):
        td='2014-09-14T0'+time[i]
        h=aia[0].header;d=aia[0].data 
        h['TELESCOP']='sdo';h['INSTRUME']='';h['WAVELNTH']=94;h['WAVEUNIT']='angstrom';h['CDELT2']=50.0;h['DATE-OBS']=td;h['T_REC']=td;h['DATE']=td;h['T_OBS']=td;h['ISPPKTIM']=td
        h['CDELT1']=50.0;h['CRPIX1']=49;h['CRPIX2']=49
        if((i==0) & (j==0)):
            h.remove('OSCNMEAN');h.remove('OSCNRMS')
        fits.writeto('/media/rohit/MWA/20140914/mwa_maps/mwa_'+str(freq[j])+'_'+str(time[i])+'.fits',Tb[i],h)
#mm=Map('test.fits');am=Map('/media/rohit/MWA/20140914/EUV/aia.lev1.171A_2014-09-14T03_06_11.34Z.image_lev1.fits')
#mm.plot();plt.show()

mwa108list=sorted(glob.glob('/media/rohit/MWA/20140914/mwa_maps/mwa_108*.fits'))
mwa161list=sorted(glob.glob('/media/rohit/MWA/20140914/mwa_maps/mwa_161*.fits'))
mwa240list=sorted(glob.glob('/media/rohit/MWA/20140914/mwa_maps/mwa_240*.fits'))
mwatsec=[0]*len(mwa240list)
for i in range(len(mwa240list)):
    mapp=Map(mwa240list[i])
    mwatsec[i]=produce_tstring(mapp)
mwatsec=np.array(mwatsec)

for i in range(200,950):
    #f,ax=plt.subplots(2,1,figsize=(8,15));ax0=ax[0];ax1=ax[1] # 171 is the base
    ii="%04d"%i;mwaidx=ut.find_predecessor(mwatsec,tidx171[i])[0]
    f,ax=plt.subplots(2,2,figsize=(14,14));ax00=ax[0,0];ax01=ax[0,1];ax10=ax[1,0];ax11=ax[1,1]
    map171=pickle.load(open(list171fits[i],'rb'));map193b=pickle.load(open(list193bfits[i],'rb'))
    map171b=pickle.load(open(list171bfits[i],'rb'))
    map171r=pickle.load(open(list171rfits[i-1],'rb'))
    map171.plot(axes=ax00);map171b.plot(axes=ax01,vmin=0,vmax=50);map171r.plot(axes=ax10,vmin=0,vmax=50);map193b.plot(axes=ax11,vmin=0,vmax=80)
    mwa240=Map(mwa240list[mwaidx]);mwa161=Map(mwa161list[mwaidx]);mwa108=Map(mwa108list[mwaidx])
    lev2=np.array([80,90])*u.percent;ax01.set_title('MWA: '+mwa108list[mwaidx].split('_')[3].split('.fits')[0])
    #dcp_data=mapvla_v[i][0].data*100/mapvla_i[i].data*(mapvla_v[i][0].data/np.nanmax(mapvla_v[i][0].data));dcp_data[np.isnan(dcp_data)]=0;dcp_data[np.where(dcp_data<0)]=0;dcp_data[np.where(dcp_data>100)]=0
    xlvla=mwa240.center.Tx.value-50.0*int(mwa240.data.shape[0]/2);xrvla=mwa240.center.Tx.value+50.0*int(mwa240.data.shape[0]/2);ylvla=mwa240.center.Ty.value-50.0*int(mwa240.data.shape[1]/2);yrvla=mwa240.center.Ty.value+50.0*int(mwa240.data.shape[0]/2)
    if(mwa240.data.max()!=0):
        mwa240.draw_contours(axes=ax00,levels=lev2,origin='lower',colors='red',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mwa240.draw_contours(axes=ax01,levels=lev2,origin='lower',colors='red',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mwa240.draw_contours(axes=ax10,levels=lev2,origin='lower',colors='red',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mwa240.draw_contours(axes=ax11,levels=lev2,origin='lower',colors='red',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    if(mwa161.data.max()!=0):
        mwa161.draw_contours(axes=ax00,levels=lev2,origin='lower',colors='green',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mwa161.draw_contours(axes=ax10,levels=lev2,origin='lower',colors='green',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mwa161.draw_contours(axes=ax01,levels=lev2,origin='lower',colors='green',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mwa161.draw_contours(axes=ax11,levels=lev2,origin='lower',colors='green',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    if(mwa108.data.max()!=0):
        mwa108.draw_contours(axes=ax00,levels=lev2,origin='lower',colors='blue',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mwa108.draw_contours(axes=ax01,levels=lev2,origin='lower',colors='blue',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mwa108.draw_contours(axes=ax10,levels=lev2,origin='lower',colors='blue',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
        mwa108.draw_contours(axes=ax11,levels=lev2,origin='lower',colors='blue',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    ax00.set_xlim([400,1200]);ax00.set_ylim([-800,0]);ax01.set_xlim([400,1200]);ax01.set_ylim([-800,0]);ax10.set_xlim([400,1200]);ax10.set_ylim([-800,0]);ax11.set_xlim([400,1200]);ax11.set_ylim([-800,0])
    plt.savefig('/media/rohit/MWA/20140914/pngs/img_'+str(ii)+'.png')
    #ax0.set_xlim([-850,-700]);ax0.set_ylim([200,350])
    #ax0.text(-840,210,str(np.round(dd0.meta['crval3']/1.e9,4))+' GHz',color='blue')
    plt.close()

list171fits=sorted(glob.glob('/media/rohit/MWA/20140914/EUV/fits/*171*.fits'))
m0=Map(list171fits[0])
for i in range(200,len(list171fits)):
    m=Map(list171fits[i])#;m.plot(vmin=0,vmax=200)
    ii="%04d"%i#plt.xlim([400,1200]);plt.ylim([-800,0])
    plt.imshow(m.data-m0.data,aspect='auto',origin=0,cmap='sdoaia171',vmin=0,vmax=50)
    plt.savefig('/media/rohit/MWA/20140914/pngs/aia171_'+str(ii)+'.png')
    plt.close()

