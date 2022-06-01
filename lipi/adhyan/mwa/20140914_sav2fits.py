import matplotlib as mpl
mpl.use('TkAgg') 
import matplotlib.pyplot as plt
from astropy.io import fits
from sunpy.map import Map
import glob
import pickle
from astropy.coordinates import SkyCoord
import astropy.units as u
from surya.utils import main as ut
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import os
from scipy.io import readsav
from astropy.time import Time
import sunpy
from sunpy import sun
from dateutil import parser
import numpy as np
from reproject import reproject_interp


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
    header['HGLT_OBS'] = sunpy.coordinates.get_horizons_coord('SUN',time=header['DATE-OBS']).lat.value#sunpy.coordinates.get_sun_B0(header['DATE-OBS']).value
    header['RSUN_REF'] = sunpy.sun.constants.radius.value
    header['RSUN_OBS'] = 958.11#sunpy.coordinates.sun.angular_radius(header['DATE-OBS']).value
    header['DSUN_OBS'] = 147.86e9#sunpy.coordinates.get_sunearth_distance(header['DATE-OBS']).to(u.meter).value
    d=reduce_size(dat,n)
    smp = smap.Map(d, header)
    return smp


def produce_tstring(date):
    #date=mapp.date
    #hhmmss=' '+str(date.hour)+':'+str(date.minute)+':'+str(date.second)+'.'+str(date.microsecond/1.e6).split('.')[1]
    hhmmss=' '+str(date.split('T')[1].split(':')[0])+':'+str(date.split('T')[1].split(':')[1])+':'+str(date.split('T')[1].split(':')[2])
    sec=ut.hms2sec_c(hhmmss)
    return sec

def get_sunpy_maps_rot(f,m,wave,inst,filename):
    print('Reading...'+f[0])
    n=len(f)
    for i in range(n):
        ii="%04d"%i
        if(inst=='HMI'):
            maplist=idl2sunpy_hmi(f[i],m)
        if(inst=='AIA'):
            maplist=idl2sunpy_sdo(f[i],m,wave,inst)
        tt=maplist.meta['date-obs'].split('T')[1];ti=produce_tstring(maplist.meta['date-obs']);
        pickle.dump(maplist,open(str(filename)+'_'+str(tt.replace(':','-'))+'_'+str(int(ti))+'.sunpy','wb'))

def get_sunpy_maps(f):
    print('Reading...'+f[0])
    n=len(f);maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(n):
        maplist[i]=Map(f[i])
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time

def get_sunpy_basediff_maps(f,m,wave,inst,filename):
    print('Reading...'+f[0])
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
    print('Reading...'+f[0])
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

dump_submaps=1
if(dump_submaps):
    list193=sorted(glob.glob('/media/rohit/MWA/20140914/EUV/fits/*193*rot.sav'))[0:5]
    list171=sorted(glob.glob('/media/rohit/MWA/20140914/EUV/fits/*171*rot.sav'))[0:5]
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

list171fits=sorted(glob.glob('/media/rohit/MWA/20140914/sunpy_maps/pickle/map*171*.p'))
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
    tidx171[i]=int(list171fits[i].split('_')[4][0:2])*3600+int(list171fits[i].split('_')[4][2:4])*60+int(list171fits[i].split('_')[4][4:6]);tidx171_str[i]=list171fits[i].split('_')[4]
for i in range(len(list193fits)):
    tidx193[i]=int(list193fits[i].split('_')[3].split('.')[0]);tidx193_str[i]=list193fits[0].split('_')[2]
for i in range(len(listhmifits)):
    tidxhmi[i]=int(listhmifits[i].split('_')[3].split('.')[0]);tidxhmi_str[i]=listhmifits[0].split('_')[2]
tidx171=np.array(tidx171);tidx193=np.array(tidx193);tidxhmi=np.array(tidxhmi)

######################## MWA #################################3

mwa108list=sorted(glob.glob('/media/rohit/MWA/20140914/mwa_maps/mwa_108*.fits'))
mwa161list=sorted(glob.glob('/media/rohit/MWA/20140914/mwa_maps/mwa_161*.fits'))
mwa240list=sorted(glob.glob('/media/rohit/MWA/20140914/mwa_maps/mwa_240*.fits'))
mwatsec=[0]*len(mwa240list)
for i in range(len(mwa240list)):
    mapp=Map(mwa240list[i]);dd=mapp.date.value
    mwatsec[i]=produce_tstring(dd)
mwatsec=np.array(mwatsec)
listTb=sorted(glob.glob('/media/rohit/MWA/20140914/Tb*.p'))
list_171=sorted(glob.glob('/media/rohit/MWA/20140914/sunpy_maps/pickle/map*171*.p'))

sys.exit()
####################### Extrapolation #######################
    
def get_fieldlines(x,y,z,out_wcs,t):
    dd1=hdu_submap.pixel_to_world(x[1:]*u.pix,y[1:]*u.pix);dd=dd1.transform_to("heliocentric")
    xkm=0;ykm=0;zkm=0;t='01:54:00'
    sc = SkyCoord(xkm*u.km, ykm*u.km, zkm*u.km,obstime="2014/09/14T"+str(t), observer="earth", frame="heliocentric")
    b_hp=sc.transform_to(frames.Helioprojective(observer='earth'))
    b_proj = utils.skycoord_to_pixel(b_hp, out_wcs)
    return x,y,z,bx,by,bz,b_hp,b_proj


from astropy import wcs
w = wcs.WCS(naxis=2)
w.wcs.crpix = [-234.75, 8.3393]
w.wcs.cdelt = np.array([-0.066667, 0.066667])
w.wcs.crval = [0, -90]
w.wcs.ctype = ["CRLN-CEA", "CRLN-CEA"]

from sunpy.coordinates import get_body_heliographic_stonyhurst
pix=[0,0]#np.meshgrid(np.arange(334),np.arange(334))
hpc_arr=hdu_submap.wcs.pixel_to_world(pix[0],pix[1])
hgcarr_arr=hpc_arr.transform_to(frames.HeliographicCarrington(obstime='2014-09-14T01:56:00',observer='earth'))
#hcc_arr=hpc_arr.transform_to(frames.Heliocentric)
#hpc_arr=hgcarr_arr.wcs.pixel_to_world(pix[0],pix[1])

hgcarr_world=hgcarr_map.wcs.pixel_to_world(xceB[0],yceB[0])
hgcarr_world.to_pixel(hgcarr_map.wcs) # Check 

from astropy.time import Time
out_shape=(334,335)
obstime=Time('2014-09-14T01:56:00')
header = make_fitswcs_header(data, frame_out,scale=[360 / shape[1],180 / np.pi * (2 / shape[0]),] * u.deg / u.pix,projection_code="CEA")


out_shape=(334,335)
out_header = sunpy.map.make_fitswcs_header(out_shape,hgcarr_arr,scale=[0.6,0.6]*u.arcsec/u.pix,reference_pixel=[0,0]*u.pix,projection_code="CEA")
out_wcs = WCS(out_header)
hgcarr_map=Map(bz,out_header)
carr_coord=hgcarr_map.wcs.pixel_to_world(xceB[0],yceB[0])
carr_coord_m=SkyCoord(carr_coord.lon,carr_coord.lat, 0*(6.9e5+zceB[0]*726)*u.km, observer='earth',
        obstime=hgcarr_map.observer_coordinate.obstime.value,frame=frames.HeliographicCarrington)
hpc_coord_m=carr_coord_m.transform_to(frames.Helioprojective)

earth = get_body_heliographic_stonyhurst('earth', hdu_submap.date)
out_wcs.heliographic_observer = earth
output, footprint = reproject_interp(hgcarr_map, hdu_submap.wcs, out_shape)
outmap = sunpy.map.Map((output, out_header))




listeB=sorted(glob.glob('/sdata/20140914_hmi/hmi/hmi*0.vtk.csv'))
exmap=readsav('/sdata/20140914_hmi/hmi/hmi_20140914_015615_magnetogram.fitsrot.sav.cube.sav');bz=exmap['bout'][-1][0]
hdubz = fits.PrimaryHDU(bz)


eBid=x[0]*len(listeB);xceB=[0]*len(listeB);yceB=[0]*len(listeB);zceB=[0]*len(listeB);eB=[0]*len(listeB)
for i in range(len(listeB)):
    eBidx[i]=int(listeB[i].split('_')[3][0:2])*3600+int(listeB[i].split('_')[3][2:4])*60+int(listeB[i].split('_')[3][4:6])
    eBcsv=pandas.read_csv(listeB[i]);eB[i]=pandas.DataFrame(eBcsv['Magnitude']).to_numpy().flatten()
    xceB[i]=pandas.DataFrame(eBcsv['Points:0']).to_numpy().flatten()
    yceB[i]=pandas.DataFrame(eBcsv['Points:1']).to_numpy().flatten()
    zceB[i]=pandas.DataFrame(eBcsv['Points:2']).to_numpy().flatten()
bb=pickle.load(open(list_171[0],'rb'))
hmi=sunpy.map.Map('/sdata/20140914_hmi/hmi/hmi_20140914_015445_magnetogram.fits')
hdu=sunpy.map.Map(bb[0])
hhdu=hdu.meta
hp_coord=hdu.reference_coordinate.transform_to(frames.Helioprojective(observer="earth"))
hp_hcc=hdu.reference_coordinate.transform_to(frames.Heliocentric(observer="earth"))
#### PIXEL TO HCC #######

top_right = SkyCoord(800 * u.arcsec, -200 * u.arcsec, frame=hdu.coordinate_frame)
bottom_left = SkyCoord(600 * u.arcsec, -400 * u.arcsec, frame=hdu.coordinate_frame)
hmi_submap = hmi.submap(bottom_left=bottom_left, top_right=top_right)

t='01:54:00'


from astropy.wcs.utils import wcs_to_celestial_frame
from astropy.coordinates.matrix_utilities import matrix_product, matrix_transpose, rotation_matrix
from astropy.coordinates.transformations import _matrix_hcc_to_hpc


dd1=hdu_submap.pixel_to_world(xceB[i]*u.pix,yceB[i]*u.pix);dd=dd1.transform_to("heliocentric")
hcc=sunpy.coordinates.frames.Heliocentric()
sunpy.coordinates.transformations.hpc_to_hcc(dd1,frames.Heliocentric(observer="earth",obstime="2014/09/14T"+str(t)))
sc = SkyCoord(xkm*u.km, ykm*u.km, zkm*u.km,obstime="2014/09/14T"+str(t), observer="earth", frame="heliocentric")
out_header = sunpy.map.make_fitswcs_header((4096,4096),hp_coord)
out_wcs = WCS(out_header)

#### HCC TO HPC ########


#### HPC TO PIXEL ######



out_header = sunpy.map.make_fitswcs_header((4096,4096),hp_coord)
out_wcs = WCS(out_header)
top_right = SkyCoord(800 * u.arcsec, -200 * u.arcsec, frame=hdu.coordinate_frame)
bottom_left = SkyCoord(600 * u.arcsec, -400 * u.arcsec, frame=hdu.coordinate_frame)
hdu_submap = hdu.submap(bottom_left=bottom_left, top_right=top_right)
b_proj=[0]*15;b_proj_pix=[0]*15;res=726.0
#xceB=[20+np.arange(5)]*15;yceB=[20+np.arange(5)]*15;zceB=[np.arange(5)]*15
for i in range(15):
    t=listeB[i].split('_')[3][0:2]+':'+listeB[i].split('_')[3][2:4]+':'+listeB[i].split('_')[3][4:6]
    #x,y,z,bx,by,bz,b_hp,b_proj[i]=get_fieldlines(xceB*res,yceB*res,zceB*res,out_wcs,t)

    dd1=hdu_submap.pixel_to_world(xceB[i]*u.pix,yceB[i]*u.pix);dd=dd1.transform_to("heliocentric")
    xkm=dd.x.value/1000.;ykm=dd.y.value/1000.;zkm=0#dd.z.value/1000.+res*zceB[i]
    sc = SkyCoord(xkm*u.km, ykm*u.km, zkm*u.km,obstime="2014/09/14T"+str(t), observer="earth", frame="heliocentric")
    b_hp=sc.transform_to(frames.Helioprojective(observer='earth'))
    b_proj[i] = utils.skycoord_to_pixel(b_hp, hdu_submap.wcs)

hdu_submap.plot()
plt.plot(b_proj[0][0],b_proj[0][1],'.',markersize=0.5)
plt.show()
#xrange=[600,800]
#yrange=[-400,-200]


# xrange=[705,855]; yrange=[-77,-277]

out_shape = (300, 300)
earth = get_body_heliographic_stonyhurst('earth', mymap.date)
out_wcs.heliographic_observer = earth
output, footprint = reproject_interp(mymap, out_wcs, out_shape)
outmap = sunpy.map.Map((output, out_header))
outmap.plot_settings = mymap.plot_settings
aiamap=Map('/media/rohit/VLA/20160409_EUV/171/ssw_cutout_20160409_184434_AIA_171_.fts')

sys.exit()
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


for i in range(40,950):
    ii="%04d"%i;mwaidx=ut.find_predecessor(mwatsec,tidx171[i])[0]
    map171,map171b,map171r=pickle.load(open(list_171[i],'rb'))
    mwa240=Map(mwa240list[mwaidx]);mwa161=Map(mwa161list[mwaidx]);mwa108=Map(mwa108list[mwaidx])
    f,ax=plt.subplots(1,1,figsize=(14,14),subplot_kw={'projection': map171})
    map171r.plot(axes=ax,vmin=-50,vmax=50,cmap='coolwarm')#;map193b.plot(axes=ax11,vmin=0,vmax=80)
    lev2=np.array([40,70])*u.percent;ax01.set_title('MWA: '+mwa108list[mwaidx].split('_')[3].split('.fits')[0])
    xlvla=mwa240.center.Tx.value-50.0*int(mwa240.data.shape[0]/2);xrvla=mwa240.center.Tx.value+50.0*int(mwa240.data.shape[0]/2);ylvla=mwa240.center.Ty.value-50.0*int(mwa240.data.shape[1]/2);yrvla=mwa240.center.Ty.value+50.0*int(mwa240.data.shape[0]/2)
    if(mwa240.data.max()!=0):
        mwa240.draw_contours(axes=ax,levels=lev2,origin='lower',colors='red',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    if(mwa161.data.max()!=0):
        mwa161.draw_contours(axes=ax,levels=lev2,origin='lower',colors='green',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    if(mwa108.data.max()!=0):
        mwa108.draw_contours(axes=ax,levels=lev2,origin='lower',colors='blue',linewidths=4,extent=[xlvla,xrvla,ylvla,yrvla])#[xlpix,xrpix,ylpix,yrpix])
    plt.xlim([2048,4098]);plt.ylim([0,2048])
    plt.savefig('/media/rohit/MWA/20140914/pngs/subr_'+str(ii)+'.png')
    plt.show()
    #ax0.set_xlim([-850,-700]);ax0.set_ylim([200,350])
    #ax0.text(-840,210,str(np.round(dd0.meta['crval3']/1.e9,4))+' GHz',color='blue')
    plt.close()

for i in range(40,950):
    #f,ax=plt.subplots(2,1,figsize=(8,15));ax0=ax[0];ax1=ax[1] # 171 is the base
    ii="%04d"%i;mwaidx=ut.find_predecessor(mwatsec,tidx171[i])[0]
    map171,map171b,map171r=pickle.load(open(list_171[i],'rb'))
    f,ax=plt.subplots(2,2,figsize=(14,14),subplot_kw={'projection': map171});ax00=ax[0,0];ax01=ax[0,1];ax10=ax[1,0];ax11=ax[1,1]
    map171.plot(axes=ax00);map171b.plot(axes=ax01,vmin=0,vmax=50);map171r.plot(axes=ax10,vmin=0,vmax=50)#;map193b.plot(axes=ax11,vmin=0,vmax=80)
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
    #ax00.set_xlim([400,1200]);ax00.set_ylim([-800,0]);ax01.set_xlim([400,1200]);ax01.set_ylim([-800,0]);ax10.set_xlim([400,1200]);ax10.set_ylim([-800,0]);ax11.set_xlim([400,1200]);ax11.set_ylim([-800,0])
    ax00.set_xlim([2048,1200]);ax00.set_ylim([0,2048])
    plt.savefig('/media/rohit/MWA/20140914/pngs/sub_'+str(ii)+'.png')
    plt.show()
    #ax0.set_xlim([-850,-700]);ax0.set_ylim([200,350])
    #ax0.text(-840,210,str(np.round(dd0.meta['crval3']/1.e9,4))+' GHz',color='blue')
    plt.close()

def mwa_pixel2_aia_pixel(x,y):
    
    return




expt=[0]*len(list_171)
for i in range(len(list_171)):
    map171,map171b,map171r=pickle.load(open(list_171[i],'rb'))
    expt[i]=map171b.exposure_time.value


for i in range(20,1063):
    ii="%04d"%i;mwaidx=ut.find_predecessor(mwatsec,tidx171[i])[0]
    map171,map171b,map171r=pickle.load(open(list_171[i],'rb'))
    mwa240=Map(mwa240list[mwaidx]);mwa161=Map(mwa161list[mwaidx]);mwa108=Map(mwa108list[mwaidx])
    #if(int(map171b.exposure_time.value)==2):
    f,ax=plt.subplots(1,1,figsize=(14,14),subplot_kw={'projection': map171})
    p=map171r.plot(axes=ax,vmin=-10,vmax=10,cmap='spring')
    s=ax.contour(mwa240.data/np.max(mwa240.data),axes=ax,levels=np.array([80,90])/100.,colors='blue',origin='lower',extent=[-2118.67,6214.67,-2118.67,6214.67],linewidths=4)
    s=ax.contour(mwa161.data/np.max(mwa161.data),axes=ax,levels=np.array([80,90])/100.,colors='green',origin='lower',extent=[-2118.67,6214.67,-2118.67,6214.67],linewidths=4)
    s=ax.contour(mwa108.data/np.max(mwa108.data),axes=ax,levels=np.array([80,90])/100.,colors='red',origin='lower',extent=[-2118.67,6214.67,-2118.67,6214.67],linewidths=4)
    ax.set_xlim([2048,4096]);ax.set_ylim([0,2048])
    plt.text(0.1, 0.05, '108 MHz', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes,color='r',fontsize=20)
    plt.text(0.1, 0.1, '161 MHz', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes,color='g',fontsize=20)
    plt.text(0.1, 0.15, '240 MHz', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes,color='b',fontsize=20)
    plt.savefig('/media/rohit/MWA/20140914/pngs/sub_'+str(ii)+'.png')
    plt.close()

list171fits=sorted(glob.glob('/media/rohit/MWA/20140914/EUV/fits/*171*.fits'))
m0=Map(list171fits[0])
for i in range(200,len(list171fits)):
    m=Map(list171fits[i])#;m.plot(vmin=0,vmax=200)
    ii="%04d"%i#plt.xlim([400,1200]);plt.ylim([-800,0])
    plt.imshow(m.data-m0.data,aspect='auto',origin=0,cmap='sdoaia171',vmin=0,vmax=50)
    plt.savefig('/media/rohit/MWA/20140914/pngs/aia171_'+str(ii)+'.png')
    plt.close()

