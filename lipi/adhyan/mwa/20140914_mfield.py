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
aiamap=Map('/media/rohit/MWA/20140914/EUV/fits/aia.lev1.171A_2014-09-14T01_57_35.34Z.image_lev1.fits')

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

####### Radio Contours #############
aa=pickle.load(open('/media/rohit/MWA/20140914/Tb_20140914_187-188_sub_test.p','rb'),encoding='bytes');Tb240=np.array(aa[0]);maxX=[0]*len(Tb240);maxY=[0]*len(Tb240);xc90_240=[0]*len(Tb240);yc90_240=[0]*len(Tb240)
aa=pickle.load(open('/media/rohit/MWA/20140914/Tb_20140914_125-126_sub_test.p','rb'),encoding='bytes');Tb160=np.array(aa[0]);xc90_160=[0]*len(Tb160);yc90_160=[0]*len(Tb160)
aa=pickle.load(open('/media/rohit/MWA/20140914/Tb_20140914_084-085_sub_test.p','rb'),encoding='bytes');Tb108=np.array(aa[0]);xc90_108=[0]*len(Tb108);yc90_108=[0]*len(Tb108)
xc90_240_,yc90_240_,xc90_160_,yc90_160_,xc90_108_,yc90_108_=0,0,0,0,0,0
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

plot_one=1
if plot_one:
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=aiamap)
    p0=aiamap.plot(axes=ax,aspect='auto')
    hp_lon0 =  bshp0x* u.arcsec;hp_lat0 = bshp0y* u.arcsec
    hp_lon1 = bshp1x * u.arcsec;hp_lat1 = bshp1y* u.arcsec
    seeds00 = SkyCoord(hp_lon0.ravel(), hp_lat0.ravel(),frame=aiamap.coordinate_frame)
    ax.plot_coord(seeds00, color='brown', marker='s', markersize=0.5,alpha=0.8,linestyle='None')
    seeds01 = SkyCoord(hp_lon1.ravel(), hp_lat1.ravel(),frame=aiamap.coordinate_frame)
    ax.plot_coord(seeds01, color='black', marker='s', markersize=0.5,alpha=0.8,linestyle='None')
    seeds1 = SkyCoord(xc90_240*u.arcsec, yc90_240*u.arcsec,frame=aiamap.coordinate_frame)
    ax.plot_coord(seeds1, color='tab:red', marker='s', markersize=10.0,alpha=0.6,linestyle='None')
    seeds2 = SkyCoord(xc90_160*u.arcsec, yc90_160*u.arcsec,frame=aiamap.coordinate_frame)
    ax.plot_coord(seeds2, color='tab:green', marker='s', markersize=10.0,alpha=0.6,linestyle='None')
    seeds3 = SkyCoord(xc90_108.mean()*u.arcsec, yc90_108.mean()*u.arcsec,frame=aiamap.coordinate_frame)
    ax.plot_coord(seeds3, color='tab:blue', marker='s', markersize=10.0,alpha=1.0,linestyle='None')
    seeds4 = SkyCoord(710*u.arcsec, -314*u.arcsec,frame=hmimap.coordinate_frame)
    ax.plot_coord(seeds4, color='tab:cyan', marker='o', markersize=10.0,alpha=1.0,linestyle='None')
    #bl=SkyCoord(-850*u.arcsec,300*u.arcsec,frame=aiamap.coordinate_frame)
    #tr=SkyCoord(-550*u.arcsec,600*u.arcsec,frame=aiamap.coordinate_frame)
    #blx=aiamap.world_to_pixel(bla)[0].value;bly=aiamap.world_to_pixel(bla)[1].value
    #trx=aiamap.world_to_pixel(tra)[0].value;rty=aiamap.world_to_pixel(tra)[1].value
    #ax.set_xlim([blx,trx]);ax.set_ylim([bly,rty])
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
    ax0.plot(tmin,xc90_240,'o-',label='240 MHz');ax0.plot(tmin,xc90_160,'o-',label='160 MHz');ax0.plot(tmin,xc90_108,'o-',label='108 MHz')
    ax1.plot(tmin,yc90_240,'o-',label='240 MHz');ax1.plot(tmin,yc90_160,'o-',label='160 MHz');ax1.plot(tmin,yc90_108,'o-',label='108 MHz')
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
