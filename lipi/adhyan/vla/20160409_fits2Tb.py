import numpy as np
from astropy.io import fits
from surya.radio import get_maps as mp
from sunpy import sun
import astropy.units as u
import sys
import glob
from sunpy.map import Map
from astropy.coordinates import SkyCoord
import astropy.units as u
from surya.utils import main as ut
import matplotlib as mpl
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt


def produce_tstring(mapp):
    date=mapp.date
    hhmmss=' '+str(date.hour)+':'+str(date.minute)+':'+str(date.second)
    sec=ut.hms2sec_c(hhmmss)
    return sec

def get_sunpy_maps(f):
    print 'Reading...'+f[0]
    n=len(f);maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(n):
        maplist[i]=Map(f[i])
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time

def get_evla_submap(f,xbl,ybl,xtr,ytr):
    print 'Reading...'+f[0]
    #xcen=(xbl+xtr)*0.5;ycen=(ybl+ytr)*0.5
    n=len(f);maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(n):
        g=fits.open(f[i])
        g[0].header['CRVAL1']=0;g[0].header['CRVAL2']=0
        map_=Map(g[0].data,g[0].header)
        #cent_pix=map_.world_to_pixel(SkyCoord(xcen*u.arcsec, ycen*u.arcsec, frame=map_.coordinate_frame)) 
        bl = SkyCoord(xbl*u.arcsec, ybl*u.arcsec, frame=map_.coordinate_frame)
        tr = SkyCoord(xtr*u.arcsec, ytr*u.arcsec, frame=map_.coordinate_frame)
        maplist[i]=map_.submap(bl,tr)
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time

#modify_fits=0
#if(modify_fits):

def modify_fits(f):
    data=fits.open(f,mode='update')
    h=data[0].header;d=data[0].data[0][0]
    # 1sec for Tb submap (141,609) pixel: (-582",500")
    p0=-26.22
    nsize=1024
    nd=np.zeros((nsize,nsize))
    #nd[440:444,120:124]=1 # 442-yaxis, 122 x-axis
    nd[314:570,0:252]=d[:,4:]
    rnd=mp.rotateimage(nd, p0, nsize/2, nsize/2)
    h['ctype1'] = 'HPLN-TAN';h['ctype2'] = 'HPLT-TAN'
    h['cdelt1']=2.0;h['cdelt2']=2.0
    h['cunit1']='arcsec';h['cunit2']='arcsec'
    h['CRPIX1']=nsize/2;h['CRPIX2']=nsize/2;h['CRVAL1']=0;h['CRVAL2']=0
    h.append(('p_angle', p0))
    h.append(('dsun_obs', sun.sunearth_distance(h['date-obs']).to(u.meter).value))
    h.append(('rsun_obs', sun.solar_semidiameter_angular_size(h['date-obs']).value))
    h.append(('rsun_ref', sun.constants.radius.value))
    h.append(('hgln_obs', 0.))
    h.append(('hglt_obs', sun.heliographic_solar_center(h['date-obs'])[1].value))
    h['BUNIT'] = 'K';h['BTYPE'] = 'Brightness Temperature'
    bmaj0=np.radians(h['bmaj']);bmin0=np.radians(h['bmin'])
    #h['bmaj']=bmaj0;h['bmin']=bmin0
    beam_area = bmaj0 * bmin0 * np.pi / (4. * np.log(2.))
    k_b = 1.38e-23;c_l = 3.e8;nu=h['crval3']
    factor = 2. * k_b * nu ** 2 / c_l ** 2  # SI unit
    jy_to_si = 1e-26;factor2=100
    Tbrnd=rnd*jy_to_si / beam_area / factor * factor2
    data[0].data=Tbrnd.reshape(1,1,nsize,nsize)
    data.flush()
    data.close()

filelist=sorted(glob.glob('*.FITS'))
for i in range(len(filelist)):
    modify_fits(filelist[i])


sys.exit()
list1600=sorted(glob.glob('/media/rohit/VLA/20160409_EUV/1600/s*fts'))
map1600,data1600,time1600=get_sunpy_maps(list1600)
mapvla,datavla,timevla=get_evla_submap(['test1.FITS'],-1000,100,-600,400)
y1=map1600[15].plot()
bl = SkyCoord(-1000*u.arcsec, 100*u.arcsec, frame=mapvla[0].coordinate_frame)
tr = SkyCoord(-600*u.arcsec, 400*u.arcsec, frame=mapvla[0].coordinate_frame)
x1=mapvla[0].submap(bl,tr)
x1.draw_contours(levels=[50,60,70,80,90]*u.percent,colors='r',linewidths=3,extent=y1.get_extent())
plt.show()

