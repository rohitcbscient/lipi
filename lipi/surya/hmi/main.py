import numpy as np
from surya.utils import main as ut
from sunpy.map import Map
import astropy.units as u
from astropy.coordinates import SkyCoord

def get_submap(f,xl,xr,yl,yr):
    h,d=ut.read_fits(f)
    d=d[::-1,::-1] # HMI images are inverted in both the axis
    cenpx=int(h[48])
    cenpy=int(h[49])
    del_=h[52]
    lx,ly=d.shape[0],d.shape[1]
    lxarc,lyarc=(np.linspace(0,lx,lx)-cenpx)*del_,(np.linspace(0,ly,ly)-cenpy)*del_
    hmxl=ut.find_nearest(lxarc,xl)[0]
    hmxr=ut.find_nearest(lxarc,xr)[0]
    hmyl=ut.find_nearest(lyarc,yl)[0]
    hmyr=ut.find_nearest(lyarc,yr)[0]
    d_=d[hmyl:hmyr,hmxl:hmxr]
    print hmyl,hmyr,hmxl,hmxr
    return d_

def get_submap_(f,xl,xr,yl,yr):
    h,d=ut.read_fits(f)
    mymap=Map(f)
    xcor=[xl,xr]*u.arcsec
    ycor=[yl,yr]*u.arcsec
    bl=SkyCoord(xl*u.arcsec, yl*u.arcsec, frame=mymap.coordinate_frame)
    tr=SkyCoord(xr*u.arcsec, yr*u.arcsec, frame=mymap.coordinate_frame)
    submap=mymap.submap(bl,tr)
    return submap

