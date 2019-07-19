from surya.utils import main as ut
from sunpy.map import Map
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel as sp
import astropy.units as u
from astropy.coordinates import SkyCoord
import sys

def get_submap(f,xl,xr,yl,yr):
    h,d=ut.read_fits(f)
    w=WCS(f)
    mymap_=Map(f)
    mymap=mymap_.rotate()
    xcor=[xl,xr]*u.arcsec
    ycor=[yl,yr]*u.arcsec
    bl=SkyCoord(xl*u.arcsec, yl*u.arcsec, frame=mymap.coordinate_frame)
    tr=SkyCoord(xr*u.arcsec, yr*u.arcsec, frame=mymap.coordinate_frame)
    submap=mymap.submap(bl,tr)
    xlpix,ylpix=sp(bl,w)
    xrpix,yrpix=sp(tr,w)
    return submap,xlpix,xrpix,ylpix,yrpix


def get_submap_hmi(f,xl,xr,yl,yr,ins):
    if(ins!='hmi'):
        print 'Use get_submap for non-hmi maps... Existing..'
        sys.exit()
    h,d=ut.read_fits(f)
    w=WCS(f)
    mymap_=Map(f)
    ang=h['CROTA2']*u.deg
    mymap=mymap_.rotate(angle=ang)
    xcor=[xl,xr]*u.arcsec
    ycor=[yl,yr]*u.arcsec
    bl=SkyCoord(xl*u.arcsec, yl*u.arcsec, frame=mymap.coordinate_frame)
    tr=SkyCoord(xr*u.arcsec, yr*u.arcsec, frame=mymap.coordinate_frame)
    submap=mymap.submap(bl,tr)
    #inv=submap.data[::-1,::-1]
    #setattr(submap, 'data_rot', inv)
    xlpix,ylpix=sp(bl,w)
    xrpix,yrpix=sp(tr,w)
    return submap,xlpix,xrpix,ylpix,yrpix


