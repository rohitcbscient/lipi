from surya.utils import main as ut
from sunpy.map import Map
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel as sp
import astropy.units as u
from astropy.coordinates import SkyCoord
import sys

def get_submap(f,xl,xr,yl,yr):
    '''
    Inputs: fits file, xl, xr, yl, yr (coordinates for cropping)
    Outputs: submap, xl, xr, yl, yr (in pixels)
    '''
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

def get_pixels(x,y,coord_frame,w):
    '''
    Inputs: x,y (wcs coordinate in arcsec),coordinate frame (e.g. sunpy.coordinates), WCS keywords
    Outputs: submap, xl, xr, yl, yr (in pixels)
    '''
    xcor=x*u.arcsec
    ycor=y*u.arcsec
    bl=SkyCoord(xcor, ycor, frame=coord_frame)
    xlpix,ylpix=sp(bl,w)
    return xlpix,ylpix

def get_submap_hmi(f,xl,xr,yl,yr,ins):
    '''
    ONLY FOR HMI
    Inputs: fits file, xl, xr, yl, yr (coordinates for cropping)
    Outputs: submap, xl, xr, yl, yr (in pixels)
    '''
    if(ins!='hmi'):
        print('Use get_submap for non-hmi maps... Existing..')
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


def read_gxs_magstr(box):
    '''
    Reads the gx simulator magnetic extrapolated structures
    Inputs:
    Outputs:
    '''
    

