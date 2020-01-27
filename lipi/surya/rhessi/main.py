import numpy as np
from surya.utils import main as ut
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel as sp
from astropy.wcs import WCS
from sunpy.map import Map

def get_submap(f,tstring,xrh,yrh,rhsize,rhres,xl,xr,yl,yr):
    data_r,hdur=fits.getdata(f,0,header=True)
    xlrh_,ylrh_=xrh-1.0*rhsize/2,yrh-1.0*rhsize/2
    xrrh_,yrrh_=xrh+1.0*rhsize/2,yrh+1.0*rhsize/2
    xarrayrh=np.linspace(xlrh_,xrrh_,rhsize)
    yarrayrh=np.linspace(ylrh_,yrrh_,rhsize)
    xlrh,xrrh=ut.find_nearest(xarrayrh,xl)[0],ut.find_nearest(xarrayrh,xr)[0]
    ylrh,yrrh=ut.find_nearest(yarrayrh,yl)[0],ut.find_nearest(yarrayrh,yr)[0]
    data_low=data_r[:,0,ylrh:yrrh,xlrh:xrrh]
    data_high=data_r[:,1,ylrh:yrrh,xlrh:xrrh]
    tsec=[0]*len(tstring)
    for i in range(len(tstring)):
        tsec[i]=ut.hms2sec_c(' '+tstring[i])
    tsec=np.array(tsec)
    return data_low,data_high,tsec


def get_submap_(rh,xl,xr,yl,yr,reffile):
    xcor=[xl,xr]*u.arcsec
    ycor=[yl,yr]*u.arcsec
    bl=SkyCoord(xl*u.arcsec, yl*u.arcsec, frame=rh.coordinate_frame)
    tr=SkyCoord(xr*u.arcsec, yr*u.arcsec, frame=rh.coordinate_frame)
    submap=rh.submap(bl,tr)
    return submap

def get_pix(reffile,xl,xr,yl,yr):
    refmap=Map(reffile)
    w=WCS(reffile)
    bl_ref=SkyCoord(xl*u.arcsec, yl*u.arcsec, frame=refmap.coordinate_frame)
    tr_ref=SkyCoord(xr*u.arcsec, yr*u.arcsec, frame=refmap.coordinate_frame)
    xlpix,ylpix=sp(bl_ref,w)
    xrpix,yrpix=sp(tr_ref,w)
    return xlpix,xrpix,ylpix,yrpix



