# Extrapolation
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from sunpy.map import Map
from reproject import reproject_interp
import sunpy
from sunpy.coordinates import get_body_heliographic_stonyhurst
from scipy.io import readsav
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from surya.plot import progressbar
from astropy import units as u
from astropy.wcs import utils


def get_gxs_sav2carr_map(index,dat,date,num,hres):
    '''
    Input: 
    index of the sav file from gx_simulator file
    dat: data to be transformed
    date: str E.g. 2016-04-09T18:34:17.00
    num: Number of the layer
    hres: height resolution in km
    Output: Sunpy Map in the Carrington Coordinates
    '''
    crval1=index['crval1'];crval2=index['crval2'];crpix1=index['crpix1'];crpix2=index['crpix2'];ctype1=index['ctype1'];ctype2=index['ctype2'];cdelt1=index['cdelt1'];cdelt2=index['cdelt2'];cunit1=index['cunit1'];cunit2=index['cunit2']
    hdu = fits.PrimaryHDU(dat);list_all=list(index.dtype.names);list_all.remove('COMMENT');list_all.remove('HISTORY');list_all.remove('SIMPLE');list_all.remove('BITPIX');list_all.remove('NAXIS');list_all.remove('DATE_D$OBS')
    index['WCSNAME'],index['CTYPE1'],index['CUNIT1'],index['CTYPE2'],index['CUNIT2']=['Carrington-Heliographic'],['CRLN-CEA'],['deg'],['CRLT-CEA'],['deg']
    index['DATE_OBS']=[date]
    ii=0
    for idx in list_all:
        hdu.header.update({str(idx):index[list_all[ii]][0]})
        ii=ii+1
    source_height=695700000+num*hres*1000.
    if 'RSUN_REF' in hdu.header:
        hdu.header.remove('RSUN_REF')
        hdu.header.append(('RSUN_REF',source_height))
    else:
        hdu.header.append(('RSUN_REF',source_height))
    hhdu=hdu.header
    #hdul = fits.HDUList([hdu])
    mymap=Map(dat,hhdu);hp_coord=mymap.reference_coordinate.transform_to(frames.Helioprojective(observer="earth"))
    out_shape = (mymap.data.shape[0],mymap.data.shape[1]);out_header = sunpy.map.make_fitswcs_header(mymap.data,hp_coord)
    out_wcs = WCS(out_header)
    output, footprint = reproject_interp(mymap, out_wcs, out_shape)
    outmap = sunpy.map.Map((output, out_header));outmap.plot_settings = mymap.plot_settings
    return outmap,mymap

def get_gxs_sav2hpp(exfile,date):
    '''
    Input:
    filename of gx_simulator file
    date: str E.g. 2016-04-09T18:34:17.00
    dr: height resolution in km
    Output:
    List of the sunpy maps
    '''
    mapex=readsav(exfile);expolB=mapex['box']
    bx,by,bz=expolB['bx'][0],expolB['by'][0],expolB['bz'][0];babs=np.sqrt(bx*bx+by*by+bz*bz)
    index=expolB['index'][0];hnum=bx.shape[0];sunpylist=[0]*hnum;dr=expolB['dr'][0][0]*696340.;carrmap=[0]*hnum
    for i in range(hnum):    
        sunpylist[i],carrmap[i]=get_gxs_sav2carr_map(index,babs[i],date,i,dr)
        Bar=progressbar.ProgressBar(hnum,i);progressbar.ShowBar(Bar)
    return sunpylist,carrmap

def get_fieldlines(file_name):
    file_ = open(file_name)
    with open(file_name) as f:
        lines=f.readlines()
    numline = len(lines)-1
    bx=[0]*numline;by=[0]*numline;bz=[0]*numline
    x=[0]*numline;y=[0]*numline;z=[0]*numline
    vtkidx=[0]*numline;i=0
    for row in lines[1:]:
        row=row.split(',')
        vtkidx[i]=float(row[1])
        bx[i]=float(row[1]);by[i]=float(row[2]);bz[i]=float(row[3])
        x[i]=float(row[5]);y[i]=float(row[6]);z[i]=float(row[7])
        i=i+1
    vtkidx=np.array(vtkidx)
    x=np.array(x);y=np.array(y);z=np.array(z)
    bx=np.array(bx);by=np.array(by);bz=np.array(bz)
    return x,y,z,bx,by,bz

def transform_fieldlines(x,y,z,bx,by,bz,date,out_wcs,mymap):
    '''
    Input:
    x,y,z,bx,by,bz
    Carrington Map
    out_wcs: Output WCS
    date: E.g. 2016/04/09T18:45:00
    Output:x,y,z,bx,by,bz,b_hp,b_proj
    '''
    dd1=mymap.pixel_to_world(x[1:]*u.pix,y[1:]*u.pix);dd=dd1.transform_to("heliocentric")
    xkm=dd.cartesian.x.value;ykm=dd.cartesian.y.value;zkm=dd.cartesian.z.value
    #source_height=695700000+num*hres*1000.
    sc = SkyCoord(xkm*u.km, ykm*u.km, zkm*u.km,obstime=date, observer="earth", frame="heliocentric")
    b_hp=sc.transform_to(frames.Helioprojective(observer='earth'))
    b_hp_pix = utils.skycoord_to_pixel(b_hp, out_wcs)
    b_carr=sc.transform_to(frames.HeliographicCarrington(observer='earth'))
    b_carr_pix = utils.skycoord_to_pixel(b_carr, mymap.wcs)
    return x,y,z,bx,by,bz,b_hp,b_hp_pix,b_carr,b_carr_pix



def transform_fieldlines_2SO(x,y,z,bx,by,bz,date,out_wcs,mymap,euimap):
    '''
    Input:
    x,y,z,bx,by,bz
    Carrington Map
    out_wcs: Output WCS
    date: E.g. 2016/04/09T18:45:00
    Output:x,y,z,bx,by,bz,b_hp,b_proj
    '''
    dd1=mymap.pixel_to_world(x[1:]*u.pix,y[1:]*u.pix);dd=dd1.transform_to("heliocentric")
    xkm=dd.cartesian.x.value;ykm=dd.cartesian.y.value;zkm=dd.cartesian.z.value
    #source_height=695700000+num*hres*1000.
    sc = SkyCoord(xkm*u.km, ykm*u.km, zkm*u.km,obstime=date, observer="earth", frame="heliocentric")
    b_carr=sc.transform_to(frames.HeliographicCarrington(observer='earth'))
    b_carr_pix = utils.skycoord_to_pixel(b_carr, mymap.wcs)
    ####
    b_hp=sc.transform_to(euimap.coordinate_frame)
    b_hp_pix = utils.skycoord_to_pixel(b_hp, out_wcs)
    return x,y,z,bx,by,bz,b_hp,b_hp_pix,b_carr,b_carr_pix

