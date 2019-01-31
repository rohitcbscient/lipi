import numpy as np
from surya.utils import main as ut

def get_map(f):
    h,d=ut.read_fits(f)
    d=d[::-1,:] # RA:Y DEC:X
    return d,h

def get_map_cube(f):
    h,d=ut.read_fits(f)
    d=d[:,::-1,:] # RA:Y DEC:X
    return d,h

def get_submap_cube(f,xl,xr,yl,yr,res,dec):
    h,d=ut.read_fits(f)
    xl_,xr_,yl_,yr_=get_crop_idx(d[0],xl,xr,yl,yr,res,dec)
    d_=d[:,yl_:yr_,xl_:xr_]
    d_=d_[:,::-1,:] # RA:Y DEC:X
    return d_,h

def get_submap(f,xl,xr,yl,yr,res,dec):
    h,d=ut.read_fits(f)
    xl_,xr_,yl_,yr_=get_crop_idx(d,xl,xr,yl,yr,res,dec)
    d_=d[yl_:yr_,xl_:xr_]
    d_=d_[::-1,:] # RA:Y DEC:X
    return d_,h
    
def get_crop_idx(d,xl,xr,yl,yr,res,dec):
    x_offset,y_offset=get_offset(dec)
    xy=res*np.linspace(-1*d.shape[-2]/2,d.shape[-1]/2,d.shape[0]+1)
    xl_,xr_=ut.find_nearest(xy*x_offset,xl),ut.find_nearest(xy*x_offset,xr)
    yl_,yr_=ut.find_nearest(xy*y_offset,yl),ut.find_nearest(xy*y_offset,yr)
    return xl_[0],xr_[0],yl_[0],yr_[0]

def get_offset(dec):
    x_offset=np.cos(dec*np.pi/180)
    y_offset=(1/np.cos(dec*np.pi/180))
    return x_offset,y_offset


