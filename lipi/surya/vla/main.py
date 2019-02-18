import numpy as np
from surya.utils import main as ut

def get_map(f):
    h,d=ut.read_fits(f)
    #d=d[::-1,:] # RA:Y DEC:X
    return d,h

def get_map_cube(f):
    h,d=ut.read_fits(f)
    #d=d[:,::-1,:] # RA:Y DEC:X
    return d,h

def get_submap_cube(f,xl,xr,yl,yr,res,dec):
    h,d=ut.read_fits(f)
    xl_,xr_,yl_,yr_=get_crop_idx(d[0],xl,xr,yl,yr,res,dec)
    d_=d[:,yl_:yr_,xl_:xr_]
    #d_=d_[:,::-1,:] # RA:Y DEC:X
    return d_,h

def get_submap(f,xl,xr,yl,yr,res,dec):
    h,d=ut.read_fits(f)
    xl_,xr_,yl_,yr_=get_crop_idx(d,xl,xr,yl,yr,res,dec)
    d_=d[yl_:yr_,xl_:xr_]
    #d_=d_[::-1,:] # RA:Y DEC:X
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

def fit_centroid(d,l):
    bimage=ut.get_bimage(d,l/np.max(d))
    xc,yc,w,l,ang=ut.fitEllipse(bimage)
    return xc,yc,w,l,ang

def fit_centroid_max(d,ll):
    bimage=ut.get_bimage(d,ll)
    xc,yc,w,l,ang=ut.fitEllipse(bimage)
    return xc,yc,w,l,ang

def fit_centroid_map(d,l,xl,yl,pix,dec):
    xc=[0]*len(l)
    yc=[0]*len(l)
    l_=[0]*len(l)
    w=[0]*len(l)
    ang=[0]*len(l)
    if(len(d.shape)==3):
        ds=d.shape[0]
    if(len(d.shape)==2):
        ds=1
    j=0
    for ll in l:
        xc[j]=[0]*ds
        yc[j]=[0]*ds
        l_[j]=[0]*ds
        w[j]=[0]*ds
        ang[j]=[0]*ds
        for i in range(ds):
            xc[j][i],yc[j][i],w[j][i],l_[j][i],ang[j][i]=fit_centroid_max(d[i],ll)
        j=j+1
    xoff,yoff=get_offset(dec)
    xc=xl+np.array(xc)*pix*xoff
    yc=yl+np.array(yc)*pix*yoff
    l=np.array(l_)
    w=np.array(w)
    ang=np.array(ang)
    return xc,yc,l,w,ang



