import numpy as np
from surya.utils import main as ut
from astropy.io import fits

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
    

