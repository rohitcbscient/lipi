from surya.utils import main as ut
from scipy.io import readsav
import numpy as np

def read_submap(f,ff):
    data_=readsav(f)
    data=data_['submap'+str(ff)]
    n=len(data)
    ts=[0]*n
    map_=[0]*n
    ti=[0]*n
    for i in range(n):
        map_[i]=data[i][0]
        time=data[i][5]
        ts[i]=ut.hms2sec_c(time)
        ti[i]=time.split(' ')[1]
    return map_,ts,ti,data[0][1],data[0][2]

def get_submap(aiafile,res,w):
    aiamap,aiats,aiatime,xaia,yaia=read_submap(aiafile,w)
    xlaia_,ylaia_=xaia-res*aiamap[0].shape[0]/2,yaia-res*aiamap[0].shape[1]/2
    xraia_,yraia_=xaia+res*aiamap[0].shape[0]/2,yaia+res*aiamap[0].shape[1]/2
    xarrayaia=np.linspace(xlaia_,xraia_,aiamap[0].shape[0])
    yarrayaia=np.linspace(ylaia_,yraia_,aiamap[0].shape[1])
    return aiamap,xarrayaia,yarrayaia

def get_submap_all(aialist,wav):
    aiam=[0]*len(wav)
    i=0
    for w in wav:
        aiam[i]=read_submap(aialist[i],w)[0]
        i=i+1
    return aiam


def get_nearest(xl,xr,yl,yr,xarrayaia,yarrayaia):
    xlaia,xraia=ut.find_nearest(xarrayaia,xl)[0],ut.find_nearest(xarrayaia,xr)[0]
    ylaia,yraia=ut.find_nearest(yarrayaia,yl)[0],ut.find_nearest(yarrayaia,yr)[0]
    return xlaia,xraia,ylaia,yraia

def get_submap_crop(aiafile,res,w,xl,xr,yl,yr):
    aiamap,xarrayaia,yarrayaia=get_submap(aiafile,res,w)
    xlaia,xraia,ylaia,yraia=get_nearest(xl,xr,yl,yr,xarrayaia,yarrayaia)
    n=len(aiamap)
    cmap=[0]*n
    for i in range(n):
        cmap[i]=aiamap[i][ylaia:yraia,xlaia:xraia]
    return cmap

