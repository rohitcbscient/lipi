import matplotlib.pyplot as plt
import numpy as np
from surya.plot import main as pl

plt.style.use('/home/i4ds1807205/scripts/general/plt_style.py')


def euv_vla(cmap,vlasubmap,xl,xr,yl,yr):
    lev=np.linspace(0.35,20.95,10)*7
    plt.imshow(cmap[40],origin=True,extent=[xl,xr,yl,yr],cmap='hot')
    plt.contour(vlasubmap[::-1,:],extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='white')
    plt.xlabel('arcsec')
    plt.ylabel('arcsec')
    #plt.title(t)
    plt.grid(True)
    plt.show()

def euv_vla_qs_centroids(ccmap,map_qs,freq,xc,yc,xl,xr,yl,yr,bmaj,bmin,angle,t):
    lev=np.linspace(map_qs.min(),map_qs.max(),10)
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(ccmap,origin=True,extent=[xl,xr,yl,yr],cmap='BuGn')
    #plt.contour(map_qs,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='white')
    ss=ax.scatter(xc,yc,c=freq,s=80,cmap='spring')
    plt.colorbar(ss,label='Frequency (MHz)')
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    pl.add_beam(ax,xl+10,yl+10,bmaj, bmin,angle)
    plt.title(t)
    ax.grid(True)

def euv_vla_rhessi_qs_centroids(ccmap,map_qs,rhmap_l,rhmap_h,lev,freq,xc,yc,xl,xr,yl,yr,bmaj,bmin,angle,t):
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(ccmap,origin=True,extent=[xl,xr,yl,yr],cmap='BuGn')
    plt.contour(rhmap_l,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='blue')
    plt.contour(rhmap_h,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='green')
    ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    plt.colorbar(ss,label='Frequency (MHz)')
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    pl.add_beam(ax,xl+10,yl+10,bmaj, bmin,angle)
    plt.title(t)
    ax.grid(True)

if __name__=='__main__':
    main();
else:
    print 'plotting module for 20120225....'


