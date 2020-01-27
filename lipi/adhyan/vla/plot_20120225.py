import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from surya.plot import main as pl
from surya.utils import main as ut
from sunpy.map import Map
from surya.rhessi import main as rh
import astropy.units as u

plt.style.use('/home/i4ds1807205/scripts/general/plt_style.py')

def plot_hmi_rhessi(hmimap,rhmap,xl,xr,yl,yr):
    hmi_=Map(hmimap)
    ax = plt.subplot(projection=hmi_)
    hmi_.plot(axes=ax,aspect='auto',cmap='gray',vmin=-100,vmax=100)
    rhxlpix,rhxrpix,rhylpix,rhyrpix=rh.get_pix(hmimap,xl,xr,yl,yr)
    rhmap.draw_contours(levels=np.linspace(0,90,10)*u.percent,extent=[rhxlpix,rhxrpix,rhylpix,rhyrpix],axes=ax)
    plt.xlim([int(rhxlpix),int(rhxrpix)])
    plt.ylim([int(rhylpix),int(rhyrpix)])
    plt.show()

def hmi_euv_map_paper(hmi,ccmap,ccmap1,xc,yc,lev_1,lev_2,xl,xr,yl,yr):
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(hmi,origin=True,extent=[xl,xr,yl,yr],cmap='gray',interpolation='none',vmin=-400,vmax=400)
    #plt.contourf(ccmap1/np.max(ccmap1),extent=[xl,xr,yl,yr],levels=lev_2,cmap='Oranges',alpha=0.1)
    plt.contour(ccmap/np.max(ccmap),extent=[xl,xr,yl,yr],linewidth=12,levels=lev_1,colors='red')
    #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    plt.colorbar(im,label='B (Gauss)')
    plt.legend()
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    #plt.title(t)
    ax.grid(True)
    plt.show()

def euv_map_paper(ccmap,ccmap1,xc,yc,lev_1,lev_2,xl,xr,yl,yr):
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(ccmap/np.max(ccmap),origin=True,extent=[xl,xr,yl,yr],cmap='sdoaia94',interpolation='none',vmin=0.002,vmax=0.5)
    #plt.contourf(ccmap1/np.max(ccmap1),extent=[xl,xr,yl,yr],levels=lev_2,cmap='Oranges',alpha=0.1)
    #plt.contour(ccmap1/np.max(ccmap1),extent=[xl,xr,yl,yr],linewidths=8,levels=lev_2,colors='k')
    #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    #plt.colorbar(ss,label='Frequency (MHz)')
    plt.legend()
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    #plt.title(t)
    ax.grid(True)
    plt.show()

def hmi_map_inv_lines(hmi,ccmap,xc,yc,lev_1,xl,xr,yl,yr):
    a,b,c,d,e,f=68,77,79,103,125,159
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(hmi,origin=True,extent=[xl,xr,yl,yr],cmap='gray',interpolation='none',vmin=-400,vmax=400)
    plt.contour(hmi,levels=[0.0],extent=[xl,xr,yl,yr],colors='blue',alpha=0.4)
    #plt.contour(euv94/np.max(euv94),extent=[xl,xr,yl,yr],levels=np.linspace(0.4,0.8,2),colors='blue')
    #plt.contour(euv131/np.max(euv131),extent=[xl,xr,yl,yr],levels=lev_1,colors='green')
    #plt.contour(euv171/np.max(euv171),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.contour(euv304/np.max(euv304),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.text(xl+3, yl+5, '94 $\AA$', color='blue',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+10, '131 $\AA$', color='green',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+15, '171 $\AA$', color='red',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+15, '304 $\AA$', color='red',bbox=dict(facecolor='white', alpha=0.1))
    #plt.contour(euv193/np.max(euv193),extent=[xl,xr,yl,yr],levels=lev_1,colors='orange')
    #plt.contour(euv211/np.max(euv211),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.contour(euv1600/np.max(euv1600),extent=[xl,xr,yl,yr],levels=lev_1,colors='black')
    #plt.contour(euv1700/np.max(euv1700),extent=[xl,xr,yl,yr],levels=lev_1,colors='blue')
    plt.contourf(ccmap/np.max(ccmap),extent=[xl,xr,yl,yr],alpha=0.6,cmap='YlGn')
    #plt.contour(rhmap_l,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='blue')
    #plt.contour(rhmap_h,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='green')
    #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    #plt.colorbar(ss,label='Frequency (MHz)')
    plt.legend()
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    #plt.title(t)
    ax.grid(True)
    plt.show()

def centroid_map(hmi,ccmap,xc,yc,lev_1,xl,xr,yl,yr):
    a,b,c,d,e,f=68,77,79,103,125,159
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(hmi,origin=True,extent=[xl,xr,yl,yr],cmap='gray',interpolation='none',vmin=-400,vmax=400)
    #plt.contour(euv94/np.max(euv94),extent=[xl,xr,yl,yr],levels=np.linspace(0.4,0.8,2),colors='blue')
    #plt.contour(euv131/np.max(euv131),extent=[xl,xr,yl,yr],levels=lev_1,colors='green')
    #plt.contour(euv171/np.max(euv171),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.contour(euv304/np.max(euv304),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.text(xl+3, yl+5, '94 $\AA$', color='blue',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+10, '131 $\AA$', color='green',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+15, '171 $\AA$', color='red',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+15, '304 $\AA$', color='red',bbox=dict(facecolor='white', alpha=0.1))
    #plt.contour(euv193/np.max(euv193),extent=[xl,xr,yl,yr],levels=lev_1,colors='orange')
    #plt.contour(euv211/np.max(euv211),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.contour(euv1600/np.max(euv1600),extent=[xl,xr,yl,yr],levels=lev_1,colors='black')
    #plt.contour(euv1700/np.max(euv1700),extent=[xl,xr,yl,yr],levels=lev_1,colors='blue')
    plt.contourf(ccmap/np.max(ccmap),extent=[xl,xr,yl,yr],alpha=0.4,cmap='YlGn')
    #plt.contour(rhmap_l,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='blue')
    #plt.contour(rhmap_h,extent=[xl,xr,yl,yr],levels=lev,linewidths=2,colors='green')
    #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    #plt.colorbar(ss,label='Frequency (MHz)')
    plt.errorbar(xc[a].mean(),yc[a].mean(),xerr=xc[a].std(),yerr=yc[a].std(),elinewidth=5,c='b',label='A')
    plt.errorbar(xc[b].mean(),yc[b].mean(),xerr=xc[b].std(),yerr=yc[b].std(),elinewidth=5,c='c',label='B')
    plt.errorbar(xc[c].mean(),yc[c].mean(),xerr=xc[c].std(),yerr=yc[c].std(),elinewidth=5,c='y',label='C')
    plt.errorbar(xc[d].mean(),yc[d].mean(),xerr=xc[d].std(),yerr=yc[d].std(),elinewidth=5,c='m',label='D')
    plt.errorbar(xc[e].mean(),yc[e].mean(),xerr=xc[e].std(),yerr=yc[e].std(),elinewidth=5,c='g',label='E')
    plt.errorbar(xc[f].mean(),yc[f].mean(),xerr=xc[f].std(),yerr=yc[f].std(),elinewidth=5,c='r',label='F')
    plt.legend()
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    #plt.title(t)
    ax.grid(True)
    plt.show()

if __name__=='__main__':
    main();
else:
    print 'plotting module for 20120225....'


