#import mayavi.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from surya.plot import main as pl
from surya.utils import main as ut
from sunpy.map import Map
from surya.rhessi import main as rh
import astropy.units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches

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

def babs_paper(hmi,hm1,hm2,hm3,ccmap,ccmap1,lev_1,lev_2,xl,xr,yl,yr,xl1,xr1,yl1,yr1):
    fig=plt.figure(figsize=(15,15))
    ax1=fig.add_subplot(221,aspect='auto')
    im1=ax1.imshow(hmi,origin=True,extent=[xl1,xr1,yl1,yr1],cmap='jet',interpolation='none',vmin=0,vmax=200)
    #plt.contourf(ccmap1/np.max(ccmap1),extent=[xl,xr,yl,yr],levels=lev_2,cmap='Oranges',alpha=0.1)
    ax1.contour(ccmap/np.max(ccmap),extent=[xl,xr,yl,yr], linestyles='--',linewidth=14,levels=lev_1,colors='k')
    #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    ax1.legend()
    ax1.set_xlabel('Solar X (arcsec)')
    ax1.set_ylabel('Solar Y (arcsec)')
    ax1.set_ylim(yl,yr)
    ax1.set_xlim(xl,xr)
    #plt.title(t)
    rect = patches.Rectangle((485,351),10,7,linewidth=3,edgecolor='k',facecolor='none')
    ax1.add_patch(rect)
    ax1.text(0.7, 0.1, '(C) 734 km', color='k',transform=ax1.transAxes,size=13, backgroundcolor='white', weight='bold')
    ax2=fig.add_subplot(222,aspect='auto')
    im2=ax2.imshow(hm1,origin=True,extent=[xl1,xr1,yl1,yr1],cmap='jet',interpolation='none',vmin=0,vmax=200)
    #plt.contourf(ccmap1/np.max(ccmap1),extent=[xl,xr,yl,yr],levels=lev_2,cmap='Oranges',alpha=0.1)
    ax2.contour(ccmap/np.max(ccmap),extent=[xl,xr,yl,yr],linestyles='--',linewidth=14,levels=lev_1,colors='k')
    #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    ax2.set_xlabel('Solar X (arcsec)')
    ax2.set_ylabel('Solar Y (arcsec)')
    ax2.set_ylim(yl,yr)
    ax2.set_xlim(xl,xr)
    #plt.title(t)
    rect = patches.Rectangle((485,351),10,7,linewidth=3,edgecolor='k',facecolor='none')
    ax2.add_patch(rect)
    ax2.text(0.7, 0.1, '(D) 4400 km', color='k',transform=ax2.transAxes,size=13,backgroundcolor='white', weight='bold')
    ax3=fig.add_subplot(223,aspect='auto')
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    im3=ax3.imshow(hm2,origin=True,extent=[xl1,xr1,yl1,yr1],cmap='jet',interpolation='none',vmin=0,vmax=200)
    #plt.contourf(ccmap1/np.max(ccmap1),extent=[xl,xr,yl,yr],levels=lev_2,cmap='Oranges',alpha=0.1)
    fig.colorbar(im3, label='B (Gauss)', cax=cax, orientation='vertical')
    ax3.contour(ccmap/np.max(ccmap),extent=[xl,xr,yl,yr],linestyles='--',linewidth=14,levels=lev_1,colors='k')
    #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    ax3.set_xlabel('Solar X (arcsec)')
    ax3.set_ylabel('Solar Y(arcsec)')
    ax3.set_ylim(yl,yr)
    ax3.set_xlim(xl,xr)
    #plt.title(t)
    rect = patches.Rectangle((485,351),10,7,linewidth=3,edgecolor='k',facecolor='none')
    ax3.text(0.7, 0.1, '(E) 7340 km', transform=ax3.transAxes,size=13, backgroundcolor='white',weight='bold')
    ax3.add_patch(rect)
    ax4=fig.add_subplot(224,aspect='auto')
    ax4.plot(np.arange(len(hm3))*734,hm3,'o-')
    ax4.text(0.8, 0.1, '(B)', transform=ax4.transAxes,size=15, weight='bold')
    ax4.set_xlabel('Radial height (km)')
    ax4.set_ylabel('B (Gauss)')
    ax4.semilogx()
    ax4.axhline(y=159,linestyle='--',color='k',label='B = 159 G')
    ax4.axvline(x=2700,linestyle='-',color='k',label='r = 2700 km')
    ax4.legend()
    plt.show()

def hmi_euv_map_paper(hmi,ccmap,ccmap1,xc,yc,lev_1,lev_2,xl,xr,yl,yr):
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(hmi,origin=True,extent=[xl,xr,yl,yr],cmap='gray',label='A',interpolation='none',vmin=-400,vmax=400)
    #plt.contourf(ccmap1/np.max(ccmap1),extent=[xl,xr,yl,yr],levels=lev_2,cmap='Oranges',alpha=0.1)
    plt.contour(ccmap/np.max(ccmap),extent=[xl,xr,yl,yr],linewidth=14,levels=lev_1,colors='blue')
    ccmap[0:40,30:80]=0
    plt.contour(ccmap/np.max(ccmap),extent=[xl,xr,yl,yr],linewidth=14,levels=lev_1,colors='red')
    plt.text(xl+5, yl+5, 'Southern Ribbon', color='blue')#,bbox=dict(facecolor='white', alpha=0.1))
    plt.text(xr-35, yr-35, 'Northern Ribbon', color='red')#,bbox=dict(facecolor='white', alpha=0.1))
    #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    plt.colorbar(im,label='B (Gauss)')
    plt.legend()
    ax.set_xlabel('Solar X (arcsec)')
    ax.set_ylabel('Solar Y (arcsec)')
    ax.annotate('N',xy=(xl+2,yl+58),color='white')
    ax.annotate('W',xy=(xl+7,yl+52),color='white')
    ax.arrow(xl+2,yl+52,4.5,0,head_width=1.0,fc="white", ec="white", head_length=1)
    ax.arrow(xl+2,yl+52,0,4.5,head_width=1.0,fc="white", ec="white", head_length=1)
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    #plt.title(t)
    ax.grid(True)
    plt.show()

def euv_map_paper(ccmap,ccmap1,rhmap_l,rhmap_h,xc,yc,lev_1,lev_2,lev,xl,xr,yl,yr):
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(ccmap/np.max(ccmap),origin=True,extent=[xl,xr,yl,yr],cmap='sdoaia94',interpolation='none',vmin=0.002,vmax=0.5)
    #plt.contourf(ccmap1/np.max(ccmap1),extent=[xl,xr,yl,yr],levels=lev_2,cmap='Oranges',alpha=0.1)
    #plt.contour(ccmap1/np.max(ccmap1),extent=[xl,xr,yl,yr],linewidths=8,levels=lev_2,colors='k')
    #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
    #plt.colorbar(ss,label='Frequency (MHz)')
    ax.contour(rhmap_l,extent=[xl,xr,yl,yr],levels=lev,linewidths=3,colors='yellow')
    ax.contour(rhmap_h,extent=[xl,xr,yl,yr],levels=lev,linewidths=3,colors='magenta')
    ax.text(xl+5, yr-30, '6-10 keV', color='magenta',bbox=dict(facecolor='white', alpha=0.))
    ax.text(xl+5, yr-35, '10-18 keV', color='yellow',bbox=dict(facecolor='white', alpha=0.))
    plt.legend()
    ax.set_xlabel('Solar X (arcsec)')
    ax.set_ylabel('Solar Y (arcsec)')
    ax.annotate('N',xy=(xl+2,yl+58),color='white')
    ax.annotate('W',xy=(xl+7,yl+52),color='white')
    ax.arrow(xl+2,yl+52,4.5,0,head_width=1.0,fc="white", ec="white", head_length=1)
    ax.arrow(xl+2,yl+52,0,4.5,head_width=1.0,fc="white", ec="white", head_length=1)
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
    ax.set_xlabel('Solar X (arcsec)')
    ax.set_ylabel('Solar Y (arcsec)')
    ax.annotate('N',xy=(xl+2,yl+58),color='white')
    ax.annotate('W',xy=(xl+7,yl+52),color='white')
    ax.arrow(xl+2,yl+52,4.5,0,head_width=1.0,fc="white", ec="white", head_length=1)
    ax.arrow(xl+2,yl+52,0,4.5,head_width=1.0,fc="white", ec="white", head_length=1)
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    #plt.title(t)
    ax.grid(True)
    plt.show()

def centroid_map(hmi1,ccmap,xc,yc,lev_1,xl,xr,yl,yr):
    a,b,c,d,e,f=68,77,79,103,125,159
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(111,aspect='auto')
    im=ax.imshow(hmi1,origin=True,extent=[xl,xr,yl,yr],cmap='gray',interpolation='none',vmin=-400,vmax=400)
    #plt.contour(euv94/np.max(euv94),extent=[xl,xr,yl,yr],levels=np.linspace(0.4,0.8,2),colors='blue')
    #plt.contour(euv131/np.max(euv131),extent=[xl,xr,yl,yr],levels=lev_1,colors='green')
    #plt.contour(euv171/np.max(euv171),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.contour(euv304/np.max(euv304),extent=[xl,xr,yl,yr],levels=lev_1,colors='red')
    #plt.text(xl+3, yl+5, '94 $\AA$', color='blue',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+10, '131 $\AA$', color='green',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+15, '171 $\AA$', color='red',bbox=dict(facecolor='white', alpha=0.1))
    #plt.text(xl+3, yl+15, '304 $\AA$', color='red',bbox=dict(facecolor='white', alpha=0.1))
    plt.text(xr-11,yl+5,'(A)',color='k',backgroundcolor='white')
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
    #plt.errorbar(xc[f].mean(),yc[f].mean(),xerr=xc[f].std(),yerr=yc[f].std(),elinewidth=5,c='r',label='F')
    plt.legend(loc=2)
    left, bottom, width, height = [0.55, 0.57, 0.22, 0.22]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.imshow(hmi1,origin=True,extent=[xl,xr,yl,yr],cmap='gray',interpolation='none',vmin=-400,vmax=400)
    ax2.contourf(ccmap/np.max(ccmap),extent=[xl,xr,yl,yr],alpha=0.4,cmap='YlGn')
    ax2.set_ylim(352,362)
    ax2.set_xlim(474,484)
    ax2.errorbar(xc[a].mean(),yc[a].mean(),xerr=xc[a].std(),yerr=yc[a].std(),elinewidth=5,c='b',label='A')
    ax2.errorbar(xc[b].mean(),yc[b].mean(),xerr=xc[b].std(),yerr=yc[b].std(),elinewidth=5,c='c',label='B')
    ax2.errorbar(xc[c].mean(),yc[c].mean(),xerr=xc[c].std(),yerr=yc[c].std(),elinewidth=5,c='y',label='C')
    ax2.errorbar(xc[d].mean(),yc[d].mean(),xerr=xc[d].std(),yerr=yc[d].std(),elinewidth=5,c='m',label='D')
    ax2.errorbar(xc[e].mean(),yc[e].mean(),xerr=xc[e].std(),yerr=yc[e].std(),elinewidth=5,c='g',label='E')
    #ax2.errorbar(xc[f].mean(),yc[f].mean(),xerr=xc[f].std(),yerr=yc[f].std(),elinewidth=5,c='r',label='F')
    ax.set_xlabel('Solar X (arcsec)')
    ax.set_ylabel('Solar Y (arcsec)')
    ax.set_ylim(yl,yr)
    ax.set_xlim(xl,xr)
    #plt.title(t)
    ax.grid(True)
    plt.show()



def plot_3d_B(d,bxc,byc,im,rd_x,rd_y):
    #d=d.swapaxes(0,2)
    #bmin,bmax=130,880
    #d[np.where((d<bmin))]=0 
    #d[np.where((d>bmax))]=0
    fig = mlab.figure()
    xs,ys,zs=d.shape[0],d.shape[1],d.shape[2]
    xc_min,xc_max=bxc-1.01679*xs/2,bxc+1.01679*xs/2
    yc_min,yc_max=byc-1.01679*ys/2,byc+1.01679*ys/2
    X, Y, Z = np.mgrid[yc_min:yc_max:ys*1j, xc_min:xc_max:zs*1j,0:50:xs*1j]
    mlab.contour3d(X,Y,Z,d.swapaxes(0,2).swapaxes(0,1), colormap="YlGnBu",contours=90,vmin=1,vmax=200,opacity=0.5)
    mlab.colorbar(title='B (Gauss)')
    mlab.outline(color=(0, 0, 0))
    axes = mlab.axes(color=(0, 0, 0), nb_labels=5)
    Xc,Yc=np.mgrid[yc_min:yc_max:ys*1j, xc_min:xc_max:xs*1j]
    im_=im/im.max()
    im_[np.where(im_<0.8)]=np.nan
    obj=mlab.contour_surf(x, y, im_,extent=[yc_min,yc_max,xc_min,xc_max,5,6],contours=5,color=(1,0,0))
    mlab.text3d(460, 320, 27, 'Southern Ribbon', scale=(3, 3, 3),color=(1,1,0))
    mlab.text3d(480, 350, 27, 'Northern Ribbon', scale=(3, 3, 3),color=(1,1,0))
    #obj=mlab.imshow(im,extent=[rd_x[0],rd_x[-1],rd_y[0],rd_y[-1],0,0.1],opacity=0.5)
    #obj=mlab.imshow(im,opacity=0.6,colormap='YlOrRd',vmin=5.e6,vmax=9.e6)
    #obj.actor.orientation = [0, 0, 0] # the required orientation (deg)
    #obj.actor.position = [rd_x[int(rd_x.shape[0]/2)], rd_y[int(rd_y.shape[0]/2)], 0] # the required position 
    #obj.actor.scale = [2.5, 2.5, 2.5] # the required scale
    mlab.show()


if __name__=='__main__':
    main();
else:
    print('plotting module for 20120225....')


