import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
from astropy.io import fits
from sunpy.map import Map
from astropy.coordinates import SkyCoord
import astropy.units as u
from surya.utils import main as ut
import matplotlib as mpl
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable


allmaps=pickle.load(open('/media/rohit/VLA/20160409/20160409_submap_50ms.p','rb'))
timevla=allmaps['vla']['timevla']
time94=allmaps['aia94']['time94'];time131=allmaps['aia131']['time131'];time335=allmaps['aia335']['time335']
time1600=allmaps['aia1600']['time1600'];time1700=allmaps['aia1700']['time1700'];time171=allmaps['aia171']['time171']
freq=np.linspace(0.997,1.245,32)


file_='/media/rohit/VLA/20160409/sun_L_20160409T1844-1846UT.50ms.cal.pol.LL.spw.0_16-31.time.18:45:10.-18:45:10.05.fits'
file_='/media/rohit/VLA/20160409/images_12s/image_184500_184600.FITS'
g=fits.open(file_)
g[0].header['CRVAL1']=-731.12;g[0].header['CRVAL2']=243.5
g.writeto('/media/rohit/VLA/20160409/sun_L_20160409T1844-1846UT.50ms.cal.pol.LL.spw.0_16-31.time.18:45:10.-18:45:10.05.fits')
map_=Map(file_)
map_.plot()
plt.show()

compmap=Map(allmaps['aia171']['map171'][0],allmaps['vla']['mapvla'][1400],composite=True)
compmap.set_levels(1,[30,40,50,60,70,80,90],percent=True)
compmap.plot()
#plt.xlim([400,900]);plt.ylim([50,550])
plt.show()


cc=allmaps['aia1600']['map1600'][0]
dd=allmaps['vla']['mapvla'][1400];dd.data[np.isnan(dd.data)]=0
p=cc.plot()
lev1=(1.5e7/dd.data.max())*np.array([40,50,60,70,80,90])*u.percent
xl = SkyCoord(dd.center.Tx.value-0.5*(2.0/0.6)*dd.data.shape[0],dd.center.Ty.value-0.5*(2.0/0.6)*dd.data.shape[1], frame=cc.coordinate_frame, unit=(u.arcsec, u.arcsec))
xlpix=cc.world_to_pixel(xl).x.value;ylpix=cc.world_to_pixel(xl).y.value
yr = SkyCoord(dd.center.Tx.value+0.5*(2.0/0.6)*dd.data.shape[0],dd.center.Ty.value+0.5*(2.0/0.6)*dd.data.shape[1], frame=cc.coordinate_frame, unit=(u.arcsec, u.arcsec))
xrpix=cc.world_to_pixel(yr).x.value;yrpix=cc.world_to_pixel(yr).y.value
dd.draw_contours(levels=lev1,colors='r',linewidths=3,extent=[xlpix,xrpix,ylpix,yrpix])
plt.text(401,250,str(np.round(dd.meta['crval3']/1.e9,4))+' GHz',color='r')
plt.text(401,50,'Contours: 0.6,0.75,0.9,1.0,1.2,1.3 MK',color='yellow')
plt.xlim([400,900]);plt.ylim([50,550])
plt.show()

