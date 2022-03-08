import numpy as np
import matplotlib.pyplot as plt
from surya.aia import main as aia
from surya.radio import get_maps
from astropy.io import fits

aiafile='/home/i4ds1807205/20141118/saia_00193_fd_20141118_223342.fts'
mwafile240='/home/i4ds1807205/20141118/240MHz.fits'
mwafile218='/home/i4ds1807205/20141118/218MHz.fits'
mwafile108='/home/i4ds1807205/20141118/108MHz.fits'

mwatime='06:55:20'
# RA and DEC 15 35 19.72 -19 17 18.8
mway,mwax=1182,462
pa=20.3406
mwa240=fits.open(mwafile240)
mwadata240=mwa240[0].data[0][0]
map240_=mwadata240[mwax-100:mwax+100,mway-100:mway+100]
map240=get_maps.rotateimage(map240_,pa*1,100,100)
aia=fits.open(aiafile)

aiamap=aia[0].data

mwa218=fits.open(mwafile218)
mwadata218=mwa218[0].data[0][0]
map218_=mwadata218[456-100:456+100,1184-100:1184+100]
map218=get_maps.rotateimage(map218_,pa*1,100,100)

mwa108=fits.open(mwafile108)
mwadata108=mwa108[0].data[0][0]
map108_=mwadata108[457-100:457+100,1142-100:1142+100]
map108=get_maps.rotateimage(map108_,pa*1,100,100)
aialim=(1024*2.4/2)/3600.
mwalim=(200*90./2)/3600.

plt.imshow(aiamap,origin='lower',extent=[-aialim,aialim,-aialim,aialim],vmin=0,vmax=1200,cmap='binary') # 1229=512*2.4
plt.contour(map108/np.max(map108),origin='lower',levels=[0.8,0.9],colors='red',linewidths=2,extent=[-mwalim,mwalim,-mwalim,mwalim]) # 40*100
plt.contour(map108/np.max(map108),origin='lower',levels=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1],colors='blue',linewidths=2,extent=[-mwalim,mwalim,-mwalim,mwalim]) # 40*100
plt.contour(map108/np.max(map108),origin='lower',levels=[0.001,0.0025,0.005],colors='green',linewidths=2,extent=[-mwalim,mwalim,-mwalim,mwalim]) # 40*100
plt.contour(map108/np.max(map108),origin='lower',levels=[0.00005],colors='gray',linewidths=2,extent=[-mwalim,mwalim,-mwalim,mwalim]) # 40*100
plt.xlabel('X (degrees)')
plt.ylabel('Y (degrees)')
#plt.xlim([-1200,1200])
#plt.ylim([-1200,1200])
plt.show()

