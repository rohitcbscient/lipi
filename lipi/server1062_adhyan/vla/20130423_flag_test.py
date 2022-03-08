import glob
import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt

path='/nas08-data02/vladata/20130423/L-Band/rohit_analysis/flag_tests/'
imgL=sorted(glob.glob(path+'images_LL/*.FITS'),key=os.path.getmtime)
imgR=sorted(glob.glob(path+'images_RR/*.FITS'),key=os.path.getmtime)
imgLf=sorted(glob.glob(path+'images_LL_flag/*.FITS'),key=os.path.getmtime)
imgRf=sorted(glob.glob(path+'images_RR_flag/*.FITS'),key=os.path.getmtime)

TbmaxL=np.zeros(len(imgL));TbmaxR=np.zeros(len(imgR));TbmaxLf=np.zeros(len(imgLf));TbmaxRf=np.zeros(len(imgRf))

for i in range(128):
    fl=fits.open(imgL[i]);dl=fl[0].data
    fr=fits.open(imgR[i]);dr=fr[0].data
    flf=fits.open(imgLf[i]);dlf=flf[0].data
    frf=fits.open(imgRf[i]);drf=frf[0].data
    TbmaxL[i]=np.nanmax(dl);TbmaxR[i]=np.nanmax(dr)
    TbmaxLf[i]=np.nanmax(dlf);TbmaxRf[i]=np.nanmax(drf)


freq=np.linspace(994,2006,128)
plt.plot(freq,TbmaxR/1.e6,'o-',color='red',markersize=5,label='RR')
plt.plot(freq,TbmaxRf/1.e6,'o--',color='red',markersize=5,label='RR (flag)')
plt.plot(freq,TbmaxL/1.e6,'o-',color='blue',markersize=5,label='LL')
plt.plot(freq,TbmaxLf/1.e6,'o--',color='blue',markersize=5,label='LL (flag)')
plt.legend()
plt.xlabel('Frequency (MHz)')
plt.ylabel('T$_B$ (MK)')
plt.show()

