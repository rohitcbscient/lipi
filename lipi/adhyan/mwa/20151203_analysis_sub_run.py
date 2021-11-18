import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
import pickle

filelist=sorted(glob.glob('*.FITS'))
data_max=[0]*len(filelist)
data_min=[0]*len(filelist)
data_mean=[0]*len(filelist)

i=0
for f in filelist:
    h=fits.open(f)
    data=h[0].data[0][0]
    data_max[i]=data.max()
    data_min[i]=data.min()
    data_mean[i]=data.mean()
    i=i+1

plt.plot(data_mean,'o-')
plt.xlabel('Time (half sec)')
plt.ylabel('Amp')
plt.ylim([-3e-4,7e-4])
plt.show()


