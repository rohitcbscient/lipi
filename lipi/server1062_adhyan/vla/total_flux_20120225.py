import numpy as np
from astropy.io import fits
import pickle
import glob

filelist=sorted(glob.glob('*.FITS') ,key=lambda f: int(filter(str.isdigit, f)))
sumall=[0]*len(filelist)
freq=[0]*len(filelist)
for i in range(len(filelist)):
    d=fits.open(filelist[i])
    header=d[0].header
    data=d[0].data[0][0][123:149,106:142]
    sumall[i]=np.nanmax(data)
    freq[i]=header['crval3']

sumall=np.array(sumall)
freq=np.array(freq)
sumall=sumall.reshape((96,832))
freq=freq.reshape((96,832))
