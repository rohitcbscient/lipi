import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

d7=fits.open('2050.1s.cal.spw.7.time.20:47:09~20:47:10.FITS')
d6=fits.open('2050.1s.cal.spw.6.time.20:47:09~20:47:10.FITS')
d5=fits.open('2050.1s.cal.spw.5.time.20:47:09~20:47:10.FITS')

d=np.concatenate((d5[0].data[0],d6[0].data[0],d7[0].data[0]))

rms_=np.std(d[:,0:500,:],axis=(1,2))
max_=np.max(d,axis=(1,2))

