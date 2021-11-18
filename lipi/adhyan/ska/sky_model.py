from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np

aaM=fits.open('/home/i4ds1807205/skymodels/lambda_mollweide_haslam408_nofilt.fits')
hM,dM=aaM[1].header,aaM[1].data;wcsM=WCS(hM)
aaH=fits.open('/home/i4ds1807205/skymodels/lambda_haslam408_nofilt.fits')
hH,dH=aaH[1].header,aaH[1].data;wcsH=WCS(hH)
aaZ=fits.open('/home/i4ds1807205/skymodels/lambda_zea_haslam408_nofilt.fits')
hZ,dZ=aaZ[1].header,aaZ[1].data;wcsZ=WCS(hZ)


plt.figure()
ax0=plt.subplot(211,projection=wcsZ)
ax1=plt.subplot(212,projection=wcsM)
im0=ax0.imshow(np.log10(dZ),origin='lower');ax0.set_title('ZEA Projection')
im1=ax1.imshow(np.log10(dM),origin='lower');ax1.set_title('Mollweide Projection')
plt.colorbar(im0,ax=ax0,label='Log10(Temperature (MK))')
plt.colorbar(im1,ax=ax1,label='Log10(Temperature (MK))')
plt.show()

