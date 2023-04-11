# Script to write the cst files
import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits
from astropy.wcs import WCS

gauss_fits=fits.open('/home/rohit/karabo/karabo-pipeline/karabo/test/cstfile_test/gauss.fits')
aperture_fits=fits.open('/home/rohit/karabo/karabo-pipeline/karabo/test/cstfile_test/beam_S0000_TIME_SEP_CHAN_SEP_AMP_XX.fits')
gauss=gauss_fits[0].data[0][0];aperture=aperture_fits[0].data[0][0]
gauss_head=gauss_fits[0].header;aperture_head=aperture_fits[0].header
gx,gy=np.where(gauss==np.nanmax(gauss))[0];apx,apy=np.where(aperture==np.nanmax(aperture))
dist=np.sqrt((gx-apx)**2 + (gy-apy)**2)
wg = WCS(gauss_head);ag=WCS(aperture_head)
wg.wcs_pix2world(gx, gy, 1)

f,ax=plt.subplots(2,1);ax0=ax[0];ax1=ax[1]
ax0.imshow(gauss,origin='lower')
ax1.imshow(aperture,origin='lower')
plt.show()


