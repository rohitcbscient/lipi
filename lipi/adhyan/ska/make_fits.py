from astropy.wcs import WCS
from astropy import wcs
from astropy.io import fits
import numpy as np

f=fits.open('/home/rohit/simulations/gleam/GLEAM_EGC_v2.fits')
data=f[1].data
ra=data['RAJ2000'];dec=data['DEJ2000']

#### Construct Fits file ######
w = wcs.WCS(naxis=2)
w.wcs.crpix = [0, 0]
w.wcs.cdelt = np.array([-0.066667, 0.066667])
w.wcs.crval = [0, -90]
w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
w.wcs.set_pv([(2, 1, 45.0)])
px, py = w.wcs_world2pix(ra, dec, 1)


