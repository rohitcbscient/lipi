import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

aa=fits.open('sky_eval.fits')
data_sky=aa[0].data
header_sky=aa[0].header

bb=fits.open('cont_eval.fits')
data_cont=bb[0].data
header_cont=bb[0].header
