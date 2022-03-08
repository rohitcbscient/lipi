import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

conpath='/home/rohit/rascil_data/rascil_config/'
lowska=np.loadtxt(conpath+'LOW_SKA-TEL-SKO-0000422_Rev3.txt',delimiter=',')
midska=np.loadtxt(conpath+'MID_SKA-TEL-INSA-0000537_Rev05.txt')

mid_beam_real=fits.open('/home/rohit/rascil_data/rascil-master-data-models/MID_FEKO_VP_Ku_45_12179_real.fits');hreal=mid_beam_real[0].header
mid_beam_imag=fits.open('/home/rohit/rascil_data/rascil-master-data-models/MID_FEKO_VP_Ku_45_12179_imag.fits');himag=mid_beam_imag[0].header

low_ska_beam=fits.open('/home/rohit/rascil_data/rascil-master-data-models/SKA1_LOW_beam.fits');ska_low_head=low_ska_beam[0].header;ska_low_data=low_ska_beam[0].data

