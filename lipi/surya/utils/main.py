import numpy as np
from scipy.io import readsav
from astropy.io import fits

def hms2sec_c(time):
    sec=int(time.split(' ')[1].split(':')[0])*3600+int(time.split(' ')[1].split(':')[1])*60+float(time.split(' ')[1].split(':')[2])
    return sec

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]

def read_fits(f):
        data_,h = fits.getdata(str(f),0,header=True)
        return h,data_

def radec_array(header,n,data):
        ra0=header['CRVAL1']
        dec0=header['CRVAL2']
        ra=data*0
        dec=data*0
        dra=np.linspace(-n/2,n/2,n)*header['CDELT1']
        ddec=np.linspace(-n/2,n/2,n)*header['CDELT2']
        ra=ra0+dra
        dec=dec0+ddec
        ra_dec=np.meshgrid(ra,dec)
        return ra_dec



