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

def makeGaussian(size, fwhmx, fwhmy , center=None):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    return np.exp(-1*((x-x0)**2)/fwhmx**2) * np.exp(-1*((y-y0)**2) / fwhmy**2)

def flux2Tb(flux,bmaj,bmin,nu):
    '''
    Input: flux (SFU/beam), beam size (bmaj in arcsec), bmin (arcsec), frequency (GHz)
    Output: Brightness temperature (K)
    cmd: flux2Tb(flux,bmaj,bmin,nu)
    Ref: https://science.nrao.edu/facilities/vla/proposing/TBconv
    '''
    Tb=1222*(flux*1.e7)/(nu*nu*bmaj*bmin)
    return Tb


def Tb2flux(Tb,bmaj,bmin,nu):
    '''
    Input:  Brightness temperature (K), beam size (bmaj in arcsec), bmin (arcsec), frequency (GHz)
    Output: Flux (Jy/beam)
    cmd: flux2Tb(flux,bmaj,bmin,nu)
    Ref: https://science.nrao.edu/facilities/vla/proposing/TBconv
    '''
    f=Tb*(nu*nu*bmaj*bmin)/1222
    f=f/1.e7 # mJy ---> SFU
    return f

