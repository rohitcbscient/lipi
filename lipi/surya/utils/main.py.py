import numpy as np
from scipy.io import readsav
from astropy.io import fits
import cv2

def rad2asec(rad):
    '''
    INPUT
    In radians
    OUTPUT
    arcsec
    '''
    return rad*(180.0*3600/(np.pi))

def hms2sec_c(time):
    sec=int(time.split(' ')[1].split(':')[0])*3600+int(time.split(' ')[1].split(':')[1])*60+float(time.split(' ')[1].split(':')[2])
    return sec

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]

def find_predecessor(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if(array[idx]-value>0):
        idx=idx-1
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
    Output: Flux (SFU/beam)
    cmd: flux2Tb(flux,bmaj,bmin,nu)
    Ref: https://science.nrao.edu/facilities/vla/proposing/TBconv
    '''
    f=Tb*(nu*nu*bmaj*bmin)/1222
    f=f/1.e7 # mJy ---> SFU
    return f


def fitEllipse(data): # Data = binarized array
    #Akshay#
    '''
    Sources: <http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1167130>
    '''
    m = cv2.moments(data, 2)
    xc = m['m10'] / m['m00']    # x-coordinate of centre of ellipse        
    yc = m['m01'] / m['m00']    # y-coordinate of centre of ellipse
    a = (m['m20'] / m['m00']) - (xc**2)
    b = 2 * ((m['m11'] / m['m00']) - (xc * yc))
    c = (m['m02'] / m['m00']) - (yc**2)
    theta = .5 * (np.arctan2(b, (a - c)))        # Is in radians
    w = np.sqrt(6 * (a + c - np.sqrt(b**2 + (a-c)**2)))    # Length of minor axis of ellipse
    l = np.sqrt(6 * (a + c + np.sqrt(b**2 + (a-c)**2))) # Length of major axis of ellipse
    if (theta < 0):
        angle = theta*180/np.pi + 180        # Angle made by major axis of ellipse with x-axis (now converted to degrees)
    else:
            angle = theta*180/np.pi
    return xc,yc,w,l,angle

def get_bimage(data,lev):
        ndata=data/np.min(abs(data[data!=0]))
        ndata1=cv2.cvtColor(ndata.astype(np.float32),cv2.COLOR_GRAY2BGR)
        ndata2=cv2.cvtColor(ndata1, cv2.COLOR_BGR2GRAY)
        ret,bimage = cv2.threshold(ndata2,np.max(ndata2)*lev,np.max(ndata2),cv2.THRESH_BINARY)
        bimage=bimage/np.max(bimage)
        return bimage

def fit_1d(freq,Tb):
	x=np.log10(np.array(freq)*1.e6)
	y=np.log10(Tb)
	fit = np.polyfit(x, y, 1,cov=True)
	return fit[0][0],fit[0][1],np.sqrt(np.diag(fit[1]))[0],np.sqrt(np.diag(fit[1]))[1]


