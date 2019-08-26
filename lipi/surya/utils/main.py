import numpy as np
from scipy.io import readsav
from astropy.io import fits
import cv2
import abel


def cart2polar(data):
    '''
    Convert cartisean to polar array
    Input: Cart Array
    Output: Polar Array, r- array, theta- array
    '''
    p,r,th=abel.tools.polar.reproject_image_into_polar(data)
    return p,r,th


def rad2asec(rad):
    '''
    INPUT
    In radians
    OUTPUT
    arcsec
    '''
    return rad*(180.0*3600/(np.pi))

def hms2sec_c(time):
    '''
    Input: time (HH:MM:SS)
    Output: time in sec
    '''
    sec=int(time.split(' ')[1].split(':')[0])*3600+int(time.split(' ')[1].split(':')[1])*60+float(time.split(' ')[1].split(':')[2])
    return sec

def sec2hms_c(time0,tres,i):
    '''
    Output time = Input time+ tres*i
    Input:
    time0,tres,i
    time0 in HH:MM:SS
    tres: time resolution
    i: Increment step
    Output:
    new time in HH:MM:SS
    '''
    hour=float(time0.split(':')[0])
    mini=float(time0.split(':')[1])
    sec=float(time0.split(':')[2])+tres*i
    isec=int(sec/60.)
    if(isec>0):
            sec=sec-60*isec
            mini=mini+isec
            if(mini>=60):
                    mini=mini-60
                    hour=hour+1
    time=str(int(hour))+':'+str(int(mini))+':'+str(sec)
    return time


def find_nearest(array, value):
    '''
    Input: 
    Array
    Input Value
    Output:
    index of closest value
    Array value of the index
    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]

def find_predecessor(array, value):
    '''
    Input:
    array
    value
    Output:
    index of closest predecessor
    Array value of closest predecessor
    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if(array[idx]-value>0):
        idx=idx-1
    return idx,array[idx]

def read_fits(f):
    '''
    Input: fits file path
    Output: header, data
    '''
    data_,h = fits.getdata(str(f),0,header=True)
    return h,data_

def radec_array(header,n,data):
    '''
    Input:
    Header, Array dimensions, data
    Output:
    RA-DEC array from simple grid formation
    Use it with caution
    '''
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


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def make2DGaussian(size, fwhmx, fwhmy , center=None):
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
    Input:
    Input Binary data
    Output:
    xc, yc, w, l, angle
    x-coordinate, y-coordinate, major axis, minor axis, angle made by major axis of ellipse with x-axis (now converted to degrees)
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
    '''
    Input:
    Input Data, level w.r.t max
    Output:
    Binary image
    '''
    ndata=data/np.min(abs(data[data!=0]))
    ndata1=cv2.cvtColor(ndata.astype(np.float32),cv2.COLOR_GRAY2BGR)
    ndata2=cv2.cvtColor(ndata1, cv2.COLOR_BGR2GRAY)
    ret,bimage = cv2.threshold(ndata2,np.max(ndata2)*lev,np.max(ndata2),cv2.THRESH_BINARY)
    bimage=bimage/np.max(bimage)
    return bimage

def fit_1d_Tb_freq(freq,Tb):
    '''
    Fit straight line to Tb vs freq
    Input:
    frequency, Tb array
    Output:
    slope fit
    intercept fit
    Uncertainity in slope fit
    Uncertainity in intercept fit
    '''
    x=np.log10(np.array(freq)*1.e6)
    y=np.log10(Tb)
    fit = np.polyfit(x, y, 1,cov=True)
    return fit[0][0],fit[0][1],np.sqrt(np.diag(fit[1]))[0],np.sqrt(np.diag(fit[1]))[1]

def fit_1d(freq,Tb):
    '''
    Fit straight line to y vs x
    Input:
    y, x
    Output:
    slope fit
    intercept fit
    Uncertainity in slope fit
    Uncertainity in intercept fit
    '''
    x=np.log10(np.array(freq))
    y=np.log10(Tb)
    fit = np.polyfit(x, y, 1,cov=True)
    return fit[0][0],fit[0][1],np.sqrt(np.diag(fit[1]))[0],np.sqrt(np.diag(fit[1]))[1]

def dec_dms2deg(dec_d,dec_m,dec_s):
    '''
    Input:
    DEC in degree, minutes and seconds format
    Output:
    DEC in degrees
    '''
    sign=dec_d/abs(dec_d)
    dec_=dec_d*1.0+sign*dec_m*(1/60.)+sign*dec_s*(1/3600.) # In degrees
    return dec_

def ra_hms2deg(ra_h,ra_m,ra_s):
    '''
    Input:
    RA in hours, minutes and seconds format
    Output:
    RA in degrees
    '''
    sign=ra_h/abs(ra_h)
    ra_=ra_h*15.0+sign*ra_m*(15/60.)+sign*ra_s*(15/3600.) # In degrees
    return ra_


