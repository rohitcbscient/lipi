import numpy as np
from astropy.io import fits
from surya.utils import main as ut
import pickle
from scipy import interpolate
from astropy import wcs
from scipy import ndimage

#def get_Tbmaps(f):

def solar_center_pixel(fitsfile,file_):
    '''
    Input:
    fitsfile: MWA Radio FITS file
    file_: NASA Horizon time file
    Output:
    X-coordinate
    Y-coordinate
    time: Time of the FITS file header
    '''
    hdulist = fits.open(fitsfile)
    time=hdulist[0].header['DATE-OBS'].split('T')[1]
    ra,dec=interpolate_ra_dec(file_,time)
    w = wcs.WCS(hdulist[0].header)
    pix = w.wcs_world2pix(np.array([[ra,dec,0,0]]), 1)
    return int(np.round(pix[0][0])),int(np.round(pix[0][1])),time

def rotateimage(data,angle,xc,yc):
    '''
    Input: 
    data,angle,xc,yc
    data: 2-D input array
    angle: rotation angle
            Negative angle imples anti-clockwise rotation and vice versa
    xc: X-coordinate for rotation axis
    yc: Y-coordinate for rotation axis
    '''
    padX=[data.shape[1]-xc,xc]
    padY=[data.shape[0]-yc,yc]
    imgP=np.pad(data,[padY,padX],'constant')
    imgR=ndimage.rotate(imgP,angle*-1,reshape=False)
    return imgR[padY[0]:-padY[1],padX[0]:-padX[1]]


def time_array(d):
    '''
    Input: Full NASA Horizon array which include time in HHMMSS format
    Output: Array of time in sec
    '''
    t_=np.zeros(d[:,0].shape[0])
    for i in range(d[:,0].shape[0]):
        t_[i]=float(list(str(int(d[:,0][i])))[-2]+list(str(int(d[:,0][i])))[-1])+float(list(str(int(d[:,0][i])))[-4]+list(str(int(d[:,0][i])))[-3])*60+float(list(str(int(d[:,0][i])))[-5])*3600
    return t_

def interpolate_ra(file_,time_):
    '''
    Input: 
    NASA Horizon file path
    Time for which you seek RA
    Output: 
    RA for the input time
    '''
    d=np.loadtxt(file_)
    tarray=time_array(d)
    ra=ut.ra_hms2deg(d[:,1],d[:,2],d[:,3])
    fra = interpolate.interp1d(tarray,ra)
    return fra(time_)


def interpolate_dec(file_,time_):
    '''
    Input: 
    NASA Horizon file path
    Time for which you seek DEC
    Output: 
    DEC for the input time
    '''
    d=np.loadtxt(file_)
    tarray=time_array(d)
    dec=ut.dec_dms2deg(d[:,4],d[:,5],d[:,6])
    fdec = interpolate.interp1d(tarray,dec)
    return fdec(time_)

def interpolate_ra_dec(file_,time):
    '''
    Input: 
    NASA Horizon file path
    Time for which you seek RA-DEC
    Output: 
    RA-DEC for the input time
    '''
    nt=float(time.split(':')[0])*3600+float(time.split(':')[1])*60+float(time.split(':')[2])
    dec=interpolate_dec(file_,nt)
    ra=interpolate_ra(file_,nt)
    return ra,dec

def mean_flux(file_,f,baseline_filelist,res):
    '''
    Input:
    file_,f,baseline_filelist,res
    NASA Horizon file path
    MWA frequency band label
    MWA Baseline List
    Time Resolution 
    Output: 
    flux array in MWA format of DS
    time string array
    time second array
    '''
    bb=0
    flux=[0]*len(baseline_filelist)
    std_flux=[0]*len(baseline_filelist)
    for b in baseline_filelist:
        aa=pickle.load(open(str(file_)+str(f)+'_T'+str(b)+'.p','r'))
        aa[17][3][0][np.isnan(aa[17][3][0])]=0
        flux[bb]=aa[17][3][0]
        bb=bb+1
    time=[0]*flux[0].shape[1]
    timesec=[0]*flux[0].shape[1]
    for i in range(flux[0].shape[1]):
        t=aa[14].split(' ')[1]
        time[i]=ut.sec2hms_c(t,res,i)
        timesec[i]=ut.hms2sec_c(' '+time[i])
    flux=np.array(flux)
    return flux,time,timesec



def compute_Tb(f,xc,yc,del_,angle,res,freq,n,S_sun_t_):
    '''
    INPUT:
    (f,xc,yc,del_,angle,res,freq,n,S_sun_t_)
    f:FITS file, xc: central position of x, yc: central position of y
    del_:number of solar pixels or size of the submap
    angle: Position angle of the Sun
    res: pixel size of the image
    freq: Frequency
    n: Cutoff sigma
    S_sun_t_: Flux from the dynamic spectrum
    OUTPUT:
    (Tb,flux_,mean,std,bmaj,bmin,bpa,ndata)
    Tb: Brightness temperature map
    flux_: Flux density map
    mean: mean of noise
    std: STD of the noise
    bmaj: Major axis (in arcmin)
    bmin: Minor axis (in arcmin)
    bpa: Beam angle (in degree)
    ndata: Jy/beam just rotated
    '''
    fit=fits.open(f)
    bmin,bmaj,bpa=fit[0].header['BMIN'],fit[0].header['BMAJ'],fit[0].header['BPA']
    data_=fit[0].data[0][0]
    ra_,dec_=ut.radec_array(fit[0].header,data_.shape[0],data_)
    ra=ra_[yc-del_:yc+del_,xc-del_:xc+del_]
    dec=dec_[yc-del_:yc+del_,xc-del_:xc+del_]
    data=data_[yc-del_:yc+del_,xc-del_:xc+del_]
    data_rot=rotateimage(data,angle,del_-1,del_-1)
    noise=data_[0:xc-del_*5,:]
    mean=np.mean(noise)
    std=np.std(noise)
    max_=np.max(data_rot)
    lev=(mean+n*std)/max_
    bimage=ut.get_bimage(data_rot,lev)
    omega=np.pi*res*res*(np.pi/(3600*180.0))**2
    omega_e=np.pi*bmaj*bmin*(np.pi/180.)**2
    beam_pix=omega_e/omega
    img_pix=np.sum(bimage)
    ndata=data_rot*bimage
    flux_jy=ndata*(img_pix/beam_pix) # flux(Jy/beam)* number of beam
    flux_=(flux_jy)*S_sun_t_/np.sum(flux_jy) # Scaling fluxes
    Tb=((flux_*1.e-22)*(3.e8)**2)/(2*(1.38e-23)*freq*freq*omega) # In K
    return Tb,flux_,mean,std,bmaj,bmin,bpa,ndata

def get_residuals(f,xc,yc,del_,angle):
    '''
    Input:
    Residual fitsfile, centre x-coordinates, centre y-coordinates, half-dimensions of the map
    Output:
    residual map at Sun, total residual maps
    '''
    fit=fits.open(f)
    bmin,bmaj,bpa=fit[0].header['BMIN'],fit[0].header['BMAJ'],fit[0].header['BPA']
    data_=fit[0].data[0][0]
    ra_,dec_=ut.radec_array(fit[0].header,data_.shape[0],data_)
    ra=ra_[yc-del_:yc+del_,xc-del_:xc+del_]
    dec=dec_[yc-del_:yc+del_,xc-del_:xc+del_]
    data=data_[yc-del_:yc+del_,xc-del_:xc+del_]
    data_rot=rotateimage(data,angle,del_-1,del_-1)
    return data_rot,data_

def scale_residuals(Tb,flux,ndata,res):
    '''
    Scale the residuals to Brightness temp and flux
    Inputs:
    Tb,flux,ndata,res
    Outputs:
    residual in K, SFU, Tb factor and flux factor
    '''
    Tb_fac=np.sum(Tb)/np.sum(ndata)
    flux_fac=np.sum(flux)/np.sum(ndata)
    flux_res=res*flux_fac
    Tb_res=res*Tb_fac
    return Tb_res,flux_res,Tb_fac,flux_fac


