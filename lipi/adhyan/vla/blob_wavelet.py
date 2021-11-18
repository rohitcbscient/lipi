import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import argrelmax
import pickle
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
from astropy.io import fits
import glob
import get_ellipse_gaussian as eg
from sunpy.map import Map
from astropy.coordinates import SkyCoord
import astropy.units as u
from surya.utils import main as ut


# 2D Gaussian model
def funce(xy,x0,y0,sigmax0,sigmay0,rot0,H0,x1,y1,sigmax1,sigmay1,rot1,H1):
    x, y = xy
    fxy0 = eg.make_gauss (xy, [x0,y0], [sigmax0, sigmay0], np.deg2rad(rot0))*H0
    fxy1 = eg.make_gauss (xy, [x1,y1], [sigmax1, sigmay1], np.deg2rad(rot1))*H1
    I = fxy0 + fxy1
    return I
def func(xy,x0,y0,sigmax0,sigmay0,H0,x1,y1,sigmax1,sigmay1,H1):
    x, y = xy
    Ax0 = 1.0 / (2 * sigmax0**2);Ax1 = 1.0 / (2 * sigmax1**2);Ay0 = 1.0 / (2 * sigmay0**2);Ay1 = 1.0 / (2 * sigmay1**2)
    I = H0 * np.exp((-Ax0 * (x - x0)**2) -Ay0*((y - y0)**2)) + H1 * np.exp((-Ax1 * (x - x1)**2) -Ay1*((y - y1)**2))
    return I
def funcs(xy,x0,y0,sigmax0,sigmay0,H0,sigmax1,sigmay1,H1):
    x, y = xy
    Ax0 = 1.0 / (2 * sigmax0**2);Ax1 = 1.0 / (2 * sigmax1**2);Ay0 = 1.0 / (2 * sigmay0**2);Ay1 = 1.0 / (2 * sigmay1**2)
    I = H0 * np.exp((-Ax0 * (x - x0)**2) -Ay0*((y - y0)**2)) + H1 * np.exp((-Ax1 * (x - 32.8)**2) -Ay1*((y - 35)**2))
    return I

# Generate 2D gaussian
def generate(x0, y0, sigma, H):
    x = np.arange(0, max(x0, y0) * 2 + sigma, 1)
    y = np.arange(0, max(x0, y0) * 2 + sigma, 1)
    xx, yy = np.meshgrid(x, y)
    I = func((xx, yy), x0=x0, y0=y0, sigma=sigma, H=H)
    return xx, yy, I

def fit(image,inin,inis):
    #x0,y0,sigmax0,sigmay0,H0=inin;sigmax1,sigmay1,H1=inis
    x0,y0,sigmax0,sigmay0,rot0,H0=inin;x1,y1,sigmax1,sigmay1,rot1,H1=inis
    # Prepare fitting
    x = np.arange(0, image.shape[1], 1)
    y = np.arange(0, image.shape[0], 1)
    xx, yy = np.meshgrid(x, y)
    lower = [30, 20, 0, 0, 0, 1.e6, 20, 20, 0, 0, 0, 1.e6]
    #upper = [image.shape[0]-30, image.shape[1]-20, 11, max(*image.shape), 180, np.nanmax(image)*100,image.shape[0]-30, image.shape[1]-20, 11, max(*image.shape), 180, np.nanmax(image)*100]
    upper = [image.shape[0]-30, image.shape[1]-20, 11, max(*image.shape), 180, 1.e11,image.shape[0]-30, image.shape[1]-20, 11, max(*image.shape), 180, 1.e11]
    bounds = [lower, upper]
    initial_guess = [x0, y0, sigmax0,sigmay0, rot0, H0, x1, y1, sigmax1, sigmay1, rot1, H1]
    pred_params, uncert_cov = curve_fit(funce, (xx.ravel(), yy.ravel()), image.ravel() ,p0=initial_guess, bounds=bounds, maxfev=10000)
    nimg=funce((xx.ravel(), yy.ravel()),pred_params[0],pred_params[1],pred_params[2],pred_params[3],pred_params[4],pred_params[5],pred_params[6],pred_params[7],pred_params[8],pred_params[9],pred_params[10],pred_params[11])
    # Get residual
    predictions = funce((xx, yy), *pred_params)
    rms = np.sqrt(np.mean((image.ravel() - predictions.ravel())**2))
    print(initial_guess,xx.ravel().shape,yy.ravel().shape,bounds)
    print("Initial Guess : ", initial_guess)
    print("Predicted params : ", pred_params)
    print("Residual : ", rms)
    return nimg,pred_params,rms,uncert_cov

def plot(image,nimg,params,g1,g2):
    fig, (ax0,ax1) = plt.subplots(2,1)
    ax0.imshow(image, cmap=plt.cm.BrBG, interpolation='nearest', origin='lower')
    ax1.imshow(nimg, cmap=plt.cm.BrBG, interpolation='nearest', origin='lower')
    ax0.scatter(params[0], params[1], s=10, c="red", marker="x")
    ax0.scatter(params[6], params[7], s=10, c="red", marker="x")
    ax0.scatter(g1[0], g1[1], s=20, c="white", marker="o")
    ax0.scatter(g2[0], g2[1], s=20, c="white", marker="o")
    plt.show()

def produce_tstring(mapp):
    date=mapp.date.ymdhms
    hhmmss=' '+str(date[3])+':'+str(date[4])+':'+str(date[5])
    sec=ut.hms2sec_c(hhmmss)
    return sec


def get_evla_submap(f,xbl,ybl,xtr,ytr):
    print ('Reading...'+f[0])
    #xcen=(xbl+xtr)*0.5;ycen=(ybl+ytr)*0.5
    n=len(f);maplist=[0]*n;datalist=[0]*n;time=[0]*n
    for i in range(n):
        g=fits.open(f[i])
        #g[0].header['CRVAL1']=0;g[0].header['CRVAL2']=0
        #g[0].header['CRVAL1']=-770;g[0].header['CRVAL2']=220
        #g[0].header['CRVAL1']=-731.12;g[0].header['CRVAL2']=243.5
        #g[0].data[np.isnan(g[0].data)]=0
        map_=Map(g[0].data[0][0],g[0].header)
        #cent_pix=map_.world_to_pixel(SkyCoord(xcen*u.arcsec, ycen*u.arcsec, frame=map_.coordinate_frame)) 
        bl = SkyCoord(xbl*u.arcsec, ybl*u.arcsec, frame=map_.coordinate_frame)
        tr = SkyCoord(xtr*u.arcsec, ytr*u.arcsec, frame=map_.coordinate_frame)
        maplist[i]=map_
        datalist[i]=maplist[i].data
        time[i]=produce_tstring(maplist[i])
    time=np.array(time)
    return maplist,datalist,time

aa=fits.open(sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_0/*spw.0*.FITS'))[421])
m=90
inin=[125-m,142-m,4,4,30,7.e6];inis=[122.8-m,125-m,4,4,30,8.e7]
data=aa[0].data[0][0];data=data[m:-m,m:-m];data[np.isnan(data)]=0
nimg_,params,rms,uncert_cov = fit(data, inin, inis);nimg=nimg_.reshape((2*(128-m),2*(128-m)))
plot(data,nimg, params,inin,inis)


dump=1
if(dump):
    spc=['0-15','16-31','32-47','48-63'];x0=[0]*8;x1=[0]*8;y0=[0]*8;y1=[0]*8;sigx0=[0]*8;sigy0=[0]*8;sigx1=[0]*8;sigy1=[0]*8;rot0=[0]*8;rot1=[0]*8;tb1=[0]*8;tb0=[0]*8
    ex0=[0]*8;ex1=[0]*8;ey0=[0]*8;ey1=[0]*8;esigx0=[0]*8;esigy0=[0]*8;esigx1=[0]*8;esigy1=[0]*8;etb1=[0]*8;etb0=[0]*8
    for i in range(8):
        x0[i]=[0]*4;x1[i]=[0]*4;y0[i]=[0]*4;y1[i]=[0]*4;sigx0[i]=[0]*4;sigy0[i]=[0]*4;sigx1[i]=[0]*4;sigy1[i]=[0]*4;rot0[i]=[0]*4;rot1[i]=[0]*4;tb1[i]=[0]*4;tb0[i]=[0]*4
        for j in range(4):
            listvla_r=sorted(glob.glob('/media/rohit/VLA/20160409/images_50ms_RR/spw_'+str(i)+'/*spw.*'+spc[j]+'*.FITS'))[0:2000]
            x0[i][j]=[0]*2000;x1[i][j]=[0]*2000;y0[i][j]=[0]*2000;y1[i][j]=[0]*2000;sigx0[i][j]=[0]*2000;sigy0[i][j]=[0]*2000;sigx1[i][j]=[0]*2000;sigy1[i][j]=[0]*2000;rot0[i][j]=[0]*2000;rot1[i][j]=[0]*2000;tb1[i][j]=[0]*2000;tb0[i][j]=[0]*2000
            for k in range(2000):
                mr,dr,tr=get_evla_submap([listvla_r[k]],0,-1,0,-1)
                hr=mr[0];data=hr.data;data=data[m:-m,m:-m];data[np.isnan(data)]=0
                nimg_,params,rms,uncert_cov = fit(data, inin, inis);nimg=nimg_.reshape((2*(128-m),2*(128-m)))
                x0[i][j][k]=params[0];y0[i][j][k]=params[1];sigx0[i][j][k]=params[2];sigy0[i][j][k]=params[3];rot0[i][j][k]=params[4];tb0[i][j][k]=params[5];x1[i][j][k]=params[6];y1[i][j][k]=params[7];sigx1[i][j][k]=params[8];sigy1[i][j][k]=params[9];rot1[i][j][k]=params[10];tb1[i][j][k]=params[11]
    x0=np.array(x0).reshape(32,2000)
    y0=np.array(y0).reshape(32,2000)
    x1=np.array(x1).reshape(32,2000)
    y1=np.array(y1).reshape(32,2000)
    sigx0=np.array(sigx0).reshape(32,2000)
    sigy0=np.array(sigy0).reshape(32,2000)
    sigx1=np.array(sigx1).reshape(32,2000)
    sigy1=np.array(sigy1).reshape(32,2000)
    rot0=np.array(rot0).reshape(32,2000)
    rot1=np.array(rot1).reshape(32,2000)
    tb0=np.array(tb0).reshape(32,2000)
    tb1=np.array(tb1).reshape(32,2000)
    pickle.dump([[x0,y0,sigx0,sigy0,rot0,tb0],[x1,y1,sigx1,sigy1,rot1,tb1]],open('/media/rohit/VLA/20160409/blob_param.p','wb'))



# The fit performs well without bounds

