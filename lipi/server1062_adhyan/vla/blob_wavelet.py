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
import os
from multiprocessing import Process
import concurrent.futures as cf
import lmfit


# 2D Gaussian model
def funce(xy,x0,y0,sigmax0,sigmay0,rot0,x1,y1,sigmax1,sigmay1,rot1,w0,w1):
    x, y = xy
    fxy0 = eg.make_gauss (xy, [x0,y0], [sigmax0, sigmay0], np.deg2rad(rot0));fxy0=fxy0/np.sum(fxy0)
    fxy1 = eg.make_gauss (xy, [x1,y1], [sigmax1, sigmay1], np.deg2rad(rot1));fxy1=fxy1/np.sum(fxy1)
    I = (w0*fxy0 + w1*fxy1)#; I = I/np.sum(I)
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

def residue(param,data,eps):
    x = np.arange(0, data.shape[1], 1)
    y = np.arange(0, data.shape[0], 1)
    xy = np.meshgrid(x, y)
    model=funce(xy,param['x0'],param['y0'],param['sigmax0'],param['sigmay0'],param['rot0'],param['x1'],param['y1'],param['sigmax1'],param['sigmay1'],param['rot1'],param['w0'],param['w1'])
    return abs(model-data)/(data)

# Generate 2D gaussian
def generate(x0, y0, sigma, H):
    x = np.arange(0, max(x0, y0) * 2 + sigma, 1)
    y = np.arange(0, max(x0, y0) * 2 + sigma, 1)
    xx, yy = np.meshgrid(x, y)
    I = func((xx, yy), x0=x0, y0=y0, sigma=sigma, H=H)
    return xx, yy, I


def plot(image,nimg,params,g1,g2):
    fig, (ax0,ax1) = plt.subplots(2,1,sharex=True)
    im0=ax0.imshow(image, cmap=plt.cm.BrBG, interpolation='nearest', origin='lower',vmin=1.e6,vmax=1.e8)
    ax1.imshow(nimg, cmap=plt.cm.BrBG, interpolation='nearest', origin='lower',vmin=1.e6,vmax=1.e8)
    #im0=ax0.imshow(image, cmap=plt.cm.BrBG, interpolation='nearest', origin='lower',vmin=1.e-5,vmax=0.001)
    #ax1.imshow(nimg, cmap=plt.cm.BrBG, interpolation='nearest', origin='lower',vmin=1.e-5,vmax=0.001)
    ax0.scatter(params[0], params[1], s=10, c="red", marker="x")
    ax0.scatter(params[5], params[6], s=10, c="red", marker="x")
    ax0.scatter(g1[0], g1[1], s=20, c="white", marker="o")
    ax0.scatter(g2[0], g2[1], s=20, c="white", marker="o")
    ax1.scatter(params[0], params[1], s=10, c="red", marker="x")
    ax1.scatter(params[5], params[6], s=10, c="red", marker="x")
    ax1.scatter(g1[0], g1[1], s=20, c="white", marker="o")
    ax1.scatter(g2[0], g2[1], s=20, c="white", marker="o")
    fig.colorbar(im0)
    plt.show()



def get_evla_submap(f,xbl,ybl,xtr,ytr):
    #print ('Reading...'+f[0])
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
    return maplist,datalist

def fit(image,inin,inis):
    #x0,y0,sigmax0,sigmay0,H0=inin;sigmax1,sigmay1,H1=inis
    x0,y0,sigmax0,sigmay0,rot0=inin;x1,y1,sigmax1,sigmay1,rot1=inis
    # Prepare fitting
    x = np.arange(0, image.shape[1], 1)
    y = np.arange(0, image.shape[0], 1)
    xx, yy = np.meshgrid(x, y);xy=xx,yy
    #lower = [b, b, 1, 1, 0, 1.e6, b, b, 1, 1, 0, 1.e6]
    ##upper = [image.shape[0]-30, image.shape[1]-20, 11, max(*image.shape), 180, np.nanmax(image)*100,image.shape[0]-30, image.shape[1]-20, 11, max(*image.shape), 180, np.nanmax(image)*100]
    #upper = [image.shape[0]-b, image.shape[1]-b, 11, 11, 180, 1.e8,image.shape[0]-b, image.shape[1]-b, 11, 11, 180, 1.e8]
    #bounds = [lower, upper]
    #initial_guess = [x0, y0, sigmax0,sigmay0, rot0, H0, x1, y1, sigmax1, sigmay1, rot1, H1];eps=0.1*image
    #pred_params, uncert_cov = curve_fit(funce, (xx.ravel(), yy.ravel()), image.ravel() ,p0=initial_guess, bounds=bounds, method='dogbox',maxfev=50000)
    #nimg=funce((xx.ravel(), yy.ravel()),pred_params[0],pred_params[1],pred_params[2],pred_params[3],pred_params[4],pred_params[5],pred_params[6],pred_params[7],pred_params[8],pred_params[9],pred_params[10],pred_params[11])
    #LMFIT
    ylim=36;params=lmfit.Parameters();y0_,x0_=np.where(image==np.nanmax(image))
    if(y0_>ylim):
        x1=x0_;y1=y0_;y1min=y0_-1;y1max=y0_+1;x1min=x0_-1;x1max=x0_+1;x0=x0_;x0min=x0_-8;x0max=x0_+8;y0=y0_;y0min=y0_-16;y0max=y0_-8;w0=0.3;w0min=0.1;w0max=0.5;w1=0.7;w1min=0.5;w1max=0.9
    if(y0_<=ylim):
        x0,y0=x0_,y0_;y0min=y0_-1;y0max=y0_+1;x0min=x0_-1;x0max=x0_+1;x1=x0_;x1min=x0_-8;x1max=x0_+8;y1=y0_;y1min=y0_+8;y1max=y0_+16;w1=0.3;w1min=0.1;w1max=0.5;w0=0.7;w0min=0.5;w0max=0.9
    params.add(name='x0', value=x0, min=x0min,max=x0max);params.add(name='y0', value=y0, min=y0min,max=y0max);params.add(name='sigmax0', value=sigmax0, min=3,max=9);params.add(name='sigmay0', value=sigmay0, min=3,max=9);params.add(name='rot0', value=rot0, min=0,max=180)
    params.add(name='x1', value=x1, min=x1min,max=x1max);params.add(name='y1', value=y1, min=y1min,max=y1max);params.add(name='sigmax1', value=sigmax1, min=3,max=9);params.add(name='sigmay1', value=sigmay1, min=3,max=9);params.add(name='rot1', value=rot1, min=0,max=180);eps=1
    params.add(name='w0',value=w0, min=w0min,max=w0max);params.add(name='w1',value=w1, min=w1min,max=w1max)
    out=lmfit.minimize(residue, params, args=(image/np.nansum(image), eps), method='least_squares')
    # Get residual
    #predictions = funce((xx, yy), *pred_params)
    #rms = np.sum((image.ravel() - predictions.ravel())**2)/(image.flatten().shape[0]**2-12)
    #print(initial_guess,xx.ravel().shape,yy.ravel().shape,bounds)
    #print("Initial Guess : ", initial_guess)
    #print("Predicted params : ", pred_params)
    #print("Residual : ", rms)
    i=0;oname=[0]*len(out.params.items());ovalue=[0]*len(out.params.items());ostd=[0]*len(out.params.items())
    for name,param in out.params.items():
        oname[i]=name;ovalue[i]=param.value;ostd[i]=param.stderr
        i=i+1
    #return nimg,pred_params,rms,uncert_cov
    nimg=funce((xx, yy), ovalue[0],ovalue[1],ovalue[2],ovalue[3],ovalue[4],ovalue[5],ovalue[6],ovalue[7],ovalue[8],ovalue[9],ovalue[10],ovalue[11])*np.nansum(image)
    return nimg,oname,ovalue,ostd

check=0
if(check):
    aa=fits.open(sorted(glob.glob('/nas08-data02/rohit/20160409_analysis/images_50ms_RR/spw_0/*spw.0*.FITS'))[1375])
    #aa=fits.open(sorted(glob.glob('/nas08-data02/rohit/20160409_analysis/1844/fits_0_16-31/*spw.0*FITS'))[1195])
    mm=100
    inin=[125-mm,142-mm,4,4,30];inis=[122.8-mm,125-mm,4,4,30]
    data=aa[0].data[0][0];data=data[mm:-mm,mm:-mm];data[np.isnan(data)]=0;data1=data/np.sum(data)
    #nimg_,params,rms,uncert_cov = fit(data, inin, inis);nimg=nimg_.reshape((2*(128-mm),2*(128-mm)))
    res = fit(data, inin, inis);nimg=res[0];params=res[2]
    plot(data,nimg, params,inin,inis)
    sys.exit()


def dump_pickle(idx_):
    i,spc1=idx_[0],idx_[1]
    #i=0;spc1='0_0-15';data=[0]*2000;y=[0]*2000
    listvla_r=sorted(glob.glob('/nas08-data02/rohit/20160409_analysis/images_50ms_RR/spw_'+str(i)+'/*spw.*'+str(spc1)+'*.FITS'))[0:2000]#[850:870]
    for k in range(len(listvla_r)):
        #mr,dr=get_evla_submap([listvla_r[k]],0,-1,0,-1);data=dr[0];data[np.isnan(data)]=0;hr=mr[0];ymax,x=np.where(hr.data==np.nanmax(hr.data));y[k]=hr.reference_coordinate.Ty.value+(ymax-(hr.reference_pixel.y.value-1))*hr.scale.axis1.value
        kk="%04d" % k
        if(os.path.isfile('/nas08-data02/rohit/20160409_analysis/b_new/blob_param_spw_'+str(i)+'_'+str(spc1)+'_'+str(kk)+'.p')==False):
            try:
                mr,dr=get_evla_submap([listvla_r[k]],0,-1,0,-1)
                #mm=90;inin=[125-mm,142-mm,4,4,30,1.e7];inis=[122.8-mm,125-mm,4,4,30,8.e7]
                m=100;inin=[125-m,142-m,4,4,30];inis=[122.8-m,125-m,4,4,30]
                hr=mr[0];data=hr.data;data=data[m:-m,m:-m];data[np.isnan(data)]=0
                #nimg_,params,rms,uncert_cov = fit(data, inin, inis);nimg=nimg_.reshape((2*(128-m),2*(128-m)))
                print listvla_r[k]
                res = fit(data, inin, inis);nimg=res[0];params=res[2];err=res[3]
                x00=hr.reference_coordinate.Tx.value+((params[0]+m)-(hr.reference_pixel.x.value-1))*hr.scale.axis1.value;x10=hr.reference_coordinate.Tx.value+((params[5]+m)-(hr.reference_pixel.x.value-1))*hr.scale.axis1.value
                y00=hr.reference_coordinate.Ty.value+((params[1]+m)-(hr.reference_pixel.y.value-1))*hr.scale.axis1.value;y10=hr.reference_coordinate.Ty.value+((params[6]+m)-(hr.reference_pixel.y.value-1))*hr.scale.axis1.value # 0 --> N 1--> S
                x0=x00;y0=y00;sigx0=params[2]*hr.scale.axis1.value;sigy0=params[3]*hr.scale.axis1.value;rot0=params[4];x1=x10;y1=y10;sigx1=params[7]*hr.scale.axis1.value;sigy1=params[8]*hr.scale.axis1.value;rot1=params[9];tb0=nimg[int(params[1]),int(params[0])];tb1=nimg[int(params[6]),int(params[5])];w0=params[10];w1=params[11]
                pickle.dump([[x0,y0,sigx0,sigy0,rot0,tb0],[x1,y1,sigx1,sigy1,rot1,tb1],w0,w1,err,data,nimg,data-nimg],open('/nas08-data02/rohit/20160409_analysis/b_new/blob_param_spw_'+str(i)+'_'+str(spc1)+'_'+str(kk)+'.p','wb'),protocol=2)
            except:
                pass
    return x0


spc=['0-15','16-31','32-47','48-63'];n=0
spc=['48-63'];n=0
spw=[6]
#spw=[0,1,2,3,4,5,6,7]
#spw=[0,1];spc=['0-15','16-31'];n=0
allidx=[0]*len(spw)*len(spc)
for l in range(len(spw)):
    for m in range(len(spc)):
        allidx[n]=[spw[l],spc[m]]
        n=n+1
with cf.ProcessPoolExecutor(max_workers=1) as exe:
    #results=[exe.submit(dump_pickle,ii) for ii in allidx] # 1st way
    results = exe.map(dump_pickle,allidx) # 2nd way

## Threads: Good for I/O bounds
## Processes: Good for CPU bound process

sys.exit()
do_multiprocessing=0
if(do_multiprocessing):
    processes=[]
    num_processes=8
    for i in range(num_processes):
        p=Process(target=dump_pickle,args=allidx[i])
        processes.append(p)
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    print('####End###')

dump=0
if(dump):
    listqs=sorted(glob.glob('/nas08-data02/rohit/20160409_analysis/1844/fits_0_16-31/*spw.0*FITS'))
    qx0=[0]*len(listqs);qx1=[0]*len(listqs);qy0=[0]*len(listqs);qy1=[0]*len(listqs);qsigx0=[0]*len(listqs);qsigy0=[0]*len(listqs);qsigx1=[0]*len(listqs);qsigy1=[0]*len(listqs);qrot0=[0]*len(listqs);qrot1=[0]*len(listqs);qtb1=[0]*len(listqs);qtb0=[0]*len(listqs)
    for k in range(3981,len(listqs)):
        kk="%04d" % k
        if(os.path.isfile('/nas08-data02/rohit/20160409_analysis/qb/qblob_param_'+str(kk)+'.p')==False):
            m,d=get_evla_submap([listqs[k]],0,-1,0,-1)
            hl=m[0];data=hl.data;data=data[mm:-mm,mm:-mm];data[np.isnan(data)]=0
            nimg_,params,rms,uncert_cov = fit(data, inin, inis);nimg=nimg_.reshape((2*(128-mm),2*(128-mm)))
            x00=hl.reference_coordinate.Tx.value+((params[0]+90)-(hl.reference_pixel.x.value-1))*hl.scale.axis1.value;x10=hl.reference_coordinate.Tx.value+((params[6]+90)-(hl.reference_pixel.x.value-1))*hl.scale.axis1.value
            y00=hl.reference_coordinate.Ty.value+((params[1]+90)-(hl.reference_pixel.y.value-1))*hl.scale.axis1.value;y10=hl.reference_coordinate.Ty.value+((params[7]+90)-(hl.reference_pixel.y.value-1))*hl.scale.axis1.value # 0 --> N 1--> S
            qx0[k]=x00;qy0[k]=y00;qsigx0[k]=params[2]*hl.scale.axis1.value;qsigy0[k]=params[3]*hl.scale.axis1.value;qrot0[k]=params[4];qtb0[k]=params[5];qx1[k]=x10;qy1[k]=y10;qsigx1[k]=params[8]*hl.scale.axis1.value;qsigy1[k]=params[9]*hl.scale.axis1.value;qrot1[k]=params[10];qtb1[k]=params[11]
            pickle.dump([[qx0[k],qy0[k],qsigx0[k],qsigy0[k],qrot0[k],qtb0[k]],[qx1[k],qy1[k],qsigx1[k],qsigy1[k],qrot1[k],qtb1[k]]],open('/nas08-data02/rohit/20160409_analysis/qb/qblob_param_'+str(kk)+'.p','wb'))
    qx0=np.array(qx0)
    qy0=np.array(qy0)
    qx1=np.array(qx1)
    qy1=np.array(qy1)
    qsigx0=np.array(qsigx0)
    qsigy0=np.array(qsigy0)
    qsigx1=np.array(qsigx1)
    qsigy1=np.array(qsigy1)
    qrot0=np.array(qrot0)
    qrot1=np.array(qrot1)
    qtb0=np.array(qtb0)
    qtb1=np.array(qtb1)






# The fit performs well without bounds

