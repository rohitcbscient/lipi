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
def funce(xy,x0,y0,sigmax0,sigmay0,rot0):
    x, y = xy
    fxy0 = eg.make_gauss(xy, [x0,y0], [sigmax0, sigmay0], np.deg2rad(rot0))
    I = fxy0 ; I = I/np.nansum(I)
    return I

def residue(param,data,eps):
    x = np.arange(0, data.shape[1], 1)
    y = np.arange(0, data.shape[0], 1)
    xy = np.meshgrid(x, y)
    model=funce(xy,param['x0'],param['y0'],param['sigmax0'],param['sigmay0'],param['rot0'])
    return abs(data-model)/(data+1)


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
    return maplist,datalist

def fit(image,inin):
    #x0,y0,sigmax0,sigmay0,H0=inin;sigmax1,sigmay1,H1=inis
    im=image/np.sum(image)
    x0,y0,sigmax0,sigmay0,rot0=inin
    # Prepare fitting
    x = np.arange(0, image.shape[1], 1)
    y = np.arange(0, image.shape[0], 1)
    xx, yy = np.meshgrid(x, y);xy=xx,yy
    #LMFIT
    params=lmfit.Parameters()
    params.add(name='x0', value=55.0, min=50,max=60);params.add(name='y0', value=83.0, min=75,max=90);params.add(name='sigmax0', value=4, min=3,max=10);params.add(name='sigmay0', value=4, min=3,max=10);params.add(name='rot0', value=rot0, min=0,max=180);eps=0.1
    out=lmfit.minimize(residue, params, args=(im, eps), method='leastsq')
    # Get residual
    i=0;oname=[0]*len(out.params.items());ovalue=[0]*len(out.params.items());ostd=[0]*len(out.params.items())
    for name,param in out.params.items():
        oname[i]=name;ovalue[i]=param.value;ostd[i]=param.stderr
        i=i+1
    #return nimg,pred_params,rms,uncert_cov
    nimg=funce((xx, yy), ovalue[0],ovalue[1],ovalue[2],ovalue[3],ovalue[4])*np.nansum(image)
    return nimg,oname,ovalue,ostd

check=0
if(check):
    f=sorted(glob.glob('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/images_LL/2013*.FITS'))[860]
    mr,dr=get_evla_submap([f],0,-1,0,-1);mm=30
    hr=mr[0];data=hr.data;data=data[mm:-mm,mm:-mm];data[np.isnan(data)]=0;y0max,x0max=np.where(data==np.nanmax(data))
    inin=[x0max[0],y0max[0],4,4,30]
    #nimg_,params,rms,uncert_cov = fit(data, inin, inis);nimg=nimg_.reshape((2*(128-mm),2*(128-mm)))
    res = fit(data, inin);nimg=res[0];params=res[2]
    x00=hr.reference_coordinate.Tx.value+((params[0]+mm)-(hr.reference_pixel.x.value-1))*hr.scale.axis1.value
    y00=hr.reference_coordinate.Ty.value+((params[1]+mm)-(hr.reference_pixel.y.value-1))*hr.scale.axis1.value
    x0=x00;y0=y00;sigx0=params[2]*hr.scale.axis1.value;sigy0=params[3]*hr.scale.axis1.value;rot0=params[4];tb0=nimg[int(params[0]),int(params[1])]
    plot(data,nimg, params,inin,inis)
    sys.exit()


def dump_pickle(idx_):
    i,spc1=idx_[0],idx_[1]
    listvla_r=sorted(glob.glob('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/images_LL/2013*'+'.spw.'+str(i)+'_'+str(spc1)+'-'+str(spc1)+'*.FITS'))
    for k in range(len(listvla_r)):
        kk="%04d" % k
        if(os.path.isfile('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/blob_fit/blob_param_spw_'+str(i)+'_'+str(spc1)+'_'+str(kk)+'.p')==False):
            try:
                print listvla_r[k]
                mr,dr=get_evla_submap([listvla_r[k]],0,-1,0,-1);m=40
                hr=mr[0];data=hr.data;data=data[m:-m,m:-m];data[np.isnan(data)]=0;y0max,x0max=np.where(data==np.nanmax(data));inin=[x0max[0],y0max[0],4,4,30]
                #nimg_,params,rms,uncert_cov = fit(data, inin, inis);nimg=nimg_.reshape((2*(128-m),2*(128-m)))
                res = fit(data, inin);nimg=res[0];params=res[2];err=res[3]
                x00=hr.reference_coordinate.Tx.value+((params[0]+m)-(hr.reference_pixel.x.value-1))*hr.scale.axis1.value
                y00=hr.reference_coordinate.Ty.value+((params[1]+m)-(hr.reference_pixel.y.value-1))*hr.scale.axis1.value
                x0=x00;y0=y00;sigx0=params[2]*hr.scale.axis1.value;sigy0=params[3]*hr.scale.axis1.value;rot0=params[4];tb0=nimg[int(params[1]),int(params[0])]
                pickle.dump([x0,y0,sigx0,sigy0,rot0,tb0],open('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/blob_fit/blob_param_spw_'+str(i)+'_'+str(spc1)+'_'+str(kk)+'.p','wb'))
            except:
                pass
    return x0


#spc=['0-15','16-31','32-47','48-63'];n=0
spc=list(np.arange(64));n=0
spw=[0,1]
#spw=[0,5]
allidx=[0]*len(spw)*len(spc)
for l in range(len(spw)):
    for m in range(len(spc)):
        allidx[n]=[spw[l],spc[m]]
        n=n+1
with cf.ProcessPoolExecutor(max_workers=4) as exe:
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

