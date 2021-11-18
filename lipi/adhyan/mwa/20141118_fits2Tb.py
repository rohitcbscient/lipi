import numpy as np
import pickle
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from surya.radio import get_maps as maps
from surya.utils import main as ut
import sys


# start_time: 06:55:14, end_time: 07:00:10, obstime:06:55:19
#np.std(np.concatenate((aa[17][3][0][7:26:,10],aa[17][3][0][36:56:,10])))
baseline_filelist=['000-008','000-009','000-010','008-009','008-010','009-010']
#freq_filelist=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
freq_filelist=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
baseline_level=['Tile011-Tile021','Tile011-Tile022','Tile011-Tile023','Tile021-Tile022','Tile021-Tile023','Tile022-Tile023']
#freq=[109.0,121.0,134.0,147.0,162.0,180.0,198.0,218.0,241.0]
freq=[108,120,133,145,161,179,197,218,240]
pickle_path='/media/rohit/MWA/20141118/pickle/'
img_path='/media/rohit/MWA/20141118/images/'
ids='20141118_'
imgfiles=sorted(glob.glob(img_path+'*I.fits'))
del_=100
angle=20.34 #Negative angle imples clockwise rotation and vice versa
n=0.0
res=99	
tres=0.5
levi=[0.35]
levf=[0.35]
S_sun_fbmean=[0]*len(freq_filelist)
S_sun_bmean=[0]*len(freq_filelist)
ff=0
for f in freq_filelist:
    bb=0
    file_=[0]*len(baseline_filelist)
    for b in baseline_filelist:
        file_[bb]=pickle_path+'flux_V1_1100328928_'+f+'.ms_T'+b+'.p'
        bb=bb+1
    S_sun,Tb_sun,time,timesec=maps.mean_flux_gen(file_,f,baseline_filelist,0.5,0) # READ FLUXES
    S_sun_fbmean[ff]=np.average(S_sun,axis=(0,1))
    S_sun_bmean[ff]=np.average(S_sun,axis=0)
    ff=ff+1

S_sun_bmean=np.array(S_sun_bmean)
S_sun_fbmean=np.array(S_sun_fbmean)

imgfiles.remove(imgfiles[0])
imgfiles.remove(imgfiles[2])
img_freq=[108,120,179,197,218,240]
Tb=[0]*len(imgfiles)
flux=[0]*len(imgfiles)
noise_mean=[0]*len(imgfiles)
noise_std=[0]*len(imgfiles)
Tb_nsig=[0]*len(imgfiles)
Tb_res=[0]*len(imgfiles)
flux_res=[0]*len(imgfiles)
ndata=[0]*len(imgfiles)
nmax=[0]*len(imgfiles)
resi=[0]*len(imgfiles)
i=0
centre=[0]*len(imgfiles)
imgsec0=timesec[0]
tmpsec=[0]*(len(imgfiles)+1)
imgsec=[0]*len(imgfiles)
imgtime=[0]*len(imgfiles)
s_sun=[0]*len(imgfiles)
bimage=[0]*len(imgfiles)
bmaj=[0]*len(imgfiles)
bmin=[0]*len(imgfiles)
bpa=[0]*len(imgfiles)
polTb=[0]*len(imgfiles)
tmpsec[0]=0
for img in imgfiles:
    j=freq.index(img_freq[i])
    h,data_=ut.read_fits(img)
    #xc,yc,tt=maps.solar_center_pixel(img,file_solar_center)
    data=data_[0][0]
    data=np.vstack((np.hstack((data,1.e-4*np.ones((1280,100)))),1.e-4*np.ones((100,1380))))
    bimage[i]=ut.get_bimage(data,0.05)
    xc,yc,w,l,ang=ut.fitEllipse(bimage[i])
    centre[i]=np.array([xc,yc])
    print img,imgtime[i],xc,yc,np.std(data[0:int(xc)-del_*5,:]),np.max(data)
    s_sun[i]=S_sun_fbmean[j,11]
    Tb[i],flu_,noise_mean[i],noise_std[i],bmaj[i],bmin[i],bpa[i],ndata[i]=maps.compute_Tb(img,int(xc),int(yc),del_,angle,res,freq[j]*1.e6,n,s_sun[i])
    polTb[i],r,th=ut.cart2polar(Tb[i])
    nmax[i]=np.max(ndata[i])
    i=i+1

pickle.dump([Tb,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles],open('Tb_'+str(ids)+f+'.p','wb'))

