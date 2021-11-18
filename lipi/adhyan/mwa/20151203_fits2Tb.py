import numpy as np
import pickle
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from surya.radio import get_maps as maps
from surya.utils import main as ut
import sys


baseline_filelist=['000-008','000-009','000-010','008-009','008-010','009-010']
freq_filelist=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
baseline_level=['Tile011-Tile021','Tile011-Tile022','Tile011-Tile023','Tile021-Tile022','Tile021-Tile023','Tile022-Tile023']
freq=[109.0,121.0,134.0,147.0,162.0,180.0,198.0,218.0,241.0]
j=7
file_='/nas08-data02/rohit/20151203_MWA/new_pickle/'
res_path='/nas08-data02/rohit/20151203_MWA/new_ms/fits_residual/'
img_path='/nas08-data02/rohit/20151203_MWA/new_ms/fits/'
ids='20151203_'
f=freq_filelist[j]#'187-188'
freq=freq[j]*1.e6#240*1.e6
imgfiles=sorted(glob.glob(img_path+'*'+str(f)+'*image.FITS'))
resfiles=sorted(glob.glob(res_path+'*'+str(f)+'*residual.FITS'))
del_=50
angle=-15.4 #Negative angle imples clockwise rotation and vice versa
n=5.0
res=50	
file_solar_center='/nas08-data02/rohit/20151203_MWA/20151203_nasa_horizon.dat'
tres=0.5
levi=[0.35]
levf=[0.35]
S_sun,Tb_sun,time,timesec=maps.mean_flux_pfence(file_+ids,f,baseline_filelist,0.5,0) # READ FLUXES
S_sun_t=np.average(S_sun,axis=(0,1))
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
tmpsec[0]=0
for img in imgfiles:
        tmpsec[i+1]=float(img.split('.')[1])
        if(tmpsec[i+1]<tmpsec[i]):
            imgsec0=imgsec0+296
        imgsec[i]=imgsec0+tmpsec[i+1]*0.5
        imgtime[i]=time[ut.find_nearest(timesec,imgsec[i])[0]]
	h,data_=ut.read_fits(img)
        resfile=img.split('fits/')[1].split('.ima')[0]+'.residual.FITS'
        hres,dres=ut.read_fits(res_path+resfile)
	xc,yc,tt=maps.solar_center_pixel(img,file_solar_center)
	centre[i]=np.array([xc,yc])
	print img,imgtime[i],xc,yc
        s_sun[i]=S_sun_t[ut.find_nearest(timesec,imgsec[i])[0]]
	Tb[i],flu_,noise_mean[i],noise_std[i],bmaj,bmin,bpa,ndata[i]=maps.compute_Tb(img,int(xc),int(yc),del_,angle,res,freq,n,s_sun[i])
        resi[i]=maps.get_residuals(res_path+resfile,int(xc)+250,int(yc)+250,del_,angle)[0]
        Tb_res[i],flux_res[i],tb_frac,flux_frac=maps.scale_residuals(Tb[i],flu_,ndata[i],resi[i])
	nmax[i]=np.max(ndata[i])
	i=i+1

pickle.dump([Tb,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles],open('Tb_'+str(ids)+f+'.p','wb'))
pickle.dump([resi,Tb_res,flux_res],open('res_'+str(ids)+f+'.p','wb'))

