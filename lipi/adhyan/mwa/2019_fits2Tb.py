import numpy as np
import pickle
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from surya.radio import get_maps as maps
from surya.utils import main as ut
import sys


baseline_filelist=['000-007']
freq_filelist=['084-084','089-089','128-128','150-150','188-188']
baseline_level=['Tile011-Tile021','Tile011-Tile022','Tile011-Tile023','Tile021-Tile022','Tile021-Tile023','Tile022-Tile023']
freq=[107.5,113.9,163.8,192.0,240.6]
j=0
file_='/media/rohit/MWA/20190405_PSP_MWA/sun_1238472032/all_pickle/'
img_path='/media/rohit/MWA/20190405_PSP_MWA/sun_1238472032/all_images/'
ids='flux_V1_1238472032_'
f=freq_filelist[j]#'187-188'
freq=freq[j]*1.e6#240*1.e6
imgfiles=sorted(glob.glob(img_path+'*'+str(f)+'*image.FITS'))
del_=125
angle=27.0 #Negative angle imples clockwise rotation and vice versa
n=5.0
res=20	
file_solar_center='/nas08-data02/rohit/20151203_MWA/20151203_nasa_horizon.dat'
tres=0.5
levi=[0.35]
levf=[0.35]
centre=[[666,222],[669,223],[],[],[],[]]
S_sun,Tb_sun,time,timesec=maps.mean_flux_pfence(file_+ids,f,baseline_filelist,10,1) # READ FLUXES
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
        imgsec[i]=imgsec0+10
        imgtime[i]=time[ut.find_nearest(timesec,imgsec[i])[0]]
	h,data_=ut.read_fits(img)
        if(j==0 or j==1):
            bimage=ut.get_bimage(data_[0][0],0.35)
            xc,yc,w,l,ang=ut.fitEllipse(bimage)
            #xc=centre[j][0];yc=centre[j][1]
            centre[i]=np.array([xc,yc])
        else:
            yc=np.where(data_[0][0]==np.nanmax(data_[0][0]))[0]-4
            xc=np.where(data_[0][0]==np.nanmax(data_[0][0]))[1]-20
            centre[i]=np.array([xc,yc])
	print img,imgtime[i],xc,yc
        s_sun[i]=S_sun_t[i]
	Tb[i],flu_,noise_mean[i],noise_std[i],bmaj,bmin,bpa,ndata[i]=maps.compute_Tb(img,int(xc),int(yc),del_,angle,res,freq,n,s_sun[i])
	nmax[i]=np.max(ndata[i])
	i=i+1

pickle.dump([Tb,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles],open('Tb_'+str(ids)+f+'.p','wb'))

