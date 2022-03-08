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
freq=[108.0,120.0,133.0,145.0,161.0,179.0,197.0,217.0,240.0]
j=3
#Tb_path='/nas08-data02/rohit/20151203_MWA/Tb_new/'
Tb_path='/nas08-data02/rohit/20151203_MWA/new_ms/fits/'
#ref_path='/nas08-data02/rohit/20151203_sub/ref_images/'
#img_path='/nas08-data02/rohit/20151203_sub/179MHz/images/fits/'
#Tb_write_path='/nas08-data02/rohit/20151203_sub_run_mean_subamp/197MHz/images_run_mean_real/fits/'
Tb_write_path='/nas08-data02/rohit/20151203_sub_run_mean_subamp/running_median/'
ref_path='/nas08-data02/rohit/20151203_sub_old/ref_images/'
#img_path='/nas08-data02/rohit/20151203_sub_run_mean_subamp/197MHz/images_run_mean_real/fits/'
img_path='/nas08-data02/rohit/20151203_sub_run_mean_subamp/running_median/'+str(int(freq[j]))+'MHz/image/fits/'
ids='20151203_'
f=freq_filelist[j]#'187-188'
freq=freq[j]*1.e6#240*1.e6
#imgfiles=sorted(glob.glob(img_path+'*'+str(f)+'*image.FITS'))
#imgfiles=sorted(glob.glob(img_path+'*image.FITS'))
imgfiles=sorted(glob.glob(img_path+'*image.FITS'),key=lambda x:float(x.split('~')[0].split(':')[-2])*60+float(x.split('~')[0].split(':')[-1]))
del_=50
angle=-15.4 #Negative angle imples clockwise rotation and vice versa
n=5.0
res=50	
Tb,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(Tb_path+'/Tb_'+str(ids)+f+'.p','rb'))
ref_h,ref_d=ut.read_fits(ref_path+'ref_'+f+'.image.FITS')
ref_d=ref_d[0][0]
print Tb[9].max(),ndata[9].max(),np.nanmax(ref_d)
fact=Tb[9].max()/ndata[9].max()#*np.nanmax(ref_d)
ref_d[np.isnan(ref_d)]=0
Tb_sub=[0]*len(imgfiles)
Tb_sub_std=[0]*len(imgfiles)
xc=[0]*len(imgfiles)
yc=[0]*len(imgfiles)
time_sec=[0]*len(imgfiles)
time_string=[0]*len(imgfiles)
i=0
for img in imgfiles:
    img_h,img_d=ut.read_fits(img)
    img_d=img_d[0][0]
    #time_string[i]=img.split('time.')[1].split('.')[0].split('~')[0]
    #time_sec[i]=ut.hms2sec_c(' '+time_string[i])
    time_string[i]=img_h['DATE-OBS'].split('T')[1]
    time_sec[i]=ut.hms2sec_c(' '+time_string[i])
    Tb_sub_=img_d*fact
    print img,fact,np.nanmax(Tb_sub_)
    Tb_sub_std[i]=np.nanstd(Tb_sub_[0:1000,0:1000])
    Tb_sub_[np.isnan(Tb_sub_)]=0
    bimage=ut.get_bimage(ref_d,lev=0.2)
    #xc[i],yc[i],w,l,ang=ut.fitEllipse(bimage)
    xc[i],yc[i]=centre[0][0],centre[0][1]
    Tb_sub_=maps.rotateimage(Tb_sub_, angle, int(xc[i]), int(yc[i]))
    Tb_sub[i]=Tb_sub_[int(yc[i])-100:int(yc[i])+100,int(xc[i])-100:int(xc[i])+100]
    ndata[i]=img_d[int(yc[i])-100:int(yc[i])+100,int(xc[i])-100:int(xc[i])+100]
    #Tb_sub[i]=Tb_sub_[int(yc[i])-100:int(yc[i])+100,int(xc[i])-100:int(xc[i])+100]
    i=i+1

Tb_sub=np.array(Tb_sub)
ndata=np.array(ndata)
Tb_sub_std=np.array(Tb_sub_std)
xc=np.array(xc)
yc=np.array(yc)
time_sec=np.array(time_sec)

pickle.dump([Tb_sub,Tb_sub_std,xc,yc,time_string,time_sec,bmaj,bmin,bpa,ndata],open(Tb_write_path+'/Tb_'+str(ids)+f+'_sub_test.p','wb'))

data5s=Tb_sub*1.0;ndata5s=Tb_sub*1.0
for i in range(data5s.shape[0]):
    s=np.nanstd(data5s[i])
    data5s[i][np.where(data5s[i]<5*s)]=0
    ndata5s[i][np.where(ndata5s[i]>-5*s)]=0

idx=np.where((Tb_sub_std<700) & (Tb_sub_std>660))


