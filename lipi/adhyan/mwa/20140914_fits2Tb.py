import numpy as np
import pickle
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from surya.radio import get_maps as maps
from surya.utils import main as ut
import sys;import itertools


baseline_filelist=['000-008','000-009','000-010','008-009','008-010','009-010']
freq_filelist=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
baseline_level=['Tile011-Tile021','Tile011-Tile022','Tile011-Tile023','Tile021-Tile022','Tile021-Tile023','Tile022-Tile023']
freq=[109.0,121.0,134.0,147.0,162.0,180.0,198.0,218.0,241.0]
j=8
file_='/nas08-data02/rohit/20140914/pickle/'
#res_path='/nas08-data02/rohit/20151203_MWA/new_ms/fits_residual/'
def arr0_reshape(arr):
    if(len(np.array(arr).shape)==4):
        rr=np.array(arr).swapaxes(0,2).reshape(2,4,25*53)
    if(len(np.array(arr).shape)==5):
        rr=np.array(arr).swapaxes(0,3).reshape(2,4,25*53).swapaxes(0,1)
    return rr

def arr1_reshape(arr):
    return np.array(arr).swapaxes(0,2).reshape(4,4,25*53).swapaxes(0,1)

combine_pickle=0
if(combine_pickle):
    for j in range(9):
        for b in range(6):
            plist=sorted(glob.glob('/nas08-data02/rohit/20140914/pickle/flux_V1_*_'+str(freq_filelist[j])+'.ms_T'+str(baseline_filelist[b])+'.p'))
            bb=[0]*18;bb[17]=[0]*10;b2=[0]*len(plist);b3=[0]*len(plist);b4=[0]*len(plist);b5=[0]*len(plist);b6=[0]*len(plist);b7=[0]*len(plist);b8=[0]*len(plist)
            b9=[0]*len(plist);b10=[0]*len(plist);b11=[0]*len(plist);b12=[0]*len(plist);b13=[0]*len(plist);b14=[0]*len(plist);b15=[0]*len(plist);b16=[0]*len(plist)
            a2=[0]*len(plist);a3=[0]*len(plist);a4=[0]*len(plist);a5=[0]*len(plist);a6=[0]*len(plist);a7=[0]*len(plist);a8=[0]*len(plist);a9=[0]*len(plist)
            for i in range(len(plist)):
                aa=pickle.load(open(plist[i],'rb'));b2[i]=aa[2];b3[i]=aa[3];b4[i]=aa[4];b5[i]=aa[5];b6[i]=aa[6];b7[i]=aa[7];b8[i]=aa[8];b9[i]=aa[9]
                b10[i]=aa[10];b11[i]=aa[11];b12[i]=aa[12];b13[i]=aa[13];b14[i]=aa[14];b15[i]=aa[15];b16[i]=aa[16]
                a2[i]=aa[17][2];a3[i]=aa[17][3];a4[i]=aa[17][4];a5[i]=aa[17][5];a6[i]=aa[17][6];a7[i]=aa[17][7];a8[i]=aa[17][8];a9[i]=aa[17][9]
            bb[0]=aa[0];bb[1]=aa[1];bb[2]=arr0_reshape(b2);bb[3]=arr0_reshape(b3);bb[4]=arr1_reshape(b4);bb[5]=arr1_reshape(b5);bb[6]=arr1_reshape(b6)
            bb[7]=np.array(b7).flatten();bb[8]=np.array(b8).flatten();bb[9]=np.array(b9).flatten();bb[10]=np.array(b10);bb[11]=np.array(b11);bb[12]=np.array(b12);bb[13]=np.array(b13)
            #bb[14]=list(itertools.chain(*b14));bb[15]=list(itertools.chain(*b15));bb[16]=list(itertools.chain(*b16))
            bb[14]=b14;bb[15]=b15;bb[16]=b16
            bb[17][0]=aa[17][0];bb[17][1]=aa[17][1];bb[17][2]=arr1_reshape(a2);bb[17][3]=arr1_reshape(a3);bb[17][4]=arr1_reshape(a4);bb[17][5]=arr1_reshape(a5)
            bb[17][6]=arr1_reshape(a6);bb[17][7]=arr1_reshape(a7);bb[17][8]=arr1_reshape(a8);bb[17][9]=arr1_reshape(a9)
            pickle.dump(bb,open('/nas08-data02/rohit/20140914/pickle/20140914_'+str(freq_filelist[j])+'.ms_T'+str(baseline_filelist[b])+'.p','wb'))
        
img_path='/nas08-data02/rohit/20140914/'+freq_filelist[j]+'/cal_data/images_full/fits/'
ids='20140914_'
f=freq_filelist[j]#'187-188'
freq=freq[j]*1.e6#240*1.e6
imgfiles=sorted(glob.glob(img_path+'*'+str(f)+'*image.FITS'))
#resfiles=sorted(glob.glob(res_path+'*'+str(f)+'*residual.FITS'))
del_=50
angle=-23.86 #Negative angle imples clockwise rotation and vice versa
n=8.0
res=50	
file_solar_center='/nas08-data02/rohit/20140914/20140914_nasa_horizon.dat'
tres=0.5
levi=[0.35]
levf=[0.35]
S_sun,Tb_sun,time,timesec=maps.mean_flux_fave(file_+ids,f,baseline_filelist,12,0) # READ FLUXES
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
        timg=img.split('.')[4].split('~')[0]
        imgsec[i]=int(timg.split(':')[0])*3600+int(timg.split(':')[1])*60+int(timg.split(':')[2])
        imgtime[i]=time[ut.find_nearest(timesec,imgsec[i])[0]]
	h,data_=ut.read_fits(img)
        #resfile=img.split('fits/')[1].split('.ima')[0]+'.residual.FITS'
        #hres,dres=ut.read_fits(res_path+resfile)
	xc,yc,tt=maps.solar_center_pixel(img,file_solar_center)
	centre[i]=np.array([xc,yc])
	print img,imgtime[i],xc,yc
        s_sun[i]=S_sun_t[ut.find_nearest(timesec,imgsec[i])[0]]
        if(data_.max()!=0):
	    Tb[i],flu_,noise_mean[i],noise_std[i],bmaj,bmin,bpa,ndata[i]=maps.compute_Tb(img,int(xc),int(yc),del_,angle,res,freq,n,s_sun[i])
            #resi[i]=maps.get_residuals(res_path+resfile,int(xc)+250,int(yc)+250,del_,angle)[0]
            #Tb_res[i],flux_res[i],tb_frac,flux_frac=maps.scale_residuals(Tb[i],flu_,ndata[i],resi[i])
            nmax[i]=np.max(ndata[i])
        else:
            Tb[i]=np.zeros((100,100))
	i=i+1
pickle.dump([Tb,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles],open('/nas08-data02/rohit/20140914/pickle/Tb_'+str(ids)+f+'.p','wb'))
#pickle.dump([resi,Tb_res,flux_res],open('res_'+str(ids)+f+'.p','wb'))

