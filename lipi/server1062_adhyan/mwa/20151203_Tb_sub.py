import pickle
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from astropy.io import fits
import matplotlib.cm as cm

freq=[108.0,120.0,133.0,145.0,161.0,179.0,197.0,217.0,240.0]

filepath='/nas08-data02/rohit/20151203_MWA/Tb_new/'
Tb240,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_187-188.p','rb'))
Tb217,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_169-170.p','rb'))
Tb179,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_139-140.p','rb'))
Tb197,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_153-154.p','rb'))
Tb145,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_113-114.p','rb'))
Tb133,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_103-104.p','rb'))
Tb108,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_084-085.p','rb'))
Tb120,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles_Tb=pickle.load(open(filepath+'/Tb_20151203_093-094.p','rb'))
Tb240=np.array(Tb240);Tb217=np.array(Tb217);Tb179=np.array(Tb179);Tb145=np.array(Tb145);Tb197=np.array(Tb197);Tb108=np.array(Tb108[0:4010]);Tb120=np.array(Tb120);Tb133=np.array(Tb133)
subfilepath240='/nas08-data02/rohit/20151203_sub_run_mean_subamp/Tb/Tb_20151203_187-188_sub.p'
subfilepath217='/nas08-data02/rohit/20151203_sub_run_mean_subamp/Tb/Tb_20151203_169-170_sub.p'
subfilepath197='/nas08-data02/rohit/20151203_sub_run_mean_subamp/Tb/Tb_20151203_153-154_sub.p'
subfilepath179='/nas08-data02/rohit/20151203_sub_run_mean_subamp/Tb/Tb_20151203_139-140_sub.p'
subfilepath145='/nas08-data02/rohit/20151203_sub_run_mean_subamp/Tb/Tb_20151203_113-114_sub.p'
subfilepath133='/nas08-data02/rohit/20151203_sub_run_mean_subamp/Tb/Tb_20151203_103-104_sub.p'
subfilepath120='/nas08-data02/rohit/20151203_sub_run_mean_subamp/Tb/Tb_20151203_093-094_sub.p'
subfilepath108='/nas08-data02/rohit/20151203_sub_run_mean_subamp/Tb/Tb_20151203_084-085_sub.p'
Tb_sub240,Tb_sub_std240,xc,yc,time_string,time_sec,bmaj240,bmin240,bpa240,ndata=pickle.load(open(subfilepath240,'rb'))
Tb_sub217,Tb_sub_std217,xc,yc,time_string,time_sec,bmaj217,bmin217,bpa217,ndata=pickle.load(open(subfilepath217,'rb'))
Tb_sub197,Tb_sub_std197,xc,yc,time_string,time_sec,bmaj197,bmin197,bpa197,ndata=pickle.load(open(subfilepath197,'rb'))
Tb_sub179,Tb_sub_std179,xc,yc,time_string,time_sec,bmaj179,bmin179,bpa179,ndata=pickle.load(open(subfilepath179,'rb'))
Tb_sub145,Tb_sub_std145,xc,yc,time_string,time_sec,bmaj145,bmin145,bpa145,ndata=pickle.load(open(subfilepath145,'rb'))
Tb_sub133,Tb_sub_std133,xc,yc,time_string,time_sec,bmaj133,bmin133,bpa133,ndata=pickle.load(open(subfilepath133,'rb'))
Tb_sub120,Tb_sub_std120,xc,yc,time_string,time_sec,bmaj120,bmin120,bpa120,ndata=pickle.load(open(subfilepath120,'rb'))
Tb_sub108,Tb_sub_std108,xc,yc,time_string,time_sec,bmaj108,bmin108,bpa108,ndata=pickle.load(open(subfilepath108,'rb'))
#Tb_sub240=Tb_sub240[:,50:-50,50:-50];Tb_sub217=Tb_sub217[:,50:-50,50:-50];Tb_sub179=Tb_sub179[:,50:-50,50:-50];Tb_sub145=Tb_sub145[:,50:-50,50:-50];Tb_sub197=Tb_sub197[:,50:-50,50:-50]
#Tb_sub108=Tb_sub108[:,50:-50,50:-50];Tb_sub120=Tb_sub120[:,50:-50,50:-50]
#Tb_sub133=Tb_sub133[:,50:-50,50:-50]
[Tb_sub108_cor,Tb_sub120_cor,Tb_sub133_cor,Tb_sub145_cor,Tb_sub179_cor,Tb_sub197_cor,Tb_sub217_cor,Tb_sub240_cor]=pickle.load(open(filepath+'Tb_all_cor.p','rb'))

def get_subTb(f):
    syr=[0]*200
    for i in range(200):
        syr[i]=[0]*200
        for j in range(200):
            syr[i][j]=[0]*6
            Tbs=f[:,i,j].reshape((6,592))
            for l in range(6):
                y0=Tbs[l];x0=np.arange(592)
                yr0=y0[np.where((y0>(np.nanmean(y0)-2.2*np.nanstd(y0))) & (y0<(np.nanmean(y0)+2.2*np.nanstd(y0))))]
                xr0=x0[np.where((y0>(np.nanmean(y0)-2.2*np.nanstd(y0))) & (y0<(np.nanmean(y0)+2.2*np.nanstd(y0))))]
                p=np.polyfit(xr0,yr0,10)
                yrnew=np.polyval(p,x0)
                syr[i][j][l]=y0-yrnew
    syr=np.array(syr).reshape(200,200,592*6).swapaxes(0,2).swapaxes(1,2)
    return syr

Tb_sub108_cor=get_subTb(Tb_sub108)
Tb_sub120_cor=get_subTb(Tb_sub120)
Tb_sub133_cor=get_subTb(Tb_sub133)
Tb_sub145_cor=get_subTb(Tb_sub145)
Tb_sub179_cor=get_subTb(Tb_sub179)
Tb_sub197_cor=get_subTb(Tb_sub197)
Tb_sub217_cor=get_subTb(Tb_sub217)
Tb_sub240_cor=get_subTb(Tb_sub240)
pickle.dump([Tb_sub108_cor,Tb_sub120_cor,Tb_sub133_cor,Tb_sub145_cor,Tb_sub179_cor,Tb_sub197_cor,Tb_sub217_cor,Tb_sub240_cor],open(filepath+'Tb_all_cor.p','wb'))
