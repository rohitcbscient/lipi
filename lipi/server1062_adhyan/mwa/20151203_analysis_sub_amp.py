import pickle
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from astropy.io import fits
import matplotlib.cm as cm
from scipy.io import readsav

freq=[108.0,120.0,133.0,145.0,161.0,179.0,197.0,217.0,240.0]
forward_100MHz=readsav('/nas08-data02/rohit/20151203_forward/20151203_100MHz_psimas.sav')

pres=forward_100MHz['modsolstruct'][0][0]
dens=forward_100MHz['modsolstruct'][0][1]
temp=forward_100MHz['modsolstruct'][0][2]
br=forward_100MHz['modsolstruct'][0][3]
bth=forward_100MHz['modsolstruct'][0][4]
bph=forward_100MHz['modsolstruct'][0][5]
vr=forward_100MHz['modsolstruct'][0][6]
vth=forward_100MHz['modsolstruct'][0][7]
vph=forward_100MHz['modsolstruct'][0][8]
babs=np.sqrt(br*br+bth*bth+bph*bph)
vabs=np.sqrt(vr*vr+vth*vth+vph*vph)


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
Tb_sub240=Tb_sub240[:,50:-50,50:-50];Tb_sub217=Tb_sub217[:,50:-50,50:-50];Tb_sub179=Tb_sub179[:,50:-50,50:-50];Tb_sub145=Tb_sub145[:,50:-50,50:-50];Tb_sub197=Tb_sub197[:,50:-50,50:-50]
#Tb_sub108=Tb_sub108[:,50:-50,50:-50];Tb_sub120=Tb_sub120[:,50:-50,50:-50]
Tb_sub133=Tb_sub133[:,50:-50,50:-50]
#ndata=ndata[:,50:-50,50:-50]
Tb_240=np.vstack((np.nanmean(Tb240[0:1776,:,:].reshape(2,888,100,100),axis=0),np.nanmean(Tb240[-1777:-1,:,:].reshape(2,888,100,100),axis=0)))
Tb_217=np.vstack((np.nanmean(Tb217[0:1776,:,:].reshape(2,888,100,100),axis=0),np.nanmean(Tb217[-1777:-1,:,:].reshape(2,888,100,100),axis=0)))
Tb_197=np.vstack((np.nanmean(Tb197[0:1776,:,:].reshape(2,888,100,100),axis=0),np.nanmean(Tb197[-1777:-1,:,:].reshape(2,888,100,100),axis=0)))
Tb_179=np.vstack((np.nanmean(Tb179[0:1776,:,:].reshape(2,888,100,100),axis=0),np.nanmean(Tb179[-1777:-1,:,:].reshape(2,888,100,100),axis=0)))
Tb_145=np.vstack((np.nanmean(Tb145[0:1776,:,:].reshape(2,888,100,100),axis=0),np.nanmean(Tb145[-1777:-1,:,:].reshape(2,888,100,100),axis=0)))
Tb_133=np.vstack((np.nanmean(Tb133[0:1776,:,:].reshape(2,888,100,100),axis=0),np.nanmean(Tb133[-1777:-1,:,:].reshape(2,888,100,100),axis=0)))
Tb_120=np.vstack((np.nanmean(Tb120[0:1776,:,:].reshape(2,888,100,100),axis=0),np.nanmean(Tb120[-1777:-1,:,:].reshape(2,888,100,100),axis=0)))
Tb_108=np.vstack((np.nanmean(Tb108[0:1776,:,:].reshape(2,888,100,100),axis=0),np.nanmean(Tb108[-1777:-1,:,:].reshape(2,888,100,100),axis=0)))
[Tb_sub108_cor,Tb_sub120_cor,Tb_sub133_cor,Tb_sub145_cor,Tb_sub179_cor,Tb_sub197_cor,Tb_sub217_cor,Tb_sub240_cor]=pickle.load(open(filepath+'Tb_all_cor.p','rb'))

Tb_sub108_cor_sel=np.vstack((Tb_sub108_cor[10:400],Tb_sub108_cor[650:1500],Tb_sub108_cor[2400:2900]))
Tb_sub120_cor_sel=np.vstack((Tb_sub120_cor[10:400],Tb_sub120_cor[650:1500],Tb_sub120_cor[1800:2300],Tb_sub120_cor[2400:2900],Tb_sub120_cor[3000:3500]))
Tb_sub133_cor_sel=np.vstack((Tb_sub133_cor[10:400],Tb_sub133_cor[650:1500],Tb_sub133_cor[1800:2300],Tb_sub133_cor[2400:2900],Tb_sub133_cor[3000:3500]))
Tb_sub145_cor_sel=np.vstack((Tb_sub145_cor[10:400],Tb_sub145_cor[650:1500],Tb_sub145_cor[1800:2300],Tb_sub145_cor[2400:2900],Tb_sub145_cor[3000:3500]))
Tb_sub179_cor_sel=np.vstack((Tb_sub179_cor[10:400],Tb_sub179_cor[650:1500],Tb_sub179_cor[1800:2300],Tb_sub179_cor[2400:2900],Tb_sub179_cor[3000:3500]))
Tb_sub197_cor_sel=np.vstack((Tb_sub197_cor[10:400],Tb_sub197_cor[650:1500],Tb_sub197_cor[1800:2300],Tb_sub197_cor[3000:3500]))
Tb_sub217_cor_sel=np.vstack((Tb_sub217_cor[10:400],Tb_sub217_cor[650:1500],Tb_sub217_cor[1800:2300],Tb_sub217_cor[2400:2900],Tb_sub217_cor[3000:3500]))
Tb_sub240_cor_sel=np.vstack((Tb_sub240_cor[10:400],Tb_sub240_cor[650:1500],Tb_sub240_cor[1800:2300],Tb_sub240_cor[2400:2900],Tb_sub240_cor[3000:3500]))

Tblim=[4000,4000,3000,2000,2000,2000,2000,2000]
Tb_sub108_cor_sel_lim=Tb_sub108_cor_sel*1.0;Tb_sub108_cor_sel_lim[Tb_sub108_cor_sel_lim<Tblim[0]]=np.nan
Tb_sub120_cor_sel_lim=Tb_sub120_cor_sel*1.0;Tb_sub120_cor_sel_lim[Tb_sub120_cor_sel_lim<Tblim[1]]=np.nan
Tb_sub133_cor_sel_lim=Tb_sub133_cor_sel*1.0;Tb_sub133_cor_sel_lim[Tb_sub133_cor_sel_lim<Tblim[2]]=np.nan
Tb_sub145_cor_sel_lim=Tb_sub145_cor_sel*1.0;Tb_sub145_cor_sel_lim[Tb_sub145_cor_sel_lim<Tblim[3]]=np.nan
Tb_sub179_cor_sel_lim=Tb_sub179_cor_sel*1.0;Tb_sub179_cor_sel_lim[Tb_sub179_cor_sel_lim<Tblim[4]]=np.nan
Tb_sub197_cor_sel_lim=Tb_sub197_cor_sel*1.0;Tb_sub197_cor_sel_lim[Tb_sub197_cor_sel_lim<Tblim[5]]=np.nan
Tb_sub217_cor_sel_lim=Tb_sub217_cor_sel*1.0;Tb_sub217_cor_sel_lim[Tb_sub217_cor_sel_lim<Tblim[6]]=np.nan
Tb_sub240_cor_sel_lim=Tb_sub240_cor_sel*1.0;Tb_sub240_cor_sel_lim[Tb_sub240_cor_sel_lim<Tblim[7]]=np.nan

Tb_sub108_cor_sel_lim_mean=np.nanmean(Tb_sub108_cor_sel_lim,axis=0);Tb_sub108_cor_sel_lim_mean[Tb_sub108_cor_sel_lim_mean<Tblim[0]]=np.nan
Tb_sub120_cor_sel_lim_mean=np.nanmean(Tb_sub120_cor_sel_lim,axis=0);Tb_sub120_cor_sel_lim_mean[Tb_sub120_cor_sel_lim_mean<Tblim[0]]=np.nan
Tb_sub133_cor_sel_lim_mean=np.nanmean(Tb_sub133_cor_sel_lim,axis=0);Tb_sub133_cor_sel_lim_mean[Tb_sub133_cor_sel_lim_mean<Tblim[2]]=np.nan
Tb_sub145_cor_sel_lim_mean=np.nanmean(Tb_sub145_cor_sel_lim,axis=0);Tb_sub145_cor_sel_lim_mean[Tb_sub145_cor_sel_lim_mean<Tblim[3]]=np.nan
Tb_sub179_cor_sel_lim_mean=np.nanmean(Tb_sub179_cor_sel_lim,axis=0);Tb_sub179_cor_sel_lim_mean[Tb_sub179_cor_sel_lim_mean<Tblim[4]]=np.nan
Tb_sub197_cor_sel_lim_mean=np.nanmean(Tb_sub197_cor_sel_lim,axis=0);Tb_sub197_cor_sel_lim_mean[Tb_sub197_cor_sel_lim_mean<Tblim[5]]=np.nan
Tb_sub217_cor_sel_lim_mean=np.nanmean(Tb_sub217_cor_sel_lim,axis=0);Tb_sub217_cor_sel_lim_mean[Tb_sub217_cor_sel_lim_mean<Tblim[6]]=np.nan
Tb_sub240_cor_sel_lim_mean=np.nanmean(Tb_sub240_cor_sel_lim,axis=0);Tb_sub240_cor_sel_lim_mean[Tb_sub240_cor_sel_lim_mean<Tblim[7]]=np.nan

Tb_sub108_cor_mean=np.nanmean(Tb_sub108_cor_sel,axis=0)
Tb_sub120_cor_mean=np.nanmean(Tb_sub120_cor_sel,axis=0)
Tb_sub133_cor_mean=np.nanmean(Tb_sub133_cor_sel,axis=0)
Tb_sub145_cor_mean=np.nanmean(Tb_sub145_cor_sel,axis=0)
Tb_sub179_cor_mean=np.nanmean(Tb_sub179_cor_sel,axis=0)
Tb_sub197_cor_mean=np.nanmean(Tb_sub197_cor_sel,axis=0)
Tb_sub217_cor_mean=np.nanmean(Tb_sub217_cor_sel,axis=0)
Tb_sub240_cor_mean=np.nanmean(Tb_sub240_cor_sel,axis=0)

if(0):
    for i in range(3552):
        Tb_sub240[i][np.where(Tb_sub240[i]<Tb_sub_std240[i])]=0
        Tb_sub217[i][np.where(Tb_sub217[i]<Tb_sub_std217[i])]=0
        Tb_sub197[i][np.where(Tb_sub197[i]<Tb_sub_std197[i])]=0
        Tb_sub179[i][np.where(Tb_sub179[i]<Tb_sub_std179[i])]=0
        Tb_sub145[i][np.where(Tb_sub145[i]<Tb_sub_std145[i])]=0
        Tb_sub133[i][np.where(Tb_sub133[i]<Tb_sub_std133[i])]=0
        Tb_sub120[i][np.where(Tb_sub120[i]<Tb_sub_std120[i])]=0
        Tb_sub108[i][np.where(Tb_sub108[i]<Tb_sub_std108[i])]=0

Tbmax240=Tb_sub240_cor.max(axis=(1,2))
Tbmax217=Tb_sub217_cor.max(axis=(1,2))
Tbmax197=Tb_sub197_cor.max(axis=(1,2))
Tbmax179=Tb_sub179_cor.max(axis=(1,2))
Tbmax145=Tb_sub145_cor.max(axis=(1,2))
Tbmax133=Tb_sub133_cor.max(axis=(1,2))
Tbmax120=Tb_sub120_cor.max(axis=(1,2))
Tbmax108=Tb_sub108_cor.max(axis=(1,2))


DR=[0]*4011;mx=[0]*4011;
for i in range(4011):
    mx[i]=ndata[i].max()
    DR[i]=ndata[i].max()/noise_std[i]

from surya.utils import main as ut

def get_position(idx,Tb_sub):
    cx=[0]*len(idx);cy=[0]*len(idx)
    for i in range(len(idx)):
        if(Tb_sub[idx[i]].max()!=0):
            bi=ut.get_bimage(Tb_sub[idx[i]],0.95)
            x,y=ut.fitEllipse(bi)[0:2]
            cx[i]=(x-100)*50;cy[i]=(y-100)*50
        else:
            cx[i]=0;cy[i]=0
        #cx[i]=((np.where(Tb_sub[idx[i]]==np.nanmax(Tb_sub[idx[i]]))[0]-50)[0]-100)*50
        #cy[i]=((np.where(Tb_sub[idx[i]]==np.nanmax(Tb_sub[idx[i]]))[1]-50)[0]-100)*50
    return np.array(cx),np.array(cy)


aiafile='/nas08-data02/rohit/20151203_EUV/aia.lev1.193A_2015-12-03T03_20_17.84Z.image_lev1.fits'
af=fits.open(aiafile)
h=af[0].header;d=af[0].data

lascoc2_file='/nas08-data02/rohit/20151203_LASCO/22572578_1.fts'
c2=fits.open(lascoc2_file)
hc=c2[0].header;dc=c2[0].data

yy=Tb_179;cy=cy179;cx=cx179;Tbmax=Tbmax179;idx=idx179
yy=Tb_120;cy=cy120;cx=cx120;Tbmax=Tbmax120;idx=idx120
levels=[0.3,0.5,0.7,0.9]
levels_sub=[0.5,0.6,0.7,0.8,0.9]
#plt.imshow(yy[0]/1.e6,aspect='equal',interpolation='none',extent=[-2500,2500,-2500,2500],origin='lower',vmin=0,vmax=0.5,cmap='YlGnBu')
plt.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
plt.contour(yy[0]/np.nanmax(yy[0]), levels,extent=[-2500,2500,-2500,2500],  colors='black',linewidths=1)
ax0=plt.scatter(cy,cx,c=Tbmax[idx],s=7,cmap='jet',vmin=1000,vmax=8000)
plt.colorbar(ax0,label='$T_B$(K)')
plt.xlabel('Solar X (arcsec)')
plt.ylabel('Solar Y (arcsec)')
plt.title('Location (T$_B$>5000 K)')
plt.xlim(-1600,1600);plt.ylim(-1600,1600)
plt.show()

Tb_108=Tb108;Tb_120=Tb120;Tb_133=Tb133;Tb_145=Tb145;Tb_179=Tb179;Tb_197=Tb197;Tb_217=Tb217;Tb_240=Tb240
Tbmax108[Tbmax108<1.e3]=np.nan;Tbmax120[Tbmax120<1.e3]=np.nan;Tbmax133[Tbmax133<1.e3]=np.nan;Tbmax145[Tbmax145<1.e3]=np.nan
Tbmax217[Tbmax217<1.e3]=np.nan;Tbmax240[Tbmax240<1.e3]=np.nan;Tbmax179[Tbmax179<1.e3]=np.nan;Tbmax197[Tbmax197<1.e3]=np.nan
Tb_sub_std108[np.where(Tb_sub_std108<100)]=np.nan;Tb_sub_std120[np.where(Tb_sub_std120<100)]=np.nan;Tb_sub_std133[np.where(Tb_sub_std133<100)]=np.nan;Tb_sub_std145[np.where(Tb_sub_std145<100)]=np.nan
Tb_sub_std179[np.where(Tb_sub_std179<100)]=np.nan;Tb_sub_std197[np.where(Tb_sub_std197<100)]=np.nan;Tb_sub_std217[np.where(Tb_sub_std217<100)]=np.nan;Tb_sub_std240[np.where(Tb_sub_std240<100)]=np.nan

t=np.arange(3552)*0.5
f,ax=plt.subplots(2,4,figsize=(12, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.plot(t,Tbmax108,'o-',markersize=4,label='(a) 108 MHz',color='k')
ax00.plot(t,Tb_sub_std108,'o-',markersize=4,color='r');ax00.set_ylabel('T$_{B}$ (K)');ax00.set_ylim(1.e2,5.e6);ax00.set_yscale('log')
ax01.plot(t,Tbmax120,'o-',markersize=4,label='(b) 120 MHz',color='k')
ax01.plot(t,Tb_sub_std120,'o-',markersize=4,color='r');ax01.set_ylim(1.e2,5.e6);ax01.set_yscale('log')
ax02.plot(t,Tbmax133,'o-',markersize=4,label='(c) 133 MHz',color='k')
ax02.plot(t,Tb_sub_std133,'o-',markersize=4,color='r');ax02.set_ylim(1.e2,5.e6);ax02.set_yscale('log')
ax03.plot(t,Tbmax145,'o-',markersize=4,label='(d) 145 MHz',color='k')
ax03.plot(t,Tb_sub_std145,'o-',markersize=4,color='r');ax03.set_ylim(1.e2,5.e6);ax03.set_yscale('log')
ax10.plot(t,Tbmax179,'o-',markersize=4,label='(e) 179 MHz',color='k')
ax10.plot(t,Tb_sub_std179,'o-',markersize=4,color='r');ax10.set_ylabel('T$_{B}$ (K)');ax10.set_ylim(1.e2,1.e6);ax10.set_yscale('log');ax10.set_xlabel('Time (s)')
ax11.plot(t,Tbmax197,'o-',markersize=4,label='(f) 197 MHz',color='k')
ax11.plot(t,Tb_sub_std197,'o-',markersize=4,color='r');ax11.set_ylim(1.e2,5.e6);ax11.set_yscale('log');ax11.set_xlabel('Time (s)')
ax12.plot(t,Tbmax217,'o-',markersize=4,label='(g) 217 MHz',color='k')
ax12.plot(t,Tb_sub_std217,'o-',markersize=4,color='r');ax12.set_ylim(1.e2,5.e6);ax12.set_yscale('log');ax12.set_xlabel('Time (s)')
ax13.plot(t,Tbmax240,'o-',markersize=4,label='(h) 240 MHz',color='k')
ax13.plot(t,Tb_sub_std240,'o-',markersize=4,color='r');ax13.set_ylim(1.e2,5.e6);ax13.set_yscale('log');ax13.set_xlabel('Time (s)')
ax00.legend();ax01.legend();ax02.legend();ax03.legend();ax10.legend();ax11.legend();ax12.legend();ax13.legend()
plt.show()

plt.plot(Tbmax240,'-')
plt.plot(Tbmax217,'-')
plt.plot(Tbmax179,'-')
plt.plot(Tbmax120,'-')
plt.plot(Tbmax108,'-')
plt.xlabel('Time (sec)')
plt.ylabel('T$_B$ (K)')
plt.show()
#Tblim=1000
#Tb_sub108[np.where(Tb_sub108<Tblim)]=0;Tb_sub120[np.where(Tb_sub120<Tblim)]=0;Tb_sub133[np.where(Tb_sub133<Tblim)]=0;Tb_sub145[np.where(Tb_sub145<Tblim)]=0;Tb_sub217[np.where(Tb_sub217<Tblim)]=0;Tb_sub240[np.where(Tb_sub240<Tblim)]=0
#Tb_sub179[np.where(Tb_sub179<Tblim)]=0;Tb_sub197[np.where(Tb_sub197<Tblim)]=0

l1=np.linspace(0,200,51)#;l2=np.linspace(20,80,13)[1:]
xxx,yyy=np.meshgrid(l1, l1)
Tb_region=[0]*8;Tb_region_mean=[0]*8;Tb_region1=[0]*8;Tb_region2=[0]*8;Tb_region3=[0]*8;Tb_region4=[0]*8;Tb_region5=[0]*8;Tb_region6=[0]*8
Tb_cregion=[0]*8;Tb_cregion_mean=[0]*8;Tb_cregion1=[0]*8;Tb_cregion2=[0]*8;Tb_cregion3=[0]*8;Tb_cregion4=[0]*8;Tb_cregion5=[0]*8;Tb_cregion6=[0]*8
ts1=[0]*8;ts2=[0]*8;ts3=[0]*8;ts4=[0]*8;ts5=[0]*8;ts6=[0]*8;ts7=[0]*8;ts8=[0]*8
#Tb_suball=[Tb_sub108,Tb_sub120,Tb_sub133,Tb_sub145,Tb_sub179,Tb_sub197,Tb_sub217,Tb_sub240]
#Tb_suball=[Tb_sub108_cor[:,50:-50,50:-50],Tb_sub120_cor[:,50:-50,50:-50],Tb_sub133_cor[:,50:-50,50:-50],Tb_sub145_cor[:,50:-50,50:-50],Tb_sub179_cor[:,50:-50,50:-50],Tb_sub197_cor[:,50:-50,50:-50],Tb_sub217_cor[:,50:-50,50:-50],Tb_sub240_cor[:,50:-50,50:-50]]
#Tb_suball=[Tb_sub108_cor_sel[:,50:-50,50:-50],Tb_sub120_cor_sel[:,50:-50,50:-50],Tb_sub133_cor_sel[:,50:-50,50:-50],Tb_sub145_cor_sel[:,50:-50,50:-50],Tb_sub179_cor_sel[:,50:-50,50:-50],Tb_sub197_cor_sel[:,50:-50,50:-50],Tb_sub217_cor_sel[:,50:-50,50:-50],Tb_sub240_cor_sel[:,50:-50,50:-50]]
Tb_suball=[Tb_sub108_cor_sel,Tb_sub120_cor_sel,Tb_sub133_cor_sel,Tb_sub145_cor_sel,Tb_sub179_cor_sel,Tb_sub197_cor_sel,Tb_sub217_cor_sel,Tb_sub240_cor_sel]
idx108=np.where(Tb_sub108_cor_sel.max(axis=(1,2))>Tblim[0])[0];cx108,cy108=get_position(np.arange(Tb_sub108_cor_sel.shape[0]),Tb_sub108_cor_sel)
idx120=np.where(Tb_sub120_cor_sel.max(axis=(1,2))>Tblim[1])[0];cx120,cy120=get_position(np.arange(Tb_sub120_cor_sel.shape[0]),Tb_sub120_cor_sel)
idx133=np.where(Tb_sub133_cor_sel.max(axis=(1,2))>Tblim[2])[0];cx133,cy133=get_position(np.arange(Tb_sub133_cor_sel.shape[0]),Tb_sub133_cor_sel)
idx145=np.where(Tb_sub145_cor_sel.max(axis=(1,2))>Tblim[3])[0];cx145,cy145=get_position(np.arange(Tb_sub145_cor_sel.shape[0]),Tb_sub145_cor_sel)
idx179=np.where(Tb_sub179_cor_sel.max(axis=(1,2))>Tblim[4])[0];cx179,cy179=get_position(np.arange(Tb_sub179_cor_sel.shape[0]),Tb_sub179_cor_sel)
idx197=np.where(Tb_sub197_cor_sel.max(axis=(1,2))>Tblim[5])[0];cx197,cy197=get_position(np.arange(Tb_sub197_cor_sel.shape[0]),Tb_sub197_cor_sel)
idx217=np.where(Tb_sub217_cor_sel.max(axis=(1,2))>Tblim[6])[0];cx217,cy217=get_position(np.arange(Tb_sub217_cor_sel.shape[0]),Tb_sub217_cor_sel)
idx240=np.where(Tb_sub240_cor_sel.max(axis=(1,2))>Tblim[7])[0];cx240,cy240=get_position(np.arange(Tb_sub240_cor_sel.shape[0]),Tb_sub240_cor_sel)
idxr1=[0]*8;idxr2=[0]*8;idxr3=[0]*8;idxr4=[0]*8;idxr5=[0]*8;idxr6=[0]*8
for j in range(8):
    Tb_region[j]=[0]*50
    for i in range(50):
        Tb_region[j][i]=[0]*50
        for k in range(50):
            Tb_region[j][i][k]=np.nanmean(Tb_suball[j][:,int(xxx[k+1][i+1]-3):int(xxx[k+1][i+1])+3,int(yyy[k+1][i+1]-3):int(yyy[k+1][i+1]+3)],axis=(1,2))
    Tb_region_mean[j]=np.nanmean(Tb_region[j],axis=2)
    Tb_region1[j]=np.array(Tb_region[j])[26,28];Tb_region2[j]=np.array(Tb_region[j])[24,18];Tb_region3[j]=np.array(Tb_region[j])[26,23]
    Tb_region4[j]=np.array(Tb_region[j])[21,24];Tb_region5[j]=np.array(Tb_region[j])[20,16];Tb_region6[j]=np.array(Tb_region[j])[16,24]
    Tb_cregion1[j]=Tb_region1[j]*1.0;Tb_cregion2[j]=Tb_region2[j]*1.0;Tb_cregion3[j]=Tb_region3[j]*1.0;Tb_cregion4[j]=Tb_region4[j]*1.0;Tb_cregion5[j]=Tb_region5[j]*1.0;Tb_cregion6[j]=Tb_region6[j]*1.0
    Tb_cregion1[j][Tb_cregion1[j]<Tblim[j]]=0
    Tb_cregion2[j][Tb_cregion2[j]<Tblim[j]]=0
    Tb_cregion3[j][Tb_cregion3[j]<Tblim[j]]=0
    Tb_cregion4[j][Tb_cregion4[j]<Tblim[j]]=0
    Tb_cregion5[j][Tb_cregion5[j]<Tblim[j]]=0
    Tb_cregion6[j][Tb_cregion6[j]<Tblim[j]]=0
    idx1=np.where(Tb_cregion1[j]!=0)[0];diffs = np.diff(idx1) != 1;indexes = np.nonzero(diffs)[0] + 1;ts1[j] = np.split(idx1, indexes)
    idx2=np.where(Tb_cregion2[j]!=0)[0];diffs = np.diff(idx2) != 1;indexes = np.nonzero(diffs)[0] + 1;ts2[j] = np.split(idx2, indexes)
    idx3=np.where(Tb_cregion3[j]!=0)[0];diffs = np.diff(idx3) != 1;indexes = np.nonzero(diffs)[0] + 1;ts3[j] = np.split(idx3, indexes)
    idx4=np.where(Tb_cregion4[j]!=0)[0];diffs = np.diff(idx4) != 1;indexes = np.nonzero(diffs)[0] + 1;ts4[j] = np.split(idx4, indexes)
    idx5=np.where(Tb_cregion5[j]!=0)[0];diffs = np.diff(idx5) != 1;indexes = np.nonzero(diffs)[0] + 1;ts5[j] = np.split(idx5, indexes)
    idx6=np.where(Tb_cregion6[j]!=0)[0];diffs = np.diff(idx6) != 1;indexes = np.nonzero(diffs)[0] + 1;ts6[j] = np.split(idx6, indexes)
    idxr1[j]=idx1;idxr2[j]=idx2;idxr3[j]=idx3;idxr4[j]=idx4;idxr5[j]=idx5;idxr6[j]=idx6


## HISTOGRAM
nts1=[0]*8;tts1=[0]*8;nhist1=[0]*8
nts2=[0]*8;tts2=[0]*8;nhist2=[0]*8
nts3=[0]*8;tts3=[0]*8;nhist3=[0]*8
nts4=[0]*8;tts4=[0]*8;nhist4=[0]*8
nts5=[0]*8;tts5=[0]*8;nhist5=[0]*8
nts6=[0]*8;tts6=[0]*8;nhist6=[0]*8
bins=np.arange(1,10)*0.5
for i in range(8):
    nts1[i]=len(ts1[i]);tts1[i]=[0]*(len(ts1[i])-1)
    for j in range(len(ts1[i])-1):
        tts1[i][j]=len(ts1[i][j])*0.5
    nhist1[i]=np.histogram(tts1[i],bins=bins)[0]
    nts2[i]=len(ts2[i]);tts2[i]=[0]*(len(ts2[i])-1)
    for j in range(len(ts2[i])-1):
        tts2[i][j]=len(ts2[i][j])*0.5
    nhist2[i]=np.histogram(tts2[i],bins=bins)[0]
    nts3[i]=len(ts3[i]);tts3[i]=[0]*(len(ts3[i])-1)
    for j in range(len(ts3[i])-1):
        tts3[i][j]=len(ts3[i][j])*0.5
    nhist3[i]=np.histogram(tts3[i],bins=bins)[0]
    nts4[i]=len(ts4[i]);tts4[i]=[0]*(len(ts4[i])-1)
    for j in range(len(ts4[i])-1):
        tts4[i][j]=len(ts4[i][j])*0.5
    nhist4[i]=np.histogram(tts4[i],bins=bins)[0]
    nts5[i]=len(ts5[i]);tts5[i]=[0]*(len(ts5[i])-1)
    for j in range(len(ts5[i])-1):
        tts5[i][j]=len(ts5[i][j])*0.5
    nhist5[i]=np.histogram(tts5[i],bins=bins)[0]
    nts6[i]=len(ts6[i]);tts6[i]=[0]*(len(ts6[i])-1)
    for j in range(len(ts6[i])-1):
        tts6[i][j]=len(ts6[i][j])*0.5
    nhist6[i]=np.histogram(tts6[i],bins=bins)[0]
nhist1=np.array(nhist1);nhist2=np.array(nhist2);nhist3=np.array(nhist3);nhist4=np.array(nhist4);nhist5=np.array(nhist5);nhist6=np.array(nhist6)


from mpl_toolkits.axes_grid1 import make_axes_locatable

f,ax=plt.subplots(2,3)
ax[0,0].imshow(nhist1,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[0,0].set_yticks(np.linspace(0,7,8));ax[0,0].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[0,0].set_xticks(np.linspace(0,8,9));ax[0,0].set_xticklabels(bins)
ax[0,0].text(6, 1, 'Region 1', horizontalalignment='center', verticalalignment='center',fontsize=15)
ax[0,1].imshow(nhist2,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[0,1].set_yticks(np.linspace(0,7,8));ax[0,1].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[0,1].set_xticks(np.linspace(0,8,9));ax[0,1].set_xticklabels(bins)
ax[0,1].text(6, 1, 'Region 2', horizontalalignment='center', verticalalignment='center',fontsize=15)
im1=ax[0,2].imshow(nhist3,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[0,2].set_yticks(np.linspace(0,7,8));ax[0,2].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[0,2].set_xticks(np.linspace(0,8,9));ax[0,2].set_xticklabels(bins)
ax[0,2].text(6, 1, 'Region 3', horizontalalignment='center', verticalalignment='center',fontsize=15)
ax[1,0].imshow(nhist4,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[1,0].set_yticks(np.linspace(0,7,8));ax[1,0].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[1,0].set_xticks(np.linspace(0,8,9));ax[1,0].set_xticklabels(bins)
ax[1,0].text(6, 1, 'Region 4', horizontalalignment='center', verticalalignment='center',fontsize=15)
ax[1,1].imshow(nhist5,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[1,1].set_yticks(np.linspace(0,7,8));ax[1,1].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[1,1].set_xticks(np.linspace(0,8,9));ax[1,1].set_xticklabels(bins)
ax[1,1].text(6, 1, 'Region 5', horizontalalignment='center', verticalalignment='center',fontsize=15)
ax[1,2].imshow(nhist6,aspect='auto',origin=0,vmin=5,vmax=300,cmap='jet',norm=matplotlib.colors.LogNorm())
ax[1,2].set_yticks(np.linspace(0,7,8));ax[1,2].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[1,2].set_xticks(np.linspace(0,8,9));ax[1,2].set_xticklabels(bins)
ax[1,2].text(6, 1, 'Region 6', horizontalalignment='center', verticalalignment='center',fontsize=15)
divider = make_axes_locatable(ax[0,2])
cax = divider.append_axes('right', size='5%', pad=0.05)
f.colorbar(im1, cax=cax, orientation='vertical')
ax[0,0].set_ylabel('Frequency (MHz)');ax[1,1].set_xlabel('Time bins (sec)')
plt.show()





ds_region1=np.zeros((8,2740));ds_cregion1=np.zeros((8,2740))
ds_region2=np.zeros((8,2740));ds_cregion2=np.zeros((8,2740))
ds_region3=np.zeros((8,2740));ds_cregion3=np.zeros((8,2740))
ds_region4=np.zeros((8,2740));ds_cregion4=np.zeros((8,2740))
ds_region5=np.zeros((8,2740));ds_cregion5=np.zeros((8,2740))
ds_region6=np.zeros((8,2740));ds_cregion6=np.zeros((8,2740))
for i in range(8):
    if(i==0):
        ds_region1[i][0:1240]=Tb_region1[i][0:1240];ds_cregion1[i][0:1240]=Tb_cregion1[i][0:1240];ds_region1[i][1740:2240]=Tb_region1[i][1240:];ds_cregion1[i][1740:2240]=Tb_cregion1[i][1240:]
        ds_region2[i][0:1240]=Tb_region2[i][0:1240];ds_cregion2[i][0:1240]=Tb_cregion2[i][0:1240];ds_region2[i][1740:2240]=Tb_region2[i][1240:];ds_cregion2[i][1740:2240]=Tb_cregion2[i][1240:]
        ds_region3[i][0:1240]=Tb_region3[i][0:1240];ds_cregion3[i][0:1240]=Tb_cregion3[i][0:1240];ds_region3[i][1740:2240]=Tb_region3[i][1240:];ds_cregion3[i][1740:2240]=Tb_cregion3[i][1240:]
        ds_region4[i][0:1240]=Tb_region4[i][0:1240];ds_cregion4[i][0:1240]=Tb_cregion4[i][0:1240];ds_region4[i][1740:2240]=Tb_region4[i][1240:];ds_cregion4[i][1740:2240]=Tb_cregion4[i][1240:]
        ds_region5[i][0:1240]=Tb_region5[i][0:1240];ds_cregion5[i][0:1240]=Tb_cregion5[i][0:1240];ds_region5[i][1740:2240]=Tb_region5[i][1240:];ds_cregion5[i][1740:2240]=Tb_cregion5[i][1240:]
        ds_region6[i][0:1240]=Tb_region6[i][0:1240];ds_cregion6[i][0:1240]=Tb_cregion6[i][0:1240];ds_region6[i][1740:2240]=Tb_region6[i][1240:];ds_cregion6[i][1740:2240]=Tb_cregion6[i][1240:]
    if(i==5):
        ds_region1[i][0:1740]=Tb_region1[i][0:1740];ds_cregion1[i][0:1740]=Tb_cregion1[i][0:1740];ds_region1[i][2240:]=Tb_region1[i][1740:];ds_cregion1[i][2240:]=Tb_cregion1[i][1740:]
        ds_region2[i][0:1740]=Tb_region2[i][0:1740];ds_cregion2[i][0:1740]=Tb_cregion2[i][0:1740];ds_region2[i][2240:]=Tb_region2[i][1740:];ds_cregion2[i][2240:]=Tb_cregion2[i][1740:]
        ds_region3[i][0:1740]=Tb_region3[i][0:1740];ds_cregion3[i][0:1740]=Tb_cregion3[i][0:1740];ds_region3[i][2240:]=Tb_region3[i][1740:];ds_cregion3[i][2240:]=Tb_cregion3[i][1740:]
        ds_region4[i][0:1740]=Tb_region4[i][0:1740];ds_cregion4[i][0:1740]=Tb_cregion4[i][0:1740];ds_region4[i][2240:]=Tb_region4[i][1740:];ds_cregion4[i][2240:]=Tb_cregion4[i][1740:]
        ds_region5[i][0:1740]=Tb_region5[i][0:1740];ds_cregion5[i][0:1740]=Tb_cregion5[i][0:1740];ds_region5[i][2240:]=Tb_region5[i][1740:];ds_cregion5[i][2240:]=Tb_cregion5[i][1740:]
        ds_region6[i][0:1740]=Tb_region6[i][0:1740];ds_cregion6[i][0:1740]=Tb_cregion6[i][0:1740];ds_region6[i][2240:]=Tb_region6[i][1740:];ds_cregion6[i][2240:]=Tb_cregion6[i][1740:]
    else:
        sh=Tb_region1[i].shape[0]
        ds_region1[i][0:sh]=Tb_region1[i];ds_cregion1[i][0:sh]=Tb_cregion1[i]
        ds_region2[i][0:sh]=Tb_region2[i];ds_cregion2[i][0:sh]=Tb_cregion2[i]
        ds_region3[i][0:sh]=Tb_region3[i];ds_cregion3[i][0:sh]=Tb_cregion3[i]
        ds_region4[i][0:sh]=Tb_region4[i];ds_cregion4[i][0:sh]=Tb_cregion4[i]
        ds_region5[i][0:sh]=Tb_region5[i];ds_cregion5[i][0:sh]=Tb_cregion5[i]
        ds_region6[i][0:sh]=Tb_region6[i];ds_cregion6[i][0:sh]=Tb_cregion6[i]

#Tb_region=np.array(Tb_region);Tb_region_flat=Tb_region.reshape(8,625,3552);Tb_region_mean=np.nanmean(Tb_region,axis=3);Tb_region_std=np.nanstd(Tb_region,axis=3)
xcr_idx=xx[1:,1:].reshape(2500)-2;ycr_idx=yy[1:,1:].reshape(2500)-2;xc_region=50*(xcr_idx-100);yc_region=50*(ycr_idx-100);xxr=50*(xx[1:,1:]-100);yyr=50*(yy[1:,1:]-100)
#Tb_region1=Tb_region[:,15,17,:];Tb_region2=Tb_region[:,13,6,:];Tb_region3=Tb_region[:,11,9,:];Tb_region4=Tb_region[:,9,17,:];Tb_region5=Tb_region[:,11,14,:];Tb_region6=Tb_region[:,4,8,:];Tb_region7=Tb_region[:,13,6,:];Tb_region8=Tb_region[:,2,5,:]

f,ax=plt.subplots(2,3,sharex=True,sharey=True)
ax[0,0].imshow(np.array(ds_cregion1),aspect='auto',origin=0,vmin=2000,vmax=4000,extent=[0,2740*0.5,0,7],cmap='YlOrRd')
ax[0,1].imshow(np.array(ds_cregion2),aspect='auto',origin=0,vmin=2000,vmax=4000,extent=[0,2740*0.5,0,7],cmap='YlOrRd')
ax[0,2].imshow(np.array(ds_cregion3),aspect='auto',origin=0,vmin=2000,vmax=4000,extent=[0,2740*0.5,0,7],cmap='YlOrRd')
ax[1,0].imshow(np.array(ds_cregion4),aspect='auto',origin=0,vmin=2000,vmax=4000,extent=[0,2740*0.5,0,7],cmap='YlOrRd')
ax[1,1].imshow(np.array(ds_cregion5),aspect='auto',origin=0,vmin=2000,vmax=4000,extent=[0,2740*0.5,0,7],cmap='YlOrRd')
im=ax[1,2].imshow(np.array(ds_cregion6),aspect='auto',origin=0,vmin=2000,vmax=4000,extent=[0,2740*0.5,0,7],cmap='YlOrRd')
ax[0,0].text(1000, 5, 'Region 1', horizontalalignment='center', verticalalignment='center',fontsize=15)
ax[0,1].text(1000, 5, 'Region 2', horizontalalignment='center', verticalalignment='center',fontsize=15)
ax[0,2].text(1000, 5, 'Region 3', horizontalalignment='center', verticalalignment='center',fontsize=15)
ax[1,0].text(1000, 5, 'Region 4', horizontalalignment='center', verticalalignment='center',fontsize=15)
ax[1,1].text(1000, 5, 'Region 5', horizontalalignment='center', verticalalignment='center',fontsize=15)
ax[1,2].text(1000, 5, 'Region 6', horizontalalignment='center', verticalalignment='center',fontsize=15)
ax[0,0].set_yticks(np.linspace(0,6,8)+0.5);ax[0,0].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[1,0].set_yticks(np.linspace(0,6,8)+0.5);ax[1,0].set_yticklabels(['108','120','133','145','179','197','217','240'])
ax[1,0].set_xticks(np.linspace(180,2500,4)*0.5);ax[1,1].set_xticks(np.linspace(180,2500,4)*0.5);ax[1,2].set_xticks(np.linspace(180,2500,4)*0.5)
ax[1,1].set_xlabel('Time (sec)');ax[0,0].set_ylabel('Frequency (MHz)')
f.colorbar(im,label='$T_B$ (MK)')
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()




from matplotlib import patches
k=6
f,ax=plt.subplots(1,1)
y=np.nanmean(Tb197[0:100],axis=0)/1.e6;y[y==0]=np.nan
im=ax.imshow(y,extent=[-2500,2500,-2500,2500],origin=0,vmin=0.01,vmax=0.4)
e1=patches.Ellipse((xxr[26,28],yyr[26,28]),240,240,edgecolor='y',fill=False,linewidth=5);ax.text(xxr[26,28]-50,yyr[26,28]-50,'1')
e2=patches.Ellipse((xxr[24,18],yyr[24,18]),240,240,edgecolor='y',fill=False,linewidth=5);ax.text(xxr[24,18]-50,yyr[24,18]-50,'2')
e3=patches.Ellipse((xxr[26,23],yyr[26,23]),240,240,edgecolor='y',fill=False,linewidth=5);ax.text(xxr[26,23]-50,yyr[26,23]-50,'3')
e4=patches.Ellipse((xxr[21,24],yyr[21,24]),240,240,edgecolor='y',fill=False,linewidth=5);ax.text(xxr[21,24]-50,yyr[21,24]-50,'4')
e5=patches.Ellipse((xxr[20,16],yyr[20,16]),240,240,edgecolor='y',fill=False,linewidth=5);ax.text(xxr[20,16]-50,yyr[20,16]-50,'5')
e6=patches.Ellipse((xxr[16,24],yyr[16,24]),240,240,edgecolor='y',fill=False,linewidth=5);ax.text(xxr[16,24]-50,yyr[16,24]-50,'6')
e7=patches.Ellipse((0,0),1920,1920,linestyle='--',edgecolor='k',fill=False,linewidth=2)
#e7=patches.Ellipse((xxr[2,10],yyr[2,10]),240,240,edgecolor='y',fill=False,linewidth=5);ax.text(xxr[2,10]-50,yyr[2,10]-50,'7')
#e8=patches.Ellipse((xxr[2,5],yyr[2,5]),240,240,edgecolor='y',fill=False,linewidth=5);ax.text(xxr[2,5]-50,yyr[2,5]-50,'8')
ax.add_patch(e1);ax.add_patch(e2);ax.add_patch(e3);ax.add_patch(e4);ax.add_patch(e5);ax.add_patch(e6);ax.add_patch(e7)
ax.set_xlabel('Solar X (arcsec)');ax.set_ylabel('Solar Y (arcsec)')
f.colorbar(im,label='$T_B$ (MK)')
plt.show()


l=np.arange(1024)-512;lmesh=np.meshgrid(l,l);rmesh=np.sqrt(lmesh[0]*lmesh[0]+lmesh[1]*lmesh[1]);dc[np.where(rmesh<187)]=0;dc=dc.astype('float16');dc[np.where(dc==0)]=np.nan
levels_sub=[0.3,0.4,0.5,0.6,0.7,0.8,0.9];levels_sub1=[0.6,0.7,0.8,0.9];levels_sub2=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
levels_all=[levels_sub1,levels_sub1,levels_sub1,levels_sub,levels_sub,levels_sub1,levels_sub,levels_sub2]
Tb_sub_cor_mean_all=[Tb_sub108_cor_mean,Tb_sub120_cor_mean,Tb_sub133_cor_mean,Tb_sub145_cor_mean,Tb_sub179_cor_mean,Tb_sub197_cor_mean,Tb_sub217_cor_mean,Tb_sub240_cor_mean]
c=['tab:olive','tab:pink','tab:purple','tab:red','tab:green','tab:cyan','tab:orange','tab:blue']
ylab=['108 MHz','120 MHz','133 MHz','145 MHz','179 MHz','197 MHz','217 MHz','240 MHz']
for i in range(len(Tb_sub_cor_mean_all)):
    ii="%02d"%i 
    f,ax=plt.subplots(1,1,figsize=(10,10))
    ax.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
    ax.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
    ax.contour(Tb_sub_cor_mean_all[i]/np.nanmax(Tb_sub_cor_mean_all[i]),levels_all[i],extent=[-5000,5000,-5000,5000], colors=c[i],linewidths=2,alpha=0.9);ax.text(-3000,-1800-200*i,ylab[i],fontweight='bold',fontsize=15,color=c[i])
    plt.xlim(-3500,3500);plt.ylim(-3500,3500);ax.set_xlabel('Solar X (arcsec)');ax.set_ylabel('Solar Y (arcsec)')
    plt.savefig('/nas08-data02/rohit/20151203_sub_run_mean_subamp/pngs_contours/cont_'+str(ii)+'.png')
    plt.close()



plt.plot(np.arange(3552)*0.5,Tb_sub240[:,26,23],'o-',label='Before Correction')
plt.plot(np.arange(3552)*0.5,Tb_sub240_cor[:,26,23],'o-',label='After Correction')
plt.ylim(-2500,2500);plt.ylabel('T$_{B}$ (K)');plt.xlabel('Time (sec)');plt.legend()
plt.show()

plt.hist(Tb_sub240[:,26,23],bins=400,label='Before Correction',linewidth=4,histtype='step')
plt.hist(Tb_sub240_cor[:,26,23],bins=400,label='After Correction',linewidth=4,histtype='step')
plt.legend();plt.xlabel('T$_{B}$ (K)');plt.xlim(-3000,3000);plt.ylabel('Number')
plt.show()

#l=np.arange(1024)-512;lmesh=np.meshgrid(l,l);rmesh=np.sqrt(lmesh[0]*lmesh[0]+lmesh[1]*lmesh[1]);dc[np.where(rmesh<187)]=0;dc.astype(np.float16)[np.where(dc==0)]=np.nan
l=np.arange(1024)-512;lmesh=np.meshgrid(l,l);rmesh=np.sqrt(lmesh[0]*lmesh[0]+lmesh[1]*lmesh[1]);dc[np.where(rmesh<187)]=0;dc=dc.astype(np.float16);dc[np.where(dc==0)]=np.nan
levels_sub=[0.3,0.4,0.5,0.6,0.7,0.8,0.9];levels_sub1=[0.6,0.7,0.8,0.9];levels_sub2=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
f,ax=plt.subplots(1,1,figsize=(10,10))
ax.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
ax.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
ax.contour(Tb_sub240_cor_mean/np.nanmax(Tb_sub240_cor_mean),levels_sub2,extent=[-5000,5000,-5000,5000],  colors='tab:blue',linewidths=4,alpha=0.7);ax.text(-3000,-3200,'240 MHz',fontweight='bold',fontsize=15,color='tab:blue')
ax.contour(Tb_sub217_cor_mean/np.nanmax(Tb_sub217_cor_mean),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:orange',linewidths=4,alpha=0.7);ax.text(-3000,-3000,'217 MHz',fontweight='bold',fontsize=15,color='tab:orange')
ax.contour(Tb_sub197_cor_mean/np.nanmax(Tb_sub197_cor_mean),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:cyan',linewidths=4,alpha=0.8);ax.text(-3000,-2800,'197 MHz',fontweight='bold',fontsize=15,color='tab:cyan')
ax.contour(Tb_sub179_cor_mean/np.nanmax(Tb_sub179_cor_mean),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:green',linewidths=4);ax.text(-3000,-2600,'179 MHz',fontweight='bold',fontsize=15,color='tab:green')
ax.contour(Tb_sub145_cor_mean/np.nanmax(Tb_sub145_cor_mean),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:red',linewidths=4,alpha=0.8);ax.text(-3000,-2400,'145 MHz',fontweight='bold',fontsize=15,color='tab:red')
ax.contour(Tb_sub133_cor_mean/np.nanmax(Tb_sub133_cor_mean),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:purple',linewidths=4);ax.text(-3000,-2200,'133 MHz',fontweight='bold',fontsize=15,color='tab:purple')
ax.contour(Tb_sub120_cor_mean/np.nanmax(Tb_sub120_cor_mean),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:pink',linewidths=4);ax.text(-3000,-2000,'120 MHz',fontweight='bold',fontsize=15,color='tab:pink')
ax.contour(Tb_sub108_cor_mean/np.nanmax(Tb_sub108_cor_mean),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:olive',linewidths=4);ax.text(-3000,-1800,'108 MHz',fontweight='bold',fontsize=15,color='tab:olive')
plt.xlim(-3500,3500);plt.ylim(-3500,3500);ax.set_xlabel('Solar X (arcsec)');ax.set_ylabel('Solar Y (arcsec)')
plt.show()


l=np.arange(1024)-512;lmesh=np.meshgrid(l,l);rmesh=np.sqrt(lmesh[0]*lmesh[0]+lmesh[1]*lmesh[1]);dc[np.where(rmesh<187)]=0;dc=dc.astype(np.float16);dc[np.where(dc==0)]=np.nan
levels_sub=[0.3,0.4,0.5,0.6,0.7,0.8,0.9];levels_sub1=[0.6,0.7,0.8,0.9];levels_sub2=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
for i in range(3552):
    f,ax=plt.subplots(1,1,figsize=(10,10))
    ax.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
    ax.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
    ii="%02d"%i
    ax.contour(Tb_sub240_cor[i],[3000,4000],extent=[-5000,5000,-5000,5000],colors='tab:blue',linewidths=4);ax.text(-3000,-3200,'240 MHz',fontweight='bold',fontsize=15,color='tab:blue')
    ax.contour(Tb_sub217_cor[i],[3500,4500],extent=[-5000,5000,-5000,5000],colors='tab:orange',linewidths=4);ax.text(-3000,-3000,'217 MHz',fontweight='bold',fontsize=15,color='tab:orange')
    ax.contour(Tb_sub197_cor[i],[3500,4500],extent=[-5000,5000,-5000,5000],colors='tab:cyan',linewidths=2);ax.text(-3000,-2800,'197 MHz',fontweight='bold',fontsize=15,color='tab:cyan')
    ax.contour(Tb_sub179_cor[i],[3500,4500],extent=[-5000,5000,-5000,5000],  colors='tab:green',linewidths=2);ax.text(-3000,-2600,'179 MHz',fontweight='bold',fontsize=15,color='tab:green')
    ax.contour(Tb_sub145_cor[i],[4000,9000],extent=[-5000,5000,-5000,5000],  colors='tab:red',linewidths=2,alpha=0.8);ax.text(-3000,-2400,'145 MHz',fontweight='bold',fontsize=15,color='tab:red')
    ax.contour(Tb_sub133_cor[i],[5000,9000],extent=[-5000,5000,-5000,5000],  colors='tab:purple',linewidths=2);ax.text(-3000,-2200,'133 MHz',fontweight='bold',fontsize=15,color='tab:purple')
    ax.contour(Tb_sub120_cor[i],[8000,9000],extent=[-5000,5000,-5000,5000],  colors='tab:pink',linewidths=2);ax.text(-3000,-2000,'120 MHz',fontweight='bold',fontsize=15,color='tab:pink')
    ax.contour(Tb_sub108_cor[i],[8000,9000],extent=[-5000,5000,-5000,5000],  colors='tab:olive',linewidths=2);ax.text(-3000,-1800,'108 MHz',fontweight='bold',fontsize=15,color='tab:olive')
    plt.xlim(-4000,4000);plt.ylim(-4000,4000)
    plt.savefig('/nas08-data02/rohit/20151203_sub_run_mean_subamp/pngs_all_freq/cont_'+str(ii)+'.png')
    plt.close()

l=np.arange(1024)-512;lmesh=np.meshgrid(l,l);rmesh=np.sqrt(lmesh[0]*lmesh[0]+lmesh[1]*lmesh[1]);dc[np.where(rmesh<187)]=0;dc=dc.astype(np.float16);dc[np.where(dc==0)]=np.nan
levels_sub=[0.3,0.4,0.5,0.6,0.7,0.8,0.9];levels_sub1=[0.6,0.7,0.8,0.9];levels_sub2=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
for i in range(350,3552):
    f,ax=plt.subplots(1,1,figsize=(10,10))
    ax.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='binary')
    ax.imshow(dc.astype(np.float)[::-1,::-1],origin=0,extent=[-6093,6093,-6093,6093],vmin=3.2e3,vmax=4.5e3,cmap='binary')
    ii="%04d"%i
    #ax.contour(Tb_sub240_cor[i],[np.nanstd(Tb_sub240_cor[i])*5.0,np.nanstd(Tb_sub240_cor[i])*10.0],extent=[-5000,5000,-5000,5000],colors='tab:blue',linewidths=4);ax.text(-3000,-3200,'240 MHz',fontweight='bold',fontsize=15,color='tab:blue')
    #ax.contour(Tb_sub217_cor[i],[np.nanstd(Tb_sub240_cor[i])*5.0,np.nanstd(Tb_sub240_cor[i])*10.0],extent=[-5000,5000,-5000,5000],colors='tab:orange',linewidths=4);ax.text(-3000,-3000,'217 MHz',fontweight='bold',fontsize=15,color='tab:orange')
    #ax.contour(Tb_sub197_cor[i],[np.nanstd(Tb_sub240_cor[i])*5.0,np.nanstd(Tb_sub240_cor[i])*10.0],[3500,4500],extent=[-5000,5000,-5000,5000],colors='tab:cyan',linewidths=2);ax.text(-3000,-2800,'197 MHz',fontweight='bold',fontsize=15,color='tab:cyan')
    #ax.contour(Tb_sub179_cor[i],[np.nanstd(Tb_sub240_cor[i])*5.0,np.nanstd(Tb_sub240_cor[i])*10.0],[3500,4500],extent=[-5000,5000,-5000,5000],  colors='tab:green',linewidths=2);ax.text(-3000,-2600,'179 MHz',fontweight='bold',fontsize=15,color='tab:green')
    #ax.contour(Tb_sub145_cor[i],[np.nanstd(Tb_sub240_cor[i])*5.0,np.nanstd(Tb_sub240_cor[i])*10.0],[4000,9000],extent=[-5000,5000,-5000,5000],  colors='tab:red',linewidths=2,alpha=0.8);ax.text(-3000,-2400,'145 MHz',fontweight='bold',fontsize=15,color='tab:red')
    #ax.contour(Tb_sub133_cor[i],[np.nanstd(Tb_sub240_cor[i])*5.0,np.nanstd(Tb_sub240_cor[i])*10.0],[5000,9000],extent=[-5000,5000,-5000,5000],  colors='tab:purple',linewidths=2);ax.text(-3000,-2200,'133 MHz',fontweight='bold',fontsize=15,color='tab:purple')
    #ax.contour(Tb_sub120_cor[i],[np.nanstd(Tb_sub240_cor[i])*5.0,np.nanstd(Tb_sub240_cor[i])*10.0],[8000,9000],extent=[-5000,5000,-5000,5000],  colors='tab:pink',linewidths=2);ax.text(-3000,-2000,'120 MHz',fontweight='bold',fontsize=15,color='tab:pink')
    ax.contour(Tb_sub108_cor[i],[np.nanstd(Tb_sub108_cor[i])*5.0,np.nanstd(Tb_sub108_cor[i])*10.0],extent=[-5000,5000,-5000,5000],  colors='tab:blue',linewidths=3);ax.text(-3000,-1800,'108 MHz',fontweight='bold',fontsize=15,color='tab:olive')
    ax.contour(Tb_sub108_cor[i],[-1*np.nanstd(Tb_sub108_cor[i])*5.0,-1*np.nanstd(Tb_sub108_cor[i])*10.0],extent=[-5000,5000,-5000,5000],  colors='tab:red',linewidths=3);ax.text(-3000,-1800,'108 MHz',fontweight='bold',fontsize=15,color='tab:olive')
    ax.contour(Tb108[i],[np.nanstd(Tb108[i])*5.0,np.nanstd(Tb108[i])*10.0],extent=[-2500,2500,-2500,2500],  colors='tab:green',linewidths=3)
    plt.xlim(-6000,6000);plt.ylim(-6000,6000)
    plt.savefig('/nas08-data02/rohit/20151203_sub_run_mean_subamp/pngs_all_freq/sub_cont_108MHz_'+str(ii)+'.png')
    plt.close()

levels_sub=[0.3,0.4,0.5,0.6,0.7,0.8,0.9];levels_sub1=[0.6,0.7,0.8,0.9];levels_sub2=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
f,ax=plt.subplots(1,1,figsize=(10,10))
#im=ax.imshow(np.log10(temp),origin=0,extent=[-2500,2500,-2500,2500],vmin=6.0,vmax=6.2,cmap='YlGn')
im=ax.imshow(babs,origin=0,extent=[-2500,2500,-2500,2500],vmin=0,vmax=5,cmap='YlGn')
#f.colorbar(im,label='log($T_B$ (K))')
f.colorbar(im,label='B (G)')
ax.contour(Tb_sub240_cor_mean/np.nanmax(Tb_sub240_cor_mean),levels_sub2,extent=[-5000,5000,-5000,5000],  colors='tab:blue',linewidths=2,alpha=0.7);ax.text(-3000,-3200,'240 MHz',fontweight='bold',fontsize=15,color='tab:blue')
ax.contour(Tb_sub217_cor_mean/np.nanmax(Tb_sub217_cor_mean),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:orange',linewidths=2,alpha=0.7);ax.text(-3000,-3000,'217 MHz',fontweight='bold',fontsize=15,color='tab:orange')
ax.contour(Tb_sub197_cor_mean/np.nanmax(Tb_sub197_cor_mean),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:cyan',linewidths=2,alpha=0.8);ax.text(-3000,-2800,'197 MHz',fontweight='bold',fontsize=15,color='tab:cyan')
ax.contour(Tb_sub179_cor_mean/np.nanmax(Tb_sub179_cor_mean),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:green',linewidths=2);ax.text(-3000,-2600,'179 MHz',fontweight='bold',fontsize=15,color='tab:green')
ax.contour(Tb_sub145_cor_mean/np.nanmax(Tb_sub145_cor_mean),levels_sub,extent=[-5000,5000,-5000,5000],  colors='tab:red',linewidths=2,alpha=0.8);ax.text(-3000,-2400,'145 MHz',fontweight='bold',fontsize=15,color='tab:red')
ax.contour(Tb_sub133_cor_mean/np.nanmax(Tb_sub133_cor_mean),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:purple',linewidths=2);ax.text(-3000,-2200,'133 MHz',fontweight='bold',fontsize=15,color='tab:purple')
ax.contour(Tb_sub120_cor_mean/np.nanmax(Tb_sub120_cor_mean),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:pink',linewidths=2);ax.text(-3000,-2000,'120 MHz',fontweight='bold',fontsize=15,color='tab:pink')
ax.contour(Tb_sub108_cor_mean/np.nanmax(Tb_sub108_cor_mean),levels_sub1,extent=[-5000,5000,-5000,5000],  colors='tab:olive',linewidths=2);ax.text(-3000,-1800,'108 MHz',fontweight='bold',fontsize=15,color='tab:olive')
plt.xlim(-3500,3500);plt.ylim(-3500,3500);ax.set_xlabel('Solar X (arcsec)');ax.set_ylabel('Solar Y (arcsec)')
plt.show()

k=2
lTb=-3000;hTb=20000
t=np.arange(len(Tb_region1[k]))*0.5#;Tb_region1[k][(Tb_region1[k]<lTb) | (Tb_region1[k]>hTb)]=0;Tb_region2[k][(Tb_region2[k]<lTb) | (Tb_region2[k]>hTb)]=0
t1=np.arange(1740)*0.5
#Tb_region3[k][(Tb_region3[k]<lTb) | (Tb_region3[k]>hTb)]=0;Tb_region4[k][(Tb_region4[k]<lTb) | (Tb_region4[k]>hTb)]=0;Tb_region5[k][(Tb_region5[k]<lTb) | (Tb_region5[k]>hTb)]=0
#Tb_region6[k][(Tb_region6[k]<lTb) | (Tb_region6[k]>hTb)]=0;Tb_region7[k][(Tb_region7[k]<lTb) | (Tb_region7[k]>hTb)]=0;Tb_region8[k][(Tb_region8[k]<lTb) | (Tb_region8[k]>hTb)]=0
f,ax=plt.subplots(5,1,figsize=(20,10),sharex=True)
ax[0].plot(t,Tb_region1[7],'-',label='240 MHz');ax[0].set_ylim(-8.e3,7.0e4);ax[0].legend(loc=2,fontsize=10)
ax[0].plot(t,Tb_region1[6],'-',label='217 MHz');ax[0].set_ylim(-8.e3,7.0e4);ax[0].legend(loc=2,fontsize=10)
ax[1].plot(t,Tb_region2[7],'-',label='240 MHz');ax[1].set_ylim(-8.e3,3.0e4);ax[1].legend(loc=2,fontsize=10)
ax[1].plot(t,Tb_region2[6],'-',label='217 MHz');ax[1].set_ylim(-8.e3,3.0e4);ax[1].legend(loc=2,fontsize=10)
ax[2].plot(t,Tb_region3[2],'-',label='133 MHz');ax[2].set_ylim(-8.e3,1.e5);ax[2].legend(loc=2,fontsize=10)
ax[2].plot(t,Tb_region3[1],'-',label='120 MHz');ax[2].set_ylim(-8.e3,1.e5);ax[2].legend(loc=2,fontsize=10)
ax[2].plot(t1,Tb_region3[0],'-',label='108 MHz');ax[2].set_ylim(-8.e3,1.e5);ax[2].legend(loc=2,fontsize=10)
ax[3].plot(t1,Tb_region5[0],'-',label='108 MHz');ax[3].set_ylim(-5.e3,1.e4);ax[3].legend(loc=2,fontsize=10)
ax[3].plot(t,Tb_region5[1],'-',label='120 MHz');ax[3].set_ylim(-5.e3,1.e4);ax[3].legend(loc=2,fontsize=10)
ax[4].plot(t1,Tb_region6[0],'-',label='108 MHz');ax[4].set_ylim(-5.e3,1.e4);ax[4].legend(loc=2,fontsize=10)
ax[4].plot(t,Tb_region6[1],'-',label='120 MHz');ax[4].set_ylim(-5.e3,1.e4);ax[4].legend(loc=2,fontsize=10)
ax[0].annotate("Region 1",xy=(1300,1.e3),fontsize=15)
ax[1].annotate("Region 2",xy=(1300,1.e3),fontsize=15)
ax[2].annotate("Region 3",xy=(1300,1.e3),fontsize=15)
ax[3].annotate("Region 5",xy=(1300,1.e3),fontsize=15)
ax[4].annotate("Region 6",xy=(1300,1.e3),fontsize=15)
ax[3].set_ylabel('$T_B$ (K)');ax[4].set_xlabel('Time (sec)')
plt.show()

k=2
lTb=-3000;hTb=20000
t=np.arange(len(Tb_region1[k]))*0.5;Tb_region1[k][(Tb_region1[k]<lTb) | (Tb_region1[k]>hTb)]=0;Tb_region2[k][(Tb_region2[k]<lTb) | (Tb_region2[k]>hTb)]=0
#Tb_region3[k][(Tb_region3[k]<lTb) | (Tb_region3[k]>hTb)]=0;Tb_region4[k][(Tb_region4[k]<lTb) | (Tb_region4[k]>hTb)]=0;Tb_region5[k][(Tb_region5[k]<lTb) | (Tb_region5[k]>hTb)]=0
#Tb_region6[k][(Tb_region6[k]<lTb) | (Tb_region6[k]>hTb)]=0;Tb_region7[k][(Tb_region7[k]<lTb) | (Tb_region7[k]>hTb)]=0;Tb_region8[k][(Tb_region8[k]<lTb) | (Tb_region8[k]>hTb)]=0
f,ax=plt.subplots(6,1,figsize=(10,10),sharex=True)
ax[0].plot(t,Tb_region1[k],'-',label='Region 1');ax[0].set_ylim(-3.e3,1.0e4);ax[0].legend(loc=2,fontsize=10)
ax[1].plot(t,Tb_region2[k],'-',label='Region 2');ax[1].set_ylim(-3.e3,1.0e4);ax[1].legend(loc=2,fontsize=10)
ax[2].plot(t,Tb_region3[k],'-',label='Region 3');ax[2].set_ylim(-5.e3,1.e5);ax[2].legend(loc=2,fontsize=10)
ax[3].plot(t,Tb_region4[k],'-',label='Region 4');ax[3].set_ylim(-5.e3,1.e4);ax[3].legend(loc=2,fontsize=10)
ax[4].plot(t,Tb_region5[k],'-',label='Region 5');ax[4].set_ylim(-5.e3,1.e4);ax[4].legend(loc=2,fontsize=10)
ax[5].plot(t,Tb_region6[k],'-',label='Region 6');ax[5].set_ylim(-5.e3,1.e4);ax[5].legend(loc=2,fontsize=10)
ax[3].set_ylabel('$T_B$ (K)');ax[5].set_xlabel('Time (sec)')
plt.show()

k=7
lTb=-3000;hTb=20000
t=np.arange(len(Tb_region1[k]))*0.5;Tb_region1[k][(Tb_region1[k]<lTb) | (Tb_region1[k]>hTb)]=0;Tb_region2[k][(Tb_region2[k]<lTb) | (Tb_region2[k]>hTb)]=0
Tb_region3[k][(Tb_region3[k]<lTb) | (Tb_region3[k]>hTb)]=0;Tb_region4[k][(Tb_region4[k]<lTb) | (Tb_region4[k]>hTb)]=0;Tb_region5[k][(Tb_region5[k]<lTb) | (Tb_region5[k]>hTb)]=0
Tb_region6[k][(Tb_region6[k]<lTb) | (Tb_region6[k]>hTb)]=0;Tb_region7[k][(Tb_region7[k]<lTb) | (Tb_region7[k]>hTb)]=0;Tb_region8[k][(Tb_region8[k]<lTb) | (Tb_region8[k]>hTb)]=0
f,ax=plt.subplots(8,1,figsize=(10,10),sharex=True)
ax[0].plot(t,Tb_region1[k],'-',label='Region 1');ax[0].set_ylim(-3.e3,1.2e4);ax[0].legend(loc=2,fontsize=10)
ax[1].plot(t,Tb_region2[k],'-',label='Region 2');ax[1].set_ylim(-3.e3,1.2e4);ax[1].legend(loc=2,fontsize=10)
ax[2].plot(t,Tb_region3[k],'-',label='Region 3');ax[2].set_ylim(-3.e3,1.2e4);ax[2].legend(loc=2,fontsize=10)
ax[3].plot(t,Tb_region4[k],'-',label='Region 4');ax[3].set_ylim(-2.e3,2.e3);ax[3].legend(loc=2,fontsize=10)
ax[4].plot(t,Tb_region5[k],'-',label='Region 5');ax[4].set_ylim(-2.e3,2.e3);ax[4].legend(loc=2,fontsize=10)
ax[5].plot(t,Tb_region6[k],'-',label='Region 6');ax[5].set_ylim(-2.e3,2.e3);ax[5].legend(loc=2,fontsize=10)
ax[6].plot(t,Tb_region7[k],'-',label='Region 7');ax[6].set_ylim(-1.5e3,1.5e3);ax[6].legend(loc=2,fontsize=10)
ax[7].plot(t,Tb_region8[k],'-',label='Region 8');ax[7].set_ylim(-1.5e3,1.5e3);ax[7].legend(loc=2,fontsize=10)
ax[4].set_ylabel('$T_B$ (K)');ax[7].set_xlabel('Time (sec)')
plt.show()

k=0;bins_=10
plt.hist(Tb_region1[k],bins=bins_,label='Region 1',histtype='step',linewidth=4)
plt.hist(Tb_region2[k],bins=bins_,label='Region 2',histtype='step',linewidth=4)
plt.hist(Tb_region5[k],bins=bins_,label='Region 5',histtype='step',linewidth=4)
plt.hist(Tb_region6[k],bins=bins_,label='Region 6',histtype='step',linewidth=4)
plt.legend();plt.xlabel('$T_B$ (K)');plt.ylabel('Numbers')
plt.show()

k=0;bins_=15
plt.hist(Tb_region1[7],bins=bins_,label='Region 1, 240 MHz',normed=1,histtype='step',linewidth=4)
plt.hist(Tb_region2[6],bins=bins_,label='Region 2, 217 MHz',normed=1,histtype='step',linewidth=4)
plt.hist(Tb_region3[0],bins=bins_,label='Region 3, 108 MHz',normed=1,histtype='step',linewidth=4)
plt.hist(Tb_region5[0],bins=bins_,label='Region 5, 108 MHz',normed=1,histtype='step',linewidth=4)
plt.hist(Tb_region6[1],bins=bins_,label='Region 6, 120 MHz',normed=1,histtype='step',linewidth=4)
plt.legend();plt.xlabel('$T_B$ (K)');plt.ylabel('Numbers')
plt.show()

rate1=np.zeros(8);rate2=np.zeros(8);rate3=np.zeros(8);rate4=np.zeros(8);rate5=np.zeros(8);rate6=np.zeros(8);err=np.zeros(8)
for i in range(8):
    rate1[i]=np.sum(Tb_cregion1[i])/len(np.where(Tb_cregion1[i]!=0)[0])
    rate2[i]=np.sum(Tb_cregion2[i])/len(np.where(Tb_cregion2[i]!=0)[0])
    rate3[i]=np.sum(Tb_cregion3[i])/len(np.where(Tb_cregion3[i]!=0)[0])
    rate4[i]=np.sum(Tb_cregion4[i])/len(np.where(Tb_cregion4[i]!=0)[0])
    rate5[i]=np.sum(Tb_cregion5[i])/len(np.where(Tb_cregion5[i]!=0)[0])
    rate6[i]=np.sum(Tb_cregion6[i])/len(np.where(Tb_cregion6[i]!=0)[0])
    err[i]=np.std(Tb_region4[i])

rateall=np.array((rate1,rate2,rate3,rate4,rate5,rate6))
freq=np.array([108.0, 120.0, 133.0, 145.0, 179.0, 197.0, 217.0, 240.0])
c=['b','g','r','cyan','magenta','orange']
for i in range(6):
    plt.plot(freq,rateall[i],'o',color=c[i])
    plt.errorbar(freq,rateall[i],yerr=err[i],color=c[i],label='Region '+str(i+1))
plt.legend();plt.xlabel('Frequency (MHz)');plt.ylabel('$T_B$ (K)')
plt.show()


for i in range(8):
    print freq[i],' MHz &', np.round(rate1[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region1[i][0:500])/1.e3,2),'&', np.round(rate2[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region2[i][0:500])/1.e3,2),'&', \
            np.round(rate3[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region3[i][0:500])/1.e3,2),'&',np.round(rate4[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region4[i][0:500])/1.e3,2), '&', \
            np.round(rate5[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region5[i][0:500])/1.e3,2), '&', np.round(rate6[i]/1.e3,2),' $\pm$',np.round(np.std(Tb_region6[i][0:500])/1.e3,2),'\\\\'
for i in range(8):
    print freq[i],' MHz &', np.round(rate1[i]/1.e3,2),' $\pm$',np.round((err[i])/1.e3,2),'&', np.round(rate2[i]/1.e3,2),' $\pm$',np.round((err[i])/1.e3,2),'&', \
            np.round(rate3[i]/1.e3,2),' $\pm$',np.round((err[i])/1.e3,2),'&',np.round(rate4[i]/1.e3,2),' $\pm$',np.round((err[i])/1.e3,2), '&', \
            np.round(rate5[i]/1.e3,2),' $\pm$',np.round((err[i])/1.e3,2), '&', np.round(rate6[i]/1.e3,2),' $\pm$',np.round((err[i])/1.e3,2),'\\\\'


plt.plot(Tb_region1[k],'-')
plt.plot(Tb_region2[k],'-')
plt.plot(Tb_region3[k],'-')
plt.plot(Tb_region4[k],'-')
plt.plot(Tb_region5[k],'-')
plt.plot(Tb_region6[k],'-')
plt.plot(Tb_region7[k],'-')
plt.show()

plt.imshow(Tb_sub197.mean(axis=0),origin=0)
plt.show()


levels_sub=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax01.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax00.contour(np.nanmean(Tb108[650:1500],axis=0)/np.nanmax(np.nanmean(Tb108[650:1500],axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14)
ax01.contour(np.nanmean(Tb120,axis=0)/np.nanmax(np.nanmean(Tb120,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax01.annotate('(b) 120 MHz',xy=(-50,2000),fontsize=14)
ax02.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax02.contour(np.nanmean(Tb133,axis=0)/np.nanmax(np.nanmean(Tb133,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14)
ax03.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax03.contour(np.nanmean(Tb145,axis=0)/np.nanmax(np.nanmean(Tb145,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14)
ax10.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax10.contour(np.nanmean(Tb179,axis=0)/np.nanmax(np.nanmean(Tb179,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14)
ax11.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax11.contour(np.nanmean(Tb197[0:500],axis=0)/np.nanmax(np.nanmean(Tb197[0:500],axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14)
ax12.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax12.contour(np.nanmean(Tb217,axis=0)/np.nanmax(np.nanmean(Tb217,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14)
ax13.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax13.contour(np.nanmean(Tb240,axis=0)/np.nanmax(np.nanmean(Tb240[10:100],axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14)
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
#f.tight_layout()
#f.suptitle('$T_{B} >$ '+str(Tblim)+' K', fontsize="large")
plt.show()

f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax01.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax00.contour(np.nanmean(Tb_sub108[650:1500],axis=0)/np.nanmax(np.nanmean(Tb_sub108[650:1500],axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14)
ax01.contour(np.nanmean(Tb_sub120,axis=0)/np.nanmax(np.nanmean(Tb_sub120,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax01.annotate('(b) 120 MHz',xy=(-50,2000),fontsize=14)
ax02.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax02.contour(np.nanmean(Tb_sub133,axis=0)/np.nanmax(np.nanmean(Tb_sub133,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14)
ax03.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax03.contour(np.nanmean(Tb_sub145,axis=0)/np.nanmax(np.nanmean(Tb_sub145,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14)
ax10.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax10.contour(np.nanmean(Tb_sub179,axis=0)/np.nanmax(np.nanmean(Tb_sub179,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14)
ax11.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax11.contour(np.nanmean(Tb_sub197,axis=0)/np.nanmax(np.nanmean(Tb_sub197,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14)
ax12.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax12.contour(np.nanmean(Tb_sub217,axis=0)/np.nanmax(np.nanmean(Tb_sub217,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14)
ax13.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax13.contour(np.nanmean(Tb_sub240,axis=0)/np.nanmax(np.nanmean(Tb_sub240,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14)
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
#f.tight_layout()
#f.suptitle('$T_{B} >$ '+str(Tblim)+' K', fontsize="large")
plt.show()


levels_sub=[0.5,0.6]
f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax01.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax00.contour(np.nanmean(Tb_sub108_cor[650:1500],axis=0)/np.nanmax(np.nanmean(Tb_sub108_cor[650:1500],axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14)
ax01.contour(np.nanmean(Tb_sub120_cor,axis=0)/np.nanmax(np.nanmean(Tb_sub120_cor,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax01.annotate('(b) 120 MHz',xy=(-50,2000),fontsize=14)
ax02.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax02.contour(np.nanmean(Tb_sub133_cor[:2000],axis=0)/np.nanmax(np.nanmean(Tb_sub133_cor[:2000],axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14)
ax03.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax03.contour(np.nanmean(Tb_sub145_cor[:3500],axis=0)/np.nanmax(np.nanmean(Tb_sub145_cor[:3500],axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14)
ax10.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax10.contour(np.nanmean(Tb_sub179_cor,axis=0)/np.nanmax(np.nanmean(Tb_sub179_cor,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14)
ax11.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax11.contour(np.nanmean(Tb_sub197_cor,axis=0)/np.nanmax(np.nanmean(Tb_sub197_cor,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14)
ax12.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax12.contour(np.nanmean(Tb_sub217_cor,axis=0)/np.nanmax(np.nanmean(Tb_sub217_cor,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14)
ax13.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
ax13.contour(np.nanmean(Tb_sub240_cor,axis=0)/np.nanmax(np.nanmean(Tb_sub240_cor,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14)
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
plt.show()

f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.hist(np.nanmean(Tb_sub108_cor_sel,axis=0).flatten(),bins=50,label='108 MHz')
ax01.hist(np.nanmean(Tb_sub120_cor_sel,axis=0).flatten(),bins=50,label='120 MHz')
ax02.hist(np.nanmean(Tb_sub133_cor_sel,axis=0).flatten(),bins=50,label='133 MHz')
ax03.hist(np.nanmean(Tb_sub145_cor_sel,axis=0).flatten(),bins=50,label='145 MHz')
ax10.hist(np.nanmean(Tb_sub179_cor_sel,axis=0).flatten(),bins=50,label='179 MHz')
ax11.hist(np.nanmean(Tb_sub197_cor_sel,axis=0).flatten(),bins=50,label='197 MHz')
ax12.hist(np.nanmean(Tb_sub217_cor_sel,axis=0).flatten(),bins=50,label='217 MHz')
ax13.hist(np.nanmean(Tb_sub240_cor_sel,axis=0).flatten(),bins=50,label='240 MHz')
ax00.set_yscale('log');ax01.set_yscale('log');ax02.set_yscale('log');ax03.set_yscale('log');ax10.set_yscale('log');ax11.set_yscale('log');ax12.set_yscale('log');ax13.set_yscale('log')
ax00.set_ylim(1.e2,5.e4);ax01.set_ylim(1.e2,5.e4);ax02.set_ylim(1.e2,5.e4);ax03.set_ylim(1.e2,5.e4);ax10.set_ylim(1.e2,5.e4);ax11.set_ylim(1.e2,5.e4);ax12.set_ylim(1.e2,5.e4);ax13.set_ylim(1.e2,5.e4)
ax00.legend(fontsize=14);ax01.legend(fontsize=14);ax02.legend(fontsize=14);ax03.legend(fontsize=14);ax10.legend(fontsize=14);ax11.legend(fontsize=14);ax12.legend(fontsize=14);ax13.legend(fontsize=14)
ax10.set_xlabel('$T_B$ (K)');ax11.set_xlabel('$T_B$ (K)');ax12.set_xlabel('$T_B$ (K)');ax13.set_xlabel('$T_B$ (K)')
ax10.set_ylabel('Occupancy');ax00.set_ylabel('Occupancy')
plt.show()

print np.nanmean(Tb_sub108_cor[650:1500],axis=0).sum(),np.nanmean(Tb_sub120_cor[650:1500],axis=0).sum(),np.nanmean(Tb_sub133_cor[:2000],axis=0).sum(),np.nanmean(Tb_sub145_cor[:3500],axis=0).sum(),np.nanmean(Tb_sub179_cor,axis=0).sum(),np.nanmean(Tb_sub197_cor,axis=0).sum(),np.nanmean(Tb_sub217_cor,axis=0).sum(),np.nanmean(Tb_sub240_cor,axis=0).sum()
print np.nanmean(Tb_sub108[650:1500],axis=0).sum(),np.nanmean(Tb_sub120[650:1500],axis=0).sum(),np.nanmean(Tb_sub133[:2000],axis=0).sum(),np.nanmean(Tb_sub145[:3500],axis=0).sum(),np.nanmean(Tb_sub179,axis=0).sum(),np.nanmean(Tb_sub197,axis=0).sum(),np.nanmean(Tb_sub217,axis=0).sum(),np.nanmean(Tb_sub240,axis=0).sum()


from mpl_toolkits.axes_grid1 import make_axes_locatable
f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.imshow(np.nansum(Tb_sub108_cor_sel,axis=0)/1.e3,extent=[-5000,5000,-5000,5000],origin=0,vmin=-2.2e2,vmax=2.2e2,cmap='coolwarm')
ax00.contour(np.nanmean(Tb_sub108_cor_sel,axis=0),extent=[-5000,5000,-5000,5000],levels=[-5*np.nanmean(Tb_sub108_cor_sel,axis=0).std(),5*np.nanmean(Tb_sub108_cor_sel,axis=0).std()])
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14,color='k')
ax01.imshow(np.nansum(Tb_sub120_cor_sel,axis=0)/1.e3,extent=[-5000,5000,-5000,5000],origin=0,vmin=-2.2e2,vmax=2.2e2,cmap='coolwarm')
ax01.contour(np.nanmean(Tb_sub120_cor_sel,axis=0),extent=[-5000,5000,-5000,5000],levels=[-5*np.nanmean(Tb_sub120_cor_sel,axis=0).std(),5*np.nanmean(Tb_sub120_cor_sel,axis=0).std()])
ax01.annotate('(a) 120 MHz',xy=(-50,2000),fontsize=14,color='k')
ax02.imshow(np.nansum(Tb_sub133_cor_sel,axis=0)/1.e3,extent=[-5000,5000,-5000,5000],origin=0,vmin=-2.2e2,vmax=2.2e2,cmap='coolwarm')
ax02.contour(np.nanmean(Tb_sub133_cor_sel,axis=0),extent=[-5000,5000,-5000,5000],levels=[-5*np.nanmean(Tb_sub133_cor_sel,axis=0).std(),5*np.nanmean(Tb_sub133_cor_sel,axis=0).std()])
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14,color='k')
im03=ax03.imshow(np.nansum(Tb_sub145_cor_sel,axis=0)/1.e3,extent=[-5000,5000,-5000,5000],origin=0,vmin=-2.2e2,vmax=2.2e2,cmap='coolwarm')
ax03.contour(np.nanmean(Tb_sub145_cor_sel,axis=0),extent=[-5000,5000,-5000,5000],levels=[-5*np.nanmean(Tb_sub145_cor_sel,axis=0).std(),5*np.nanmean(Tb_sub145_cor_sel,axis=0).std()])
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14,color='k')
ax10.imshow(np.nansum(Tb_sub179_cor_sel,axis=0)/1.e3,extent=[-5000,5000,-5000,5000],origin=0,vmin=-2.2e2,vmax=2.2e2,cmap='coolwarm')
ax10.contour(np.nanmean(Tb_sub179_cor_sel,axis=0),extent=[-5000,5000,-5000,5000],levels=[-5*np.nanmean(Tb_sub179_cor_sel,axis=0).std(),7*np.nanmean(Tb_sub179_cor_sel,axis=0).std()],colors='maroon')
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14,color='k')
ax11.imshow(np.nansum(Tb_sub197_cor_sel,axis=0)/1.e3,extent=[-5000,5000,-5000,5000],origin=0,vmin=-2.2e2,vmax=2.2e2,cmap='coolwarm')
ax11.contour(np.nanmean(Tb_sub197_cor_sel,axis=0),extent=[-5000,5000,-5000,5000],levels=[-5*np.nanmean(Tb_sub197_cor_sel,axis=0).std(),7*np.nanmean(Tb_sub197_cor_sel,axis=0).std()],colors='maroon')
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14,color='k')
ax12.imshow(np.nansum(Tb_sub217_cor_sel,axis=0)/1.e3,extent=[-5000,5000,-5000,5000],origin=0,vmin=-2.2e2,vmax=2.2e2,cmap='coolwarm')
ax12.contour(np.nanmean(Tb_sub217_cor_sel,axis=0),extent=[-5000,5000,-5000,5000],levels=[-5*np.nanmean(Tb_sub217_cor_sel,axis=0).std(),7*np.nanmean(Tb_sub217_cor_sel,axis=0).std()],colors='maroon')
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14,color='k')
im13=ax13.imshow(np.nansum(Tb_sub240_cor_sel,axis=0)/1.e3,extent=[-5000,5000,-5000,5000],origin=0,vmin=-2.2e2,vmax=2.2e2,cmap='coolwarm')
ax13.contour(np.nanmean(Tb_sub240_cor_sel,axis=0),extent=[-5000,5000,-5000,5000],levels=[-5*np.nanmean(Tb_sub240_cor_sel,axis=0).std(),7*np.nanmean(Tb_sub240_cor_sel,axis=0).std()],colors='maroon')
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14,color='k')
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
ax00.set_xlim(-2500,2500);ax00.set_ylim(-2500,2500);ax01.set_xlim(-2500,2500);ax01.set_ylim(-2500,2500);ax02.set_xlim(-2500,2500);ax02.set_ylim(-2500,2500);ax03.set_xlim(-2500,2500);ax03.set_ylim(-2500,2500)
ax10.set_xlim(-2500,2500);ax10.set_ylim(-2500,2500);ax11.set_xlim(-2500,2500);ax11.set_ylim(-2500,2500);ax12.set_xlim(-2500,2500);ax12.set_ylim(-2500,2500);ax13.set_xlim(-2500,2500);ax13.set_ylim(-2500,2500)
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
c00=plt.Circle((0,0),960,fill=False,color='k');ax00.add_artist(c00);c01=plt.Circle((0,0),960,fill=False,color='k');ax01.add_artist(c01);c02=plt.Circle((0,0),960,fill=False,color='k');ax02.add_artist(c02);c03=plt.Circle((0,0),960,fill=False,color='k');ax03.add_artist(c03)
c10=plt.Circle((0,0),960,fill=False,color='k');ax10.add_artist(c10);c11=plt.Circle((0,0),960,fill=False,color='k');ax11.add_artist(c11);c12=plt.Circle((0,0),960,fill=False,color='k');ax12.add_artist(c12);c13=plt.Circle((0,0),960,fill=False,color='k');ax13.add_artist(c13)
divider = make_axes_locatable(ax03);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im03, cax=cax, label='$T_B$ (kK)', orientation='vertical')
divider = make_axes_locatable(ax13);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im13, cax=cax, label='$T_B$ (kK)', orientation='vertical')
plt.show()


f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.imshow(Tb_sub108_cor_sel_lim_mean,extent=[-5000,5000,-5000,5000],origin=0,vmin=3000,vmax=8000,cmap='YlOrRd')
#ax00.imshow(Tb_sub108_cor_sel_lim_mean/np.nanmax(Tb_sub108_cor_sel_lim_mean),extent=[-5000,5000,-5000,5000],origin=0,vmin=0.1,vmax=1,cmap='YlOrRd')
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14,color='w')
ax01.imshow(Tb_sub120_cor_sel_lim_mean,extent=[-5000,5000,-5000,5000],origin=0,vmin=3000,vmax=8000,cmap='YlOrRd')
#ax01.imshow(Tb_sub120_cor_sel_lim_mean/np.nanmax(Tb_sub120_cor_sel_lim_mean),extent=[-5000,5000,-5000,5000],origin=0,vmin=0.5,vmax=1,cmap='YlOrRd')
ax01.annotate('(a) 120 MHz',xy=(-50,2000),fontsize=14,color='w')
ax02.imshow(Tb_sub133_cor_sel_lim_mean,extent=[-5000,5000,-5000,5000],origin=0,vmin=3000,vmax=8000,cmap='YlOrRd')
#ax02.imshow(Tb_sub133_cor_sel_lim_mean/np.nanmax(Tb_sub133_cor_sel_lim_mean),extent=[-5000,5000,-5000,5000],origin=0,vmin=0.7,vmax=1,cmap='YlOrRd')
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14,color='w')
im03=ax03.imshow(Tb_sub145_cor_sel_lim_mean,extent=[-5000,5000,-5000,5000],origin=0,vmin=3000,vmax=8000,cmap='YlOrRd')
#im03=ax03.imshow(Tb_sub145_cor_sel_lim_mean/np.nanmax(Tb_sub145_cor_sel_lim_mean),extent=[-5000,5000,-5000,5000],origin=0,vmin=0.1,vmax=1,cmap='YlOrRd')
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14,color='w')
ax10.imshow(Tb_sub179_cor_sel_lim_mean,extent=[-5000,5000,-5000,5000],origin=0,vmin=3000,vmax=8000,cmap='YlOrRd')
#ax10.imshow(Tb_sub179_cor_sel_lim_mean/np.nanmax(Tb_sub179_cor_sel_lim_mean),extent=[-5000,5000,-5000,5000],origin=0,vmin=0.1,vmax=1,cmap='YlOrRd')
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14,color='w')
ax11.imshow(Tb_sub197_cor_sel_lim_mean,extent=[-5000,5000,-5000,5000],origin=0,vmin=3000,vmax=8000,cmap='YlOrRd')
#ax11.imshow(Tb_sub197_cor_sel_lim_mean/np.nanmax(Tb_sub179_cor_sel_lim_mean),extent=[-5000,5000,-5000,5000],origin=0,vmin=0.1,vmax=1,cmap='YlOrRd')
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14,color='w')
ax12.imshow(Tb_sub217_cor_sel_lim_mean,extent=[-5000,5000,-5000,5000],origin=0,vmin=3000,vmax=8000,cmap='YlOrRd')
#ax12.imshow(Tb_sub217_cor_sel_lim_mean/np.nanmax(Tb_sub217_cor_sel_lim_mean),extent=[-5000,5000,-5000,5000],origin=0,vmin=0.1,vmax=1,cmap='YlOrRd')
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14,color='w')
im13=ax13.imshow(Tb_sub240_cor_sel_lim_mean,extent=[-5000,5000,-5000,5000],origin=0,vmin=3000,vmax=8000,cmap='YlOrRd')
#im13=ax13.imshow(Tb_sub240_cor_sel_lim_mean/np.nanmax(Tb_sub240_cor_sel_lim_mean),extent=[-5000,5000,-5000,5000],origin=0,vmin=0.1,vmax=1,cmap='YlOrRd')
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14,color='w')
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
c00=plt.Circle((0,0),960,fill=False,color='k');ax00.add_artist(c00);c01=plt.Circle((0,0),960,fill=False,color='k');ax01.add_artist(c01);c02=plt.Circle((0,0),960,fill=False,color='k');ax02.add_artist(c02);c03=plt.Circle((0,0),960,fill=False,color='k');ax03.add_artist(c03)
c10=plt.Circle((0,0),960,fill=False,color='k');ax10.add_artist(c10);c11=plt.Circle((0,0),960,fill=False,color='k');ax11.add_artist(c11);c12=plt.Circle((0,0),960,fill=False,color='k');ax12.add_artist(c12);c13=plt.Circle((0,0),960,fill=False,color='k');ax13.add_artist(c13)
divider = make_axes_locatable(ax03);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im03, cax=cax, label='$T_B$ (K)', orientation='vertical')
divider = make_axes_locatable(ax13);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im13, cax=cax, label='$T_B$ (K)', orientation='vertical')
plt.show()


f,ax=plt.subplots(2,4,figsize=(16, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
ax00.imshow(np.nanmean(Tb_sub108[650:1500],axis=0),extent=[-2500,2500,-2500,2500],origin=0,vmin=-100,vmax=1500)
ax00.annotate('(a) 108 MHz',xy=(-50,2000),fontsize=14,color='w')
ax01.imshow(np.nanmean(Tb_sub120[600:],axis=0),extent=[-2500,2500,-2500,2500],origin=0,vmin=-100,vmax=1500)
ax01.annotate('(a) 120 MHz',xy=(-50,2000),fontsize=14,color='w')
ax02.imshow(np.nanmean(Tb_sub133[:2000],axis=0),extent=[-2500,2500,-2500,2500],origin=0,vmin=-100,vmax=1500)
ax02.annotate('(c) 133 MHz',xy=(-50,2000),fontsize=14,color='w')
im03=ax03.imshow(np.nanmean(Tb_sub145[:3500],axis=0),extent=[-2500,2500,-2500,2500],origin=0,vmin=-100,vmax=1500)
ax03.annotate('(d) 145 MHz',xy=(50,2000),fontsize=14,color='w')
ax10.imshow(np.nanmean(Tb_sub179,axis=0),extent=[-2500,2500,-2500,2500],origin=0,vmin=-100,vmax=1500)
ax10.annotate('(e) 179 MHz',xy=(50,2000),fontsize=14,color='w')
ax11.imshow(np.nanmean(Tb_sub197,axis=0),extent=[-2500,2500,-2500,2500],origin=0,vmin=-100,vmax=1500)
ax11.annotate('(f) 197 MHz',xy=(50,2000),fontsize=14,color='w')
ax12.imshow(np.nanmean(Tb_sub217,axis=0),extent=[-2500,2500,-2500,2500],origin=0,vmin=-100,vmax=1500)
ax12.annotate('(g) 217 MHz',xy=(50,2000),fontsize=14,color='w')
im13=ax13.imshow(np.nanmean(Tb_sub240,axis=0),extent=[-2500,2500,-2500,2500],origin=0,vmin=-100,vmax=1500)
ax13.annotate('(h) 240 MHz',xy=(50,2000),fontsize=14,color='w')
ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
ax01.set_yticks([]);ax02.set_yticks([]);ax03.set_yticks([]);ax11.set_yticks([]);ax12.set_yticks([]);ax13.set_yticks([]);ax00.set_xticks([]);ax01.set_xticks([]);ax02.set_xticks([]);ax03.set_xticks([])
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0.05)
c00=plt.Circle((0,0),960,fill=False,color='k');ax00.add_artist(c00);c01=plt.Circle((0,0),960,fill=False,color='k');ax01.add_artist(c01);c02=plt.Circle((0,0),960,fill=False,color='k');ax02.add_artist(c02);c03=plt.Circle((0,0),960,fill=False,color='k');ax03.add_artist(c03)
c10=plt.Circle((0,0),960,fill=False,color='k');ax10.add_artist(c10);c11=plt.Circle((0,0),960,fill=False,color='k');ax11.add_artist(c11);c12=plt.Circle((0,0),960,fill=False,color='k');ax12.add_artist(c12);c13=plt.Circle((0,0),960,fill=False,color='k');ax13.add_artist(c13)
divider = make_axes_locatable(ax03);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im03, cax=cax, orientation='vertical',label='T$_{B}$ (K)')
divider = make_axes_locatable(ax13);cax = divider.append_axes('right', size='5%', pad=0.05);f.colorbar(im13, cax=cax, orientation='vertical',label='T$_{B}$ (K)')
plt.show()

Tblim_list=[1000,2000,3000,4000,5000,6000,7000,10000];j=0;levels_list=[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
Tblim_list=[0,0,0,0,0,0,0,000];j=0;levels_list=[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
for Tblim in Tblim_list:
    print Tblim, 'K'
    Tb_sub108[np.where(Tb_sub108<Tblim)]=0;Tb_sub120[np.where(Tb_sub120<Tblim)]=0;Tb_sub133[np.where(Tb_sub133<Tblim)]=0;Tb_sub145[np.where(Tb_sub145<Tblim)]=0;Tb_sub217[np.where(Tb_sub217<Tblim)]=0;Tb_sub240[np.where(Tb_sub240<Tblim)]=0
    Tb_sub179[np.where(Tb_sub179<Tblim)]=0;Tb_sub197[np.where(Tb_sub197<Tblim)]=0
    if(j>4):
        levels_sub=levels_list[4:]
    else:
        levels_sub=levels_list[j:]
    f,ax=plt.subplots(2,4,figsize=(12, 8));ax00=ax[0,0];ax01=ax[0,1];ax02=ax[0,2];ax03=ax[0,3]
    ax10=ax[1,0];ax11=ax[1,1];ax12=ax[1,2];ax13=ax[1,3]
    ax00.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax01.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    if(j==0 or j==1):
        ax00.contour(np.nanmean(Tb_sub108,axis=0)/np.nanmax(np.nanmean(Tb_sub108,axis=0)), levels_sub,extent=[-5000,5000,-5000,5000],  colors='red',linewidths=1)
        ax00.annotate('(a) 108 MHz',xy=(50,4200))
        ax01.contour(np.nanmean(Tb_sub120,axis=0)/np.nanmax(np.nanmean(Tb_sub120,axis=0)), levels_sub,extent=[-5000,5000,-5000,5000],  colors='red',linewidths=1)
        ax01.annotate('(b) 120 MHz',xy=(50,4200))
    else:
        ax00.contour(np.nanmean(Tb_sub108,axis=0)[50:-50,50:-50]/np.nanmax(np.nanmean(Tb_sub108,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
        ax00.annotate('(a) 108 MHz',xy=(50,2000))
        ax01.contour(np.nanmean(Tb_sub120,axis=0)[50:-50,50:-50]/np.nanmax(np.nanmean(Tb_sub120,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
        ax01.annotate('(b) 120 MHz',xy=(50,2000))
    ax02.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax02.contour(np.nanmean(Tb_sub133,axis=0)/np.nanmax(np.nanmean(Tb_sub133,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax02.annotate('(c) 133 MHz',xy=(50,2000))
    ax03.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax03.contour(np.nanmean(Tb_sub145,axis=0)/np.nanmax(np.nanmean(Tb_sub145,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax03.annotate('(d) 145 MHz',xy=(50,2000))
    ax10.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax10.contour(np.nanmean(Tb_sub179,axis=0)/np.nanmax(np.nanmean(Tb_sub179,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax10.annotate('(e) 179 MHz',xy=(50,2000))
    ax11.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax11.contour(np.nanmean(Tb_sub197,axis=0)/np.nanmax(np.nanmean(Tb_sub197,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax11.annotate('(f) 197 MHz',xy=(50,2000))
    ax12.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax12.contour(np.nanmean(Tb_sub217,axis=0)/np.nanmax(np.nanmean(Tb_sub217,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax12.annotate('(g) 217 MHz',xy=(50,2000))
    ax13.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax13.contour(np.nanmean(Tb_sub240,axis=0)/np.nanmax(np.nanmean(Tb_sub240,axis=0)), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax13.annotate('(h) 240 MHz',xy=(50,2000))
    ax10.set_xlabel('Solar X (arcsec)');ax10.set_ylabel('Solar Y (arcsec)');ax11.set_xlabel('Solar X (arcsec)');ax12.set_xlabel('Solar X (arcsec)');ax13.set_xlabel('Solar X (arcsec)');ax00.set_ylabel('Solar Y (arcsec)')
    f.tight_layout()
    f.suptitle('$T_{B} >$ '+str(Tblim)+' K', fontsize="large")
    f.savefig('pngs_mean/Tb_'+str("%03d"%j)+'.png',dpi=100)
    plt.close()
    j=j+1




for j in range(3500):
    f,ax=plt.subplots(4,4,figsize=(12, 10))
    ax00=ax[0,0];ax10=ax[1,0];ax20=ax[2,0];ax30=ax[3,0]
    ax01=ax[0,1];ax11=ax[1,1];ax21=ax[2,1];ax31=ax[3,1]
    ax02=ax[0,2];ax12=ax[1,2];ax22=ax[2,2];ax32=ax[3,2]
    ax03=ax[0,3];ax13=ax[1,3];ax23=ax[2,3];ax33=ax[3,3]
    ax00.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax00.contour(Tb_108[j]/np.nanmax(Tb_108[j]), levels,extent=[-2500,2500,-2500,2500],  colors='black',linewidths=1);ax00.annotate('108 MHz',xy=(50,2000))
    ax00.contour(Tb_sub108[j]/np.nanmax(Tb_sub108[j]), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax01.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax01.contour(Tb_120[j]/np.nanmax(Tb_120[j]), levels,extent=[-2500,2500,-2500,2500],  colors='black',linewidths=1);ax01.annotate('120 MHz',xy=(50,2000))
    ax01.contour(Tb_sub120[j]/np.nanmax(Tb_sub120[j]), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax02.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax02.contour(Tb_133[j]/np.nanmax(Tb_133[j]), levels,extent=[-2500,2500,-2500,2500],  colors='black',linewidths=1);ax02.annotate('133 MHz',xy=(50,2000))
    ax02.contour(Tb_sub133[j]/np.nanmax(Tb_sub133[j]), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax03.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax03.contour(Tb_145[j]/np.nanmax(Tb_145[j]), levels,extent=[-2500,2500,-2500,2500],  colors='black',linewidths=1);ax03.annotate('145 MHz',xy=(50,2000))
    ax03.contour(Tb_sub145[j]/np.nanmax(Tb_sub145[j]), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax10.plot(t,Tbmax108,'o-',markersize=4,color='k');ax10.set_yscale('log');ax10.set_ylim([1.e3,1.e6]);ax10.axvline(x=t[j],color='gray')
    ax11.plot(t,Tbmax120,'o-',markersize=4,color='k');ax11.set_yscale('log');ax11.set_ylim([1.e3,1.e6]);ax11.axvline(x=t[j],color='gray')
    ax12.plot(t,Tbmax133,'o-',markersize=4,color='k');ax12.set_yscale('log');ax12.set_ylim([1.e3,1.e6]);ax12.axvline(x=t[j],color='gray')
    ax13.plot(t,Tbmax145,'o-',markersize=4,color='k');ax13.set_yscale('log');ax13.set_ylim([1.e3,1.e6]);ax13.axvline(x=t[j],color='gray')
    ax20.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax20.contour(Tb_145[j]/np.nanmax(Tb_145[j]), levels,extent=[-2500,2500,-2500,2500],  colors='black',linewidths=1);ax20.annotate('145 MHz',xy=(50,2000))
    ax20.contour(Tb_sub145[j]/np.nanmax(Tb_sub145[j]), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax21.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax21.contour(Tb_217[j]/np.nanmax(Tb_217[j]), levels,extent=[-2500,2500,-2500,2500],  colors='black',linewidths=1);ax21.annotate('217 MHz',xy=(50,2000))
    ax21.contour(Tb_sub217[j]/np.nanmax(Tb_sub217[j]), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax22.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax22.contour(Tb_217[j]/np.nanmax(Tb_217[j]), levels,extent=[-2500,2500,-2500,2500],  colors='black',linewidths=1);ax22.annotate('217 MHz',xy=(50,2000))
    ax22.contour(Tb_sub217[j]/np.nanmax(Tb_sub217[j]), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax23.imshow(d,aspect='equal',interpolation='none',extent=[-1229,1229,-1229,1229],origin='lower',vmin=0,vmax=700,cmap='YlGnBu')
    ax23.contour(Tb_240[j]/np.nanmax(Tb_240[j]), levels,extent=[-2500,2500,-2500,2500],  colors='black',linewidths=1);ax23.annotate('240 MHz',xy=(50,2000))
    ax23.contour(Tb_sub240[j]/np.nanmax(Tb_sub240[j]), levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax30.plot(t,Tbmax145,'o-',markersize=4,color='k');ax30.set_yscale('log');ax30.set_ylim([1.e3,1.e6]);ax30.axvline(x=t[j],color='gray')
    ax31.plot(t,Tbmax217,'o-',markersize=4,color='k');ax31.set_yscale('log');ax31.set_ylim([1.e3,1.e6]);ax31.axvline(x=t[j],color='gray')
    ax32.plot(t,Tbmax217,'o-',markersize=4,color='k');ax32.set_yscale('log');ax32.set_ylim([1.e3,1.e6]);ax32.axvline(x=t[j],color='gray')
    ax33.plot(t,Tbmax240,'o-',markersize=4,color='k');ax33.set_yscale('log');ax33.set_ylim([1.e3,1.e6]);ax33.axvline(x=t[j],color='gray')
    f.tight_layout()
    f.savefig('images_pngs/Tb_'+str("%03d"%j)+'.png',dpi=100)
    plt.close()

ax1.add_subplot(411)

    
sys.exit()

nc=8
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
cc=cm.rainbow(np.linspace(0, 1, nc))
idx1=[16,40,51,60,69,79,105,133]
f=plt.figure()
ax1=f.add_subplot(211)
im1=ax1.imshow(Tb_240[0]/1.e6,aspect='auto',interpolation='none',extent=[-2500,2500,-2500,2500],origin='lower',vmin=0,vmax=0.5,cmap='YlGnBu')
ax1.contour(Tb_240[0]/np.nanmax(Tb_240[0]), levels,extent=[-2500,2500,-2500,2500],  colors='white',linewidths=1)
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
ax2=f.add_subplot(212)
ax2.plot(Tbmax[1500:1650]/1.e3,'o-',color='k',markersize=5)
ax1.contour(Tb_sub[1500+idx1[1]]/np.nanmax(Tb_sub[1500+idx1[1]]),extent=[-2500,2500,-2500,2500],levels=[0.5,0.7,0.9],colors='blue')
ax1.contour(Tb_sub[1500+idx1[6]]/np.nanmax(Tb_sub[1500+idx1[6]]),extent=[-2500,2500,-2500,2500],levels=[0.5,0.7,0.9],colors='red')
i=0
for y,c in zip(Tbmax[1500:1650][idx1]/1.e3,cc):
    plt.rcParams['lines.markersize']=20
    ax2.scatter(idx1[i],y,color=c,marker='*')
    yy=(np.where(Tb_sub[1500+idx1[i]]==np.nanmax(Tb_sub[1500+idx1[i]]))[0]-50)*50
    xx=(np.where(Tb_sub[1500+idx1[i]]==np.nanmax(Tb_sub[1500+idx1[i]]))[1]-50)*50
    ax1.scatter(xx,yy,color=c,edgecolors='black')
    plt.rcParams["scatter.marker"] = '*';plt.rcParams['lines.markersize']=20
    i=i+1
ax1.set_xlim(-1400,1400);ax1.set_ylim(-1400,1400)
ax2.set_ylabel("$T_B$ (kK)")
ax1.set_ylabel("Solar Y (arcsec)")
ax1.set_xlabel("Solar X (arcsec)")
ax2.set_xlabel("Time (HH:MM:SS UT)")
ax2.set_xticks([0,50,100,150])
ax2.set_xticklabels(['03:54:45','03:55:35','03:56:25','03:57:15'])
f.tight_layout()
f.colorbar(im1, cax=cax, orientation='vertical',label="T$_{B}$ (MK)")
plt.show()



sys.exit()
for i in range(len(Tb_108)):
    levels=np.array([0.2,0.3,0.5,0.7,0.9])
    levels_sub=np.array([0.1,0.2,0.3,0.5,0.7,0.9])*1.e4
    fig=plt.figure(figsize=(8,10))
    ax=fig.add_subplot(211)
    im=ax.imshow(Tb_108[i]/1.e6,aspect='equal',interpolation='none',extent=[-2500,2500,-2500,2500],origin='lower',vmin=0,vmax=0.8,cmap='YlGnBu')
    CS1=ax.contour(Tb_108[i]/np.nanmax(Tb_108[i]), levels,extent=[-2500,2500,-2500,2500],  colors='white',linewidths=1)
    #CS=ax.contour(Tb_sub[i]/np.nanmax(Tb_sub[i]), levels,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    CS=ax.contour(Tb_sub[i], levels_sub,extent=[-2500,2500,-2500,2500],  colors='red',linewidths=1)
    ax.set_xlabel('X (arcsec)')
    ax.set_ylabel('Y (arcsec)')
    #ax.set_xlim(-3000,3000)
    #ax.set_ylim(-3000,3000)
    r1=16.*60;r2=32.*60
    circ1=plt.Circle((0.5,0.5), radius=r1, color='brown', linewidth=4,fill=False)
    circ2=plt.Circle((0.5,0.5), radius=r2, color='brown', linewidth=4,fill=False)
    ax.add_patch(circ2);ax.add_patch(circ1)
    ax.grid(True)
    ax.set_title('Contours: 1, 2, 3, 5, 7, 9 kK')
    ax=fig.add_subplot(212)
    ax.plot(Tb_sub.max(axis=(1,2)))
    ax.axvline(x=i,color='k')
    ax.set_xlabel('Time (sec)');ax.set_ylabel('$T_B$ (K)')
    fig.colorbar(im,label='(MK)')
    fig.savefig('pngs_real_108MHz/Tb_'+str("%03d"%i)+'.png',dpi=100)
    plt.close()

