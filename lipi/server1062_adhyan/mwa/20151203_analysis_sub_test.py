import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob
import sunpy
from sunpy.map import Map
from astropy.io import fits
from surya.utils import main as ut

plt.style.use('/nas08-data02/rohit/scripts/general/plt_style.py')
file_161='/nas08-data02/rohit/20151203_MWA/new_ms/fits/Tb_20151203_125-126.p'
subfile_161='/nas08-data02/rohit/20151203_sub_run_mean/Tb/Tb_20151203_125-126_sub.p'
file_240='/nas08-data02/rohit/20151203_MWA/new_ms/fits/Tb_20151203_187-188.p'
subfile_240='/nas08-data02/rohit/20151203_sub_run_mean/Tb/Tb_20151203_187-188_sub.p'
file_217='/nas08-data02/rohit/20151203_MWA/new_ms/fits/Tb_20151203_169-170.p'
subfile_217='/nas08-data02/rohit/20151203_sub_run_mean/Tb/Tb_20151203_169-170_sub.p'
file_179='/nas08-data02/rohit/20151203_MWA/new_ms/fits/Tb_20151203_139-140.p'
subfile_179='/nas08-data02/rohit/20151203_sub_run_mean/Tb/Tb_20151203_139-140_sub.p'
file_120='/nas08-data02/rohit/20151203_MWA/new_ms/fits/Tb_20151203_093-094.p'
subfile_120='/nas08-data02/rohit/20151203_sub_run_mean/Tb/Tb_20151203_093-094_sub.p'
subfile_145='/nas08-data02/rohit/20151203_sub_run_mean/Tb/Tb_20151203_113-114_sub.p'
file_145='/nas08-data02/rohit/20151203_MWA/new_ms/fits/Tb_20151203_113-114.p'
file_108='/nas08-data02/rohit/20151203_MWA/new_ms/fits/Tb_20151203_084-085.p'
subfile_108='/nas08-data02/rohit/20151203_sub_run_mean/Tb/Tb_20151203_084-085_sub.p'
imgpath='/nas08-data02/rohit/20151203_sub_run_mean/pngs/'

def get_Tb(file_,subfile_):
    print "Reading "+file_+'....'
    data=pickle.load(open(file_,'rb'))
    subdata=pickle.load(open(subfile_,'rb'))
    Tb=np.array(data[0])
    Tbsub_=subdata[0]
    #ndata=np.array(data[7])
    #ndata_sub=np.array(subdata[-1])[:,50:-50,50:-50]
    Tbsub=[0]*Tbsub_.shape[0]
    for i in range(Tbsub_.shape[0]):
        #Tbsub[i]=Tbsub_[i][50:-50,50:-50]
        Tbsub[i]=Tbsub_[i][50:-50,50:-50]*0+data[0][::2][i].max()/data[-5][::2][i].max()*subdata[-1][i][50:-50,50:-50]
    Tbsub=np.array(Tbsub)
    ts_data=data[9];ts_subdata=subdata[4]
    tssec_data=[0]*len(ts_data);tssec_subdata=[0]*len(ts_subdata)
    for i in range(len(ts_data)):
        tssec_data[i]=ut.hms2sec_c(' '+ts_data[i])
    tssec_data=np.array(tssec_data)
    for i in range(len(ts_subdata)):
        tssec_subdata[i]=ut.hms2sec_c(' '+ts_subdata[i])
    tssec_subdata=np.array(tssec_subdata)
    Tb_disk=[0]*len(ts_subdata)
    for i in range(len(ts_subdata)):
        idx=ut.find_nearest(tssec_data,tssec_subdata[i])[0]
        Tb_disk[i]=Tb[idx]
    Tb_disk=np.array(Tb_disk)
    return data,subdata,Tb_disk,Tbsub,tssec_data,tssec_subdata,ts_data,ts_subdata

data108,subdata108,Tb_disk108,Tbsub108,tssec_data108,tssec_subdata108,ts_data108,ts_subdata108=get_Tb(file_108,subfile_108)
data120,subdata120,Tb_disk120,Tbsub120,tssec_data120,tssec_subdata120,ts_data120,ts_subdata120=get_Tb(file_120,subfile_120)
data145,subdata145,Tb_disk145,Tbsub145,tssec_data145,tssec_subdata145,ts_data145,ts_subdata145=get_Tb(file_145,subfile_145)
data161,subdata161,Tb_disk161,Tbsub161,tssec_data161,tssec_subdata161,ts_data161,ts_subdata161=get_Tb(file_161,subfile_161)
data179,subdata179,Tb_disk179,Tbsub179,tssec_data179,tssec_subdata179,ts_data179,ts_subdata179=get_Tb(file_179,subfile_179)
data217,subdata217,Tb_disk217,Tbsub217,tssec_data217,tssec_subdata217,ts_data217,ts_subdata217=get_Tb(file_217,subfile_217)
data240,subdata240,Tb_disk240,Tbsub240,tssec_data240,tssec_subdata240,ts_data240,ts_subdata240=get_Tb(file_240,subfile_240)
    
r1=[44,54,36,46];Tbsub108_r1=Tbsub108[:,r1[0]:r1[1],r1[2]:r1[3]];Tbsub240_r1=Tbsub240[:,r1[0]:r1[1],r1[2]:r1[3]];Tbsub161_r1=Tbsub161[:,r1[0]:r1[1],r1[2]:r1[3]]

cont_cut=0.0
Tbsub240_pos=np.zeros(Tbsub240.shape[0]);Tbsub240_neg=np.zeros(Tbsub240.shape[0])
Tbsub240_mean=Tbsub240.mean(axis=0)
Tbsub179_pos=np.zeros(Tbsub179.shape[0]);Tbsub179_neg=np.zeros(Tbsub179.shape[0])
Tbsub179_mean=Tbsub179.mean(axis=0)
for i in range(Tbsub240.shape[0]):
    Tbsub240_pos[i]=np.mean(Tbsub240[i][Tbsub240[i]>0])
    Tbsub240_neg[i]=np.mean(Tbsub240[i][Tbsub240[i]<0])
    Tbsub179_pos[i]=np.mean(Tbsub179[i][Tbsub179[i]>0])
    Tbsub179_neg[i]=np.mean(Tbsub179[i][Tbsub179[i]<0])

sys.exit()
sub00='/nas08-data02/rohit/20151203_sub_run/161MHz/images/'
aiafile='/nas08-data02/rohit/20151203_EUV/aia.lev1.193A_2015-12-03T03_25_29.84Z.image_lev1.fits'
aia=fits.open(aiafile)
aiamap=Map(aiafile)

sys.exit()

j=1
for i in range(len(Tb_disk161)):
    plt.imshow(np.log10(aiamap.data),origin='lower',extent=[-1229,1229,-1229,1229],cmap='sdoaia193')
    plt.contour(Tbsub161[i],levels=[1100,1200,2000,2500,3000,4000,5000],origin='lower',extent=[-2500,2500,-2500,2500],colors='red')
    plt.title(ts_subdata161[i]);plt.xlim(-1500,1500);plt.ylim(-1500,1500)
    plt.savefig(imgpath+'/Tb_161MHz_'+"%04d"%i+'.png')
    plt.close()


for i in range(1200):
    plt.imshow(np.log10(aiamap.data),origin='lower',extent=[-1229,1229,-1229,1229],cmap='sdoaia193')
    #plt.contour(Tbsub161[i],levels=[1100,2000,2500,3000,4000,5000],origin='lower',extent=[-2500,2500,-2500,2500],colors='red')
    #plt.contour(Tbsub217[i],levels=[1100,2000,2500,3000,4000,5000],origin='lower',extent=[-2500,2500,-2500,2500],colors='green')
    #plt.contour(Tbsub240[i],levels=[1100,2000,2500,3000,4000,5000],origin='lower',extent=[-2500,2500,-2500,2500],colors='blue')
    plt.contour(Tbsub161[i]/Tbsub161[i].max(),levels=[0.7,0.8,0.9],origin='lower',extent=[-2500,2500,-2500,2500],colors='red')
    #plt.contour(Tbsub217[i]/Tbsub217[i].max(),levels=[0.7,0.8,0.9],origin='lower',extent=[-2500,2500,-2500,2500],colors='green')
    #plt.contour(Tbsub240[i]/Tbsub240[i].max(),levels=[0.7,0.8,0.9],origin='lower',extent=[-2500,2500,-2500,2500],colors='blue')
    plt.title(ts_subdata161[i]);plt.xlim(-1500,1500);plt.ylim(-1500,1500)
    plt.savefig(imgpath+'/Tb_161MHz_240MHz_'+"%04d"%i+'.png')
    plt.close()

plt.imshow(np.log10(aiamap.data),origin='lower',extent=[-1229,1229,-1229,1229],cmap='Greys')
#plt.contour(Tbsub240.sum(axis=0)/Tbsub240.sum(axis=0).max(),levels=[-0.3,-0.2,0,0.3,0.4,0.5,0.6,0.7,0.8,0.9],origin='lower',extent=[-2500,2500,-2500,2500],colors='red')
#plt.contour(Tbsub161.sum(axis=0)/Tbsub161.sum(axis=0).max(),levels=[-0.3,-0.2,0,0.3,0.4,0.5,0.6,0.7,0.8,0.9],origin='lower',extent=[-2500,2500,-2500,2500],colors='red')
#plt.contour(Tbsub217.sum(axis=0)/Tbsub217.sum(axis=0).max(),levels=[-0.3,-0.2,0,0.3,0.4,0.5,0.6,0.7,0.8,0.9],origin='lower',extent=[-2500,2500,-2500,2500],colors='red')
plt.contour(Tbsub108.sum(axis=0)/Tbsub108.sum(axis=0).max(),levels=[-0.3,-0.2,0,0.3,0.4,0.5,0.6,0.7,0.8,0.9],origin='lower',extent=[-2500,2500,-2500,2500],colors='red')
plt.contour(Tbsub120.sum(axis=0)/Tbsub120.sum(axis=0).max(),levels=[-0.3,-0.2,0,0.3,0.4,0.5,0.6,0.7,0.8,0.9],origin='lower',extent=[-2500,2500,-2500,2500],colors='blue')
#plt.contour(Tbsub145.sum(axis=0)/Tbsub145.sum(axis=0).max(),levels=[-0.3,-0.2,0,0.3,0.4,0.5,0.6,0.7,0.8,0.9],origin='lower',extent=[-2500,2500,-2500,2500],colors='blue')
#plt.contour(Tbsub179.sum(axis=0)/Tbsub179.sum(axis=0).max(),levels=[-0.3,-0.2,0,0.3,0.4,0.5,0.6,0.7,0.8,0.9],origin='lower',extent=[-2500,2500,-2500,2500],colors='blue')
plt.xlim(-1600,1600);plt.ylim(-1600,1600)
plt.xlabel('Solar X (arcsec)')
plt.ylabel('Solar Y (arcsec)')
plt.show()

label_all=['108 MHz','120 MHz','145 MHz','161 MHz','179 MHz','240 MHz']
Tbsub_all=np.array([Tbsub108,Tbsub120,Tbsub145,Tbsub161,Tbsub179,Tbsub240])
cm = plt.get_cmap('gist_rainbow')
x=np.arange(1232).reshape((7,176))
for i in range(7):
    x[i]=x[i]+120*i
x=x.flatten()
colormap = plt.cm.gist_ncar
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(Tbsub_all))])
for i in range(len(label_all)):
    plt.plot(x,Tbsub_all[i].mean(axis=(1,2)),'o',label=label_all[i])
plt.xlabel('Time (sec)')
plt.ylabel('T$_B$ (K)')
plt.legend()
plt.show()

sys.exit()

dTb=Tb-np.mean(Tb,axis=0)
dndata=ndata-np.mean(ndata,axis=0)
xcmax=[0]*3500;ycmax=[0]*3500;xcmin=[0]*3500;ycmin=[0]*3500;Tbsubmax=[0]*3500;Tbsubmin=[0]*3500
for i in range(3500):
    Tbsubmax[i]=np.max(Tbsub[i])
    Tbsubmin[i]=np.max(Tbsub[i]*-1)
    ycmax[i]=(np.where(Tbsub[i]==Tbsubmax[i])[0][0]-50)*50;xcmax[i]=(np.where(Tbsub[i]==Tbsubmax[i])[1][0]-50)*50
    ycmin[i]=(np.where(Tbsub[i]*-1==Tbsubmin[i])[0][0]-50)*50;xcmin[i]=(np.where(Tbsub[i]*-1==Tbsubmin[i])[1][0]-50)*50
Tbsubmax=np.array(Tbsubmax)
Tbsubmax[np.where(Tbsubmax>4000)]=1200
#xcmax=xcmax[300:-300];xcmin=xcmin[300:-300];ycmax=ycmax[300:-300];ycmin=ycmin[300:-300];Tbsubmax=Tbsubmax[300:-300];Tbsubmin=Tbsubmin[300:-300]

f11=360;f12=380;f21=910;f22=930;f31=1040;f32=1050;f41=1600;f42=1640;f51=1850;f52=1860
f61=2460;f62=2480;f71=2850;f72=2880;f81=3470;f82=3490
Tbflare=[Tbsub[f11:f12],Tbsub[f21:f22],Tbsub[f31:f32],Tbsub[f41:f42],Tbsub[f51:f52],Tbsub[f61:f62],Tbsub[f71:f72],Tbsub[f81:f82]]

plt.plot(np.arange(20)*0.5,Tbflare[0].max(axis=(1,2)),'o-',label='Flare 1')
plt.plot(np.arange(10)*0.5,Tbflare[4].max(axis=(1,2)),'o-',label='Flare 2')
plt.plot(np.arange(20)*0.5,Tbflare[7].max(axis=(1,2)),'o-',label='Flare 3')
plt.legend()
plt.xlabel('Time (sec)')
plt.ylabel('$T_B$ (K)')
plt.show()

sys.exit()

plot_max_all=0
if(plot_max_all):
    f,ax=plt.subplots(nrows=3)
    ax[0].hist(np.array(Tbsubmax),bins=40,histtype='step',label='SUB MAX')
    ax[1].hist(np.array(Tb).max(axis=(1,2)),bins=40,histtype='step',label='MAX')
    ax[2].hist(np.array(Tb).mean(axis=(1,2)),bins=400,histtype='step',label='MEAN')
    ax[0].legend();ax[1].legend();ax[2].legend()
    plt.show()

plot_map_max_min=0
if(plot_map_max_min):
    limpl=1800
    #plt.imshow(Tb[11],aspect='equal',origin=0,cmap='binary')
    cm1 = plt.cm.get_cmap('YlOrRd')
    cm2 = plt.cm.get_cmap('winter')
    plt.contour(Tb[11],levels=np.array([0.1,0.2,0.25,0.3,0.32])*1.e6,extent=[-2500,2500,-2500,2500],colors='black')
    sc1=plt.scatter(np.array(xcmax)[np.where(np.array(Tbsubmax)>limpl)],np.array(ycmax)[np.where(np.array(Tbsubmax)>limpl)],c=np.array(Tbsubmax)[np.where(np.array(Tbsubmax)>limpl)],s=35,cmap=cm1,vmin=limpl,vmax=1.5*limpl)
    sc2=plt.scatter(np.array(xcmin)[np.where(np.array(Tbsubmin)>limpl)],np.array(ycmin)[np.where(np.array(Tbsubmin)>limpl)],c=-1*np.array(Tbsubmin)[np.where(np.array(Tbsubmin)>limpl)],s=35,cmap=cm2,vmax=-1*limpl,vmin=-1.5*limpl)
    plt.xlim(-2500,2500);plt.ylim(-2500,2500)
    plt.colorbar(sc1,orientation='vertical',label='(K)');plt.colorbar(sc2,orientation='horizontal',label='(K)')
    plt.xlabel('Solar X (arcsec)');plt.ylabel('Solar Y (arcsec)')
    plt.show()


plot_Tb=1
if(plot_Tb):
    Tb[np.where(Tb>1.e6)]=3.5e5;Tbsub[np.where(Tbsub>1.e6)]=0
    f,ax=plt.subplots(nrows=3,sharex=True)
    ax[0].plot(np.array(Tb).max(axis=(1,2)),'o-')
    ax[1].plot(np.arange(Tbsub.shape[0])+11,np.array(Tbsub).max(axis=(1,2)),'o-')
    ax[2].plot(np.arange(Tbsub.shape[0])+11,np.array(Tbsub).min(axis=(1,2)),'o-')
    ax[0].set_ylabel('Tb (K)')
    ax[1].set_ylabel('Max Tbsub (K)')
    ax[2].set_ylabel('Min Tbsub (K)')
    ax[2].set_xlabel('Time (sec)')
    ax[2].set_ylim(-3500,0);ax[1].set_ylim(0,3500)
    ax[1].axhline(y=1800,color='k')
    plt.show()
sys.exit()
for i in range(1500):
    idx="%03d"%i
    plt.imshow(Tb[i+11],aspect='equal',origin=0,cmap='binary')
    plt.contour(Tb[i+11],levels=np.array([0.1,0.2,0.25,0.3,0.32])*1.e6,colors='yellow')
    plt.contour(Tbsub[i],levels=np.array([0.15,0.20,0.25,0.3])*1.e4,colors='red')
    plt.savefig('pngs/Tb_'+idx+'.png')
    plt.close()

