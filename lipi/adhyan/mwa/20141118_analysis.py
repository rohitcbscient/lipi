import numpy as np
import matplotlib.pyplot as plt
from surya.aia import main as aia
from surya.radio import get_maps
from astropy.io import fits
from surya.radio import get_maps as maps
import glob
import pickle
import plot_20141118 as pl

aiafile='/home/i4ds1807205/20141118/saia_00193_fd_20141118_223342.fts'

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
ids='20151203_'
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

####################### IMAGING ##############################

img_freq=[108,120,179,197,218,240]
Tb,noise_mean,noise_std,bmaj,bmin,bpa,centre,ndata,s_sun,imgtime,imgsec,imgfiles=pickle.load(open('Tb_20151203_187-188.p','rb'))
Tb=np.array(Tb)
snr=np.nanmax(ndata)/np.array(noise_std)
tick=[2,5,10,15,20,25,30,35]
xlabel_list=['06:55:19','06:55:39','06:56:04','06:56:29','06:56:54','06:57:19','06:57:44','06:58:09','06:58:34','06:58:59','06:59:24','06:59:49']
xtick_list=[10,50,100,150,200,250,300,350,400,450,500,550]
ytick_list=[0,64,88,152,176,240,264,328,352,416,440,504,528,592,616,680,704,768]
ylabel_list=[107.52,108.8,119.04,120.32,131.84,133.12,144.64,145.92,160.00,161.28,177.92,179.2,195.84,197.12,216.32,217.6,239.36,240.64]
#pl.ds_plot(S_sun_bmean,ytick_list,xtick_list,ylabel_list,xlabel_list,3,40,'(SFU)',tick)


Tb_map=[0]*6
Tb_max=[0]*6
size=[0]*6
for i in range(6):
    lev=np.array([0.01,0.03,0.05,0.08,0.1,0.5,0.9])*Tb[-1].max()
    Tb_map[i]=Tb[i]
    Tb_map[i][Tb[i]<0.5e4]=0.0
    Tb_max[i]=np.mean(Tb[i][Tb[i]>lev[-2]])/1.e6
    size[i]=Tb[i][Tb[i]>lev[-2]].shape[0]*99*99/3600
    #pl.plot_map(Tb_map[i],str(img_freq[i])+'MHz',bmin[i],bmaj[i],lev)

#mwafile240='/home/i4ds1807205/20141118/240MHz.fits'
#mwafile218='/home/i4ds1807205/20141118/218MHz.fits'
#mwafile108='/home/i4ds1807205/20141118/108MHz.fits'

#mwatime='06:55:20'
# RA and DEC 15 35 19.72 -19 17 18.8
#mway,mwax=1182,462
#pa=20.3406
#mwa240=fits.open(mwafile240)
#mwadata240=mwa240[0].data[0][0]
#map240_=mwadata240[mwax-100:mwax+100,mway-100:mway+100]
#map240=get_maps.rotateimage(map240_,pa*1,100,100)
#aia=fits.open(aiafile)

#aiamap=aia[0].data

#mwa218=fits.open(mwafile218)
#mwadata218=mwa218[0].data[0][0]
#map218_=mwadata218[456-100:456+100,1184-100:1184+100]
#map218=get_maps.rotateimage(map218_,pa*1,100,100)

#mwa108=fits.open(mwafile108)
#mwadata108=mwa108[0].data[0][0]
#map108_=mwadata108[457-100:457+100,1142-100:1142+100]
#map108=get_maps.rotateimage(map108_,pa*1,100,100)
#aialim=(1024*2.4/2)/3600.
#mwalim=(200*90./2)/3600.



plot_euv_mwa=0
if(plot_euv_mwa):
    plt.imshow(aiamap,origin='lower',extent=[-aialim,aialim,-aialim,aialim],vmin=0,vmax=1200,cmap='binary') # 1229=512*2.4
    plt.contour(map108/np.max(map108),origin='lower',levels=[0.8,0.9],colors='red',linewidths=2,extent=[-mwalim,mwalim,-mwalim,mwalim]) # 40*100
    plt.contour(map108/np.max(map108),origin='lower',levels=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1],colors='blue',linewidths=2,extent=[-mwalim,mwalim,-mwalim,mwalim]) # 40*100
    plt.contour(map108/np.max(map108),origin='lower',levels=[0.001,0.0025,0.005],colors='green',linewidths=2,extent=[-mwalim,mwalim,-mwalim,mwalim]) # 40*100
    plt.contour(map108/np.max(map108),origin='lower',levels=[0.00005],colors='gray',linewidths=2,extent=[-mwalim,mwalim,-mwalim,mwalim]) # 40*100
    plt.xlabel('X (degrees)')
    plt.ylabel('Y (degrees)')
    #plt.xlim([-1200,1200])
    #plt.ylim([-1200,1200])
    plt.show()

