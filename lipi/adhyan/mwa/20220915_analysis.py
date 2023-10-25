import pickle
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from sunpy.map import Map
import itertools
import astropy.units as u
from scipy.ndimage.interpolation import rotate
from surya.utils import model as mdl
from surya.utils import main as ut
from scipy.io import readsav


baseline=['000-002','000-005','000-007','002-005','002-007','005-007']
chlist=[61,65,69,73,81,86,91,96,101,113,120,127,134,142,150,158,167,187,226]
flux=[0]*len(chlist);freq=[0]*len(chlist);dflux=[0]*len(chlist)
i=0
for ch in chlist:
    f2=pickle.load(open('/data/20220915_MWA/pickle/flux_V1_20220915_ch'+str(ch)+'.p','rb'),encoding='latin1')
    flux[i]=np.mean(f2[17][3][0],axis=0)
    freq[i]=ch*1.28
    dflux[i]=flux[i]-flux[i][0]
    i=i+1
freq=np.array(freq);flux=np.array(flux);dflux=np.array(dflux)

#------- Newkirk Model
height_Mm = [0]*len(freq);height_Mm_harmonic=[0]*len(freq)
for i in range(len(freq)):
    height_Mm_harmonic[i]=(float(mdl.nk_freq2r(freq[i]/2.0,1)[0])-1)*6.99e2
    if(freq[i]<240):
        height_Mm[i]=(float(mdl.nk_freq2r(freq[i],1)[0])-1)*6.99e2
    else:
        height_Mm[i]=15 # 15 Mm threshold


f,ax=plt.subplots(1,1)
ax.imshow(dflux,aspect='auto',origin='lower',interpolation=None,vmin=0.1,vmax=10, cmap='coolwarm')
#f.colorbar(label='NCCF')
ax.set_xticks([0,50,100,150,200,250,300,350])
ax.set_xticklabels(['03:29:58','03:38:18','03:46:38','03:54:58','04:03:18','04:11:38','04:19:58','04:28:18'])
ax.set_yticks(np.arange(len(chlist)))
ax.set_yticklabels(np.round(freq,2));ax.set_xlabel('Time (2022-09-15 HH:MM:SS UT)');ax.set_ylabel('Frequency (MHz)')
plt.show()

write_data=0
if(write_data):
    pa=24
    img=[0]*len(chlist);i=0;freq_mwa=[0]*len(chlist);Tb=[0]*len(chlist);S=[0]*len(chlist)
    bmin=[0]*len(chlist);bmaj=[0]*len(chlist)
    for c in chlist:
        print('Channel: ',c)
        imglist=sorted(glob.glob('/sdata/20220915_MWA/images/20220915_ch'+str(c)+'.pol.I.time.*.image_image.FITS'))[0:380]
        j=0;img[i]=[0]*len(imglist);Tb[i]=[0]*len(imglist);sumdata=[0]*len(imglist);S[i]=[0]*len(imglist)
        for im in imglist:
            aa=fits.open(im);data=aa[0].data[0][0][924:1124,924:1124]
            std=np.nanstd(aa[0].data[0][0]);data[np.where(data<5*std)]=np.nan
            img[i][j]=data;sumdata[j]=np.nansum(data)
            img_=img[i][j];img_[np.isnan(img_)]=0
            padX=[img_.shape[0] - 99, 99];padY=[img_.shape[1] - 99, 99]
            imgP = np.pad(img_, [padY, padX], 'constant')
            img[i][j] = rotate(imgP, angle=pa,reshape=False,mode='constant')
            img[i][j] = img[i][j][padY[0] : -padY[1], padX[0] : -padX[1]]
            j=j+1
        bmin[i]=aa[0].header['BMIN']*3600;bmaj[i]=aa[0].header['BMAJ']*3600
        sumdata=np.array(sumdata)
        minid=np.where(sumdata==np.nanmin(sumdata))[0][0]
        fact=flux[i]/np.nansum(img[i][0])*(freq[i]/108.)**2 # Flux Setting
        for k in range(len(imglist)):
            S[i][k]=img[i][k]*fact[k]
            Tb[i][k]=1224*S[i][k]*1.e7/(freq[i]**2*1.e-6*bmin[i]*bmaj[i])
            k=k+1
        freq_mwa[i] = np.round(aa[0].header['CRVAL3'] / 1.e6, 1)
        i=i+1
    Tb=np.array(Tb);S=np.array(S);img=np.array(img);Tbmax=np.nanmax(Tb,axis=(2,3))
    freq_mwa=np.array(freq_mwa);bmin=np.array(bmin);bmaj=np.array(bmaj)
    pickle.dump([freq_mwa,Tb,Tbmax,S,bmin,bmaj,img],open('/sdata/20220915_MWA/Tb_20220915.p','wb'))


#chan=np.arange(14)
freq_mwa,Tb,Tbmax,S,bmin,bmaj,img=pickle.load(open('/data/20220915_MWA/Tb_20220915.p','rb'))
Tb[:,:,0:70,:] = np.nan;Tb[:,:,140:,:]= np.nan
Tb[:,:,:,0:70] = np.nan;Tb[:,:,:,140:]= np.nan
#Tb[np.where(Tb<1.e3)]=np.nan#;Tb[np.where(Tb>1.e7)]=np.nan
#region1=img[:,:,105:130,105:130];region1_max=np.nanmax(region1,axis=(2,3))
#region1=img[:,:,105:120,110:135];region1_max=np.nanmax(region1,axis=(2,3))
region1=img[:,:,115:124,115:124];region1_max=np.nanmax(region1,axis=(2,3))
Tb_region1=Tb[:,:,115:124,115:124];Tb_region1_max=np.nanmax(Tb_region1,axis=(2,3));Tb_region1_mean=np.nanmean(Tb_region1,axis=(2,3))
region2=img[:,:,99:106,116:123];region2_max=np.nanmax(region2,axis=(2,3))
Tb_region2=Tb[:,:,99:106,116:123];Tb_region2_max=np.nanmax(Tb_region2,axis=(2,3));Tb_region2_mean=np.nanmean(Tb_region2,axis=(2,3))
Tb_region3=Tb[:,:,86:90,101:103];Tb_region3_max=np.nanmax(Tb_region3,axis=(2,3));Tb_region3_mean=np.nanmean(Tb_region3,axis=(2,3))


#region1=img;region1_max=np.nanmean(region1,axis=(2,3))
#region1_mean=np.nanmean(region1,axis=(2,3));Tb_region1_mean=np.nanmean(Tb_region1,axis=(2,3))
tmwa=12608+np.arange(Tb_region1_mean.shape[1])*10
freq_mwa=np.array(freq_mwa)
#for i in range(len(chlist)):
#    region1_mean[i]=region1_mean[i]-region1_mean[i].min()

Tb_pixel0=Tb[:,:,98,120];Tb_pixel0[Tb_pixel0<5.e3]=np.nan
f,ax=plt.subplots(1,1)
j=[34,36,37,38,40,41,42]
cl=['g','r','b','k','orange','brown','magenta'];al=[0.5,1.0,1.0,0.5,0.5,1.0,0.5]
ll=['03:35:40','03:36:00','03:36:10','03:36:20','03:36:40','03:36:50','03:37:00']
for i in range(7):
    ax.plot(freq[i],Tb_pixel0[:,j[i]]/1.e6,'o-',color=cl[i],alpha=al[i],label=ll[i])
ax.set_ylabel('$T_B$ (MK)');ax.set_xlabel('Frequency (MHz)')
ax.legend(loc=2);ax.set_ylim(0,10)
plt.show()

#---------------------------- new

#chlist_img=[61,65,69,77,81,86,91,96,101,107,113,120,127,134,142,150,167,177,187,210,226]
chlist_img=[61,65,69,81,86,91,96,101,107,120,127,134,142,150,167,177,187]
Tbmax_r1=[0]*len(chlist_img);Tbmax_r2=[0]*len(chlist_img)
Smax_r1=[0]*len(chlist_img);Smax_r2=[0]*len(chlist_img)
bmin=[0]*len(chlist_img);bmaj=[0]*len(chlist_img)
freq_mwa=[0]*len(chlist_img)
Tb_mean_region1=[0]*len(chlist_img);Tb_mean_region2=[0]*len(chlist_img)
Tb_r2_1=[0]*len(chlist_img);Tb_r1_1=[0]*len(chlist_img)
xc_max_reg1=[0]*len(chlist_img);xc_max_reg2=[0]*len(chlist_img)
yc_max_reg1=[0]*len(chlist_img);yc_max_reg2=[0]*len(chlist_img)
idx_r1=[0]*len(chlist_img);idx_r2=[0]*len(chlist_img)
reg='reg1'
for i in range(len(chlist_img)):
    c=chlist_img[i]
    print('Channel: '+str(c))
    freq_mwa[i],Tb,Tbmax_r1[i],S,bmin[i],bmaj[i],img=pickle.load(open('/data/20220915_MWA/Tb_20220915_'+str(reg)+'_ch'+str(c)+'.p','rb'),encoding='latin1')
    print(Tb.shape)
    Smax_r1[i]=np.nanmax(S,axis=(1,2))
    Tb_r1 = Tb[:, 115:125,112:122]
    #Tb_r1 = Tb[:,119,116]
    Tb_r1_1[i]=Tb_r1*0;xc_max_reg1[i]=[0]*len(Tb_r1);yc_max_reg1[i]=[0]*len(Tb_r1);idx_r1[i]=[0]*len(Tb_r1)
    for k in range(len(Tb_r1)):
        if(any(map(lambda x: x == 0, np.array(np.where(Tb_r1[k]==np.nanmax(Tb_r1[k]))).flatten()))==False):
            Tb_r1_1[i][k]=Tb_r1[k]
            xc_max_reg1[i][k],yc_max_reg1[i][k]=np.where(Tb[k][115:125,112:122]==np.nanmax(Tb[k][115:125,112:122]))
        else:
            Tbmax_r1[i][k] = np.nan
            Tb_r1_1[i][k] = np.nan
    Tbmax_r1[i][Tbmax_r1[i]>1.e10] = np.nan;Tbmax_r1[i][Tbmax_r1[i]<1.e5] = np.nan
    Tb_r1_1[i][Tb_r1_1[i]>1.e10] = np.nan; Tb_r1_1[i][Tb_r1_1[i]<1.e5] = np.nan
    idx_r1[i]=np.arange(len(Tbmax_r1[i]))[np.isfinite(Tbmax_r1[i])]
    #Tb_mean_region1[i]=np.nanmean(Tb_r1_1[i][idx_r1[i]],axis=(1,2))
    Tb_mean_region1[i]=Tb[:, 119,116]

reg='reg2'
for i in range(len(chlist_img)):
    c=chlist_img[i]
    print('Channel: '+str(c))
    freq_mwa[i],Tb,Tbmax_r2[i],S,bmin[i],bmaj[i],img=pickle.load(open('/data/20220915_MWA/Tb_20220915_'+str(reg)+'_ch'+str(c)+'.p','rb'),encoding='latin1')
    Smax_r2[i]=np.nanmax(S,axis=(1,2))
    Tb_r2=Tb[:,99:106,116:123]
    #Tb_r2=Tb[:,102,120]
    Tb_r2_1[i]=Tb_r2*0;idx_r2[i]=[0]*len(Tb_r2);xc_max_reg2[i]=[0]*len(Tb_r2);yc_max_reg2[i]=[0]*len(Tb_r2)
    for k in range(len(Tb_r2)):
        if(any(map(lambda x: x == 0, np.array(np.where(Tb_r2[k]==np.nanmax(Tb_r2[k]))).flatten()))==False):
            Tb_r2_1[i][k]=Tb_r2[k]
            xc_max_reg2[i][k],yc_max_reg2[i][k]=np.where(Tb[k][115:125,112:122]==np.nanmax(Tb[k][115:125,112:122]))
            idx_r2[i][k]=k
    Tbmax_r2[i][Tbmax_r2[i]>1.e10] = np.nan;Tbmax_r2[i][Tbmax_r2[i]<1.e5] = np.nan
    Tb_r2_1[i][Tb_r2_1[i]>1.e10] = np.nan; Tb_r2_1[i][Tb_r2_1[i]<1.e5] = np.nan
    idx_r2[i]=np.arange(len(Tbmax_r2[i]))[np.isfinite(Tbmax_r2[i])]
    #Tb_mean_region2[i]=np.nanmean(Tb_r2_1[i][idx_r2[i]],axis=(1,2))
    Tb_mean_region2[i]=Tb[:,102,120]

freq_mwa=np.array(freq_mwa)
Tbmax_r1=np.array(Tbmax_r1);Tbmax_r2=np.array(Tbmax_r2)
Smax_r1=np.array(Smax_r1);Smax_r2=np.array(Smax_r2)
bmin=np.array(bmin);bmaj=np.array(bmaj)
Tb_mean_region2=np.array(Tb_mean_region2);Tb_mean_region1=np.array(Tb_mean_region1)
Tb_r1_1=np.array(Tb_r1_1);Tb_r2_1=np.array(Tb_r2_1)
pickle.dump([freq_mwa,Tbmax_r1,Tbmax_r2,Tb_mean_region1,Tb_mean_region2,Smax_r1, \
             Smax_r2,bmin,bmaj,Tb_r1_1,Tb_r2_1,idx_r1,idx_r2],open('/data/20220915_MWA/20220915_Tbmax.p','wb'))

#freq_mwa,Tbmax_reg1,Smax_reg1,bmin,bmaj = pickle.load(open('/data/20220915_MWA/20220915_reg1_Tbmax.p','rb'),encoding='latin1')
#freq_mwa,Tbmax_reg2,Smax_reg2,bmin,bmaj = pickle.load(open('/data/20220915_MWA/20220915_reg2_Tbmax.p','rb'),encoding='latin1')
freq_mwa,Tbmax_r1,Tbmax_r2,Tb_mean_region1,Tb_mean_region2,Smax_r1,Smax_r2,bmin,bmaj,Tb_r1_1,Tb_r2_1,idx_r1,idx_r2= \
pickle.load(open('/data/20220915_MWA/20220915_Tbmax.p','rb'),encoding='latin1')

#----------------- HMI
hmifile=sorted(glob.glob('/sdata/20220915_hmi/hmi.m_45s.2022.09.12_03_01_30_TAI.magnetogram.fits'))
hmimap=Map(hmifile)
hmimap.plot(vmin=-0.1,vmax=0.1)
plt.show()

#---------------------- STIX X-ray

stix=readsav('/data/Dropbox/STIX-MWA/20220915/STIX_data/stix_lightcurves_spec_8-s_Earth-UT_20220915.sav')
stix_data=stix['stix_lcstr']
stix_data0=stix_data['data'][0]
tstix_low=stix_data['ut'][0][:,0]-stix_data['ut'][0][:,0][0]+ 11022 #Start time 03:03:42
tstix_high=stix_data['ut'][0][:,1]-stix_data['ut'][0][:,1][0]+11041 #Start time 03:04:01

#plt.imshow(region1_max,aspect='auto',origin='lower',interpolation=None)

#plt.imshow(img[0].mean(axis=0),aspect='auto',origin='lower',interpolation=None)
Tb_ds=np.mean(Tb[:,:,120:150,140:170],axis=(2,3))
for i in range(5):
    Tb_ds[i]=Tb_ds[i]#-Tb_ds[i][0]
plt.imshow(Tb_ds,origin='lower',aspect='auto',cmap='jet',vmin=1.e3,vmax=1.e5)
plt.show()


f,ax=plt.subplots(1,1)
ax1=ax.twinx()
c=['k','g','magenta','cyan','brown'];#fact=[7,1,0.5,0.5,0.4,1,1,1,1,1,1,1,1]
i=3
ax1.plot(tmwa, Tb_ds[i] * 5.e-4, 'o-', linewidth=1, markersize=1,label='MWA (' + str(freq_mwa[i]) + ' MHz)')
#ax1.plot(tmwa,region1_mean[i]*fact[i],'o-',linewidth=1,markersize=1,label='MWA ('+str(freq_mwa[i])+' MHz)',color=c[i])
#ax1.plot(tmwa,(region1_mean[2]-17)*0.2,'o-',linewidth=1,markersize=3,label='MWA ('+str(freq_mwa[2])+' MHz)',color='green')
#ax1.plot(tmwa,region1_mean[0:].mean(axis=0),'o-',linewidth=1,markersize=3,label='MWA ($\\nu$-Average)',color='k')
ax.plot(tstix_low,stix_data0[:,0],'-',label='STIX (6-7 keV)',color='r')
ax.plot(tstix_high,stix_data0[:,1],'-',label='STIX (16-22 keV)',color='b')
ax.set_xticks(tstix_low[::50])
ax.set_xticklabels(['03:03:42','03:14:12','03:21:14','03:28:16','03:36:20','03:44:58','03:55:11','04:03:29','04:14:13'])
ax.set_ylabel('Amplitude');ax.set_xlabel('Time (HH:MM:SS)')
ax.legend();ax1.legend(loc=4);ax.set_yscale('log');ax.set_ylabel('STIX Count Flux (counts/s/cm$^2$/keV)')
ax.set_xlabel('Time (HH:MM:SS UT) (Start Time: 03:03:42 UT)');ax1.set_ylabel('Amplitude')
#ax1.set_ylim(0.5,4)
plt.show()

#Tb_region1_max[:,272:273] = np.nan

mway=np.nanmax(Tb_region2_mean[0:2],axis=0)*1.e-5
f,ax=plt.subplots(1,1)
ax1=ax.twinx()
c=['k','g','magenta','cyan','brown'];fact=[7,1,0.5,0.5,0.4]
i=5
ax1.plot(tmwa,mway,'o-',linewidth=1,markersize=1,label='MWA (<'+str(freq_mwa[i])+' MHz) | Harmonic Height'+str(np.round(height_Mm_harmonic[i],0))+' Mm',color=c[0])
#ax1.plot(tmwa,(region1_mean[2]-17)*0.2,'o-',linewidth=1,markersize=3,label='MWA ('+str(freq_mwa[2])+' MHz)',color='green')
#ax1.plot(tmwa,region1_mean[0:].mean(axis=0),'o-',linewidth=1,markersize=3,label='MWA ($\\nu$-Average)',color='k')
ax.plot(tstix_low,stix_data0[:,0],'-',label='STIX (6-7 keV)',color='r')
ax.plot(tstix_high,stix_data0[:,1],'-',label='STIX (16-22 keV)',color='b')
ax.set_xticks(tstix_low[::50])
ax.set_xticklabels(['03:03:42','03:14:12','03:21:14','03:28:16','03:36:20','03:44:58','03:55:11','04:03:29','04:14:13'])
ax1.set_ylabel('T$_B$ ($\\times 10^5$ K)');ax.set_xlabel('Time (HH:MM:SS)')
ax.legend();ax1.legend(loc=2);ax.set_yscale('log');ax.set_ylabel('STIX Count Flux (counts/s/cm$^2$/keV)')
ax.set_xlabel('Time (HH:MM:SS UT) (Start Time: 03:03:42 UT)')
ax1.set_ylim(0.1,0.6) # Chan 0
plt.show()

mway=np.nanmax(Tb_region2_mean[15:],axis=0)*1.e-5
f,ax=plt.subplots(1,1)
ax1=ax.twinx()
c=['k','g','magenta','cyan','brown'];fact=[7,1,0.5,0.5,0.4]
i=15
ax1.plot(tmwa,mway,'o-',linewidth=1,markersize=1,label='MWA (>'+str(freq_mwa[i])+' MHz) | Harmonic Height'+str(np.round(height_Mm_harmonic[i],0))+' Mm',color=c[0])
#ax1.plot(tmwa,(region1_mean[2]-17)*0.2,'o-',linewidth=1,markersize=3,label='MWA ('+str(freq_mwa[2])+' MHz)',color='green')
#ax1.plot(tmwa,region1_mean[0:].mean(axis=0),'o-',linewidth=1,markersize=3,label='MWA ($\\nu$-Average)',color='k')
ax.plot(tstix_low,stix_data0[:,0],'-',label='STIX (6-7 keV)',color='r')
ax.plot(tstix_high,stix_data0[:,1],'-',label='STIX (16-22 keV)',color='b')
ax.set_xticks(tstix_low[::50])
ax.set_xticklabels(['03:03:42','03:14:12','03:21:14','03:28:16','03:36:20','03:44:58','03:55:11','04:03:29','04:14:13'])
ax1.set_ylabel('T$_B$ ($\\times 10^5$ K)');ax.set_xlabel('Time (HH:MM:SS)')
ax.legend();ax1.legend(loc=2);ax.set_yscale('log');ax.set_ylabel('STIX Count Flux (counts/s/cm$^2$/keV)')
ax.set_xlabel('Time (HH:MM:SS UT) (Start Time: 03:03:42 UT)')
ax1.set_ylim(0.25,10.5) # Chan 0
plt.show()


#r2_y=Tb_region2_max[:,14:16].mean(axis=1);r1_y=Tb_region1_max[:,258:261].mean(axis=1)
r2_y=Tb_region2_max[:,15];r1_y=Tb_region1_max[:,259]
r1_y0=Tb_region1_max[:,20];r2_y0=Tb_region2_max[:,0]
r2_y[10]=np.nan;r2_y[11]=np.nan;r2_y[-1]=np.nan;r1_y[-1]=np.nan
diff_r1=r1_y-r1_y0;diff_r2=r2_y-r2_y0

f,ax=plt.subplots(2,1);ax0=ax[0];ax1=ax[1]
ax0.plot(freq[0:5],diff_r1[0:5]/1.e3,'o-',label='Region 1 (04:16:30)')
ax1.plot(freq[13:],diff_r2[13:]/1.e6,'o-',label='Region 2 (03:32:30)')
ax0.set_xlabel('Frequency (MHz)');ax1.set_xlabel('Frequency (MHz)')
ax0.set_ylabel('$T_B$ (kK)');ax1.set_ylabel('$T_B$ (MK)')
ax0.legend();ax1.legend()
plt.show()

f,ax=plt.subplots(1,1)
ax.plot(freq,diff_r2/1.e6,'o-',label='Region 2')
ax.plot(freq,diff_r1/1.e6,'o-',label='Region 1')
ax.legend()
ax.set_xlabel('Frequency (MHz)');ax.set_ylabel('$T_B$ (MK)')
ax.set_ylim(0,0.6)
plt.show()

f,ax=plt.subplots(1,1)
ax.plot(freq,r2_y/1.e6,'o-',label='Region 2 (03:32:30)')
ax.plot(freq,r1_y/1.e6,'o-',label='Region 1 (04:16:30)')
#ax.plot(freq,r2_y0/1.e6,'o-',label='Region 2 (03:30:10)')
#ax.plot(freq,r1_y0/1.e6,'o-',label='Region 1 (03:33:20)')
ax.legend()
ax.set_xlabel('Frequency (MHz)');ax.set_ylabel('$T_B$ (MK)')
ax.set_ylim(1.e-3,10);ax.set_yscale('log')
plt.show()

plt.plot(region1_max[0],'o-',label='Channel '+str(chan[i]))
plt.legend()
plt.show()

region1_max1=region1_max*1.0
for j in range(5):
    region1_max1[j]=region1_max[j]-region1_max[j][0]

f,ax=plt.subplots(1,1)
ax.imshow(region1_max1,aspect='auto',origin='lower',cmap='jet',interpolation=None,vmin=0,vmax=3)
ax.set_xticks([0,50,100,150,200,250,300,350])
ax.set_xticklabels(['03:30:00','03:38:20','03:46:40','03:55:00','04:03:20','04:11:40','04:20:00','04:28:20'])
ax.set_yticks(np.arange(5));ax.set_yticklabels(freq[chan])
plt.show()

#euv094_list=sorted(glob.glob('/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220830/20220830_EUV/*.94A*.fits'))
#euv171_list=sorted(glob.glob('/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220830/20220830_EUV/*.171A*.fits'))

#------------------ AIA Analysis
euv094_list=sorted(glob.glob('/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220915_EUV/*.94A*.fits'))
euv171_list=sorted(glob.glob('/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220915_EUV/*.171A*.fits'))


aiamap171_diff=[0]*(len(euv171_list)-10)
for i in range(len(aiamap171_diff)):
    aia1=Map(euv171_list[i+10])
    aia_ref=Map(euv171_list[i])
    aiamap171_diff[i]=Map(aia1.data-aia_ref.data,aia1.meta)

aiamap094_diff=[0]*(len(euv094_list)-10)
for i in range(len(aiamap094_diff)):
    aia1=Map(euv094_list[i+10])
    aia_ref=Map(euv094_list[i])
    aiamap094_diff[i]=Map(aia1.data-aia_ref.data,aia1.meta)

for i in range(len(aiamap171_diff)):
    ii="%04d" % i
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection=aia1)
    p0=aiamap171_diff[i].plot(axes=ax,aspect='auto',vmin=-400,vmax=800,cmap='coolwarm')
    ax.set_ylim(2200,3000);ax.set_xlim(3200,4000)
    plt.savefig('/home/rohit/20220915/20220915_pngs/aia171_'+str(ii)+'.png',dpi=60)
    plt.close()

for i in range(250):
    ii = "%04d" % i
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection=aia1)
    p0 = aiamap094_diff[i].plot(axes=ax, aspect='auto', vmin=-100, vmax=100, cmap='coolwarm')
    ax.set_ylim(2200, 3000);
    ax.set_xlim(3200, 4000)
    plt.savefig('/home/rohit/20220915/20220915_pngs/aia094_' + str(ii) + '.png', dpi=60)
    plt.close()


plt.show()

time_mwa=np.arange(380)*10+3*3600+30*60
i=0;time_aia171=[0]*len(euv171_list)
for e in euv171_list:
    l=euv171_list[i].split('T')[2].split('Z')[0].split('_')
    time_aia171[i]=int(l[0])*3600+int(l[1])*60+float(l[2])
    i=i+1
time_aia171=np.array(time_aia171)


euv094_file=euv094_list[0]
map94=Map(euv094_file)


for j in range(30,50):
    idx = ut.find_predecessor(time_aia171, time_mwa[j])
    euv171_file = euv171_list[idx[0]]
    map171 = Map(euv171_file)
    aiahead1 = []
    aiahead1 = map171.meta.copy()
    aiahead1['naxis1'] = 200
    aiahead1['naxis2'] = 200
    aiahead1['CRPIX1'] = 99
    aiahead1['CRPIX2'] = 99
    aiahead1['CRVAL1'] = 0
    aiahead1['CRVAL2'] = 0
    aiahead1['CDELT1'] = 50
    aiahead1['CDELT2'] = 50
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection=map171)
    #lev_ct = [[1.0e5, 2.0e5, 3.0e5], [3.0e4,4.0e4], [2.0e4, 3.0e4], [3.e5, 5.e5], [2.e4, 3.e5]]
    lev_ct = [[1.e6], [1.0e5], [3.0e4,5.e4], [3.e5], [3.e5]]
    cc = ['r', 'g', 'b', 'k', 'brown']
    for k in range(5):
        mwadata=Tb[k][j]
        mwadata1=mwadata*0;mwadata1[90:160,90:160] = mwadata[90:160,90:160]
        mwadata1[mwadata1<0.5] = 0
        mwamap = Map(mwadata, aiahead1);mwamap1 = Map(mwadata1, aiahead1)
        ii="%04d" % j
        p0=map171.plot(axes=ax,aspect='auto',vmin=0,vmax=800)
        #frac_r1 = [0.05,0.1,0.12,0.2,0.3,0.4,0.5]
        for f in lev_ct[k]:
            #lev_r1 = np.nanmax(mwamap.data) * f
            c1 = mwamap1.contour(level=f * u.ct)
            if(len(c1)!=0):
                ax.plot_coord(c1[0], color=cc[k])
        ax.text(-1400, -700+k*500, 'Freq: '+str(np.round(freq[chan[k]],1))+' MHz', style='italic', bbox={'facecolor': cc[k], 'alpha': 0.5, 'pad': 10})
    ax.set_ylim(-1600,5000);ax.set_xlim(-1600,5000)
    plt.savefig('/home/rohit/20220915/20220915_pngs/img_'+str(ii)+'.png')
    plt.close()

for j in range(len(time_mwa)):
    idx = ut.find_predecessor(time_aia171, time_mwa[j])
    euv171_file = euv171_list[idx[0]]
    map171 = Map(euv171_file)
    aiahead1 = []
    aiahead1 = map171.meta.copy()
    aiahead1['naxis1'] = 200
    aiahead1['naxis2'] = 200
    aiahead1['CRPIX1'] = 99
    aiahead1['CRPIX2'] = 99
    aiahead1['CRVAL1'] = 0
    aiahead1['CRVAL2'] = 0
    aiahead1['CDELT1'] = 50
    aiahead1['CDELT2'] = 50
    fig = plt.figure(figsize=(8, 20))
    ax0 = fig.add_subplot(211, projection=map171)
    ax1 = fig.add_subplot(212)
    lev_ct = [[1.e5], [1.0e5], [3.0e4,5.e4], [3.e5], [3.e5]]
    cc = ['r', 'g', 'b', 'k', 'brown']
    for k in range(5):
        mwadata=Tb[k][j]
        mwadata1=mwadata*0;mwadata1[90:160,110:170] = mwadata[90:160,110:170]
        mwadata1[mwadata1<0.5] = 0
        mwamap = Map(mwadata, aiahead1);mwamap1 = Map(mwadata1, aiahead1)
        ii="%04d" % j
        p0=map171.plot(axes=ax0,aspect='auto',vmin=0,vmax=800)
        #frac_r1 = [0.05,0.1,0.12,0.2,0.3,0.4,0.5]
        for f in lev_ct[k]:
            #lev_r1 = np.nanmax(mwamap.data) * f
            c1 = mwamap1.contour(level=f * u.ct)
            if(len(c1)!=0):
                ax0.plot_coord(c1[0], color=cc[k])
        ax0.text(5400, -700+k*500, 'Freq: '+str(np.round(freq[chan[k]],1))+' MHz', style='italic', bbox={'facecolor': cc[k], 'alpha': 0.5, 'pad': 10})
    ax0.set_ylim(-1600,7000);ax0.set_xlim(-1600,7000)
    idx1=ut.find_predecessor(tstix_high,time_mwa[j])[0]
    ax1.plot(tstix_high[idx1-100:idx1+100]-tstix_high[0],stix_data0[:,1][idx1-100:idx1+100],'o-',label='STIX (16-22 keV)',color='b')
    #ax1.set_xticklabels(['03:03:42','03:14:12','03:21:14','03:28:16','03:36:20','03:44:58','03:55:11','04:03:29','04:14:13'])
    ax1.axvline(x=time_mwa[j]-tstix_high[0],color='k')
    ax1.set_ylabel('STIX Count Flux (counts/s/cm$^2$/keV)');ax1.set_xlabel('Time (HH:MM:SS) / Start time: 03:03:42')
    ax1.legend()
    plt.savefig('/home/rohit/20220915/20220915_pngs/img_ts_'+str(ii)+'.png')
    plt.close()



idx = ut.find_predecessor(time_aia171, time_mwa[j])
euv171_file = euv171_list[idx[0]]
map171 = Map(euv171_file)
aiahead1 = []
aiahead1 = map171.meta.copy()
aiahead1['naxis1'] = 200
aiahead1['naxis2'] = 200
aiahead1['CRPIX1'] = 99
aiahead1['CRPIX2'] = 99
aiahead1['CRVAL1'] = 0
aiahead1['CRVAL2'] = 0
aiahead1['CDELT1'] = 50
aiahead1['CDELT2'] = 50
xlmwa,xrmwa,ylmwa,yrmwa=-5000,5000,-5000,5000
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection=map171)
lev_ct = [[7.0e4, 1.0e5], [1.0e5,2.0e5,3.0e5], [1.0e5, 3.0e5], [8.e4, 3.e5], [2.e3,8.e4, 1.e5, 1.5e5,2.e5,3.e6]]
cc = ['r', 'k', 'b', 'g', 'brown']
for k in range(5):
    k1=k+13
    mwadata=np.nanmean(Tb[k1],axis=0)
    mwadata[np.isnan(mwadata)]=0
    #mwadata1=mwadata*0;mwadata1[90:160,90:160] = mwadata[90:160,90:160]
    #mwadata1[mwadata1<0.5] = 0
    mwamap = Map(mwadata, aiahead1)#;mwamap1 = Map(mwadata1, aiahead1)
    ii="%04d" % j
    p0=map171.plot(axes=ax,aspect='auto',vmin=0,vmax=800)
    lev1=np.array([20,40])*u.percent
    mwamap.draw_contours(axes=ax,levels=lev1,colors=cc[k],linewidths=3,extent=[xlmwa,xrmwa,ylmwa,yrmwa])
    #frac_r1 = [0.05,0.1,0.12,0.2,0.3,0.4,0.5]
    #for f in lev_ct[k]:
    #    #lev_r1 = np.nanmax(mwamap.data) * f
    #    c1 = mwamap.contour(level=f * u.ct)
    #    if(len(c1)!=0):
    #        ax.plot_coord(c1[0], color=cc[k])
    ax.text(-1400, -700+k*500, 'Freq: '+str(np.round(freq[k1],1))+' MHz', style='italic', bbox={'facecolor': cc[k], 'alpha': 0.5, 'pad': 10})
ax.set_ylim(-1600,7000);ax.set_xlim(-1600,7000)
plt.show()


j=36
idx = ut.find_predecessor(time_aia171, time_mwa[j])
euv171_file = euv171_list[idx[0]]
map171 = Map(euv171_file)
aiahead1 = []
aiahead1 = map171.meta.copy()
aiahead1['naxis1'] = 200
aiahead1['naxis2'] = 200
aiahead1['CRPIX1'] = 99
aiahead1['CRPIX2'] = 99
aiahead1['CRVAL1'] = 0
aiahead1['CRVAL2'] = 0
aiahead1['CDELT1'] = 50
aiahead1['CDELT2'] = 50
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection=map171)
lev_ct =[[3.e4], [3.0e4], [3.0e4,5.e4], [3.e5,8.e5], [3.e5,8.e5]]
#lev_ct = [[1.0e5, 0.5e5, 1.0e5], [2.0e4,3.0e4], [2.0e4, 3.0e4], [2.e5, 3.e5], [2.e3, 3.e5, 5.e5]]
cc = ['r', 'g', 'b', 'k', 'brown']
for k in range(5):
    mwadata=np.nanmean(Tb[k][j-2:j+2],axis=0)
    mwadata1=mwadata*0;mwadata1[90:160,90:160] = mwadata[90:160,90:160]
    mwadata1[mwadata1<0.5] = 0
    mwamap = Map(mwadata, aiahead1);mwamap1 = Map(mwadata1, aiahead1)
    ii="%04d" % j
    p0=map171.plot(axes=ax,aspect='auto',vmin=0,vmax=800)
    #frac_r1 = [0.05,0.1,0.12,0.2,0.3,0.4,0.5]
    for f in lev_ct[k]:
        #lev_r1 = np.nanmax(mwamap.data) * f
        c1 = mwamap1.contour(level=f * u.ct)
        if(len(c1)!=0):
            ax.plot_coord(c1[0], color=cc[k])
    ax.text(-1400, -700+k*500, 'Freq: '+str(np.round(freq[chan[k]],1))+' MHz', style='italic', bbox={'facecolor': cc[k], 'alpha': 0.5, 'pad': 10})
ax.set_ylim(-1600,7000);ax.set_xlim(-1600,7000)
plt.show()

level=[0.6]
cc = ['r', 'g', 'b', 'k', 'brown']
f,ax=plt.subplots(1,1)
for i in range(5):
    ax.contour(Tb[i].mean(axis=0)/np.nanmax(Tb[i].mean(axis=0)),levels=level,colors=cc[i],origin='lower')
pc=plt.Circle((130,130),19.2,fill=False)
ax.add_artist( pc )
plt.show()
