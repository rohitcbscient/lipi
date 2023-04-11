import pickle
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from sunpy.map import Map
import itertools
import astropy.units as u
from scipy.io import readsav
from surya.utils import main as ut
from scipy.ndimage.interpolation import rotate
from surya.utils import model as mdl

baseline=['000-002','000-005','000-007','002-005','002-007','005-007']
ds=[0]*len(baseline);i=0;flux=[0]*len(baseline)
for b in baseline:
    f1=pickle.load(open('/home/rohit/20220915/20220915_sun_ms_T'+b+'.DS.dat.p','rb'),encoding='latin1')
    ds[i]=f1[5].reshape(4,48,390)[:,:,0:380]
    f2=pickle.load(open('/home/rohit/20220915/flux_V1_20220915_sun_ms_T'+b+'.p','rb'),encoding='latin1')
    flux[i]=f2[17][3][0][0:380]
    i=i+1
freq=f1[0]-15.36
meands=np.nanmean(np.array(ds),axis=0);flux=np.nanmean(np.array(flux),axis=0)
#------- Newkirk Model
height_Mm = [0]*len(freq)
for i in range(len(freq)):
    height_Mm[i]=(float(mdl.nk_freq2r(freq[i],1)[0])-1)*6.99e2

chan=[0,12,22,38,46]
f,ax=plt.subplots(1,1)
ax.imshow(meands[0],aspect='auto',origin='lower',interpolation=None,vmin=0.1,vmax=1)
#f.colorbar(label='NCCF')
ax.set_xticks([0,50,100,150,200,250,300,350])
ax.set_xticklabels(['03:29:58','03:38:18','03:46:38','03:54:58','04:03:18','04:11:38','04:19:58','04:28:18'])
ax.set_yticks(np.arange(48))
ax.set_yticklabels(np.round(freq,2));ax.set_xlabel('Time (2022-09-15 HH:MM:SS UT)');ax.set_ylabel('Frequency (MHz)')
for c in chan:
    ax.axhline(y=c)
plt.show()

write_data=0
if(write_data):
    pa=-24
    img=[0]*len(chan);i=0;freq_mwa=[0]*len(chan);Tb=[0]*len(chan);S=[0]*len(chan);bmin=[0]*len(chan);bmaj=[0]*len(chan)
    for c in chan:
        print(c)
        imglist=sorted(glob.glob('/home/rohit/20220915/20220915_sun.chan.'+str(c)+'.pol.I.time.*.image_image.FITS'))[0:380]
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
    pickle.dump([freq_mwa,Tb,Tbmax,S,bmin,bmaj,img],open('/home/rohit/20220915/Tb_20220915.p','wb'))

#chan=np.arange(14)
freq_mwa,Tb,Tbmax,S,bmin,bmaj,img=pickle.load(open('/home/rohit/20220915/Tb_20220915.p','rb'))
img[np.where(img<1.e-5)]=np.nan
#region1=img[:,:,105:130,105:130];region1_max=np.nanmax(region1,axis=(2,3))
#region1=img[:,:,105:120,110:135];region1_max=np.nanmax(region1,axis=(2,3))
region1=img[:,:,120:150,140:170];region1_max=np.nanmax(region1,axis=(2,3))
Tb_region1=Tb[:,:,120:150,140:170];Tb_region1_max=np.nanmax(Tb_region1,axis=(2,3))
#region1=img;region1_max=np.nanmean(region1,axis=(2,3))
region1_mean=np.nanmean(region1,axis=(2,3));Tb_region1_mean=np.nanmean(Tb_region1,axis=(2,3))
tmwa=12608+np.arange(region1_mean.shape[1])*10
freq_mwa=np.array(freq_mwa)
for i in range(len(chan)):
    region1_mean[i]=region1_mean[i]-region1_mean[i].min()

Tb_pixel0=Tb[:,:,98,120];Tb_pixel0[Tb_pixel0<5.e3]=np.nan
f,ax=plt.subplots(1,1)
j=[34,36,37,38,40,41,42]
cl=['g','r','b','k','orange','brown','magenta'];al=[0.5,1.0,1.0,0.5,0.5,1.0,0.5]
ll=['03:35:40','03:36:00','03:36:10','03:36:20','03:36:40','03:36:50','03:37:00']
for i in range(7):
    ax.plot(freq[chan],Tb_pixel0[:,j[i]]/1.e6,'o-',color=cl[i],alpha=al[i],label=ll[i])
ax.set_ylabel('$T_B$ (MK)');ax.set_xlabel('Frequency (MHz)')
ax.legend(loc=2);ax.set_ylim(0,10)
plt.show()


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

f,ax=plt.subplots(1,1)
ax1=ax.twinx()
c=['k','g','magenta','cyan','brown'];fact=[7,1,0.5,0.5,0.4]
i=2
ax1.plot(tmwa,Tb_region1_max[i]*1.e-5,'o-',linewidth=1,markersize=1,label='MWA ('+str(freq_mwa[i])+' MHz) | Height'+str(np.round(height_Mm[i],0))+' Mm',color=c[i])
#ax1.plot(tmwa,(region1_mean[2]-17)*0.2,'o-',linewidth=1,markersize=3,label='MWA ('+str(freq_mwa[2])+' MHz)',color='green')
#ax1.plot(tmwa,region1_mean[0:].mean(axis=0),'o-',linewidth=1,markersize=3,label='MWA ($\\nu$-Average)',color='k')
ax.plot(tstix_low,stix_data0[:,0],'-',label='STIX (6-7 keV)',color='r')
ax.plot(tstix_high,stix_data0[:,1],'-',label='STIX (16-22 keV)',color='b')
ax.set_xticks(tstix_low[::50])
ax.set_xticklabels(['03:03:42','03:14:12','03:21:14','03:28:16','03:36:20','03:44:58','03:55:11','04:03:29','04:14:13'])
ax1.set_ylabel('T$_B$ ($\\times 10^5$ K)');ax.set_xlabel('Time (HH:MM:SS)')
ax.legend();ax1.legend(loc=3);ax.set_yscale('log');ax.set_ylabel('STIX Count Flux (counts/s/cm$^2$/keV)')
ax.set_xlabel('Time (HH:MM:SS UT) (Start Time: 03:03:42 UT)')
#ax1.set_ylim(0.25,0.5) # Chan 0
#ax1.set_ylim(2.6,4.0) # Chan 1
#ax1.set_ylim(1.6,2.5) # Chan 1
#ax1.set_ylim(0.8,1.3) # Chan 1
#ax1.set_ylim(3.6,5.2) # Chan 2
#ax1.set_ylim(4.5,6.5) # Chan 3
#ax1.set_ylim(5.2,7.8) # Chan 4
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


j=272
idx = ut.find_predecessor(time_aia171, time_mwa[j])
euv171_file = euv171_list[idx[0]]
map171 = Map(euv171_file)
aiahead1 = []
aiahead1 = map171.meta.copy()
aiahead1['naxis1'] = 264
aiahead1['naxis2'] = 264
aiahead1['CRPIX1'] = 132
aiahead1['CRPIX2'] = 132
aiahead1['CRVAL1'] = 0
aiahead1['CRVAL2'] = 0
aiahead1['CDELT1'] = 50
aiahead1['CDELT2'] = 50
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection=map171)
lev_ct = [[2.0e5, 3.0e5], [2.0e5,3.0e5], [2.0e5, 3.0e5], [2.e5, 3.e5], [2.e3, 3.e5]]
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
