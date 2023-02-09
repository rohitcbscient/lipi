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

#chan=[0,12,22,38,46]
chan=np.arange(14)
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
    img=[0]*len(chan);i=0;freq_mwa=[0]*len(chan);Tb=[0]*len(chan);S=[0]*len(chan);bmin=[0]*len(chan);bmaj=[0]*len(chan)
    for c in chan:
        print(c)
        imglist=sorted(glob.glob('/home/rohit/20220915/20220915_sun.chan.'+str(c)+'.pol.I.time.*.image_image.FITS'))[0:380]
        j=0;img[i]=[0]*len(imglist);Tb[i]=[0]*len(imglist);sumdata=[0]*len(imglist);S[i]=[0]*len(imglist)
        for im in imglist:
            aa=fits.open(im);data=aa[0].data[0][0][924:1124,924:1124]
            std=np.nanstd(aa[0].data[0][0]);data[np.where(data<5*std)]=np.nan
            img[i][j]=data;sumdata[j]=np.nansum(data)
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


freq_mwa,Tb,Tbmax,S,bmin,bmaj,img=pickle.load(open('/home/rohit/20220915/Tb_20220915.p','rb'))
img[np.where(img<1.e-5)]=np.nan
#region1=img[:,:,105:130,105:130];region1_max=np.nanmax(region1,axis=(2,3))
#region1=img[:,:,105:120,110:135];region1_max=np.nanmax(region1,axis=(2,3))
region1=img[:,:,90:135,95:150];region1_max=np.nanmax(region1,axis=(2,3))
Tb_region1=Tb[:,:,90:135,95:150];Tb_region1_max=np.nanmax(Tb_region1,axis=(2,3))
#region1=img;region1_max=np.nanmean(region1,axis=(2,3))
region1_mean=np.nanmean(region1,axis=(2,3));Tb_region1_mean=np.nanmean(Tb_region1,axis=(2,3))
tmwa=12608+np.arange(region1_mean.shape[1])*10
freq_mwa=np.array(freq_mwa)
for i in range(len(chan)):
    region1_mean[i]=region1_mean[i]-region1_mean[i].min()

#---------------------- STIX X-ray

stix=readsav('/data/Dropbox/STIX-MWA/20220915/STIX_data/stix_lightcurves_spec_8-s_Earth-UT_20220915.sav')
stix_data=stix['stix_lcstr']
stix_data0=stix_data['data'][0]
tstix_low=stix_data['ut'][0][:,0]-stix_data['ut'][0][:,0][0]+ 11022 #Start time 03:03:42
tstix_high=stix_data['ut'][0][:,1]-stix_data['ut'][0][:,1][0]+11041 #Start time 03:04:01

#plt.imshow(region1_max,aspect='auto',origin='lower',interpolation=None)

#plt.imshow(img[0].mean(axis=0),aspect='auto',origin='lower',interpolation=None)


f,ax=plt.subplots(1,1)
ax1=ax.twinx()
c=['k','g','magenta','cyan','brown'];fact=[7,1,0.5,0.5,0.4]
for i in range(len(chan)):
    ax1.plot(tmwa,region1_mean[i]*fact[i],'o-',linewidth=1,markersize=1,label='MWA ('+str(freq_mwa[i])+' MHz)',color=c[i])
#ax1.plot(tmwa,(region1_mean[2]-17)*0.2,'o-',linewidth=1,markersize=3,label='MWA ('+str(freq_mwa[2])+' MHz)',color='green')
#ax1.plot(tmwa,region1_mean[0:].mean(axis=0),'o-',linewidth=1,markersize=3,label='MWA ($\\nu$-Average)',color='k')
ax.plot(tstix_low,stix_data0[:,0],'-',label='STIX (6-7 keV)',color='r')
ax.plot(tstix_high,stix_data0[:,1],'-',label='STIX (16-22 keV)',color='b')
ax.set_xticks(tstix_low[::50])
ax.set_xticklabels(['03:03:42','03:14:12','03:21:14','03:28:16','03:36:20','03:44:58','03:55:11','04:03:29','04:14:13'])
ax.set_ylabel('Amplitude');ax.set_xlabel('Time (HH:MM:SS)')
ax.legend();ax1.legend(loc=4);ax.set_yscale('log');ax.set_ylabel('STIX Count Flux (counts/s/cm$^2$/keV)')
ax.set_xlabel('Time (HH:MM:SS UT) (Start Time: 03:03:42 UT)');ax1.set_ylabel('Amplitude')
ax1.set_ylim(0.5,4)
plt.show()


f,ax=plt.subplots(1,1)
ax1=ax.twinx()
c=['k','g','magenta','cyan','brown'];fact=[7,1,0.5,0.5,0.4]
i=4
ax1.plot(tmwa,region1_mean[i],'o-',linewidth=1,markersize=1,label='MWA ('+str(freq_mwa[i])+' MHz) | Height'+str(np.round(height_Mm[i],0))+' Mm',color=c[i])
#ax1.plot(tmwa,(region1_mean[2]-17)*0.2,'o-',linewidth=1,markersize=3,label='MWA ('+str(freq_mwa[2])+' MHz)',color='green')
#ax1.plot(tmwa,region1_mean[0:].mean(axis=0),'o-',linewidth=1,markersize=3,label='MWA ($\\nu$-Average)',color='k')
ax.plot(tstix_low,stix_data0[:,0],'-',label='STIX (6-7 keV)',color='r')
ax.plot(tstix_high,stix_data0[:,1],'-',label='STIX (16-22 keV)',color='b')
ax.set_xticks(tstix_low[::50])
ax.set_xticklabels(['03:03:42','03:14:12','03:21:14','03:28:16','03:36:20','03:44:58','03:55:11','04:03:29','04:14:13'])
ax.set_ylabel('Amplitude');ax.set_xlabel('Time (HH:MM:SS)')
ax.legend();ax1.legend(loc=4);ax.set_yscale('log');ax.set_ylabel('STIX Count Flux (counts/s/cm$^2$/keV)')
ax.set_xlabel('Time (HH:MM:SS UT) (Start Time: 03:03:42 UT)');ax1.set_ylabel('Amplitude')
#ax1.set_ylim(0.25,0.5) # Chan 0
#ax1.set_ylim(2.6,4.0) # Chan 1
#ax1.set_ylim(1.6,2.5) # Chan 1
#ax1.set_ylim(0.8,1.3) # Chan 1
#ax1.set_ylim(3.6,5.2) # Chan 2
#ax1.set_ylim(4.5,6.5) # Chan 3
ax1.set_ylim(5.2,7.8) # Chan 4
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

euv094_list=sorted(glob.glob('/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220915_EUV/*.94A*.fits'))
euv171_list=sorted(glob.glob('/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220915_EUV/*.171A*.fits'))


time_mwa=np.arange(380)*10+3*3600+30*60
i=0;time_aia171=[0]*len(euv171_list)
for e in euv171_list:
    l=euv171_list[i].split('T')[2].split('Z')[0].split('_')
    time_aia171[i]=int(l[0])*3600+int(l[1])*60+float(l[2])
    i=i+1
time_aia171=np.array(time_aia171)


euv094_file=euv094_list[0]
map94=Map(euv094_file)

for j in range(len(time_mwa)):
    mwadata=img[0][j]
    idx=ut.find_predecessor(time_aia171,time_mwa[j])
    euv171_file=euv171_list[idx[0]]
    map171=Map(euv171_file)
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
    mwadata1=mwadata*0;mwadata1[105:130,105:130] = mwadata[105:130,105:130]
    mwadata1[mwadata1<0.5] = 0
    mwamap = Map(mwadata, aiahead1);mwamap1 = Map(mwadata1, aiahead1)
    ii="%04d" % j
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection=map171)
    p0=map171.plot(axes=ax,aspect='auto',vmin=0,vmax=800)
    frac_r1 = [0.05,0.1,0.12,0.2,0.3,0.4,0.5];lev_ct=[0.5,1.0,2.0,3.0]
    for f in lev_ct:
        #lev_r1 = np.nanmax(mwamap.data) * f
        c1 = mwamap1.contour(level=f * u.ct)
        if(len(c1)!=0):
            ax.plot_coord(c1[0], color='green')
    ax.set_ylim(-1000,5000);ax.set_xlim(-1000,5000)
    plt.savefig('/home/rohit/20220915/20220915_pngs/img0_'+str(ii)+'.png')
    plt.close()

