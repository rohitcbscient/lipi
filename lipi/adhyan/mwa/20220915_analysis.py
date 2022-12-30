import pickle
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from sunpy.map import Map
import itertools
import astropy.units as u
from surya.utils import main as ut

baseline=['087-088','087-105','087-123','088-105','088-123','105-123']
ds=[0]*len(baseline);i=0
for b in baseline:
    f1=pickle.load(open('/home/rohit/20220915/20220915_sun.ms_T'+b+'.DS.dat.p','rb'),encoding='latin1')
    ds[i]=f1[5].reshape(4,48,390)
    i=i+1
freq=f1[0]-15.36
meands=np.nanmean(np.array(ds),axis=0)
chan=[0,12,22,38,46]

f,ax=plt.subplots(1,1)
ax.imshow(meands[0],aspect='auto',origin='lower',interpolation=None,vmin=0,vmax=0.01)
#f.colorbar(label='NCCF')
ax.set_xticks([0,50,100,150,200,250,300,350])
ax.set_xticklabels(['03:29:58','03:38:18','03:46:38','03:54:58','04:03:18','04:11:38','04:19:58','04:28:18'])
ax.set_yticks(np.arange(48))
ax.set_yticklabels(np.round(freq,2));ax.set_xlabel('Time (2022-09-15 HH:MM:SS UT)');ax.set_ylabel('Frequency (MHz)')
for c in chan:
    ax.axhline(y=c)
plt.show()

img=[0]*len(chan);i=0
for c in chan:
    print(c)
    imglist=sorted(glob.glob('/home/rohit/20220915/20220915_sun.chan.'+str(c)+'.pol.I.time.*.image_image.FITS'))[0:380]
    j=0;img[i]=[0]*len(imglist)
    for im in imglist:
        aa=fits.open(im);data=aa[0].data[0][0][924:1124,924:1124]
        img[i][j]=data
        j=j+1
    i=i+1
img=np.array(img)
region1=img[:,:,105:130,105:130];region1_max=np.nanmax(region1,axis=(2,3))
region1_mean=np.nanmean(region1,axis=(2,3))

plt.imshow(region1_max,aspect='auto',origin='lower',interpolation=None)

plt.imshow(img[0].mean(axis=0),aspect='auto',origin='lower',interpolation=None)


f,ax=plt.subplots(1,1)
ax.plot(region1_mean[3],'o-')
ax.set_xticks([0,50,100,150,200,250,300,350])
ax.set_xticklabels(['03:30:00','03:38:20','03:46:40','03:55:00','04:03:20','04:11:40','04:20:00','04:28:20'])
ax.set_ylabel('Amplitude');ax.set_xlabel('Time (HH:MM:SS)')
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

