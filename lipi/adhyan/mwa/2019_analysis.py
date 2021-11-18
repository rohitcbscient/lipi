import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob
from astropy.io import fits
from sunpy.map import Map

freq=[107.5,163.8,192.0,240.6]
list_=sorted(glob.glob('/media/rohit/MWA/20190405_PSP_MWA/sun_1238472032/all_images/*.p'))
cor=[1.06,1.07,1.2,1.77]

Tb=[0]*len(freq)
S=[0]*len(freq)
for i in range(len(freq)):
    d=pickle.load(open(list_[i],'rb'))
    Tb[i]=np.array(d[0])/cor[i];S[i]=np.array(d[8])/cor[i]

Tb=np.array(Tb)
Tbmean=np.nanmean(Tb,axis=1)

aiafile='/media/rohit/MWA/2019_EUV/JSOC_20200907_1363/aia.lev1_euv_12s.2019-04-05T040151Z.171.image.fits'
aia=fits.open(aiafile,checksum=True)
map171=Map(aiafile)
aiadata_=aia[1].data
aiahead=aia[1].header;aiaxc=aiahead['CRPIX1'];aiayc=aiahead['CRPIX2']
aiadata=aiadata_[int(aiaxc)-1600:int(aiaxc)+1600,int(aiayc)-1600:int(aiayc)+1600]

for i in range(len(freq)):
    plt.plot(S[i],'o-',label=str(freq[i])+' MHz')
plt.xlabel('Time (HH:MM:SS UT)')
plt.ylabel('Flux Density (SFU)')
plt.legend()
plt.ylim(0,40)
plt.xticks([0,5,10,15,20,25],['04:00:16','04:01:06','04:01:56','04:02:46','04:03:36','04:04:26'])
plt.show()

plot_tb=1
if(plot_tb):
    Tb[np.where(Tb==0)]=np.nan
    plt.plot(freq,np.nanmean(Tb,axis=(1,2,3))/1.e6,'o-')
    plt.errorbar(freq,np.nanmean(Tb,axis=(1,2,3))/1.e6,yerr=0.2*np.nanmean(Tb,axis=(1,2,3))/1.e6,uplims=[0,1,1,1])
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Mean T$_B$ (MK)')
    plt.xlim(100,250)
    plt.show()

plt.plot(freq,np.array(S).mean(axis=(1)),'o-')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Mean Solar Flux (SFU)')
plt.xlim(100,250)
plt.show()

for i in range(len(freq)):
    plt.plot(np.max(Tb[i],axis=(1,2))/1.e6,'o-',label=str(freq[i])+' MHz')
plt.xlabel('Time (HH:MM:SS UT)')
ax3.set_ylim(mwaxl,mwaxr);ax3.set_xlim(mwaxl,mwaxr)
ax3.set_ylim(mwaxl,mwaxr);ax3.set_xlim(mwaxl,mwaxr)
plt.ylabel('T$_{B}$ (MK)')
plt.legend()
plt.ylim(0.1,500)
plt.xticks([0,5,10,15,20,25],['04:00:16','04:01:06','04:01:56','04:02:46','04:03:36','04:04:26'])
plt.show()

f,ax=plt.subplots(2,2,figsize=(15,8));ax0=ax[0,0];ax1=ax[0,1];ax2=ax[1,0];ax3=ax[1,1]
mwaxl=-2500;mwaxr=2500
i0=ax0.imshow(Tbmean[0],extent=[mwaxl,mwaxr,mwaxl,mwaxr],origin=0)
i1=ax1.imshow(Tbmean[1],extent=[mwaxl,mwaxr,mwaxl,mwaxr],origin=0)
i2=ax2.imshow(Tbmean[2],extent=[mwaxl,mwaxr,mwaxl,mwaxr],origin=0)
i3=ax3.imshow(Tbmean[3],extent=[mwaxl,mwaxr,mwaxl,mwaxr],origin=0)
ax0.set_ylim(mwaxl,mwaxr);ax0.set_xlim(mwaxl,mwaxr)
ax1.set_ylim(mwaxl,mwaxr);ax1.set_xlim(mwaxl,mwaxr)
ax2.set_ylim(mwaxl,mwaxr);ax2.set_xlim(mwaxl,mwaxr)
ax3.set_ylim(mwaxl,mwaxr);ax3.set_xlim(mwaxl,mwaxr)
circ1=plt.Circle((0.5,0.5), radius=60*16, color='white',linestyle='--', linewidth=4,fill=False);ax0.add_patch(circ1)
circ2=plt.Circle((0.5,0.5), radius=60*16, color='white',linestyle='--', linewidth=4,fill=False);ax1.add_patch(circ2)
circ3=plt.Circle((0.5,0.5), radius=60*16, color='white',linestyle='--', linewidth=4,fill=False);ax2.add_patch(circ3)
circ4=plt.Circle((0.5,0.5), radius=60*16, color='white',linestyle='--', linewidth=4,fill=False);ax3.add_patch(circ4)
ax0.set_title('108 MHz');ax1.set_title('164 MHz');ax2.set_title('192 MHz');ax3.set_title('240 MHz')
ax2.set_xlabel('Solar-X (arcsec)');ax3.set_xlabel('Solar-X (arcsec)')
ax2.set_ylabel('Solar-Y (arcsec)');ax0.set_ylabel('Solar-Y (arcsec)')
plt.show()

plt.rcParams.update({'font.size': 22})
lev=[0.5,0.7,0.8,0.9]
lev_108=[0.7,0.8,0.9,0.99]
#plt.imshow(aiadata,extent=[-1229,1229,-1229,1229],cmap='gray')
map171.plot(vmin=100,vmax=2000)
mwaxl=-2500/0.6+int(aiaxc);mwaxr=2500/0.6+int(aiayc);mwayl=-2500/0.6+int(aiaxc)
plt.contour(Tbmean[0]/np.nanmax(Tbmean[0]),extent=[mwaxl,mwaxr,mwaxl,mwaxr],levels=lev_108,linewidths=2,colors='red')
plt.contour(Tbmean[1]/np.nanmax(Tbmean[1]),extent=[mwaxl,mwaxr,mwaxl,mwaxr],levels=lev,linewidths=2,colors='lawngreen')
plt.contour(Tbmean[2]/np.nanmax(Tbmean[2]),extent=[mwaxl,mwaxr,mwaxl,mwaxr],levels=lev,linewidths=2,colors='magenta')
plt.contour(Tbmean[3]/np.nanmax(Tbmean[3]),extent=[mwaxl,mwaxr,mwaxl,mwaxr],levels=lev,linewidths=2,colors='blue')
plt.xlim(-1500/0.6+int(aiaxc),1500/0.6+int(aiayc))
plt.ylim(-1500/0.6+int(aiaxc),1500/0.6+int(aiayc))
#plt.text(0.50,0.5,'240 MHz',color='blue')
#plt.text(0.50,-150,'192 MHz',color='magenta')
#plt.text(0.50,-300,'164 MHz',color='lawngreen')
#plt.text(0.50,-450,'108 MHz',color='red')
plt.title('AIA 171 $\AA$ and MWA Contours')
plt.text(0.5+3200,0.5+3200,'240 MHz',color='blue',fontweight='bold')
plt.text(0.5+3200,-150+3200,'192 MHz',color='magenta',fontweight='bold')
plt.text(0.5+3200,-300+3200,'164 MHz',color='lawngreen',fontweight='bold')
plt.text(0.5+3200,-450+3200,'108 MHz',color='red',fontweight='bold')
plt.show()



