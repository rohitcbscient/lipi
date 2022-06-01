import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.io import readsav
import glob
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable

filename='/media/rohit/MWA/20140914/rstn_20140914.SRD'
print('Reading '+str(filename)+'....')
ff=open(str(filename),"r")
l = ff.readlines()
n=30700
time=[0]*n
time_s=[0]*n
flux_245=[0]*n;flux_410=[0]*n;flux_610=[0]*n;flux_1415=[0]*n
flux_2695=[0]*n;flux_4995=[0]*n;flux_8800=[0]*n;flux_15400=[0]*n
for i in range(n):
        time[i]=l[i].split(' ')[0]
        time_s[i]=float(l[i].split(' ')[0][0:2])*3600+float(l[i].split(' ')[0][2:4])*60+float(l[i].split(' ')[0][4:6])
        l_=l[i].split(' ')
        if(l_[1][0]!='/'):
            flux_245[i]=float(l_[1][0]+'.'+l_[1][1]+l_[1][2])*10**int(l_[1][3])
        else:
            flux_245[i]=np.nan
        if(l_[2][0]!='/'):
            flux_410[i]=float(l_[2][0]+'.'+l_[2][1]+l_[2][2])*10**int(l_[2][3])
        else:
            flux_410[i]=np.nan
        if(l_[3][0]!='/'):
            flux_610[i]=float(l_[3][0]+'.'+l_[3][1]+l_[3][2])*10**int(l_[3][3])
        else:
            flux_610[i]=np.nan
        if(l_[4][0]!='/'):
            flux_1415[i]=float(l_[4][0]+'.'+l_[4][1]+l_[4][2])*10**int(l_[4][3])
        else:
            flux_1415[i]=np.nan
        if(l_[5][0]!='/'):
            flux_2695[i]=float(l_[5][0]+'.'+l_[5][1]+l_[5][2])*10**int(l_[5][3])
        else:
            flux_2695[i]=np.nan
        if(l_[6][0]!='/'):
            flux_4995[i]=float(l_[6][0]+'.'+l_[6][1]+l_[6][2])*10**int(l_[6][3])
        else:
            flux_4995[i]=np.nan
        if(l_[7][0]!='/'):
            flux_8800[i]=float(l_[7][0]+'.'+l_[7][1]+l_[7][2])*10**int(l_[7][3])
        else:
            flux_8800[i]=np.nan
        if(l_[8][0]!='/'):
            flux_15400[i]=float(l_[8][0]+'.'+l_[8][1]+l_[8][2])*10**int(l_[8][3])
        else:
            flux_15400[i]=np.nan

flux_245=np.array(flux_245);flux_410=np.array(flux_410);flux_610=np.array(flux_610)
flux_1415=np.array(flux_1415);flux_2695=np.array(flux_2695);flux_4995=np.array(flux_4995);flux_8800=np.array(flux_8800);flux_15400=np.array(flux_15400)
flux_all=np.array([flux_245,flux_410,flux_610,flux_1415,flux_2695,flux_4995,flux_8800,flux_15400])
flux_all_diff=flux_all-np.array(list(flux_all[:,0])*n).reshape(n,8).swapaxes(0,1)
freq_rstn=[245,410,610,1415,2695,4995,8800,15400]


f1='/data/Dropbox/20140914/wind/20140914.R1'
f2='/data/Dropbox/20140914/wind/20140914.R2'
f3='/data/Dropbox/20140914/wind/20140914.tnr'
goesfile='/data/Dropbox/20140914/20140914_idlsave_goes.sav'
freq=np.array([108, 120, 133, 145, 160, 179, 197, 217, 240])
freq_list=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
######### WIND WAVES FREQUENCY ################

# RAD1 - Radio Receiver Band- Frequency range    20 kHz - 1,040 kHz / No. channels   256 / Bandwidth 3 kHz (https://solar-radio.gsfc.nasa.gov/wind/instrument.html)
# RAD2 - Radio Receiver Band- Frequency range   1.075 MHz - 13.825 MHz / No. channels   256 / Bandwidth 20 kHz 
# Thermal Noise Receiver (TNR) / Frequency range    4 kHz - 256 kHz / No. channels  32 or 16 per band (5 bands) / Bandwidth 400 Hz - 6.4 kHz
freq_rad1=np.round(np.linspace(20,1040,256)/1.e3,2)
freq_rad2=np.round(np.linspace(1.075,13.825,256),2)
freq_tnr=np.round(np.concatenate((np.linspace(4,25.6,36),np.linspace(25.6,256,54))),2)

aa1=readsav(f1)['arrayb'];aa2=readsav(f2)['arrayb'];aa3=readsav(f3)['arrayb'] # 1 min time cadence
#flux_list=sorted(glob.glob('/media/rohit/MWA/20140914/pickle/flux_V1_1094693472_*.ms_T000-008.p'))
flux=[0]*len(freq_list)
for i in range(len(freq_list)):
    ff=sorted(glob.glob('/media/rohit/MWA/20140914/pickle/flux*'+str(freq_list[i])+'*T008-009.p'));flux[i]=[0]*len(ff)
    for j in range(len(ff)):
        pp=pickle.load(open(ff[j],'rb'),encoding='latin1')
        flux[i][j]=pp[17][3][0][0]
        #flux[i][j]=pp[5][0][0]
flux=np.array(flux)
flux=np.array(flux).reshape(9,53*25)
########## GOES ############################

goes=readsav(goesfile);goes18=goes['yclean'][0];goes0540=goes['yclean'][1];tgoes=goes['tarray']

plt.plot(flux[1].flatten())
plt.show()

f,ax=plt.subplots(1,1)
ax.imshow(flux,aspect='auto',origin='lower',vmin=0,vmax=150)
plt.show()

f,ax=plt.subplots(5,1);ax0=ax[0];ax1=ax[1];ax2=ax[2];ax3=ax[3];ax4=ax[4]
im4=ax4.imshow(aa3[:,60:351],aspect='auto',origin='lower',vmin=1,vmax=5.0)
im3=ax3.imshow(aa1[:,60:351],aspect='auto',origin='lower',vmin=1,vmax=2.0)
im2=ax2.imshow(aa2[:,60:351],aspect='auto',origin='lower',vmin=1,vmax=1.3)
im1=ax1.imshow(flux,aspect='auto',origin='lower',vmin=0,vmax=150)
divider = make_axes_locatable(ax1);cax = divider.append_axes("right", size="2%", pad=0.0);f.colorbar(im1,cax=cax,label='(SFU)')
for i in range(8):
    im0=ax0.plot(time_s[5400:21060],flux_all[i][5400:21060],label=str(np.round(freq_rstn[i]/1.e3,2))+' GHz')
ax01=ax0.twinx()
ax01.plot(tgoes[900:6930]+3600,goes18[900:6930],'--',color='k',label='GOES ($1-8 \AA$)')
ax01.plot(tgoes[900:6930]+3600,goes0540[900:6930],'-',color='k',label='GOES ($0.5-4.0 \AA$)');ax01.set_yscale('log')
ax0.legend(loc=1,fontsize=6);ax01.legend(loc=4,fontsize=6);ax01.set_ylabel('X-ray Flux ($W/m^{2}$)')
f.subplots_adjust(hspace=0.03)
ax0.set_ylabel('Flux (SFU)')
ax3.set_yticks(np.arange(256)[::50]);ax3.set_yticklabels(freq_rad1[::50])#;ax2.set_yscale('log');ax2.set_ylim(1,256)
ax2.set_yticks(np.arange(256)[::50]);ax2.set_yticklabels(freq_rad2[::50])#;ax3.set_yscale('log');ax3.set_ylim(1,256)
ax4.set_yticks(np.arange(90)[::15]);ax4.set_yticklabels(freq_tnr[::15])#;ax4.set_yscale('log');ax4.set_ylim(1,90)
ax1.set_yticks(np.arange(9));ax1.set_yticklabels(freq)
ax1.set_ylabel('Frequency (MHz)')
ax2.set_ylabel('Frequency (MHz)')
ax3.set_ylabel('Frequency (MHz)')
ax4.set_ylabel('Frequency (kHz)')
ax4.set_xticks(np.arange(290)[::50]);ax4.set_xticklabels(['01:30','02:10','03:00','03:50','04:40','05:30'])
ax4.set_xlabel('2014/09/14 Time (HH:MM) UT')
#f.colorbar(im0)
plt.show()

