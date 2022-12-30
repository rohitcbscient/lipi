import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.io import readsav
import glob
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sunpy.map import Map
from scipy.io import  readsav
import pandas as pd
import astropy.units as u
from surya.utils import main as ut
import matplotlib.cm as cm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap
from astropy.coordinates import SkyCoord

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
flux_all[:,12000:13000]=np.nan

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

bdiff_list=sorted(glob.glob('/sdata/fits/running_diff/*94*bdiff.fits'))[::2];max94=[0]*len(bdiff_list);i=0
for bb in bdiff_list:
    aa=Map(bb);max94[i]=np.nanmax(aa.data)
    i=i+1
#------------- RHESSI
corr_count=readsav("/media/rohit/MWA/20140914/rhessi/counts_obs.sav")
#corr_count.setflags(write=1)
obs_time=corr_count['obs_times']
obs_energies = corr_count['obs_energies']
corr_data=corr_count['obs_data']*1.0;corr_data[753:1137,:]=np.nan;corr_data[2300:2620,:]=np.nan;corr_data=corr_data.swapaxes(0,1)

#2014-09-14T01:30:28.000 TO 2014-09-14T04:53:20.000'

f,ax=plt.subplots(1,1)
ax.plot(obs_time,corr_data[2],'o-',markersize=2,label='12-25 keV')
ax.plot(obs_time,corr_data[3],'o-',markersize=2,label='25-50 keV')
ax.axvspan(obs_time[400],obs_time[700],color='yellow',alpha=0.2)
ax.axvspan(obs_time[753],obs_time[1137],color='gray',alpha=0.2)
ax.axvspan(obs_time[2300],obs_time[2620],color='gray',alpha=0.2,label='RHESSI Nights')
ax.set_yscale('log')
ax.set_xticks([obs_time[0],obs_time[500],obs_time[1000],obs_time[1500],obs_time[2000],obs_time[2500],obs_time[3000]])
ax.set_xticklabels(['01:30:28','02:03:48','02:37:08','03:10:28','03:43:48','04:14:08','04:47:28'])
ax.set_xlabel('Time (HH:MM:SS)');ax.set_ylabel('Corrected Counts (cm$^{-2}$ s$^{-1}$ keV$^{-1}$)')
ax.legend(loc=2);ax.set_ylim([1.e0,8.e3])
left, bottom, width, height = [0.45, 0.45, 0.4, 0.4]
ax2 = f.add_axes([left, bottom, width, height])
ax2.plot(obs_time[400:1000],corr_data[2][400:1000],'o-',markersize=2,label='')
ax2.plot(obs_time[400:1000],corr_data[3][400:1000],'o-',markersize=2,label='')
#ax2.set_xlabel('Strain (in/in)')
#ax2.set_ylabel('Stress (ksi)')
ax2.set_title('Rising Phase');ax2.set_yscale('log')
ax2.set_xticks([obs_time[400],obs_time[500],obs_time[600],obs_time[700]])
ax2.set_xticklabels(['02:03:48','02:10:28','02:17:08','02:23:48'])
#ax2.set_xlim([0,0.008])
#ax2.set_ylim([0,100])
plt.show()


jj=38
rhessi_list=sorted(glob.glob('/media/rohit/MWA/20140914/rhessi/fitspec_*.fits'))
rhessi_spec=fits.open(rhessi_list[jj])
title_string= rhessi_spec[1].header['DATE-OBS'].split('T')[0]+' '+rhessi_spec[1].header['DATE-OBS'].split('T')[1]
err_rhessi_spec_rate=rhessi_spec[1].data['STAT_ERR']; rhessi_spec_rate=rhessi_spec[1].data['RATE']
rhessi_chan = rhessi_spec[2].data['CHANNEL']
phmodel=rhessi_spec[1].data['PHMODEL'][0];confac= rhessi_spec[1].data['CONVFAC'][0]
rhessi_res=rhessi_spec[1].data['RESIDUAL'][0]
f,(ax0,ax1) = plt.subplots(2,1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
ax0.plot(rhessi_chan,phmodel,'o-',label='PHMODEL')
ax0.plot(rhessi_chan,confac,'o-',label='CONVFAC')
ax0.plot(rhessi_chan,rhessi_spec_rate[0],color='blue')
ax0.errorbar(rhessi_chan,rhessi_spec_rate[0],yerr=err_rhessi_spec_rate[0],color='blue',label='RATE')
ax0.legend();ax0.set_title(title_string);ax0.set_xscale('log');ax0.set_yscale('log')
ax1.plot(rhessi_chan,rhessi_res,'o-',label='RESIDUAL: $\chi^2$ = '+str(np.round(rhessi_spec[1].data['CHISQ'][0],2)));ax1.legend()
ax1.set_ylim(-10,10);ax1.set_xticks(rhessi_chan[::15]);ax1.set_xticks(rhessi_chan[::15])
ax1.set_xscale('log');ax0.set_ylabel('Flux (cm$^{-2}$ s$^{-1}$)');ax1.set_xlabel('Energy Bins (keV)')
ax1.legend()
plt.show()


tmwa,tsubmwa,Tbsubmax,xcsub90,ycsub90,maxsubX,maxsubY,pa_angle = pickle.load(open('/media/rohit/MWA/20140914/Tb_centroid.p','rb'))
Tb240=pickle.load(open('/media/rohit/MWA/20140914/pickle/Tb_20140914_187-188_sub_test.p','rb'),encoding='bytes')
Tbmap240=Tb240[0]
Tb161=pickle.load(open('/media/rohit/MWA/20140914/pickle/Tb_20140914_125-126_sub_test.p','rb'),encoding='bytes')
Tbmap161=Tb161[0]
Tb108=pickle.load(open('/media/rohit/MWA/20140914/pickle/Tb_20140914_084-085_sub_test.p','rb'),encoding='bytes')
Tbmap108=Tb108[0]


ddiff_list=sorted(glob.glob('/sdata/fits/running_diff/*171*ddiff.fits'))
bdiff_list=sorted(glob.glob('/sdata/fits/running_diff/*171*bdiff.fits'))
rhessi_image_list=sorted(glob.glob('/media/rohit/MWA/20140914/rhessi/fitsimage_15-25keV_*.fits'))
aiafiles171=sorted(glob.glob('/media/rohit/MWA/20140914/EUV/fits/*171A*.fits'))
aiafiles94=sorted(glob.glob('/media/rohit/MWA/20140914/EUV/fits/*94A*.fits'))



taiasec=np.zeros(len(ddiff_list));tmwaidx=[0]*len(ddiff_list);i=0
for t in range(len(ddiff_list)):
    aa=fits.open(ddiff_list[i]);aiadate=aa[0].header['DATE-OBS']
    ta = aiadate.split('T')[1].split(':')
    taiasec[i]=int(ta[0])*3600+int(ta[1])*60+float(ta[2].split('Z')[0])
    i=i+1
trhessi_sec=np.zeros(len(rhessi_image_list));i=0
for t in range(len(rhessi_image_list)):
    aa=fits.open(rhessi_image_list[i]);aiadate=aa[0].header['DATE_OBS']
    ta = aiadate.split('T')[1].split(':')
    trhessi_sec[i]=int(ta[0])*3600+int(ta[1])*60+float(ta[2].split('Z')[0])
    i=i+1

aiamap_171=[0]*len(aiafiles171);aia171sec=[0]*len(aiafiles171)
for i in range(len(aiafiles171)):
    print(i)
    aa=fits.open(aiafiles171[i])
    date171 = aa[0].header['DATE-OBS'].split('T')[1]
    aia171sec[i]=int(date171.split(':')[0])*3600+int(date171.split(':')[1])*60+float(date171.split(':')[2])
    aiamap_171[i]=np.nanmax(aa[0].data[1000:1800,3000:3800])
aiamap_171=np.array(aiamap_171);aia171sec=np.array(aia171sec)

aiamap_94=[0]*len(aiafiles94);aia94sec=[0]*len(aiafiles94)
for i in range(len(aiafiles94)):
    aa=fits.open(aiafiles94[i])
    date94 = aa[0].header['DATE-OBS'].split('T')[1]
    aia94sec[i]=int(date94.split(':')[0])*3600+int(date94.split(':')[1])*60+float(date94.split(':')[2])
    aiamap_94[i]=np.nanmax(aa[0].data[1000:1800,3000:3800])
aiamap_94=np.array(aiamap_94);aia94sec=np.array(aia94sec)

pickle.dump([aia94sec,aiamap_94,aia171sec,aiamap_171],open('/media/rohit/MWA/20140914/EUV/20140914_AIA_timeseries,p','wb'))

aia94sec,aiamap_94,aia171sec,aiamap_171=pickle.load(open('/media/rohit/MWA/20140914/EUV/20140914_AIA_timeseries,p','rb'))
rh_counts=readsav('/media/rohit/MWA/20140914/rhessi/counts_obs.sav')



i=20
for i in range(len(ddiff_list)):
    #i=400
    taiasec=np.array(taiasec)
    tmwaidx=(np.abs(tsubmwa - taiasec[i])).argmin()
    trhessi_idx=(np.abs(trhessi_sec- taiasec[i])).argmin()
    #-----------------------
    #-----------------------
    aiamap=Map(ddiff_list[i]);aiahead1=aiamap.meta
    aiahead1['naxis1']=imghead['NAXIS1'];aiahead1['naxis2']=imghead['NAXIS2']
    aiahead1['CRPIX1']=imghead['CRPIX1'];aiahead1['CRPIX2']=imghead['CRPIX2']
    aiahead1['CRVAL1']=imghead['CRVAL1'];aiahead1['CRVAL2']=imghead['CRVAL2']
    aiahead1['CDELT1']=imghead['CDELT1'];aiahead1['CDELT2']=imghead['CDELT2']
    img = fits.open(rhessi_image_list[trhessi_idx]);imgdata = img[0].data;imghead = img[0].header
    rhessi_map=Map(imgdata,aiahead1)
    #-----------
    aiamap=Map(ddiff_list[i]);aiahead2=aiamap.meta
    aiahead2['naxis1']=Tbmap240.shape[1];aiahead2['naxis2']=Tbmap240.shape[2]
    aiahead2['CRPIX1']=Tbmap240.shape[1]/2.;aiahead2['CRPIX2']=Tbmap240.shape[2]/2.
    aiahead2['CRVAL1']=0;aiahead2['CRVAL2']=0
    aiahead2['CDELT1']=50;aiahead2['CDELT2']=50
    mwamap240=Map(Tbmap240[tmwaidx],aiahead2)
    mwamap161=Map(Tbmap161[tmwaidx],aiahead2)
    mwamap108=Map(Tbmap108[tmwaidx],aiahead2)
    aiamap=Map(ddiff_list[i])
    #------------
    f=plt.figure()
    ii="%04d"%i
    ax=f.add_subplot(projection=aiamap)
    aiamap.plot(axes=ax,vmin=-800,vmax=4000,cmap='binary');frac_r1=0.9;frac_r2=0.8;frac_r3=0.85
    lev_r1=np.nanmax(rhessi_map.data)*frac_r1
    lev_r2 = np.nanmax(rhessi_map.data) * frac_r2
    lev_r3 = np.nanmax(rhessi_map.data) * frac_r3
    if(np.nanmax(rhessi_map.data)!=0):
        c1=rhessi_map.contour(level=lev_r1 * u.ct);ax.plot_coord(c1[0],color='yellow')
        c2=rhessi_map.contour(level=lev_r2* u.ct);ax.plot_coord(c2[0],color='yellow')
        c3=rhessi_map.contour(level=lev_r3 * u.ct);ax.plot_coord(c3[0],color='yellow')
    frac=0.9
    if(np.nanmax(mwamap240.data)!=0):
        lev240=np.nanmax(mwamap240.data)*frac;lev161=np.nanmax(mwamap161.data)*frac;lev108=np.nanmax(mwamap108.data)*frac
        c1=mwamap240.contour(level=lev240*u.ct);ax.plot_coord(c1[0],color='blue',linewidth=3)
        c2=mwamap161.contour(level=lev161*u.ct);ax.plot_coord(c2[0],color='green',linewidth=3)
        c3=mwamap108.contour(level=lev108*u.ct);ax.plot_coord(c3[0],color='red',linewidth=3)
    #ax.set_xlim([2048,4096]); ax.set_ylim([0,2048])
    ax.set_xlim([3000,4000]); ax.set_ylim([800,1800])
    ax.text(3750, 1080, '15-25 keV', color='yellow', bbox=dict(facecolor='k', alpha=0.8))
    ax.text(3750, 1000, '240 MHz', color='blue', bbox=dict(facecolor='white', alpha=0.8))
    ax.text(3750, 920, '161 MHz', color='green', bbox=dict(facecolor='white', alpha=0.8))
    ax.text(3750, 840, '108 MHz', color='red', bbox=dict(facecolor='white', alpha=0.8))
    f.savefig('/sdata/fits/pngs_rhessi/all_'+str(ii)+'.png',dpi=150)
    plt.close()


for i in range(len(ddiff_list)):
    taiasec=np.array(taiasec)
    tmwaidx=(np.abs(tsubmwa - taiasec[i])).argmin()
    trhessi_idx=(np.abs(trhessi_sec- taiasec[i])).argmin()
    img = fits.open(rhessi_image_list[trhessi_idx]);imgdata = img[0].data;imghead = img[0].header
    #-----------------------
    #-----------------------
    aiamap=Map(ddiff_list[i]);aiahead1=aiamap.meta
    aiahead1['naxis1']=imghead['NAXIS1'];aiahead1['naxis2']=imghead['NAXIS2']
    aiahead1['CRPIX1']=imghead['CRPIX1'];aiahead1['CRPIX2']=imghead['CRPIX2']
    aiahead1['CRVAL1']=imghead['CRVAL1'];aiahead1['CRVAL2']=imghead['CRVAL2']
    aiahead1['CDELT1']=imghead['CDELT1'];aiahead1['CDELT2']=imghead['CDELT2']
    rhessi_map=Map(imgdata,aiahead1)
    #-----------
    aiamap=Map(ddiff_list[i]);aiahead2=aiamap.meta
    aiahead2['naxis1']=Tbmap240.shape[1];aiahead2['naxis2']=Tbmap240.shape[2]
    aiahead2['CRPIX1']=Tbmap240.shape[1]/2.;aiahead2['CRPIX2']=Tbmap240.shape[2]/2.
    aiahead2['CRVAL1']=0;aiahead2['CRVAL2']=0
    aiahead2['CDELT1']=50;aiahead2['CDELT2']=50
    mwamap240=Map(Tbmap240[tmwaidx],aiahead2)
    mwamap161=Map(Tbmap161[tmwaidx],aiahead2)
    mwamap108=Map(Tbmap108[tmwaidx],aiahead2)
    aiamap=Map(ddiff_list[i])
    #------------
    f=plt.figure()
    ii="%04d"%i
    ax=f.add_subplot(projection=aiamap)
    aiamap.plot(axes=ax,vmin=-100,vmax=200,cmap='coolwarm');frac_r1=0.9;frac_r2=0.8;frac_r3=0.85
    lev_r1=np.nanmax(rhessi_map.data)*frac_r1
    lev_r2 = np.nanmax(rhessi_map.data) * frac_r2
    lev_r3 = np.nanmax(rhessi_map.data) * frac_r3
    if(np.nanmax(rhessi_map.data)!=0):
        c1=rhessi_map.contour(level=lev_r1 * u.ct);ax.plot_coord(c1[0],color='yellow')
        c2=rhessi_map.contour(level=lev_r2* u.ct);ax.plot_coord(c2[0],color='yellow')
        c3=rhessi_map.contour(level=lev_r3 * u.ct);ax.plot_coord(c3[0],color='yellow')
    #ax.set_xlim([2048,4096]); ax.set_ylim([0,2048])
    ax.set_xlim([3000,4000]); ax.set_ylim([800,1800])
    ax.text(3750, 1080, '15-25 keV', color='yellow', bbox=dict(facecolor='k', alpha=0.8))
    #ax.text(3750, 1000, '240 MHz', color='blue', bbox=dict(facecolor='white', alpha=0.8))
    #ax.text(3750, 920, '161 MHz', color='green', bbox=dict(facecolor='white', alpha=0.8))
    #ax.text(3750, 840, '108 MHz', color='red', bbox=dict(facecolor='white', alpha=0.8))
    cm = LinearSegmentedColormap.from_list('jet', ["#0000FF", "#00FF00" , "#FF0000"])
    colors = cm(np.linspace(0, 1, 9));pix=[0]*9
    #colors = iter(cm.jet(np.linspace(0, 1, 9)))
    for j in range(9):
        j=8-j
        seeds0 = SkyCoord(np.array(xcsub90[j])[tmwaidx]*u.arcsec, np.array(ycsub90[j])[tmwaidx]*u.arcsec,frame=aiamap.coordinate_frame)
        im=ax.plot_coord(seeds0, color=colors[j], marker='s', markersize=8,alpha=0.8,linestyle='None',markeredgecolor='white')
        pix_=aiamap.wcs.world_to_pixel(seeds0);pix[j]=[float(pix_[0]),float(pix_[1])]
    pix=np.array(pix);qx=pix[7:9,0].mean();qy=pix[7:9,1].mean();qx1=pix[0:2,0].mean();qy1=pix[0:2,1].mean()
    cbar=f.colorbar(ScalarMappable(cmap=cm, norm=plt.Normalize(0, 8)), ticks=np.arange(9), label='Frequency (MHz)')
    cbar.ax.set_yticklabels(freq)
    ax.quiver(qx,qy,qx1-qx,(-qy+qy1), color='k', scale=1, scale_units='xy',zorder=20)
    f.savefig('/sdata/fits/pngs_rhessi/all_centroids1_'+str(ii)+'.png',dpi=150)
    plt.close()




rh_fit_list=glob.glob('/media/rohit/MWA/20140914/rhessi/fitsummary_*.txt')
spidx1=[0]*len(rh_fit_list);spidx2=[0]*len(rh_fit_list);spidx_break=[0]*len(rh_fit_list)
espidx1=[0]*len(rh_fit_list);espidx2=[0]*len(rh_fit_list);espidx_break=[0]*len(rh_fit_list);chisq=[0]*len(rh_fit_list)
rh_time=[0]*len(rh_fit_list);rh_time_sec=[0]*len(rh_fit_list)
for i in range(len(rh_fit_list)):
    df = pd.read_csv(rh_fit_list[i],skiprows=7)
    df1 = df.to_records()
    rh_time[i]=df1.dtype.descr[2][0].split(' ')[2]
    rh_time_sec[i]=int(rh_time[i].split(':')[0])*3600+int(rh_time[i].split(':')[1])*60+int(float(rh_time[i].split(':')[2]))
    rec1=df1[3][1].split('     ');rec1=[x for x in rec1 if x]
    rec2=df1[4][1].split('     ');rec2=[x for x in rec2 if x]
    rec3=df1[5][1].split('     ');rec3=[x for x in rec3 if x]
    rec4=df1[6][1].split('     ');rec4=[x for x in rec4 if x]
    spidx1[i]=float(rec1[1]);spidx2[i]=float(rec3[1]);spidx_break[i]=float(rec2[1])
    espidx1[i]=float(rec1[2]);espidx2[i]=float(rec3[2]);espidx_break[i]=float(rec2[2])
    chisq[i]= float(df1[0][1].split('Chisq=')[1][0:4])
rh_idx=sorted(range(len(rh_time_sec)), key=lambda k: rh_time_sec[k])
rh_time1=[0]*len(rh_time)
for i in range(len(rh_time)):
    rh_time1[i]=rh_time[rh_idx[i]]
rh_time_sec=np.array(rh_time_sec)[rh_idx]
spidx1=np.array(spidx1)[rh_idx];spidx2=np.array(spidx2)[rh_idx];spidx_break=np.array(spidx_break)[rh_idx]
espidx1=np.array(espidx1)[rh_idx];espidx2=np.array(espidx2)[rh_idx];espidx_break=np.array(espidx_break)[rh_idx];chisq=np.array(chisq)




plt.plot(rh_time_sec,spidx1,'o-',color='blue',label='Low-Energy')
plt.errorbar(rh_time_sec,spidx1,linestyle='None',yerr=espidx1,color='blue')
plt.plot(rh_time_sec,spidx2,'o-',color='red',label='High-Energy')
plt.errorbar(rh_time_sec,spidx2,linestyle='None',yerr=espidx2,color='red')
plt.legend();plt.xlabel('Time (HH:MM:SS)');plt.ylabel('Spectral Index')
plt.xticks(rh_time_sec[::18],rh_time1[::18])
plt.show()

plt.plot(rh_time_sec,spidx_break,'o-',color='k',label='High-Energy')
plt.errorbar(rh_time_sec,spidx_break,linestyle='None',yerr=espidx_break,color='k')
plt.legend();plt.xlabel('Time (HH:MM:SS)');plt.ylabel('Break Energy (keV)')
plt.xticks(rh_time_sec[::18],rh_time1[::18])
plt.show()


rh_time=obs_time-obs_time[0]+5400
f,ax=plt.subplots(1,1)
rstny=(flux_all[-1][5400:21060]-flux_all[-1][5400])*150+2000
ax.plot(time_s[5400:21060],rstny,linewidth=0.4,label='RSTN '+str(np.round(freq_rstn[-1]/1.e3,2))+' GHz')
ax.plot(tgoes[900:6930]+3600,goes18[900:6930]*1.e9,'--',color='k',label='GOES ($1-8 \\AA$)')
ax.plot(tgoes[900:6930]+3600,goes0540[900:6930]*1.e9,'-',color='k',label='GOES ($0.5-4.0 \\AA$)')
ax.plot(aia94sec,(aiamap_94-aiamap_94[0])*1,'-',label='AIA 94$\\AA$')
#ax1=ax.twinx();ax2=ax.twinx()
rhy=np.nanmean(corr_count['obs_data'],axis=1)*1.e2;rhy[790:1150]=np.nan;rhy[2290:2620]=np.nan
ax.plot(rh_time,rhy,'-',label='RHESSI (8-40 keV)')
mwatb=Tbsubmax[-1]/1.e4;mwatb[627:631]=np.nan;mwatb[378]=np.nan;mwatb[mwatb==0]=np.nan;mwatb[631]=np.nan
ax.plot(tsubmwa,mwatb,'-',markersize=1,label='MWA (240 MHz)');ax.set_ylim(-2000,25000)
ax.legend();ax.set_ylabel('Amplitude (arbitrary)');ax.set_xlabel('Time (HH:MM UT)')
#ax.set_xticks(np.arange(5)*3600+5400);ax.set_xticklabels(['01:30','02:30','03:30','04:30','05:50'])
plt.show()



