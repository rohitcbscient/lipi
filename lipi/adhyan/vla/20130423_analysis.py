## Script for MARINA (20130423 EVENT)

import sunpy.map
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.io import readsav
import glob
import itertools

def hms2sec_c(time):
    '''
    Input: time (HH:MM:SS)
    Output: time in sec
    '''
    sec=int(time.split(':')[0])*3600+int(time.split(':')[1])*60+float(time.split(':')[2])
    return sec

def find_predecessor(array, value):
    '''
    Input:
    array
    value
    Output:
    index of closest predecessor
    Array value of closest predecessor
    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if(array[idx]-value>0):
        idx=idx-1
    return idx,array[idx]

########## AIA #####################

def read_submap(f,ff):
    '''
    Objective: Read AIA submaps
    Inputs: filename, aia wavelength as integer
    Outputs: submap, time (in sec), time string, center x, center y, dx, dy, x array, y array
    '''
    data_=readsav(f)
    data=data_['submap'+str(ff)]
    n=len(data)
    ts=[0]*n
    map_=[0]*n
    ti=[0]*n
    xarr=[0]*n
    yarr=[0]*n
    for i in range(n):
        map_[i]=data[i][0]
        mx=data[i][0].shape[0]
        my=data[i][0].shape[1]
        time=data[i][5]
        ts[i]=hms2sec_c(time.split(' ')[1])
        ti[i]=time.split(' ')[1]
        xarr[i]=np.linspace(data['xc'][i]-data['dx'][i]*(mx/2.),data['xc'][i]+data['dx'][i]*(mx/2.),mx)
        yarr[i]=np.linspace(data['yc'][i]-data['dy'][i]*(my/2.),data['yc'][i]+data['dy'][i]*(my/2.),my)
    return map_,ts,ti,data[0][1],data[0][2],data[0][3],data[0][4],np.array(xarr),np.array(yarr)

def get_submap(aiafile,w):
    aiamap,aiats,aiatime,xc,yc,dx,dy,xarrayaia,yarrayaia=read_submap(aiafile,w)
    return aiamap,xarrayaia,yarrayaia,xc,yc,dx,dy,aiats,aiatime

do_flux_density=1
if(do_flux_density):
    FBAND=['L','S']
    POL_LIST=['LL','RR']
    LSPW_LIST=['0','1','2','3','4','5','6','7']
    SSPW_LIST=['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']
    imgLLL=[0]*len(LSPW_LIST)
    imgSLL=[0]*len(SSPW_LIST)
    imgLRR=[0]*len(LSPW_LIST)
    imgSRR=[0]*len(SSPW_LIST)
    i=0
    for SPW in LSPW_LIST:
        imgLLL[i]=sorted(glob.glob('*.L.*.LL.spw.'+SPW+'*.FITS'),key=lambda f: int(filter(str.isdigit, f)))
        imgLRR[i]=sorted(glob.glob('*.L.*.RR.spw.'+SPW+'*.FITS'),key=lambda f: int(filter(str.isdigit, f)))
        i=i+1
    i=0
    for SPW in SSPW_LIST:
        imgSLL[i]=sorted(glob.glob('*.S.*.LL.spw.*'+SPW+'.FITS'),key=lambda f: int(filter(str.isdigit, f)))
        imgSRR[i]=sorted(glob.glob('*.S.*.RR.spw.*'+SPW+'.FITS'),key=lambda f: int(filter(str.isdigit, f)))
        i=i+1
    img_list_LL=list(itertools.chain(*imgLLL))+list(itertools.chain(*imgSLL))
    img_list_RR=list(itertools.chain(*imgLRR))+list(itertools.chain(*imgSRR))
    sum_flux_jyperbeam=[0]*len(img_list_LL)
    sum_flux=[0]*len(img_list_LL)
    bmin=[0]*len(img_list_LL)
    bmaj=[0]*len(img_list_LL)
    freq=[0]*len(img_list_LL)
    ii=0
    for img in img_list_LL:
        img_=fits.open(img)
        bmin[ii]=img_[0].header['BMIN'] # in degrees
        bmaj[ii]=img_[0].header['BMAJ']
        freq[ii]=img_[0].header['RESTFRQ']/1.e9
        Tb=img_[0].data
        #omega=np.pi*bmin[ii]*bmin[ii]*(np.pi/180)**2/(4*log(2.))
        omega=np.pi*0.5*0.5*(np.pi/180)**2/(4*log(2.))
        nbeam=(256*2/bmin[ii])*(256*2/bmaj[ii]) # 256 = number of pixels in image, 2" is cell size
        npix_per_beam=np.pi*0.25*(bmin[ii]/2.)*(bmaj[ii]/2.) # circle with diameter bmin/2 in pixels
        flux=Tb*freq[ii]*freq[ii]*bmin[ii]*bmaj[ii]*3600*3600/(1222*1.e7)
        sum_flux_jyperbeam[ii]=np.nansum(flux)/npix_per_beam *nbeam
        sum_flux[ii]=sum_flux_jyperbeam[ii]
        #sum_flux[ii]=np.nansum(Tb)*(2*1.38e-23*1.4e9*1.4e9)*1.e26*1.e-4/(3.e8*3.e8)*omega
        ii=ii+1
    sum_flux=np.array(sum_flux)
    plot_flux=1
    if(plot_flux):
        plt.plot(freq,sum_flux_jyperbeam,'o')
        plt.ylabel('Total Flux (SFU/beam)')
        #plt.plot(freq,sum_flux,'o-')
        #plt.ylabel('Total Flux (SFU)')
        plt.xlabel('Frequency (GHz)')
        plt.show()
    plot_beam=1
    if(plot_beam):
        plt.plot(freq,np.array(bmin)*3600,'o-',label='BMIN')
        plt.plot(freq,np.array(bmaj)*3600,'o-',label='BMAJ')
        plt.ylabel('BMIN/BMAJ (in arcsec)')
        plt.xlabel('Frequency (GHz)')
        plt.legend()
        plt.show()


wav=94
aiapath='/nas08-data02/aiadata/20130423/maps/'
filename=aiapath+'aia_'+str(wav)+'_sub_map.sav'
aiamap,xarraia,yarraia,xcaia,ycaia,dxaia,dyaia,aiats,aiatime=get_submap(filename,wav)

########## RSTN ####################
rstn=readsav('sagamore_hill_fluxes.sav')
rstn_times=rstn['times']-rstn['times'][0]
rstn_1415=rstn['sfu_1415mhz']
rstn_2695=rstn['sfu_2695mhz']
# BASE time: 09:50:40 UT
basetime_sec=9*3600+50*60+40
start_time=20*3600 # 20 hours
end_time=21*3600
st_idx=find_predecessor(basetime_sec+rstn_times,start_time)[0]
ed_idx=find_predecessor(basetime_sec+rstn_times,end_time)[0]
plt_x=rstn_times[st_idx:ed_idx]-rstn_times[st_idx:ed_idx][0]
plt.plot(plt_x,rstn_2695[st_idx:ed_idx],'o-')
plt.xticks([0,900,1800,2700,3600],['20:00','20:15','20:30','20:45','21:00'])
plt.xlabel('Time (UT)')
plt.ylabel('RSTN Flux (SFU)')
plt.title('Sagamore Hill Fluxes')
plt.show()

sys.exit()
########## VLA #####################
vlapath='/nas08-data02/vladata/20130423/images_203430_203431/'
prefix=['20130423T2030-2050','20130423T2030-2040']
band_list=['L','S']
pol_list=['LL','RR']
lspw_list=['0','1','2','3','4','5','6','7']
sspw_list=['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']
chan_list=['0~3','4~7','8~11','12~15','16~19','20~23','24~27','28~31','32~35','36~39','40~43','44~47','48~51','52~55','56~59','60~63']
tag=['selfcal','selfcal']

vlatime='20:34:30~20:34:31'

Tb=[0]*len(pol_list)
pp=0
for p in pol_list:
    Tb[pp]=[0]*len(band_list)
    freq=[0]*len(band_list)
    bb=0
    for b in band_list:
        if(b=='L'):
            spw_list=lspw_list
        if(b=='S'):
            spw_list=sspw_list
        Tb[pp][bb]=[0]*len(spw_list)
        freq[bb]=[0]*len(spw_list)
        ss=0
        for s in spw_list:
            Tb[pp][bb][ss]=[0]*len(chan_list)
            freq[bb][ss]=[0]*len(chan_list)
            cc=0
            for c in chan_list:
                filename=prefix[bb]+'.'+b+'.50ms.'+tag[bb]+'.pol.'+p+'.spw.'+s+':'+c+'.time.'+vlatime+'.FITS'
                fits_data=fits.open(vlapath+filename)
                Tb[pp][bb][ss][cc]=fits_data[0].data[0][0]
                header_=fits_data[0].header
                freq[bb][ss][cc]=header_['RESTFRQ']
                ndim=Tb[pp][bb][ss][cc].shape
                xarrvla=np.linspace(header_['CRVAL1']-header_['CDELT1']*(ndim[0]/2.),header_['CRVAL1']+header_['CDELT1']*(ndim[0]/2.),ndim[0])
                yarrvla=np.linspace(header_['CRVAL2']-header_['CDELT2']*(ndim[1]/2.),header_['CRVAL2']+header_['CDELT2']*(ndim[1]/2.),ndim[1])
                cc=cc+1
            ss=ss+1
        bb=bb+1
    pp=pp+1

freq=np.array(freq)
freq_all=np.concatenate((np.array(freq[0]).flatten(),np.array(freq[1]).flatten()))
Tb_all=[0]*2
Tb_all[0]=np.concatenate((np.array(Tb[0][0]).reshape(np.array(Tb[0][0]).shape[0]*np.array(Tb[0][0]).shape[1],np.array(Tb[0][0]).shape[2],np.array(Tb[0][0]).shape[3]),np.array(Tb[0][1]).reshape(np.array(Tb[0][1]).shape[0]*np.array(Tb[0][1]).shape[1],np.array(Tb[0][1]).shape[2],np.array(Tb[0][1]).shape[3])),axis=0)
Tb_all[1]=np.concatenate((np.array(Tb[1][0]).reshape(np.array(Tb[1][0]).shape[0]*np.array(Tb[1][0]).shape[1],np.array(Tb[1][0]).shape[2],np.array(Tb[1][0]).shape[3]),np.array(Tb[1][1]).reshape(np.array(Tb[1][1]).shape[0]*np.array(Tb[1][1]).shape[1],np.array(Tb[1][1]).shape[2],np.array(Tb[1][1]).shape[3])),axis=0)
Tb_all=np.array(Tb_all)
Tb_rmd_LL=[0]*24
Tb_rmd_RR=[0]*24
eTb_rmd_LL=[0]*24
eTb_rmd_RR=[0]*24
freq_rmd=[0]*24
mm=0
for i in range(24):
    Tb_rmd_LL[i]=Tb_all[0][i*16+mm:i*16+16-mm]
    Tb_rmd_RR[i]=Tb_all[1][i*16+mm:i*16+16-mm]
    eTb_rmd_LL[i]=Tb_all[0][i*16+mm:i*16+16-mm][:,110:150,0:30].std(axis=(1,2))
    eTb_rmd_RR[i]=Tb_all[1][i*16+mm:i*16+16-mm][:,110:150,0:30].std(axis=(1,2))
    freq_rmd[i]=freq_all[i*16+mm:i*16+16-mm]
nfrq=384-mm*24
Tb_rmd_LL=np.array(Tb_rmd_LL).reshape(nfrq,256,256)    
Tb_rmd_RR=np.array(Tb_rmd_RR).reshape(nfrq,256,256)    
freq_rmd=np.array(freq_rmd).flatten()
eTb_rmd_LL=np.array(eTb_rmd_LL).flatten()
eTb_rmd_RR=np.array(eTb_rmd_RR).flatten()

Tb_source_LL=[0]*nfrq
Tb_max_LL=[0]*nfrq
Tb_source_RR=[0]*nfrq
Tb_max_RR=[0]*nfrq
for i in range(nfrq):
    ttLL=Tb_rmd_LL[i]
    ttLL[np.isnan(ttLL)]=0
    ttRR=Tb_rmd_RR[i]
    ttRR[np.isnan(ttRR)]=0
    idxLL=np.where(ttLL==ttLL.max())
    idxRR=np.where(ttRR==ttRR.max())
    Tb_source_LL[i]=np.nanmean(Tb_rmd_LL[i])#[idxLL[0][0]-15:idxLL[0][0]+15,idxLL[1][0]-15:idxLL[1][0]+15])
    Tb_source_RR[i]=np.nanmean(Tb_rmd_RR[i])#[idxRR[0][0]-15:idxRR[0][0]+15,idxRR[1][0]-15:idxRR[1][0]+15])
    Tb_source_LL[i]=np.nanmean(Tb_rmd_LL[i][90:140,90:150])
    Tb_source_RR[i]=np.nanmean(Tb_rmd_RR[i][90:140,90:150])
    Tb_max_LL[i]=np.nanmax(Tb_rmd_LL[i])
    Tb_max_RR[i]=np.nanmax(Tb_rmd_RR[i])
Tb_source_LL=np.array(Tb_source_LL)
Tb_max_LL=np.array(Tb_max_LL)
Tb_source_RR=np.array(Tb_source_RR)
Tb_max_RR=np.array(Tb_max_RR)
############### Analysis ###################

vlaid=hms2sec_c(vlatime.split('~')[0])
vla_aiaidx=find_predecessor(aiats,vlaid)[0]


############## PLOT ########################

plot_max=1 # Plot maximum of each image for LL and RR
if(plot_max):
    plt.plot(freq_rmd/1.e9,Tb_max_LL/1.e6,'o-',label='LL',color='red')
    plt.errorbar(freq_rmd/1.e9,Tb_max_LL/1.e6,yerr=eTb_rmd_LL/1.e6,color='red')
    plt.plot(freq_rmd/1.e9,Tb_max_RR/1.e6,'o-',label='RR',color='blue')
    plt.errorbar(freq_rmd/1.e9,Tb_max_RR/1.e6,yerr=eTb_rmd_RR/1.e6,color='blue')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('T$_{B}$ (MK)')
    plt.legend()
    plt.show()

plot_LS=1 # Plot maximum of each image for LL and RR
if(plot_LS):
    plt.plot(freq_rmd[0:128]/1.e9,Tb_max_LL[0:128]/1.e6,'o-',label='L',color='red')
    plt.errorbar(freq_rmd[0:128]/1.e9,Tb_max_LL[0:128]/1.e6,yerr=eTb_rmd_LL[0:128]/1.e6,color='red')
    plt.plot(freq_rmd[128:]/1.e9,Tb_max_LL[128:]/1.e6,'o-',label='S',color='blue')
    plt.errorbar(freq_rmd[128:]/1.e9,Tb_max_LL[128:]/1.e6,yerr=eTb_rmd_LL[128:]/1.e6,color='blue')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('T$_{B}$ (MK)')
    plt.legend()
    plt.show()

plot_aia_vla=0
imagepath='/nas08-data02/vladata/20130423/L-Band/pngs/'
if(plot_aia_vla):
    for i in range(Tb_all.shape[1]):
    #for i in range(1):
        lev_2=[0.5,0.6,0.7,0.8,0.9]
        fig=plt.figure(figsize=(6,6))
        ax=fig.add_subplot(111,aspect='auto')
        im=ax.imshow(aiamap[vla_aiaidx]/np.max(aiamap[vla_aiaidx]),origin=True,extent=[xarraia[vla_aiaidx][0],xarraia[vla_aiaidx][-1],yarraia[vla_aiaidx][0],yarraia[vla_aiaidx][-1]],cmap='sdoaia94',interpolation='none',vmin=0.00001,vmax=0.2)
        #plt.contourf(ccmap1/np.max(ccmap1),extent=[xl,xr,yl,yr],levels=lev_2,cmap='Oranges',alpha=0.1)
        ax.contour(Tb_all[0][i]/np.nanmax(Tb_all[0][i]),extent=[xarrvla[0],xarrvla[-1],yarrvla[0],yarrvla[-1]],linewidths=2,levels=lev_2,colors='yellow',label=str(np.round(freq_all[0]/1.e9,2))+' GHz')
        ax.contour(Tb_all[1][i]/np.nanmax(Tb_all[1][i]),extent=[xarrvla[0],xarrvla[-1],yarrvla[0],yarrvla[-1]],linewidths=2,levels=lev_2,colors='magenta',label=str(np.round(freq_all[0]/1.e9,2))+' GHz')
        plt.text(0.1,0.1,'LL',color='yellow',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes)
        plt.text(0.1,0.2,'RR',color='magenta',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes)
        #ss=ax.scatter(xc,yc,c=freq,s=80,cmap='YlOrRd')
        #plt.colorbar(ss,label='Frequency (MHz)')
        plt.title('AIA: '+aiatime[vla_aiaidx]+'  VLA: '+vlatime+' FREQ: '+str(np.round(freq_all[i]/1.e9,3))+' GHz')
        ax.set_xlabel('arcsec')
        ax.set_ylabel('arcsec')
        ax.set_ylim(yarraia[vla_aiaidx][0],yarraia[vla_aiaidx][-1])
        ax.set_xlim(xarraia[vla_aiaidx][0],xarraia[vla_aiaidx][-1])
        #plt.title(t)
        ax.grid(True)
        plt.savefig(imagepath+'20130423_AIA94_LLRR_'+str(np.round(freq_all[i]/1.e9,3))+'GHz_'+vlatime+'_'+str("%03d"%i)+'.png')
        plt.close()



