import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
import os
import gc
import analysis_spikes_module as pl
from scipy.io import readsav
from surya.vla import main as vla
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from sunpy.map import Map
from astropy import units as u
from surya.utils import main as ut


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


data_dump_vla=0
if (data_dump_vla):
    spw=['0','1']
    chan=list(np.arange(64));Tbmax=[0]*len(spw)*len(chan);freq=[0]*len(spw)*len(chan)
    xc_max=[0]*len(spw)*len(chan);yc_max=[0]*len(spw)*len(chan);header=[0]*len(spw)*len(chan);time=[0]*len(spw)*len(chan)
    xf=[0]*len(spw)*len(chan);yf=[0]*len(spw)*len(chan)
    f=0
    for ss in spw:
        print 'SPW: ',ss
        for cc in chan:
            list_=sorted(glob.glob('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/images_LL/*.LL.spw.'+ss+'_'+str(cc)+'-'+str(cc)+'.time.*.FITS'))
            print 'SPW:',ss,'chan:',cc
            header[f]=[0]*len(list_);Tbmax[f]=[0]*len(list_);xc_max[f]=[0]*len(list_);yc_max[f]=[0]*len(list_);time=[0]*len(list_);time_sec=[0]*len(list_);xf[f]=[0]*len(list_);yf[f]=[0]*len(list_)
            l=0
            for ll in list_:    
                data=Map(ll)
                #imsize=data[f][l].shape[0]
                h=fits.getheader(ll,0,memmap=False)
                header[f][l]=h
                freq[f]=h['CRVAL3']/1.e9
                da=data.data;da[np.isnan(da)]=0;dab=ut.get_bimage(da,0.9);dfit=ut.fitEllipse(dab);xf[f][l]=data.pixel_to_world(dfit[0]*u.pixel,dfit[1]*u.pixel,1).Tx.value;yf[f][l]=data.pixel_to_world(dfit[0]*u.pixel,dfit[1]*u.pixel,1).Ty.value
                #xc_array=h['CRVAL1']+np.linspace(-imsize/2,imsize/2,imsize)*h['CDELT1']
                #yc_array=h['CRVAL2']+np.linspace(-imsize/2,imsize/2,imsize)*h['CDELT2']
                yc_,xc_=np.where(da==np.nanmax(da)) # Note the order is y first and x second
                xc_max[f][l],yc_max[f][l]=data.pixel_to_world(xc_*u.pixel,yc_*u.pixel,1).Tx.value[0],data.pixel_to_world(xc_*u.pixel,yc_*u.pixel,1).Ty.value[0]
                time[l]=ll.split('.time.')[1].split('.FITS')[0].split('-')[0]
                time_sec[l]=time[l].split(':')[-1]
                Tbmax[f][l]=np.nanmax(data.data)
                l=l+1
            f=f+1
    freq=np.array(freq);Tbmax=np.array(Tbmax);time_sec=np.array(time_sec)
    #max_=np.nanmax(data,axis=(2,3))
    xc_max1=np.array(xc_max);yc_max1=np.array(yc_max);xf=np.array(xf);yf=np.array(yf)
    datadump={'header':header,'Tbamax':Tbmax,'time':time,'time_sec':time_sec,'frequency':freq,'xc_max':xc_max1,'yc_max':yc_max1,'xf':xf,'yf':yf}
    pickle.dump(datadump,open('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/20130423T2030-2050.L.50ms.selfcal.sub.pol.LL.p','wb'))


list_=sorted(glob.glob('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/images_LL/*.LL.spw.0_30-30.time.*.FITS'))
d=[0]*5;c=['r','b','g','c','k'];ll=[31,101,196,257,314];xlaia,ylaia,xraia,yraia=407,5.5,919,517.6
for i in range(5):
    data=Map(list_[ll[i]]);data.data[np.isnan(data.data)]=0
    d[i]=data.data
    cs=plt.contour(d[i]/np.max(d[i]),levels=[0.8],extent=[xlaia,xraia,ylaia,yraia],aspect='auto',origin='lower',colors=c[i])
    cs.collections[0].set_label(str(np.round(ll[i]*0.05,3))+' sec')
plt.legend(loc=2);plt.show()

list_=glob.glob('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/images_LL/*.LL.spw.0*.time.20:45:25.8-20:45:25.85.FITS')+glob.glob('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/images_LL/*.LL.spw.1*.time.20:45:25.8-20:45:25.85.FITS')
d=[0]*5;c=['r','b','g','c','k'];ll=[5,42,73,92,105];xlaia,ylaia,xraia,yraia=407,5.5,919,517.6
for i in range(5):
    data=Map(list_[ll[i]]);data.data[np.isnan(data.data)]=0
    d[i]=data.data
    cs=plt.contour(d[i]/np.max(d[i]),levels=[0.95],extent=[xlaia,xraia,ylaia,yraia],aspect='auto',origin='lower',colors=c[i])
    cs.collections[0].set_label(str(np.round(freq[ll[i]],3))+' GHz')
plt.legend(loc=2);plt.show()

listr=sorted(glob.glob('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/images_LL/2013*'+'.spw.0_0-0*.FITS'))
tstr=[0]*len(listr)
for tt in range(len(listr)):
    tstr[tt]=listr[tt].split('-')[-2].split('time.')[1]

do_check=1
if(do_check):
    xf=[0]*2;yf=[0]*2;tbf=[0]*2;sigx0=[0]*2;sigy0=[0]*2
    for s in range(2):
        xf[s]=[0]*64;yf[s]=[0]*64;tbf[s]=[0]*64;sigx0[s]=[0]*64;sigy0[s]=[0]*64
        for c in range(64):
            xf[s][c]=[0]*348;yf[s][c]=[0]*348;tbf[s][c]=[0]*348;sigx0[s][c]=[0]*348;sigy0[s][c]=[0]*348
            for i in range(348):
                ii="%04d" % i
                bb=pickle.load(open('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/blob_fit/blob_param_spw_'+str(s)+'_'+str(c)+'_'+str(ii)+'.p','rb'))
                xf[s][c][i]=bb[0];yf[s][c][i]=bb[1];tbf[s][c][i]=bb[5];sigx0[s][c][i]=bb[2];sigy0[s][c][i]=bb[3]
    xf=np.array(xf);yf=np.array(yf);tbf=np.array(tbf);sigx0=np.array(sigx0);sigy0=np.array(sigy0)
    xf=xf.reshape(128,348);yf=yf.reshape(128,348);tbf=tbf.reshape(128,348);sigx0=sigx0.reshape(128,348);sigy0=sigy0.reshape(128,348)

#spiketime:list_[16],chan[28],spw=1,f=93,l=16
print 'Reading RR....'
dataRR=pickle.load(open('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/20130423T204525-204542.L.50ms.selfcal.sub.pol.RR.p','rb'))
print 'Reading LL....'
dataLL=pickle.load(open('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/20130423T2030-2050.L.50ms.selfcal.sub.pol.LL.p','rb'))
xc_max1_LL, yc_max1_LL = dataLL['xc_max'], dataLL['yc_max']

freq=dataLL['frequency']
time_sec=dataLL['time_sec']
time=dataLL['time']
max_LL=dataLL['Tbamax']
max_LL[47:49]=0
nan_idx=np.where(max_LL<1.0e6)
idx=np.where(max_LL>1.0e6)
list_t=glob.glob('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/images_LL/*.LL.spw.0_57-57*.FITS')
max_LL[nan_idx]=np.nan;xf[nan_idx]=np.nan;yf[nan_idx]=np.nan
xc_max1_LL[nan_idx]=0
yc_max1_LL[nan_idx]=0
tbf[nan_idx]=np.nan
rf=np.sqrt((xf-600)**2 + (yf-200)**2);rf[nan_idx]=0

xc_max1_RR, yc_max1_RR= dataRR['xc_max'],dataRR['yc_max']
max_RR=dataRR['Tbamax']
nan_idx=np.where(max_RR<1.0e6)
max_RR[nan_idx]=np.nan
xc_max1_RR[nan_idx]=0
yc_max1_RR[nan_idx]=0


f=plt.figure(figsize=(15,10))
ax=f.add_subplot(111)
im=ax.imshow(max_LL/1.e6,origin=0,aspect='auto',extent=[0,348,0,128])
ax.set_ylabel('Frequency (GHz)')
ax.set_xlabel('Time (HH:MM:SS UT)')
ax.set_xticks([0,50,100,150,200,250,300,348])
ax.set_xticklabels(['20:45:25.00','20:45:27.50','20:45:30.00','20:45:32.50','20:45:35.00','20:45:37.50','20:45:40.00','20:45:42.4'])
ax.set_yticks([0,20,40,60,80,100,120])
ax.set_yticklabels(['0.99','1.03','1.07','1.11','1.15','1.19','1.23'])
f.colorbar(im,label='T$_B$ (MK)')
plt.show()

sys.exit()
wav=94
#aiapath='/nas08-data02/aiadata/20130423/maps/'
#filename=aiapath+'aia_'+str(wav)+'_sub_map.sav'
#aiamap,xarraia,yarraia,xcaia,ycaia,dxaia,dyaia,aiats,aiatime=get_submap(filename,wav)
aiapath='/nas08-data02/aiadata/20130423/'
#filename=aiapath+'AIA20130423_204535_0171.fits';aiamap=Map(filename)
filename=aiapath+'AIA20130423_204525_0094.fits';aiamap=Map(filename)
#filename1=aiapath+'AIA20130423_204532_0171.fits';aiamap=Map(filename1)
hmifile='/nas08-data02/hmidata/20130423/hmi.m_45s.2013.04.23_20_35_15_TAI.magnetogram.fits'
hmimap=Map(hmifile)

yc_max1_LL[np.where((yc_max1_LL>288))]=np.nan;xc_max1_LL[np.where((xc_max1_LL>700))]=np.nan
yc_max1_LL[np.where((yc_max1_LL<200))]=np.nan;xc_max1_LL[np.where((xc_max1_LL<500))]=np.nan
yf[np.where((yf>288))]=np.nan;xf[np.where((xf>700))]=np.nan
yf[np.where((yf<200))]=np.nan;xf[np.where((xf<500))]=np.nan
xlaia,ylaia,xraia,yraia=407,5.5,919,517.6

f=plt.figure(figsize=(10, 10))
ax0 = f.add_subplot(111)
aiamap.plot(axes=ax0)
ax0.scatter(xf[:,0:80].flatten(),yf[:,0:80].flatten(),label='1st burst',alpha=0.5)
ax0.scatter(xf[:,80:145].flatten(),yf[:,80:145].flatten(),label='2nd burst',alpha=0.5)
ax0.scatter(xf[:,145:240].flatten(),yf[:,145:240].flatten(),label='3rd burst',alpha=0.5)
ax0.scatter(xf[:,240:].flatten(),yf[:,240:].flatten(),label='4th burst',alpha=0.5);ax0.legend()
ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
plt.show()

f=plt.figure(figsize=(10, 10))
ax0 = f.add_subplot(111)
aiamap.plot(axes=ax0)
ax0.errorbar(np.nanmean(xf[:,0:80].flatten()),np.nanmean(yf[:,0:80].flatten()),xerr=3*np.nanstd(xf[:,0:80].flatten()),yerr=3*np.nanstd(xf[:,0:80].flatten()),label='1st burst',alpha=0.9)
ax0.errorbar(np.nanmean(xf[:,80:145].flatten()),np.nanmean(yf[:,80:145].flatten()),xerr=3*np.nanstd(xf[:,80:145].flatten()),yerr=3*np.nanstd(xf[:,80:145].flatten()),label='2nd burst',alpha=0.9)
ax0.errorbar(np.nanmean(xf[:,145:240].flatten()),np.nanmean(yf[:,145:240].flatten()),xerr=3*np.nanstd(xf[:,145:240].flatten()),yerr=3*np.nanstd(xf[:,145:240].flatten()),label='3rd burst',alpha=0.9)
ax0.errorbar(np.nanmean(xf[:,240:].flatten()),np.nanmean(yf[:,240:].flatten()),xerr=3*np.nanstd(xf[:,240:].flatten()),yerr=3*np.nanstd(xf[:,240:].flatten()),label='4th burst',alpha=0.9)
ax0.legend()
ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
plt.show()

f=plt.figure(figsize=(10, 10))
ax0 = f.add_subplot(111)
aiamap.plot(axes=ax0,cmap='coolwarm')
ax0.errorbar(np.nanmean(xf[0:40].flatten()),np.nanmean(yf[0:40].flatten()),xerr=3*np.nanstd(xf[0:40].flatten()),yerr=3*np.nanstd(yf[0:40].flatten()),label='0.99-1.07 GHz',alpha=0.9)
ax0.errorbar(np.nanmean(xf[80:].flatten()),np.nanmean(yf[40:80].flatten()),xerr=3*np.nanstd(xf[40:80].flatten()),yerr=3*np.nanstd(yf[40:80].flatten()),label='1.07-1.15 GHz',alpha=0.9)
ax0.errorbar(np.nanmean(xf[80:].flatten()),np.nanmean(yf[80:].flatten()),xerr=3*np.nanstd(xf[80:].flatten()),yerr=3*np.nanstd(yf[80:].flatten()),label='1.15-1.25 GHz',alpha=0.9)
ax0.legend()
ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
plt.show()

f=plt.figure(figsize=(10, 10))
ax0 = f.add_subplot(111)
hmimap.plot(axes=ax0,vmin=-40,vmax=40,cmap='binary')
ax0.errorbar(np.nanmean(xf[:,0:80].flatten()),np.nanmean(yf[:,0:80].flatten()),xerr=3*np.nanstd(xf[:,0:80].flatten()),yerr=3*np.nanstd(xf[:,0:80].flatten()),label='1st burst',alpha=0.9)
ax0.errorbar(np.nanmean(xf[:,80:145].flatten()),np.nanmean(yf[:,80:145].flatten()),xerr=3*np.nanstd(xf[:,80:145].flatten()),yerr=3*np.nanstd(xf[:,80:145].flatten()),label='2nd burst',alpha=0.9)
ax0.errorbar(np.nanmean(xf[:,145:240].flatten()),np.nanmean(yf[:,145:240].flatten()),xerr=3*np.nanstd(xf[:,145:240].flatten()),yerr=3*np.nanstd(xf[:,145:240].flatten()),label='3rd burst',alpha=0.9)
ax0.errorbar(np.nanmean(xf[:,240:].flatten()),np.nanmean(yf[:,240:].flatten()),xerr=3*np.nanstd(xf[:,240:].flatten()),yerr=3*np.nanstd(xf[:,240:].flatten()),label='4th burst',alpha=0.9)
ax0.legend()
ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
plt.show()

f=plt.figure(figsize=(10, 10))
ax0 = f.add_subplot(111)
hmimap.plot(axes=ax0,vmin=-40,vmax=40,cmap='coolwarm')
ax0.scatter(xc_max1_LL[0:40].flatten(),yc_max1_LL[0:40].flatten(),label='0.99-1.07 GHz',alpha=0.5)
ax0.scatter(xc_max1_LL[40:80].flatten(),yc_max1_LL[40:80].flatten(),label='1.07-1.15 GHz',alpha=0.5)
ax0.scatter(xc_max1_LL[80:].flatten(),yc_max1_LL[80:].flatten(),label='1.15-1.25 GHz',alpha=0.5)
ax0.legend()
ax0.set_xlim([xlaia,xraia]);ax0.set_ylim([ylaia,yraia])
plt.show()

max_LL[np.isnan(max_LL)]=0
cs=plt.contour(max_LL,levels=[1],extent=[0,17.4,freq[0],freq[-1]])
csp=cs.collections[0].get_paths()
nfeatures=len(csp);cspall=[0]*nfeatures;ffreq=[0]*nfeatures;ftime=[0]*nfeatures;fdfreq=[0]*nfeatures;fdtime=[0]*nfeatures;fdr=[0]*nfeatures
i=0
for c in csp:
    cspall[i]=csp[i].vertices
    fdfreq[i]=cspall[i][:,1].max()-cspall[i][:,1].min()
    fdtime[i]=cspall[i][:,0].max()-cspall[i][:,0].min()
    ffreq[i]=np.mean(cspall[i][:,1])
    ftime[i]=np.mean(cspall[i][:,0])
    fdr[i]=fdfreq[i]/fdtime[i]
    i=i+1
fdfreq=np.array(fdfreq)
ffreq=np.array(ffreq)
ftime=np.array(ftime)
fdtime=np.array(fdtime)
fdr=np.array(fdr)

f=plt.figure(figsize=(15,10))
ax=f.add_subplot(111)
plt.contour(max_LL,levels=[1],extent=[0,348,0,128])
ax.set_ylabel('Frequency (GHz)')
ax.set_xlabel('Time (HH:MM:SS UT)')
ax.set_xticks([0,50,100,150,200,250,300,348])
ax.set_xticklabels(['20:45:25.00','20:45:27.50','20:45:30.00','20:45:32.50','20:45:35.00','20:45:37.50','20:45:40.00','20:45:42.4'])
ax.set_yticks([0,20,40,60,80,100,120])
ax.set_yticklabels(['0.99','1.03','1.07','1.11','1.15','1.19','1.23'])
plt.show()


f,ax=plt.subplots(2,2)
ax[0,0].hist(ffreq,bins=15)
ax[0,0].set_xlabel('Frequency (GHz)')
ax[0,1].hist(fdfreq,bins=15)
ax[0,1].set_xlabel('Frequency (GHz)')
ax[1,0].hist(fdtime,bins=15)
ax[1,0].set_xlabel('Time (sec)')
ax[1,1].hist(fdr[fdr<0.4],bins=15)
ax[1,1].set_xlabel('df/dt (GHz/sec)')
plt.show()



for i in range(len(xc_max1_LL)):
    fig,ax =plt.subplots()
    ax.imshow(np.log10(aiamap[152]),aspect='auto',extent=[xarraia[0][0],xarraia[0][-1],yarraia[0][0],yarraia[0][-1]],origin=0,cmap='YlOrBr')
    ax.set_title('AIA: 20:45:37 UT VLA: 20:45:25 + '+str(np.round(float(time_sec[i])-25.0,4))+' sec')
    cmap = mpl.cm.get_cmap('jet')
    norm = mpl.colors.Normalize(vmin=min(freq), vmax=max(freq),clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='jet')
    time_color = np.array([(mapper.to_rgba(v)) for v in freq])
    ss=plt.scatter(xc_max1_LL[i],yc_max1_LL[i],c=freq,s=90,cmap=cmap,alpha=0.5)
    #for x, y, color in zip(xc_max1[:,i], yc_max1[:,i], time_color):
    #ax.plot(x, y, 'o', color=color,markersize=6,alpha=0.4)
    #ax.errorbar(x, y, xerr=ex, yerr=ey, lw=1, capsize=3, color=color)
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_xlim([xarraia[0][0],xarraia[0][-1]])
    ax.set_ylim([yarraia[0][0],yarraia[0][-1]])
    plt.colorbar(ss,label='Frequency (GHz)')
    ii="%03d" %i
    plt.savefig('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/png/vla_cent_aia94_'+ii+'.png',dpi=70)
    plt.close()

for i in range(348):
    ii="%03d" %i
    f=plt.figure(figsize=(14, 14))
    ax0 = f.add_subplot(111)
    aiamap.plot(axes=ax0,cmap='coolwarm')
    cmap = mpl.cm.get_cmap('jet')
    norm = mpl.colors.Normalize(vmin=min(freq), vmax=max(freq),clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap='jet')
    time_color = np.array([(mapper.to_rgba(v)) for v in freq])
    ss=plt.scatter(xf[:,i],yf[:,i],c=freq,s=90,cmap=cmap,alpha=0.9,vmin=1.0,vmax=1.2)
    ax0.errorbar(np.nanmean(xf[:,i].flatten()),np.nanmean(yf[:,i].flatten()),xerr=3*np.nanstd(xf[:,i].flatten()),yerr=3*np.nanstd(yf[:,i].flatten()),color='k',alpha=0.9)
    #ax0.set_xlim([570,660]);ax0.set_ylim([210,310])
    ax0.set_xlim([570,600]);ax0.set_ylim([220,250])
    ax0.set_title('AIA: 20:45:37 UT VLA: 20:45:25 + '+str(np.round(float(time_sec[i])-25.0,4))+' sec')
    ax0.set_xlabel('arcsec');ax0.set_ylabel('arcsec');v1 = np.linspace(1, 1.2, 8, endpoint=True)
    cbar=f.colorbar(ss,label='Frequency (GHz)',ticks=v1)
    cbar.ax.set_yticklabels(["{:1.2f}".format(i) for i in v1])
    plt.savefig('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/png/vla_cent_aia94_'+ii+'.png',dpi=100)
    plt.close()

for i in range(348):
    ii="%03d" %i
    f=plt.figure(figsize=(14, 14))
    ax0 = f.add_subplot(111)
    aiamap.plot(axes=ax0,cmap='coolwarm')
    cmap = mpl.cm.get_cmap('jet')
    norm = mpl.colors.Normalize(vmin=min(freq), vmax=max(freq),clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap='jet')
    time_color = np.array([(mapper.to_rgba(v)) for v in freq])
    ss=plt.scatter(xf[:,i],yf[:,i],c=freq,s=90,cmap=cmap,alpha=0.9,vmin=1.0,vmax=1.2)
    ax0.errorbar(np.nanmean(xf[:,i].flatten()),np.nanmean(yf[:,i].flatten()),xerr=3*np.nanstd(xf[:,i].flatten()),yerr=3*np.nanstd(yf[:,i].flatten()),color='k',alpha=0.9)
    #ax0.set_xlim([570,660]);ax0.set_ylim([210,310])
    ax0.set_xlim([570,600]);ax0.set_ylim([220,250])
    ax0.set_title('AIA: 20:45:37 UT VLA: 20:45:25 + '+str(np.round(float(time_sec[i])-25.0,4))+' sec')
    ax0.set_xlabel('arcsec');ax0.set_ylabel('arcsec');v1 = np.linspace(1, 1.2, 8, endpoint=True)
    cbar=f.colorbar(ss,label='Frequency (GHz)',ticks=v1)
    cbar.ax.set_yticklabels(["{:1.2f}".format(i) for i in v1])
    plt.savefig('/nas08-data02/vladata/20130423/L-Band/sub_spikes_new/png/vla_cont_aia94_'+ii+'.png',dpi=100)
    plt.close()

plot_rr=1
if(plot_rr):
    maxcc,xcc,ycc=max_RR.swapaxes(0,1)/1.e6,xc_max1_RR.swapaxes(0,1),yc_max1_RR.swapaxes(0,1)
    fig =plt.figure()
    spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
    ax1 = fig.add_subplot(spec[0, :])
    ax1.set_title('DS')
    im1=ax1.imshow(maxcc,extent=[0,15,freq[0],freq[-1]],aspect='auto',interpolation=None,cmap='jet',origin='lower')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1,cax=cax,orientation='vertical',label='MK')
    ax2 = fig.add_subplot(spec[1, 0],sharex=ax1,sharey=ax1)
    ax2.set_title('XC')
    im2=ax2.imshow(xcc,extent=[0,15,freq[0],freq[-1]],aspect='auto',interpolation=None,cmap='jet',origin='lower',vmin=620,vmax=690)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im2,cax=cax,orientation='vertical',label='arcsec')
    ax3 = fig.add_subplot(spec[1, 1],sharex=ax2,sharey=ax2)
    ax3.set_title('YC')
    im3=ax3.imshow(ycc,extent=[0,15,freq[0],freq[-1]],aspect='auto',interpolation=None,cmap='jet',origin='lower',vmin=240,vmax=260)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im3,cax=cax,orientation='vertical',label='arcsec')
    ax2.set_xlabel('Time (sec)')
    ax3.set_xlabel('Time (sec)')
    ax1.set_ylabel('Frequency (MHz)')
    plt.show()


plot_xc_yc=1
if(plot_xc_yc):
    fig =plt.figure()
    spec = gridspec.GridSpec(ncols=1, nrows=3, figure=fig)
    ax1 = fig.add_subplot(spec[0, :])
    ax1.set_title('DS')
    im1=ax1.imshow(tbf/1.e6,extent=[0,15,freq[0],freq[-1]],aspect='auto',interpolation=None,cmap='coolwarm',origin='lower')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1,cax=cax,orientation='vertical',label='MK')
    ax2 = fig.add_subplot(spec[1, :],sharex=ax1,sharey=ax1)
    ax2.set_title('XC')
    im2=ax2.imshow(xf,extent=[0,15,freq[0],freq[-1]],aspect='auto',interpolation=None,cmap='coolwarm',origin='lower',vmin=580,vmax=590)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im2,cax=cax,orientation='vertical',label='arcsec')
    ax3 = fig.add_subplot(spec[2, :],sharex=ax2,sharey=ax2)
    ax3.set_title('YC')
    im3=ax3.imshow(yf,extent=[0,15,freq[0],freq[-1]],aspect='auto',interpolation=None,cmap='coolwarm',origin='lower',vmin=240,vmax=250)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im3,cax=cax,orientation='vertical',label='arcsec')
    ax2.set_xlabel('Time (sec)')
    ax3.set_xlabel('Time (sec)')
    ax1.set_ylabel('Frequency (MHz)')
    plt.show()

plot_rf=1
if(plot_rf):
    fig =plt.figure()
    spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig)
    ax1 = fig.add_subplot(spec[0, :])
    ax1.set_title('DS')
    im1=ax1.imshow(tbf/1.e6,extent=[0,15,freq[0],freq[-1]],aspect='auto',interpolation=None,cmap='coolwarm',origin='lower')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1,cax=cax,orientation='vertical',label='MK')
    ax3 = fig.add_subplot(spec[1, :],sharex=ax2,sharey=ax2)
    ax3.set_title('Radial distance from (600,200)')
    im3=ax3.imshow(rf,extent=[0,15,freq[0],freq[-1]],aspect='auto',interpolation=None,cmap='coolwarm',origin='lower',vmin=40,vmax=45)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im3,cax=cax,orientation='vertical',label='arcsec')
    ax2.set_xlabel('Time (sec)')
    ax3.set_xlabel('Time (sec)')
    ax1.set_ylabel('Frequency (MHz)')
    plt.show()

plot_hist=0
if(plot_hist):
    plt.hist(max_.flatten()/1.e6,bins=50,histtype='step')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$T_{B} (MK)$')
    plt.show()





