import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
import os
import gc
import analysis_spikes_module as pl
from surya.vla import main as vla

#spiketime:list_[16],chan[28],spw=1,f=93,l=16

spw=['0','1']
chan=list(np.arange(64))
data=[0]*len(spw)*len(chan)
freq=[0]*len(spw)*len(chan)
xc_max=[0]*len(spw)*len(chan)
yc_max=[0]*len(spw)*len(chan)
header=[0]*len(spw)*len(chan)
time=[0]*len(spw)*len(chan)

f=0
for ss in spw:
    print 'SPW: ',ss
    for cc in chan:
        list_=sorted(glob.glob('*.LL.spw.'+ss+'_'+str(cc)+'-'+str(cc)+'.time.*.FITS'))[0:300]
        print 'SPW:',ss,'chan:',cc
        header[f]=[0]*len(list_)
        data[f]=[0]*len(list_)
        xc_max[f]=[0]*len(list_)
        yc_max[f]=[0]*len(list_)
        time=[0]*len(list_)
        time_sec=[0]*len(list_)
        l=0
        for ll in list_:    
            data[f][l]=fits.getdata(ll,0,memmap=False)[0][0]
            imsize=data[f][l].shape[0]
            h=fits.getheader(ll,0,memmap=False)
            header[f][l]=h
            freq[f]=h['CRVAL3']/1.e9
            xc_array=h['CRVAL1']+np.linspace(-imsize/2,imsize/2,imsize)*h['CDELT1']
            yc_array=h['CRVAL2']+np.linspace(-imsize/2,imsize/2,imsize)*h['CDELT2']
            yc_,xc_=np.where(data[f][l]==np.nanmax(data[f][l])) # Note the order is y first and x second
            xc_max[f][l],yc_max[f][l]=xc_array[xc_[0]],yc_array[yc_[0]]
            time[l]=ll.split('.time.')[1].split('.FITS')[0].split('-')[0]
            time_sec[l]=time[l].split(':')[-1]
            l=l+1
        f=f+1

freq=np.array(freq)
data=np.array(data).swapaxes(0,1)
time_sec=np.array(time_sec)
max_=np.nanmax(data,axis=(2,3))
xc_max1=np.array(xc_max).swapaxes(0,1)
yc_max1=np.array(yc_max).swapaxes(0,1)
datadump={'header':header,'data':data,'time':time,'time_sec':time_sec,'frequency':freq,'x_array':xc_array,'y_array':yc_array,'xc_max':xc_max1,'yc_max':yc_max1}
pickle.dump(datadump,open('/nas08-data02/vladata/20130423/L-Band/sub_spikes/20130423T2030-2050.L.50ms.selfcal.sub.pol.LL.p','wb'))

sys.exit()
#result = Fido.search(a.Time('2013-04-23T20:35:00', '2013-04-23T20:50:00'),a.Instrument("aia"), a.Wavelength(94*u.angstrom))
# file_download = Fido.fetch(result, site='ROB',path='.')

import glob
from astropy.io import fits
from sunpy.map import Map

list131=sorted(glob.glob('/media/rohit/VLA/20130423_spikes/aia*0131*.fits'))
list94=sorted(glob.glob('/media/rohit/VLA/20130423_spikes/aia*0094*.fits'))
list193=sorted(glob.glob('/media/rohit/VLA/20130423_spikes/aia*0193*.fits'))
list211=sorted(glob.glob('/media/rohit/VLA/20130423_spikes/aia*0211*.fits'))
list171=sorted(glob.glob('/media/rohit/VLA/20130423_spikes/aia*0171*.fits'))
int131=[0]*len(list131);ts131=[0]*len(list131);tstring131=[0]*len(list131);tb131=[0]*len(list131)
for i in range(len(list131)):
    m=Map(list131[i]);int131[i]=m.data
    ts131[i]=int131[i][2450:2500,3000:3040].max()
    tb131[i]=int131[i][2450,3000]
    tstring131[i]=list131[i].split('_')[3]
int131=np.array(int131)

int94=[0]*len(list94);ts94=[0]*len(list94);tstring94=[0]*len(list94);tb94=[0]*len(list94)
for i in range(len(list94)):
    m=Map(list94[i]);int94[i]=m.data
    ts94[i]=int94[i][2450:2500,3000:3040].max()
    tb94[i]=int94[i][2450,3000]
    tstring94[i]=list94[i].split('_')[3]
int94 = np.array(int94)

int193=[0]*len(list193);ts193=[0]*len(list193);tstring193=[0]*len(list193);tb193=[0]*len(list193)
for i in range(len(list193)):
    m=Map(list193[i]);int193[i]=m.data
    ts193[i]=int193[i][2450:2500,3000:3040].max()
    tb193[i]=int193[i][2450,3000]
    tstring193[i]=list193[i].split('_')[3]
int193 = np.array(int193)

int211=[0]*len(list211);ts211=[0]*len(list211);tstring211=[0]*len(list211);tb211=[0]*len(list211)
for i in range(len(list211)):
    m=Map(list211[i]);int211[i]=m.data
    ts211[i]=int211[i][2450:2500,3000:3040].max()
    tb211[i]=int211[i][2450,3000]
    tstring211[i]=list211[i].split('_')[3]
int211 = np.array(int211)

int171=[0]*len(list171);ts171=[0]*len(list171);tstring171=[0]*len(list171);tb171=[0]*len(list171)
for i in range(len(list171)):
    m=Map(list171[i]);int171[i]=m.data
    ts171[i]=int171[i][2450:2500,3000:3040].max()
    tb171[i]=int171[i][2450,3000]
    tstring171[i]=list171[i].split('_')[3]
int171 = np.array(int171)

f,ax=plt.subplots(1,1)
#ax.plot(ts131,'o-',label='131 $\AA$')
#ax.plot(tb131,'o-',label='BKG: 131 $\AA$')
#ax.set_xticks(np.arange(len(list131))[::10]);ax.set_xticklabels(tstring131[::10])
#ax.plot(ts94,'o-',label='94 $\AA$')
#ax.plot(tb94,'o-',label='BKG: 94 $\AA$')
#ax.set_xticks(np.arange(len(list94))[::10]);ax.set_xticklabels(tstring94[::10])
#ax.plot(ts193,'o-',label='193 $\AA$')
#ax.plot(tb193,'o-',label='BKG: 193 $\AA$')
#ax.set_xticks(np.arange(len(list193))[::10]);ax.set_xticklabels(tstring193[::10])
#ax.plot(ts211,'o-',label='211 $\AA$')
#ax.plot(tb211,'o-',label='BKG: 211 $\AA$')
#ax.set_xticks(np.arange(len(list211))[::10]);ax.set_xticklabels(tstring211[::10])
ax.plot(ts171,'o-',label='171 $\AA$')
ax.plot(tb171,'o-',label='BKG: 171 $\AA$')
ax.set_xticks(np.arange(len(list171))[::10]);ax.set_xticklabels(tstring171[::10])
ax.legend();ax.set_xlabel(' Time (HHMMSS) UT');ax.set_ylabel('DN')
plt.show()

import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
nan_idx=np.where(max_<2.0e6)
max_[nan_idx]=np.nan
xc_max1[nan_idx]=np.nan
yc_max1[nan_idx]=np.nan
fig =plt.figure()
spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec[0, :])
ax1.set_title('DS')
im1=ax1.imshow(max_/1.e6,extent=[0,15,freq[0],freq[-1]],aspect='auto',interpolation=None,cmap='jet',origin='lower')
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im1,cax=cax,orientation='vertical',label='MK')
ax2 = fig.add_subplot(spec[1, 0],sharex=ax1,sharey=ax1)
ax2.set_title('XC')
im2=ax2.imshow(xc_max1,extent=[0,15,freq[0],freq[-1]],aspect='auto',interpolation=None,cmap='jet',origin='lower',vmin=620,vmax=690)
divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im2,cax=cax,orientation='vertical',label='arcsec')
ax3 = fig.add_subplot(spec[1, 1],sharex=ax2,sharey=ax2)
ax3.set_title('YC')
im3=ax3.imshow(yc_max1,extent=[0,15,freq[0],freq[-1]],aspect='auto',interpolation=None,cmap='jet',origin='lower',vmin=240,vmax=260)
divider = make_axes_locatable(ax3)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im3,cax=cax,orientation='vertical',label='arcsec')
ax2.set_xlabel('Time (sec)')
ax3.set_xlabel('Time (sec)')
ax1.set_ylabel('Frequency (MHz)')
plt.show()


plt.hist(max_.flatten()/1.e6,bins=50,histtype='step')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$T_{B} (MK)$')
plt.show()





