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





