import numpy as np
import matplotlib.pyplot as plt
import pickle
from astropy.io import fits

aa=fits.open('sun_20130423_L_1s.med.spec.fits')
bb=pickle.load(open('/home/i4ds1807205/Dropbox/20130423_marina/20130423_power.p','rb'))
idx=bb[0];tot=bb[1]*100*1.e-4/1388;max_=bb[2]*100*1.e-4
tidx=[0]*len(idx)
for i in range(len(idx)):
    tidx[i]=int(idx[i].split(':')[0])*3600+int(idx[i].split(':')[1])*60+int(idx[i].split(':')[2])
tidx=np.array(tidx)-tidx[0]+1883 # Start time: 20:00:00 UT 

ds=aa[2].data.mean(axis=(0,1))
freq=aa[1].data
time=aa[0].data
tsec=time-time[0] # Start time: 19:28:36.7 UT
tstring=['','','','','','','']
f2_idx=3952 # 20:34:30 ~ 20:34:31
f3_idx=4133 # 20:37:31 ~ 30:37:32
f4_idx=4336 # 20:40:54 ~ 20:40:55
f5_idx=4603 # 20:45:21 ~ 20:45:22
f6_idx=5247 # 20:56:05.5 ~ 20:56:06.5
back2_idx= 6922# 21:23:00.0~21:26:00.0

fig,ax = plt.subplots(3,1,sharex=True)
ax[0].plot(tsec,ds[-64:].mean(axis=0),'o',label='spw=7')
ax[1].plot(tidx,tot,'o',label='Total Power')
ax[2].plot(tidx,max_,'o',label='Maximum')
ax[0].axvline(x=f2_idx,color='r',label='Flare 3')
ax[0].axvline(x=f3_idx,color='g',label='Flare 4')
ax[0].axvline(x=f4_idx,color='cyan',label='Flare 5')
ax[0].axvline(x=f5_idx,color='magenta',label='Flare 6')
ax[0].axvline(x=f6_idx,color='orange',label='Flare 7')
ax[0].axvline(x=back2_idx,color='k',label='Background')
ax[0].legend();ax[1].legend();ax[2].legend()
ax[0].set_ylabel('Amplitude');ax[1].set_ylabel('SFU');ax[2].set_ylabel('SFU/Beam');ax[2].set_xlabel('Time (sec)/Start time: 19:28:37 UT')
ax[0].set_ylim([0.06,0.1]);ax[1].set_ylim([12,15]);ax[2].set_ylim([3.30,4.20]);ax[2].set_xlim([1943,8000])
plt.show()

