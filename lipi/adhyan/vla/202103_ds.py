import pickle
import numpy as np
import matplotlib.pyplot as plt
import glob

#plist=sorted(glob.glob('/media/rohit/VLA/20210321/base/'+'*.p'))
#d=pickle.load(open(plist[0],'rb'))
#amp=d['amplitude'].mean(axis=(0,2))
#time=d['axis_info']['time_axis']['MJDseconds']

#spec=np.load('/media/rohit/VLA/20210321/21A-118.sb39452573.eb39520716.59295.66536_203000_210000_1s.ms.dspec.npz')
spec=np.load('/media/rohit/VLA/20210321/21A-118.sb39452573.eb39520716.59295.66536_100ms_203700_203900.ms.dspec.npz')
s=spec['spec']

stix=np.loadtxt('/home/i4ds1807205/Dropbox/202103_VLA_STIX/20210322.txt')
tstix=np.arange(7200)*4
#tvla=15957+np.arange(1805)
1sec=0
if(1sec):
    tvla=15960+np.arange(1805)
tvla=
f=spec['freq']
#2021-03-19T16:04:03.197

vlaline=s[:,:,0:64,:].mean(axis=(1,2))+0.01

plot_1sec=0
if(plot_1sec):
    f,ax=plt.subplots(1,1)
    ax.plot(tstix,stix[:,0],'o-',label='STIX (4-10 keV)')
    ax.plot(tstix,stix[:,1],'o-',label='STIX (10-15 keV)')
    ax1=ax.twinx()
    ax1.plot(tvla,vlaline[0],'o-',color='r',label='VLA (LL)')
    ax1.plot(tvla,vlaline[1],'o-',color='yellow',label='VLA (RR)')
    ax.legend();ax1.legend(loc=4)
    ax.set_xlim(16166,17760)
    ax.set_xticks([16164,16164+480,16164+960,17760])
    ax.set_xticklabels(['20:33:19','20:41:19','20:49:19','20:59:55'])
    ax.set_xlabel('Time (UT)')
    ax.set_ylabel('Amplitude')
    plt.show()

lev=[0.05]
plt.imshow(s.mean(axis=(0,1)),origin=0,aspect='auto',extent=[0,120,f[0]/1.e9,f[-1]/1.e9],vmin=0.001,vmax=0.08,interpolation='None',cmap='YlGnBu')
#plt.contour(s.mean(axis=(0,1)),levels=lev,colors='r',extent=[0,120,f[0]/1.e9,f[-1]/1.e9])
plt.xlabel('Time (HH:MM:SS UT)'); plt.ylabel('Frequency (GHz)')
#plt.xticks([0,20,40,60,80,100,120],['20:37:00','20:37:20','20:37:40','20:38:00','20:38:20','20:38:40','20:39:00'])
plt.xticks([0,20,25,29,30,31,32,33,34,35,36,37,38,39,40,60,61,80,100,120],['20:37:00','20:37:20','20:37:25','20:37:29','20:37:30','20:37:31','20:37:32','20:37:33','20:37:34','20:37:35','20:37:36','20:37:37','20:37:38','20:37:39','20:37:40','20:38:00','20:38:01','20:38:20','20:38:40','20:39:00'])
plt.xlim([29,40])
plt.show()

