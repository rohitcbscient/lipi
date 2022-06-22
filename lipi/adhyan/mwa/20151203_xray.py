import matplotlib.pyplot as plt
from astropy.io import fits


aa=fits.open('/data/Dropbox/20151203/x-ray/go1520151128.fits')
goes=aa[2].data[0][1];time=aa[2].data[0][0]

f,ax=plt.subplots(1,1)
ax.plot(time,goes[:,0]/1.e-6,label='0.5-4.0$\AA$',color='blue')
ax1=ax.twinx()
ax1.plot(time,goes[:,1]/1.e-9,label='1.0-8.0$\AA$',color='green')
ax.legend();ax.set_ylabel('Flux ($\mu$W/m$^2$)',color='blue');ax.set_xlabel('Time (HH:MM UT)')
ax1.legend(loc=2);ax1.set_ylabel('Flux (nW/m$^2$)',color='green')
ax.axvline(x=12240,color='k');ax.axvline(x=13140,color='k');ax.axvline(x=13440,color='k');ax.axvline(x=14340,color='k')
ax.axvspan(12240, 13140, alpha=0.2);ax.axvspan(13440, 14340, alpha=0.2)
ax.set_xlim([7200,18000]);ax.set_ylim([0.16,0.22])
ax.set_xticks([7200,9000,10800,12600,14400,16200,18000])
ax.set_xticklabels(['02:00','02:30','03:00','03:30','04:00','04:30','05:00'])
ax.spines['left'].set_color('blue');ax1.spines['right'].set_color('green')
ax.tick_params(axis='y', colors='blue');ax1.tick_params(axis='y', colors='green')
plt.show()

