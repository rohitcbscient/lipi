from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt

screen=fits.open('/home/rohit/karabo/karabo-pipeline/result/test_screen_60s.fits')
screen_data=screen[0].data

ion_with=fits.open('/home/rohit/karabo/karabo-pipeline/result/test_ion.fits')
ion_with_data=ion_with[0].data[0][0]
ion_without=fits.open('/home/rohit/karabo/karabo-pipeline/result/test_ion_without.fits')
ion_without_data=ion_without[0].data[0][0]


f,(ax0,ax1,ax2)=plt.subplots(1,3,sharex=True,sharey=True)
ax0.imshow(ion_with_data,aspect='auto',origin='lower',vmin=0.1,vmax=10)
ax1.imshow(ion_without_data,aspect='auto',origin='lower',vmin=0.1,vmax=10)
ax2.imshow(ion_without_data-ion_with_data,aspect='auto',origin='lower',vmin=-1.e-8,vmax=1.e-8)
ax0.set_title('With Screen')
ax1.set_title('Without Screen')
ax2.set_title('Difference')
plt.show()

i=0
for i in range(240):
    ii="%02d" %i
    f,ax=plt.subplots(1,1)
    ax.imshow(screen_data[0][i],aspect='auto',origin='lower')
    ax.set_xlabel('XX');ax.set_ylabel('YY')
    ax.set_title('t='+ii+' mins')
    f.savefig('/home/rohit/simulations/ionosphere/ion_'+ii+'.png',dpi=100)
    plt.close()

