import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

f='/nas08-data02/vladata/202104_vlastix/dynspec/20210319_183000_190000.1s.ms.dspec.fits'
aa=fits.open(f)
h,d=aa[0].header,aa[0].data

spec=np.load('/nas08-data02/vladata/202104_vlastix/20210322/rohit/SUN.100ms.cal.npz')
ds=spec['spec'];freq=spec['freq']/1.e9

#plt.imshow(ds[0][0],aspect='auto',origin=0)

plt.plot(freq,ds[0][0][:,300:400].mean(axis=1),'o-')
for i in range(8):
    plt.axvline(x=freq[int(i*128.)-5],color='k')
plt.xlabel('Frequency (GHz)');plt.ylabel('T$_B$(K)')
plt.show()

