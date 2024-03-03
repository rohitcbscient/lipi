from scipy.io import readsav
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import numpy as np

goes=readsav('/data/Dropbox/20160409/fermi/idlsave_goes.sav')
gf0540=goes['lx'][1];gf1080=goes['lx'][0];gtime=goes['tarray']


fm=fits.open('/data/Dropbox/20160409/fermi/glg_cspec_n5_160409_v00_data.fits')
fm0=fm[0];fm1=fm[1];fm2=fm[2]
print('List of keywords: '+str(fm1.data.names))
fmrate=fm1.data['RATE'];fmtime=fm1.data['TIME']
fmtime=fmtime-fmtime[0]
plt.plot(fmtime,fmrate,'o-')
plt.show()


