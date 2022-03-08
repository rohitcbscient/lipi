from astropy.io import fits
import numpy as np
import glob
from surya.utils import main as ut
import matplotlib.pyplot as plt

flist=glob.glob('/nas08-data02/vladata/20130423/L-Band/sub_spikes/images_LL_DC/*.FITS')
d=[0]*len(flist);maxidx=[0]*len(flist);idx90=[0]*len(flist)
i=0
for f in flist:
    aa=fits.open(f)
    d[i]=aa[0].data[0][0];maxidx[i]=np.where(d[i]==np.nanmax(d[i]))
    da=d[i]*1.0;da[np.isnan(d[i])]=0;dab=ut.get_bimage(da,0.9);dfit=ut.fitEllipse(dab);idx90[i]=[dfit[0],dfit[1]]
    i=i+1
d=np.array(d);idx90=np.array(idx90)
plt.plot(np.arange(len(flist))*0.05,(idx90[:,0]-idx90[:,0][0])*2,'o-');plt.ylabel('X-position (arcsec)');plt.xlabel('Time (s)');plt.title('SPW: 1:56~56')
plt.show()
