from sunpy.map import Map
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav

hmifile='/media/rohit/VLA/20160409_EUV/hmi.m_45s.2016.04.09_18_45_00_TAI.magnetogram.fits'
hmimap=Map(hmifile)
hmid=hmimap.data[::-1,::-1]
hmid[np.where(hmid<-5000)]=0
hmimap=Map(hmid,hmimap.meta)


expolB_=readsav('/media/rohit/VLA/20160409_EUV/hmi.M_720s.20160409_183417.E18N10CR.CEA.NAS.sav')
expolB=expolB_['box']
bx,by,bz=expolB['bx'][0],expolB['by'][0],expolB['bz'][0]
babs=np.sqrt(bx*bx+by*by+bz*bz)
dr=expolB['dr'] # 3 values of extrapolated box in units Rsun
index=expolB['index'][0]
crval1=index['crval1']
crval2=index['crval2']
crpix1=index['crpix1']
crpix2=index['crpix2']
ctype1=index['ctype1']
ctype2=index['ctype2']
cdelt1=index['cdelt1']
cdelt2=index['cdelt2']
cunit1=index['cunit1']
cunit2=index['cunit2']
hdu = fits.PrimaryHDU(babs[0])
list_all=list(index.dtype.names)
list_all.remove('COMMENT')
list_all.remove('HISTORY')
list_all.remove('SIMPLE')
list_all.remove('BITPIX')
list_all.remove('NAXIS')
list_all.remove('DATE_D$OBS')
index['WCSNAME'],index['CTYPE1'],index['CUNIT1'],index['CTYPE2'],index['CUNIT2']=['Carrington-Heliographic'],['CRLN-CEA'],['deg'],['CRLT-CEA'],['deg']
index['DATE_OBS']=['2012-02-25T20:34:13.300']
ii=0
for idx in list_all:
        hdu.header.update({str(idx):index[list_all[ii]][0]})
        ii=ii+1

hdu.data=babs[0]
hhdu=hdu.header
hdul = fits.HDUList([hdu])
mymap=Map(babs[0],hhdu)

nrho = 60
rss = 2.5
input = pfsspy.Input(mymap, nrho, rss)
output = pfsspy.pfss(input)

