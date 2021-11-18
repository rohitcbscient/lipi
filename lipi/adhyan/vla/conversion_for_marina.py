from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.map import Map
import numpy as np
from sunpy.coordinates import frames
from reproject import reproject_interp
from scipy.io import readsav
from astropy.io import fits
import sunpy
from astropy.wcs import WCS


expolB_=readsav('/home/i4ds1807205/Dropbox/20120225_VLA_work/gregory_gyrosynchroton/hmi.M_720s.20120225_203413.W133N16CR.CEA.NAS.sav')
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
    #print idx
    #hdu.header.append((idx,index[list_all[ii]][0],[]))
    hdu.header.update({str(idx):index[list_all[ii]][0]})
    ii=ii+1

hdu.data=babs[0]
hhdu=hdu.header
hdul = fits.HDUList([hdu])
mymap=Map(babs[0],hhdu)
hp_coord=mymap.reference_coordinate.transform_to(frames.Helioprojective(observer="earth"))
###################################

out_shape = (140, 200)
out_header = sunpy.map.make_fitswcs_header(
            out_shape,
            hp_coord)
out_wcs = WCS(out_header)
earth = get_body_heliographic_stonyhurst('earth', mymap.date)
out_wcs.heliographic_observer = earth
output, footprint = reproject_interp(mymap, out_wcs, out_shape)
outmap = sunpy.map.Map((output, out_header))
outmap.plot_settings = mymap.plot_settings

