# Extrapolation
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from sunpy.map import Map
from reproject import reproject_interp
import sunpy
from sunpy.coordinates import get_body_heliographic_stonyhurst
from scipy.io import readsav
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

exfile='/media/rohit/VLA/20160409_EUV/hmi.M_720s.20160409_183417.E18N10CR.CEA.NAS.sav'
mapex=readsav(exfile);expolB=mapex['box'];bx,by,bz=expolB['bx'][0],expolB['by'][0],expolB['bz'][0];babs=np.sqrt(bx*bx+by*by+bz*bz);dr=expolB['dr']
index=expolB['index'][0];crval1=index['crval1'];crval2=index['crval2'];crpix1=index['crpix1'];crpix2=index['crpix2'];ctype1=index['ctype1'];ctype2=index['ctype2'];cdelt1=index['cdelt1'];cdelt2=index['cdelt2'];cunit1=index['cunit1'];cunit2=index['cunit2']
hdu = fits.PrimaryHDU(babs[0]);list_all=list(index.dtype.names);list_all.remove('COMMENT');list_all.remove('HISTORY');list_all.remove('SIMPLE');list_all.remove('BITPIX');list_all.remove('NAXIS');list_all.remove('DATE_D$OBS')
index['WCSNAME'],index['CTYPE1'],index['CUNIT1'],index['CTYPE2'],index['CUNIT2']=['Carrington-Heliographic'],['CRLN-CEA'],['deg'],['CRLT-CEA'],['deg']
index['DATE_OBS']=['2016-04-09T18:34:17.00']
ii=0
for idx in list_all:
    hdu.header.update({str(idx):index[list_all[ii]][0]})
    ii=ii+1

outlist=[0]*200
for i in range(200):
    hdu.data=bx[0];hhdu=hdu.header;hdul = fits.HDUList([hdu]);mymap=Map(bx[i],hhdu);hp_coord=mymap.reference_coordinate.transform_to(frames.Helioprojective(observer="earth"))
    mapex_center=[-800,220]
    out_shape = (300, 300)
    out_header = sunpy.map.make_fitswcs_header(out_shape,hp_coord)
    out_wcs = WCS(out_header)
    earth = get_body_heliographic_stonyhurst('earth', mymap.date)
    out_wcs.heliographic_observer = earth
    output, footprint = reproject_interp(mymap, out_wcs, out_shape)
    outmap = sunpy.map.Map((output, out_header))
    outmap.plot_settings = mymap.plot_settings
    outlist[i]=outmap
pickle.dump(outlist,open('/media/rohit/VLA/20160409_EUV/20160409_bx.p','wb'),protocol=2)

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection=mymap)
mymap.plot(axes=ax1)
outmap.draw_grid(color='w')
ax2 = fig.add_subplot(1, 2, 2, projection=outmap)
outmap.plot(axes=ax2)
outmap.draw_grid(color='w')
plt.show()


