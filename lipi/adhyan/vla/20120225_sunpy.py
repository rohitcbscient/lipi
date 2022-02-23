import matplotlib.pyplot as plt
from reproject import reproject_interp
from sunpy.coordinates import frames
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import numpy as np
import sunpy.map 
from sunpy.map import Map
from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.net import Fido
from sunpy.net import attrs as a
from scipy.io import readsav
from astropy.io import fits
from surya.aia import main as aia
import glob
from surya.utils import maps as umaps
import plot_20120225 as pl

plt.style.use('/home/i4ds1807205/scripts/general/plt_style.py')
wavel=94
aiapath='/home/i4ds1807205/vla_data/20120225/aia/'
aiafile=aiapath+'aia_'+str(wavel)+'_sub_map.sav'
aiafile_all=sorted(glob.glob(aiapath+'aia*.sav'))
aiafile_all.sort(key=lambda ss: len(ss))
wavelength=wavel
wav=[94,131,171,193,211,304,335,1600,1700]
xl,xr,yl,yr=450,520,320,390
#xl,xr,yl,yr=460,500,330,370
#xl,xr,yl,yr=420,540,280,420 # For Quiet Sun plotting
res=0.6
aiaxarray=np.linspace(xl,xr,100)
aiayarray=np.linspace(yl,yr,100)
#cmap=aia.get_submap(aiafile,res,wavelength)
cmap,aiatsec,aiatstring=aia.get_submap_crop(aiafile,res,wavelength,xl,xr,yl,yr)
cmap_all,aiatsec_all,aiatstring_all=aia.get_submap_all(aiafile_all,wav)
startt=' 20:46:00'
endt=' 20:50:00'
idxs,idxe=aia.get_time_idx(startt,endt,aiatsec)
ccmap=cmap[idxs:idxe]
aiactime=aiatstring[idxs:idxe]

aiamax=[0]*len(ccmap)
aiamax_xid=[0]*len(ccmap)
aiamax_yid=[0]*len(ccmap)
for j in range(len(ccmap)):
    aiamax[j]=np.max(ccmap[j])
    aiamax_xid[j]=aiaxarray[np.where(ccmap[j]==np.max(ccmap[j]))[0][0]]
    aiamax_yid[j]=aiayarray[np.where(ccmap[j]==np.max(ccmap[j]))[1][0]]

aiadata304=umaps.get_submap('/media/rohit/VLA/20120225_aia/AIA20120225_205132_0304.fits',xl,xr,yl,yr)[0]
aiadata171=umaps.get_submap('/media/rohit/VLA/20120225_aia/AIA20120225_204600_0171.fits',xl,xr,yl,yr)[0]
aiadata94=umaps.get_submap('/media/rohit/VLA/20120225_aia/AIA20120225_204902_0094.fits',xl,xr,yl,yr)[0]

#######################################


expolB_=readsav('/data/Dropbox/20120225_VLA_work/gregory_gyrosynchroton/hmi.M_720s.20120225_203413.W133N16CR.CEA.NAS.sav')
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
out_header = sunpy.map.make_fitswcs_header(out_shape,hp_coord)
out_wcs = WCS(out_header)
earth = get_body_heliographic_stonyhurst('earth', mymap.date)
out_wcs.heliographic_observer = earth
output, footprint = reproject_interp(mymap, out_wcs, out_shape)
outmap = sunpy.map.Map((output, out_header))
outmap.plot_settings = mymap.plot_settings

fig = plt.figure()

ax1 = fig.add_subplot(1, 2, 1, projection=mymap)
mymap.plot(axes=ax1)
outmap.draw_grid(color='w')

ax2 = fig.add_subplot(1, 2, 2, projection=outmap)
outmap.plot(axes=ax2)
outmap.draw_grid(color='w')
plt.close()
bradial=babs[:,61:67,106:116].max(axis=(1,2))
lev_1=np.linspace(0.18,0.180001,2) #for 304
lev_2=np.linspace(0.2,0.2001,3) #for 171
xl1,xr1,yl1,yr1=out_header['crval1']-1.0*100.5,out_header['crval1']+1.0*100.5,out_header['crval2']-1.0*70.5,out_header['crval2']+1.0*70.5
pl.babs_paper(babs[0],babs[6],babs[10],bradial,aiadata304.data*1.0,aiadata171.data*1.0, lev_1,lev_2,xl,xr,yl,yr,xl1,xr1,yl1,yr1)

