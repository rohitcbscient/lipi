from surya.utils import main as ut
from surya.aia import main as aia
from surya.vla import main as vla
import numpy as np
import glob
import matplotlib.pyplot as plt
import plot_20120225 as pl

aiapath='/home/i4ds1807205/vla_data/20120225/aia/'
aiafile=aiapath+'aia_171_sub_map.sav'
aiafile_all=sorted(glob.glob(aiapath+'aia*.sav'))
aiafile_all.sort(key=lambda ss: len(ss))
wavelength=171
wav=[94,131,171,193,211,304,335,1600,1700]
#xl,xr,yl,yr=450,510,320,380
xl,xr,yl,yr=420,540,280,420 # For Quiet Sun plotting
res=0.6
#cmap=aia.get_submap(aiafile,res,wavelength)
cmap=aia.get_submap_crop(aiafile,res,wavelength,xl,xr,yl,yr)
cmap_all=aia.get_submap_all(aiafile_all,wav)

t='20:47:15~20:47:16'
vlapath='/media/rohit/VLA/20120225_sub/fits/'
vlafile_fl=vlapath+'hcs_2050.1s.cal.spw.7.time.'+t+'.FITS'
vlafile_qs=vlapath+'hcs_2050.1s.cal.spw.7.time.20:46:02~20:46:03.FITS'
vlafile_qs='/media/rohit/VLA/selfcal/selfcal_pf_qs/hcs_2050.1s.cal.spw.7.time.20:46:00~20:46:01.FITS'
vlafile_qs0='/media/rohit/VLA/yingie_data/20120225_images/time_studies/spw_7/flare_time/hcs_2050.1s.cal.spw.7:0~3.time.20:46:01~20:46:02.FITS'
vlafile_all=sorted(glob.glob(vlapath+'hcs*.FITS'))
vlafile_all.sort(key=lambda ss: len(ss))
pix=2.5
dec=-9.13056
vlamap,vlah=vla.get_map(vlafile_qs0)
#vlamap,vlah=vla.get_map_cube(vlafile_qs)
vlasubmap,vlah=vla.get_submap(vlafile_qs0,xl,xr,yl,yr,pix,dec)
#vlasubmap,vlah=vla.get_submap_cube(vlafile_fl,xl,xr,yl,yr,pix,dec)

## Plotting
pl.euv_vla(cmap,vlasubmap,xl,xr,yl,yr)



