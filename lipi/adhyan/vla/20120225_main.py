from surya.utils import main as ut
from surya.aia import main as aia
import numpy as np
import glob
import 20120225_plots as pl

aiapath='/home/i4ds1807205/vla_data/20120225/aia/'
aiafile=aiapath+'aia_171_sub_map.sav'
aiafile_all=sorted(glob.glob(aiapath+'aia*.sav'))
aiafile_all.sort(key=lambda ss: len(ss))
wavelength=171
wav=[94,131,171,193,211,304,335,1600,1700]
xl,xr,yl,yr=450,510,320,380
res=0.6
cmap=aia.get_submap_crop(aiafile,res,wavelength,xl,xr,yl,yr)
cmap_all=aia.get_submap_all(aiafile_all,wav)

