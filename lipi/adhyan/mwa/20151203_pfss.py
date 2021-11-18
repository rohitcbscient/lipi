from datetime import datetime
import os

import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.io.fits
import pickle
import pfsspy
import pfsspy.tracing as tracing


gongfile='/media/rohit/VLA/pfss/mrzqs151203t0304c2171_235.fits.gz'
gongfile='/media/rohit/VLA/pfss/mrzqs200901t1304c2234_022.fits'
gong_map = sunpy.map.Map(gongfile)
gong_map = sunpy.map.Map(gong_map.data - np.mean(gong_map.data), gong_map.meta)
nrho = 35
rss = 2
input = pfsspy.Input(gong_map, nrho, rss)
#output = pfsspy.pfss(input)

flines=pickle.load(open('/media/rohit/VLA/pfss/flines.p','rb'))
aiafile='/media/rohit/MWA/20151203_EUV/aia.lev1.193A_2015-12-03T03_17_05.84Z.image_lev1.fits'
aiamap=sunpy.map.Map(aiafile)


fig = plt.figure()
ax = plt.subplot(1, 1, 1, projection=aiamap)
aiamap.plot(ax)
for fline in flines.closed_field_lines:
    ax.plot_coord(fline.coords, alpha=0.8, linewidth=1, color='yellow')
    # ax.set_xlim(500, 900)
    # ax.set_ylim(400, 800)
for fline in flines.open_field_lines:
    ax.plot_coord(fline.coords, alpha=0.8, linewidth=1, color='red')
plt.show()


sys.exit()


