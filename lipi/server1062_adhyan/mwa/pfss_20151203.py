import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.io.fits

import pfsspy
import pfsspy.tracing as tracing


gongfile='/nas08-data02/rohit/pfss_20151203/mrzqs151203t0304c2171_235.fits.gz'
#gongfile='/nas08-data02/rohit/pfss_20151203/mrzqs200901t1304c2234_022.fits'
gong_map = sunpy.map.Map(gongfile)
gong_map = sunpy.map.Map(gong_map.data - np.mean(gong_map.data), gong_map.meta)
nrho = 35
rss = 2
input = pfsspy.Input(gong_map, nrho, rss)
output = pfsspy.pfss(input)
# Create 5 points spaced between sin(lat)={0.35, 0.55}
s = np.linspace(-1, 1, 30)
# Create 5 points spaced between long={60, 100} degrees
phi = np.linspace(0, 350, 30)
print(f's = {s}')
print(f'phi = {phi}')
# Make a 2D grid from these 1D points
s, phi = np.meshgrid(s, phi)

plot_cmap=0
if(plot_cmap):
    m = input.map
    fig = plt.figure()
    ax = plt.subplot(projection=m)
    m.plot()
    plt.colorbar()
    ax.set_title('Input field')
    plt.show()



# Now convert the points to a coordinate object
lat = np.arcsin(s) * u.rad
lon = phi * u.deg
seeds = SkyCoord(lon.ravel(), lat.ravel(), 1.01 * const.R_sun,
                         frame=gong_map.coordinate_frame)

tracer = tracing.PythonTracer()
flines = tracer.trace(seeds, output)

aiafile='/nas08-data02/rohit/pfss_20151203/aia.lev1.193A_2015-12-03T03_17_05.84Z.image_lev1.fits'
aiamap=sunpy.map.Map(aiafile)
fig = plt.figure()
ax = plt.subplot(1, 1, 1, projection=aiamap)
aiamap.plot(ax)
for fline in flines1.closed_field_lines:
    ax.plot_coord(fline.coords, alpha=0.8, linewidth=1, color='yellow')
    # ax.set_xlim(500, 900)
    # ax.set_ylim(400, 800)
plt.show()



sys.exit()

