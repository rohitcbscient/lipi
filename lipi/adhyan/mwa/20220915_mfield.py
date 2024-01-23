import os

import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from sunpy.map import Map
from sunpy.net import Fido
from sunpy.net import attrs as a
from astropy.utils.data import download_file
from sunpy.map.header_helper import make_fitswcs_header
from sunpy.coordinates import frames
import pfsspy
from astropy.coordinates import SkyCoord

filename = download_file(
    "http://jsoc.stanford.edu/data/hmi/synoptic/hmi.Synoptic_Mr.2191.fits", cache=True
)
syn_map = Map(filename)
syn_map_data = syn_map.data
syn_map_data[np.isnan(syn_map_data)] = 0
syn_map1 = Map(syn_map_data, syn_map.meta)


time = a.Time("2022/09/12", "2022/09/16")
series = a.jsoc.Series("hmi.synoptic_mr_polfil_720s")
crot = a.jsoc.PrimeKey("CAR_ROT", 2210)


# result = Fido.search(time, series, crot, a.jsoc.Notify("rohitcbscient@gmail.com"))
# files = Fido.fetch(result)

hmi_map = Map("/sdata/20220915_hmi/hmi.m_45s.2022.09.12_03_01_30_TAI.magnetogram.fits")
new_dimensions = [128, 128] * u.pixel
hmi_map_resampled = hmi_map.resample(new_dimensions)

shape = (720 * 20, 1440 * 20)
carr_data = np.ones(shape)
# carr_coord=SkyCoord(0*u.arcsec, 0*u.arcsec, observer = 'earth', obstime=hmi_map.meta['t_obs'], frame=frames.HeliographicCarrington)
carr_coord = hmi_map_resampled.center.transform_to(frames.HeliographicCarrington)
carr_header = make_fitswcs_header(carr_data, carr_coord)
# carr_map=Map(carr_data,carr_header)
carr_map = hmi_map_resampled.reproject_to(carr_header)
carr_map.plot()
plt.show()


f, ax = plt.subplots(1, 2)
ax0 = ax[0]
ax1 = ax[1]
ax0.imshow(
    hmi_map.data,
    interpolation=None,
    vmin=-40,
    vmax=40,
    aspect="auto",
    origin="lower",
    cmap="jet",
)
ax1.imshow(
    carr_map.data,
    interpolation=None,
    vmin=-40,
    vmax=40,
    aspect="auto",
    origin="lower",
    cmap="jet",
)
plt.show()

plt.imshow(
    carr_map.data,
    interpolation=None,
    vmin=-40,
    vmax=40,
    aspect="auto",
    origin="lower",
    cmap="jet",
)
plt.colorbar()
plt.show()

nrho = 35
rss = 2.5
pfss_in = pfsspy.Input(carr_map, nrho, rss)
pfss_out = pfsspy.pfss(pfss_in)

ss_br = pfss_out.source_surface_br
# Create the figure and axes
fig = plt.figure()
ax = plt.subplot(projection=ss_br)

# Plot the source surface map
ss_br.plot()
# Plot the polarity inversion line
ax.plot_coord(pfss_out.source_surface_pils[0])
# Plot formatting
plt.colorbar()
ax.set_title("Source surface magnetic field")

plt.show()
