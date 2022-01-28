from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
from oskar.sky import Sky
from oskar.measurement_set import MeasurementSet as ms
import oskar
import pandas


tel = oskar.Telescope()
tel.load('/home/rohit/ska_pipeline/SKA-main/data/telescope.tm')
alma_conf=np.loadtxt('/home/rohit/simulations/alma/alma.cycle5.3.edited.cfg')
# UTM Zone for ALMA is 19
p = Proj(proj='utm',zone=19,ellps='WGS84', preserve_units=False);p(-120.108, 34.3)


ra_deg=0;dec_deg=0;I=1
sky_data=np.array([[20.0, -30.0, 1],[20.0, -30.5, 3]])
sources_x=np.random.normal(-40,40,1000);sources_y=np.random.normal(-40,40,1000)
s1=Sky.from_array(sky_data)

filename = "test_zenith.ms"
# Define data dimensions.
num_pols = 4
num_channels = 2
num_stations = 3
num_times = 4
num_baselines = num_stations * (num_stations - 1) // 2
ref_freq_hz = 100e6
freq_inc_hz = 100e3
exposure_sec = 1.0
interval_sec = 1.0

# Data to write are stored as numpy arrays.
uu = np.zeros([num_baselines])
vv = np.zeros_like(uu)
ww = np.zeros_like(uu)
vis = np.zeros([num_times, num_channels,num_baselines, num_pols], dtype='c8')

# Create the empty Measurement Set.
ms_ = ms.create(filename, num_stations,num_channels, num_pols,ref_freq_hz, freq_inc_hz)

# Set phase centre.
ra_rad = np.pi / 4
dec_rad = -np.pi / 4
ms_.set_phase_centre(ra_rad, dec_rad)

# Write data one block at a time.
for t in range(num_times):
    # Dummy data to write.
    time_stamp = 51544.5 * 86400.0 + t
    uu[:] = 1e0 * t + 1
    vv[:] = 1e1 * t + 2
    ww[:] = 1e2 * t + 3
    for c in range(num_channels):
        for b in range(num_baselines):
            vis[t, c, b, :] = (t * 10 + b) + 1j * (c + 1)

# Write coordinates and visibilities.
start_row = t * num_baselines
time_stamp=10
ms_.write_coords(start_row, num_baselines, uu, vv, ww,exposure_sec, interval_sec, time_stamp)




##########################


