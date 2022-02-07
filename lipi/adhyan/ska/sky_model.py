from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
from oskar.sky import Sky
from oskar.measurement_set import MeasurementSet as ms
from oskar.telescope import Telescope as tel
import oskar

ra_deg=0;dec_deg=0;I=1
sky_data=np.array([[20.0, -30.0, 1],[20.0, -30.5, 3]])
s1=Sky.from_array(sky_data)
ss=oskar.Sky();ss.load('sky.osm')
tel = oskar.Telescope()
tel.load('/home/rohit/simulations/alma/alma1.tm')

filename = "test_zenith.ms"
# Define data dimensions.
num_pols = 4
num_channels = 2
num_stations = 3
num_times = 400
num_baselines = num_stations * (num_stations - 1) // 2
ref_freq_hz = 230e9 # Observational Frequency 
freq_inc_hz = 7.8e6 #
exposure_sec = 10.0
interval_sec = 10.0

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
ms_.close()





##########################

aaM=fits.open('/home/i4ds1807205/skymodels/lambda_mollweide_haslam408_nofilt.fits')
hM,dM=aaM[1].header,aaM[1].data;wcsM=WCS(hM)
aaH=fits.open('/home/i4ds1807205/skymodels/lambda_haslam408_nofilt.fits')
hH,dH=aaH[1].header,aaH[1].data;wcsH=WCS(hH)
aaZ=fits.open('/home/i4ds1807205/skymodels/lambda_zea_haslam408_nofilt.fits')
hZ,dZ=aaZ[1].header,aaZ[1].data;wcsZ=WCS(hZ)


plt.figure()
ax0=plt.subplot(211,projection=wcsZ)
ax1=plt.subplot(212,projection=wcsM)
im0=ax0.imshow(np.log10(dZ),origin='lower');ax0.set_title('ZEA Projection')
im1=ax1.imshow(np.log10(dM),origin='lower');ax1.set_title('Mollweide Projection')
plt.colorbar(im0,ax=ax0,label='Log10(Temperature (MK))')
plt.colorbar(im1,ax=ax1,label='Log10(Temperature (MK))')
plt.show()

