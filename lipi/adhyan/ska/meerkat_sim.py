from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
from oskar.sky import Sky
from oskar.measurement_set import MeasurementSet as ms
from oskar.telescope import Telescope as tel
import oskar
import time
import sys
import os
import concurrent

precision='single'
ra_deg=45.0;dec_deg=45.0;I=1
sky = oskar.Sky.generate_grid(ra_deg, dec_deg,25, 3.0, precision=precision)
sky_data=np.array([[0.0, 0.0, 1]])
#ss=oskar.Sky()
#sky.load('sky.osm')
tel = oskar.Telescope(precision)
tel.load('/home/rohit/simulations/meerKat/meerkat.tm')
#tel.load('/home/rohit/simulations/alma/telescope.tm')
tel.set_channel_bandwidth(1.0e3)
tel.set_time_average(10.0)
tel.set_pol_mode('Scalar')
# Set station properties after stations have been defined.
tel.set_phase_centre(ra_deg, dec_deg)
tel.set_station_type('Gaussian beam')
tel.set_gaussian_station_beam_width(5.0, 100e6)

simulator = oskar.Interferometer(precision)
#simulator.set_settings_path(os.path.abspath('.'))
simulator.set_max_sources_per_chunk(sky.num_sources+1)
simulator.set_sky_model(sky)
simulator.set_telescope_model(tel)
simulator.set_observation_frequency(start_frequency_hz=950.0e6,inc_hz=5.e6,num_channels=32)
#simulator.set_observation_time(start_time_mjd_utc=51545.0, length_sec=8.0, num_time_steps=80000)
simulator.set_observation_time(start_time_mjd_utc=0, length_sec=8.0, num_time_steps=80)
simulator.set_output_measurement_set('run1.ms')
simulator.set_output_vis_file('run1.vis')
simulator.set_gpus(None)
print('Running interferometer simulator...')
simulator.run()


sys.exit()

(header, handle) = oskar.VisHeader.read('run1.vis')
#(header, handle) = oskar.VisHeader.read('example.vis')
block = oskar.VisBlock.create_from_header(header)
tasks_read = []
executor = concurrent.futures.ThreadPoolExecutor(2)
for i_block in range(header.num_blocks):
        tasks_read.append(executor.submit(block.read, header, handle,i_block))
vis = block.cross_correlations()
um=block.baseline_uu_metres();vm=block.baseline_vv_metres()
print(vis.shape)


settings = oskar.SettingsTree('oskar_sim_interferometer')
python_dict = {
   'simulator': {
       'double_precision': 'true',
       'use_gpus': 'true',
       'max_sources_per_chunk': 23000
   },
   'observation' : {
       'length': 14400.0,
       'start_frequency_hz': 132e6,
       'frequency_inc_hz': 100e3,
       'num_channels': 160,
       'num_time_steps': 240
   },
   'telescope': {
       'input_directory': '/path/to/telescope/model',
       'pol_mode': 'Scalar'
   },
   'interferometer': {
       'channel_bandwidth_hz': 100e3,
       'time_average_sec': 1.0,
       'max_time_samples_per_block': 4
   }
}
settings.from_dict(python_dict)

sys.exit()

filename = "test_zenith.ms"
# Define data dimensions.
num_pols = 4
num_channels = 2
num_stations = 64
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

