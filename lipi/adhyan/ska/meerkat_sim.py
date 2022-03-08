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
from cora.foreground import galaxy, pointsource

nside = 64
nfreq=32
fa = np.linspace(950.0, 1110.0, nfreq)*1.e6
cr = pointsource.CombinedPointSources()
cr.nside=nside
cr.frequencies = fa
unpol = cr.getsky()
pol = cr.getpolsky()
unpol=unpol.reshape(nfreq,12,nside,nside)
plt.imshow(np.log10(unpol.max(axis=(0,1))),aspect='auto')
plt.colorbar(label='Log(T (K))')
plt.show()

hZ,dZ=aaZ[1].header,aaZ[1].data;wcsZ=WCS(hZ)
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

sys.exit()

ms.open('run1.ms')
amp=ms.getdata('amplitude');u=ms.getdata('u')['u'];u=ms.getdata('v')['v']
f,ax=plt.subplots(1,2);ax0=ax[0];ax1=ax[1]
ax0.plot(u,v,'o')
ax1.plot(u,v,'o')
plt.show()

plt.figure()
ax0=plt.subplot(211,projection=wcsZ)
ax1=plt.subplot(212,projection=wcsM)
im0=ax0.imshow(np.log10(dZ),origin='lower');ax0.set_title('ZEA Projection')
im1=ax1.imshow(np.log10(dM),origin='lower');ax1.set_title('Mollweide Projection')
plt.colorbar(im0,ax=ax0,label='Log10(Temperature (MK))')
plt.colorbar(im1,ax=ax1,label='Log10(Temperature (MK))')
plt.show()

