import sys
import oskar
import matplotlib
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
import numpy as np
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
import os

plt.style.use(astropy_mpl_style)

start_frequency_hz = 2.e8
phase_centre_ra_deg = 300.0
phase_centre_dec_deg = -40.0
output_root = 'meerkat_visibilities'
os.system("rm -rf "+output_root+'.vis')
#telescope_conf='/home/rohit/simulations/meerKat/meerkat.tm'
telescope_conf='/home/rohit/simulations/alma/telescope.tm'
gleam_sky=0
if(gleam_sky):
    f=fits.open('../gleam/GLEAM_EGC_v2.fits')
    data=f[1].data
    ra=data['RAJ2000']
    dec=data['DEJ2000']
    #fpw = data['Fpwide']
    sky_array = np.column_stack((ra, dec,  np.zeros(ra.shape[0]), np.zeros(ra.shape[0]),np.zeros(ra.shape[0]), [start_frequency_hz]*ra.shape[0]))
    precision = "single"
    sky = oskar.Sky.from_array(sky_array, precision) # 307455 sources
    sky.filter_by_radius(.4, .55, phase_centre_ra_deg, phase_centre_dec_deg)
    data = sky.to_array()
    ra_filtered = data[:,0]
    dec_filtered = data[:,1]
    stokes_i_flux = data[:,2]


# Set up the sky model.
sky_grid=1
if(sky_grid):
    precision = "double"
    sky = oskar.Sky.generate_grid(phase_centre_ra_deg, phase_centre_dec_deg,5, 10, precision=precision)
    sky.append_sources(phase_centre_ra_deg, phase_centre_dec_deg, 1.0)
    skyarr=sky.to_array()


params = {
    "simulator": {
        "use_gpus": False
    },
    "observation" : {
        "num_channels": 32,
        "start_frequency_hz": start_frequency_hz,
        "frequency_inc_hz": 1.0e6,
        "phase_centre_ra_deg": phase_centre_ra_deg,
        "phase_centre_dec_deg": phase_centre_dec_deg,
        "num_time_steps": 12, # No 
        "start_time_utc": "01-01-2000 6:00:00.000",
        "length": "2:00:00.000" # No
    },
    "telescope": {
        "input_directory": telescope_conf
    },
    "interferometer": {
       "ms_filename": output_root+".ms",
       "oskar_vis_filename": output_root+".vis",
        "channel_bandwidth_hz": 1e6,
        "time_average_sec": 20 # No
    }
}
settings = oskar.SettingsTree("oskar_sim_interferometer")
settings.from_dict(params)

if precision == "single":
    settings["simulator/double_precision"] = False

# Set the sky model and run the simulation.
sim = oskar.Interferometer(settings=settings)
sky_data = np.array([[20.0, -30.0, 1, 0, 0, 0, 100.0e6, -0.7, 0.0, 0,   0,   0],
    [20.0, -30.5, 3, 2, 2, 0, 100.0e6, -0.7, 0.0, 600, 50,  45],
    [20.5, -30.5, 3, 0, 0, 2, 100.0e6, -0.7, 0.0, 700, 10, -10]])
sky = oskar.Sky.from_array(sky_data, precision)
sim.set_sky_model(sky)
sim.run()

(header, handle) = oskar.VisHeader.read(output_root+'.vis')
#(header, handle) = oskar.VisHeader.read('example.vis')
block = oskar.VisBlock.create_from_header(header)
tasks_read = []
import concurrent
executor = concurrent.futures.ThreadPoolExecutor(2)
for i_block in range(header.num_blocks):
    tasks_read.append(executor.submit(block.read, header, handle,i_block))
vis = block.cross_correlations()
print(vis.shape)

