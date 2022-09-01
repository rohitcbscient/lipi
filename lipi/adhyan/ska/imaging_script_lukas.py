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

plt.style.use(astropy_mpl_style)

f=fits.open('../gleam/GLEAM_EGC_v2.fits')
data=f[1].data
ra=data['RAJ2000']
dec=data['DEJ2000']
#fpw = data['Fpwide']

start_frequency_hz = 76e6
phase_centre_ra_deg = 250
phase_centre_dec_deg = -80

precision = "single"
sky_array = np.column_stack((ra, dec,  np.zeros(ra.shape[0]), np.zeros(ra.shape[0]), 
                             np.zeros(ra.shape[0]), [start_frequency_hz]*ra.shape[0]))
sky = oskar.Sky.from_array(sky_array, precision) # 307455 sources
sky.filter_by_radius(.5, .59, phase_centre_ra_deg, phase_centre_dec_deg)
data = sky.to_array()
ra_filtered = data[:,0]
dec_filtered = data[:,1]
stokes_i_flux = data[:,2]

params = {
    "simulator": {
        "use_gpus": False
    },
    "observation" : {
        #"num_channels": 64,
        "num_channels": 16,
        "start_frequency_hz": start_frequency_hz,
        "frequency_inc_hz": 20e6,
        "phase_centre_ra_deg": phase_centre_ra_deg,
        "phase_centre_dec_deg": phase_centre_dec_deg,
        "num_time_steps": 24,
        "start_time_utc": "01-01-2000 12:00:00.000",
        "length": "12:00:00.000"
    },
    "telescope": {
        "input_directory": "../alma/telescope.tm"
        #"input_directory": "meerkat.tm"
    },
    "interferometer": {
        "ms_filename": "visibilities_gleam_meerkat.ms",
        "oskar_vis_filename": "visibilities_gleam.vis",
        "channel_bandwidth_hz": 1e6,
        "time_average_sec": 10,
        "noise/enable":True,
        "noise/seed":'time',
        "noise/freq/start":1.e9,
        "noise/freq/inc":1.e8,
        "noise/freq/number":10,
        "noise/rms":"Range",
        "noise/rms/start":500,
        "noise/rms/end":1000
    }
}
settings = oskar.SettingsTree("oskar_sim_interferometer")
settings.from_dict(params)

if precision == "single":
    settings["simulator/double_precision"] = False

# Set the sky model and run the simulation.
sim = oskar.Interferometer(settings=settings)

sim.set_sky_model(sky)
sim.run()

sys.exit()
(header, handle) = oskar.VisHeader.read('visibilities_gleam.vis')
block = oskar.VisBlock.create_from_header(header)
tasks_read = []
import concurrent
executor = concurrent.futures.ThreadPoolExecutor(2)
for i_block in range(header.num_blocks):
    tasks_read.append(executor.submit(block.read, header, handle,i_block))
vis = block.cross_correlations()


from rascil.apps import rascil_imager
from rascil.processing_components.util.performance import (
    performance_store_dict,
    performance_environment,
)
    
def start_imager(rawargs):
    parser = rascil_imager.cli_parser()
    args = parser.parse_args(rawargs)
    performance_environment(args.performance_file, mode="w")
    performance_store_dict(args.performance_file, "cli_args", vars(args), mode="a")
    image_name = rascil_imager.imager(args)

start_imager(
    [
        '--ingest_msname','visibilities_gleam.ms',
        '--ingest_dd', '0', 
        #'--ingest_vis_nchan', '64',
        '--ingest_vis_nchan', '16',
        '--ingest_chan_per_blockvis', '4' ,
        '--ingest_average_blockvis', 'True',
        '--imaging_npixel', '2048', 
        '--imaging_cellsize', '3.878509448876288e-05',
        '--imaging_weighting', 'robust',
        '--imaging_robustness', '-0.5',
        '--clean_nmajor', '1' ,
        '--clean_algorithm', 'mmclean',
        '--clean_scales', '0', '6', '10', '30', '60',
        '--clean_fractional_threshold', '0.3',
        '--clean_threshold', '0.12e-3',
        '--clean_nmoment' ,'5',
        '--clean_psf_support', '640',
        '--clean_restored_output', 'integrated'
    ])
