# MeerKat Benchmark Run - 28th April 2022 with Sara, ETHz
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
import healpy as hp
import h5py

plt.style.use(astropy_mpl_style)

start_frequency_hz = 8.e8
phase_centre_ra_deg = 0.0;phase_centre_dec_deg = -30.721
output_root = 'meerkat_core_visibilities_test2'
os.system("rm -rf "+output_root+'.vis')
telescope_conf='/home/rohit/simulations/meerKat/meerkat_core.tm'
telconfig=np.loadtxt(telescope_conf+'/layout.txt')
position=np.loadtxt(telescope_conf+'/position.txt')

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
        "num_time_steps": 1, # No 
        "start_time_utc": "21-03-2022 10:39:00.000",
        "length": "0:00:08.000", # No
    },
    "telescope": {
        "input_directory": telescope_conf,
        "gaussian_beam/fwhm_deg": 10.0,
        "gaussian_beam/ref_freq_hz":8.e8
    },
    "interferometer": {
       "ms_filename": output_root+".ms",
       "oskar_vis_filename": output_root+".vis",
        "channel_bandwidth_hz": 1e6,
        "time_average_sec": 1 # No
    }
}

settings = oskar.SettingsTree("oskar_sim_interferometer")
settings.from_dict(params)

if precision == "single":
    settings["simulator/double_precision"] = False

# Set the sky model and run the simulation.
sim = oskar.Interferometer(settings=settings)
sky_data = np.array(([phase_centre_ra_deg, phase_centre_dec_deg, 1, 0, 0, 0, start_frequency_hz, -0.7, 0.0, 0.1,   0.1,   0]))
# RA,DEC,I,Q,U,V,f,alpha,RM, major FWHM, minor FWHM, PA
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

lam=3.e8/start_frequency_hz
uum=block.baseline_uu_metres()[0];vvm=block.baseline_vv_metres()[0];ulam=uum/lam;vlam=vvm/lam
uvdistm=np.sqrt(uum**2 +vvm**2);uvdistlam=uvdistm/lam

amp=np.sqrt(np.real(vis[0,0,:,0])**2 + np.imag(vis[0,0,:,0])**2)
amp=np.sqrt(np.real(vis[0,0,:,0])**2 + np.imag(vis[0,0,:,0])**2)
phase=np.angle(vis[0,0,:,0], deg=True)

####
hirax_file='/data/Dropbox/RadioSims/Benchmarking_karabo&hirax/files/vis_test_c0.h5'
f = h5py.File(hirax_file)
vish=np.array(f.get('vis'))
amph=np.sqrt(vish[0].flatten().real**2 + vish[0].flatten().imag**2);phaseh=np.angle(vish[0], deg=True)

f,ax=plt.subplots(2,2);ax00=ax[0,0];ax01=ax[0,1];ax10=ax[1,0];ax11=ax[1,1]
ax00.plot(uum[uvdistm<70],vvm[uvdistm<70],'o');ax00.set_xlabel('U (m)');ax00.set_ylabel('V (m)')
#ax10.plot(uvdistlam/1.e3,amp,'o');ax10.set_ylabel('Amplitude (Jy)');ax10.set_xlabel('UV distance ($k\lambda$)')
#ax11.plot(uvdistlam/1.e3,phase,'o');ax11.set_ylabel('Phase');ax11.set_xlabel('UV distance ($k\lambda$)')
ax10.plot(uvdistlam[uvdistm<70]/1.e3,amp[uvdistm<70],'o');ax10.set_ylabel('Amplitude (Jy)');ax10.set_xlabel('UV distance ($k\lambda$)')
ax11.plot(uvdistlam[uvdistm<70]/1.e3,phase[uvdistm<70],'o');ax11.set_ylabel('Phase');ax11.set_xlabel('UV distance ($k\lambda$)')
ax01.plot(telconfig[:,0],telconfig[:,1],'o');ax01.set_xlabel('X (m)');ax01.set_ylabel('Y (m)')
plt.show()


