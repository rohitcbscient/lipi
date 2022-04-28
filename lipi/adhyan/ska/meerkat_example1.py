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
phase_centre_ra_deg = 300.0;phase_centre_dec_deg = -40.0
#phase_centre_ra_deg = 45.0;phase_centre_dec_deg = 45.0
output_root = 'meerkat_core_visibilities_test2'
os.system("rm -rf "+output_root+'.vis')
telescope_conf='/home/rohit/simulations/meerKat/meerkat_core.tm'
telconfig=np.loadtxt(telescope_conf+'/layout.txt')
position=np.loadtxt(telescope_conf+'/position.txt')
#telescope_conf='/home/rohit/simulations/alma/telescope.tm'
#telescope_conf='/home/rohit/simulations/meerKat/vla.c.tm'
#telescope_conf='/home/rohit/simulations/configuration/oskar_configurations/vla.c.tm'
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

make_healpix=1
if(make_healpix):
    nside=12;npix=128
    radius = 1.0 # in degrees
    theta = 0.5 * np.pi - np.deg2rad(phase_centre_dec_deg)
    phi = np.deg2rad(phase_centre_ra_deg)
    radius = np.deg2rad(radius)
    xyz = hp.ang2vec(theta, phi) # Cartisian Coordinates
    ipix_disc = hp.query_disc(nside, xyz, radius)
    healpixsky=np.zeros(nside*npix*npix)
    healpixsky[ipix_disc]=1.0


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
        "start_time_utc": "01-01-2022 13:21:00.000",
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

#tel.set_phase_centre(phase_centre_ra_deg, phase_centre_dec_deg)
#tel.set_station_type('Gaussian beam')
#tel.set_gaussian_station_beam_width(5.0, 100e6)
# python_dict = settings.to_dict(include_defaults=True)

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
print(vis)

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
ax00.plot(uum,vvm,'o');ax00.set_xlabel('U (m)');ax00.set_ylabel('V (m)')
ax10.plot(uvdistlam/1.e3,amp,'o');ax10.set_ylabel('Amplitude (Jy)');ax10.set_xlabel('UV distance ($k\lambda$)')
ax11.plot(uvdistlam/1.e3,phase,'o');ax11.set_ylabel('Phase');ax11.set_xlabel('UV distance ($k\lambda$)')
ax01.plot(telconfig[:,0],telconfig[:,1],'o');ax01.set_xlabel('X (m)');ax01.set_ylabel('Y (m)')
plt.show()
####

f1=fits.open('/home/rohit/simulations/meerKat/meerkat_visibilities_test1_I.fits');fd=f[0].data
fd=f1[0].data

f,ax=plt.subplots(2,1);ax0=ax[0];ax1=ax[1]
im0=ax0.imshow(fd[0],aspect='auto',vmin=0.01,vmax=1,cmap='jet')
ax1.plot(fd[0][256])
plt.colorbar(im0)
plt.title('MeerKat Test Benchmark Image')
plt.show()


