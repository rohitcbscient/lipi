from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation, ObservationLong
from karabo.simulation.visibility import Visibility

from oskar.imager import Imager
import numpy as np, os, sys, pandas as pd
import tools21cm as t2c
import matplotlib.pyplot as plt

from datetime import datetime, timedelta

import astropy.units as u
import astropy.constants as cst
from astropy.io import fits
from astropy.wcs import WCS
from astropy.cosmology import Planck18
from astropy.coordinates import Angle

from smoothing import smoothing_grid

idx_f = 750

# move to scratch
os.chdir('/scratch/snx3000/rsharma/')

# --- Sky model ---
path_in = '/scratch/snx3000/mibianco/output_sdc3/dataLC_256_train_090523/data/'
path_out = '/scratch/snx3000/rsharma/subvis_tests/'
path_telescope = '/store/ska/sk014/dataset_sdc3/inputs/telescope.tm'

os.chdir(path_out)

root_name = 'coevalLC_256_train_190922_i0_dT'
rseed = 2023

with fits.open(path_in+root_name+'.fits', mode="readonly", memmap=True) as hdulist:
    hdr = hdulist[0].header
    Nx, Ny = hdr['NAXIS1'], hdr['NAXIS2']
    FoV = (abs(hdr['CDELT2']*Nx) * u.deg).to('rad')
    beam_sim = (FoV.value / Nx)**2     # from grid to point source
    RA, DEC = hdr['CRVAL1']*u.deg, hdr['CRVAL2']*u.deg
    
    freqs = hdr['CRVAL3']+np.arange(hdr['NAXIS3']) * hdr['CDELT3'] # Hz
    z = t2c.nu_to_z(freqs[idx_f]*1e-6)

    # add galactic foreground
    root_name += 'gf'
    zmax = t2c.nu_to_z(106.900)
    fov_mpc = (Planck18.comoving_transverse_distance(zmax).value * FoV.value)
    #data_gf = t2c.galactic_synch_fg(z=[z], ncells=256, boxsize=fov_mpc, rseed=2023)*1e-3 #* u.K
    #data = (hdulist[0].data[idx_f] + data_gf) * beam_sim
    data = hdulist[0].data[idx_f] * beam_sim
    #data = hdulist[0].data[idx_f] * beam_sim
    w = WCS(hdr).celestial

# get primary beam and 
path_beam = '/store/ska/sk01/sdc3-old-data-please-dont-use-it/Image/station_beam.fits'
beam_data = fits.getdata(path_beam)

beam_small = smoothing_grid(arr=beam_data[idx_f], noc=Nx, opt='mean')
print(beam_small.shape)

# get coordinates
idx_ra, idx_dec = np.arange(0, Nx).reshape(Nx, 1), np.arange(0, Ny).reshape(1, Ny)
lon, lat = w.celestial.all_pix2world(idx_ra, idx_dec, 1)
sky_grid =  np.vstack((lon[np.newaxis, ...], lat[np.newaxis, ...])).reshape(2,lon.shape[0]*lon.shape[1]).T

freq = freqs[idx_f]

print('FoV = %.2f %s' %(FoV.to('deg').value, FoV.to('deg').unit))
print('RA, DEC = (%f, %f) %s' %(RA.value, DEC.value, DEC.unit))
print('nu = %.5e Hz (idx_f = %d)' %(freq, idx_f))

# convert K to Jy
Jy2kel = (u.Jy * cst.c * cst.c / (2 * cst.k_B * (freq*u.Hz)**2)).cgs
data *= (u.K / Jy2kel).cgs.value

# use only non-zero grid points
idx_nozero = np.nonzero(data.flatten())[0]

# create sky model with columns:
# RA [deg], Dec [deg], I [Jy], Q [Jy], U [Jy], V [Jy], ref freq [Hz], alpha, rot, maj ax [arcsec], min ax [arcsec], pos angle [deg], object ID
sky_data = np.zeros((idx_nozero.size, 13))
sky_data[:,:3] = np.hstack((sky_grid[idx_nozero, :], data.flatten()[idx_nozero, np.newaxis]))
sky_data[:,9] += 2*FoV.to('arcsec').value/Nx
sky_data[:,10] += 2*FoV.to('arcsec').value/Nx
#sky = SkyModel(sky_data)

# Apply primary beam
path_beam = '/store/ska/sk01/sdc3-old-data-please-dont-use-it/Image/station_beam.fits'
beam_data = fits.getdata(path_beam)
#data *= smooth_primary_beam

# add GLEAM point sources foreground
root_name += 'gleam'
path_point = '/store/ska/sk014/dataset_sdc3/inputs/dataLC_256_train_090523_test/frg/exgf/'
gleam_data = np.loadtxt(path_point+'sdc3cat_skymodel_4deg.txt')
#gleam_data = np.loadtxt(path_point+'gleamcat_skymodel_4deg.txt')

# create inner & outter sky
ra_wrap, dec_wrap = Angle([gleam_data[:,0], gleam_data[:,1]], unit='degree').wrap_at(180 * u.deg).deg

# select inner sky sources
inner_mask = np.sqrt(ra_wrap**2+(dec_wrap+30)**2) <= 2
inner_sky = gleam_data[inner_mask]

# select outter sky sources
outter_mask = np.sqrt(ra_wrap**2+(dec_wrap+30)**2) > 2
outter_sky = gleam_data[outter_mask]
outter_sky[:,2] *= 1e-3

sky = SkyModel()
sky.add_point_sources(inner_sky)
sky.add_point_sources(outter_sky)


telescope = Telescope.read_from_file(path_telescope)

t_start = datetime(2021, 9, 21, 14, 12, 40, 0) # HA between -2h to +2h, obs start at '2021-09-21 14:12:40.1'
#t_obs = timedelta(hours=4, minutes=0, seconds=0, milliseconds=0)
t_obs = timedelta(hours=0, minutes=1, seconds=0, milliseconds=0)
t_day = t_obs

t_int = timedelta(seconds=10)
nr_tsteps = int(t_day.total_seconds() / t_int.total_seconds())
nr_days_obs = int(t_obs.total_seconds() / t_day.total_seconds())

print(' Simulating %d days observation\n time steps: %d' %(nr_days_obs, nr_tsteps))

# name of the outputfile to change manually to the time selected
run_name = '%s_test_ch%d_4h1d_%d' %(path_out+root_name, idx_f, data.shape[0])

observation_settings = Observation(phase_centre_ra_deg=RA.value,
                                        phase_centre_dec_deg=DEC.value,
                                        start_date_and_time=t_start,
                                        start_frequency_hz=freq,
                                        number_of_channels=1,
                                        number_of_time_steps=nr_tsteps,
                                        length=t_day)

simulation = InterferometerSimulation(ms_file_path=run_name+'.MS',
                                        vis_path=run_name+'.vis',
                                        use_gpus=True, use_dask=False,
                                        channel_bandwidth_hz=1e5,
                                        noise_enable=False)

visibilities = simulation.run_simulation(telescope, sky, observation_settings)

# Make Image
imager = Imager('single')
imager.set(fov_deg=FoV.to('deg').value, image_size=2048)
imager.set(input_file=run_name+'.vis', output_root=run_name)
output = imager.run(return_images=1)
image = output["images"][0]
new_w = WCS(naxis=2)
new_w.wcs.crpix = [1025, 1025]
new_w.wcs.cdelt = np.array([-0.0044444444, 0.0044444444])
new_w.wcs.crval = [0, -30]
new_w.wcs.ctype = ["RA---SIN", "DEC--SIN"]






