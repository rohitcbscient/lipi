import casacore.tables as ct
import numpy as np
import sys
import os
import oskar
from astropy.coordinates import Angle
import astropy.units as u
from karabo.imaging.imager import Imager
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from datetime import datetime, timedelta

i = int(os.environ['SLURM_ARRAY_TASK_ID'])+750
#i=0
def get_gleam(root_name):
	# add GLEAM point sources foreground
	root_name += 'gleam'
	path_point = '/store/ska/sk014/dataset_sdc3/inputs/dataLC_256_train_090523_test/frg/exgf/'
	gleam_data = np.loadtxt(path_point+'rohit_sdc3cat_skymodel_4deg.txt')
	gleam_data[:,2]=gleam_data[:,2]
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
	return inner_sky,outter_sky


def write_ms(filename,ms_new,skadc_uvw,n):
    ms_new=ms_new.reshape(1440,n)
    uv_new=skadc_uvw.reshape(1440,n,3)
    num_stations=512;num_channels=1;num_pols=1;ref_freq_hz=1.8e8;freq_inc_hz=1.e5
    num_baselines=n#int(num_stations*(num_stations-1)*0.5)
    num_times=1440 # 1440 number of time channels and 10 sec integration
    num_pols=1
    num_channels=1
    ms = oskar.MeasurementSet.create(filename, num_stations,num_channels, num_pols,ref_freq_hz, freq_inc_hz)
    ra_rad = 0.0
    dec_rad = -30.0
    exposure_sec = 1.0
    interval_sec = 1.0
    ms.set_phase_centre(ra_rad, dec_rad)
    mjd_20210921=59478.592071770836
    #UTC 2021-09-21T14:12:35.1
    #Modified Julian Day (MJD)   59478.592071770836
    #Julian Day (JD) 2459479.0920717707
    uu = np.zeros([num_baselines])
    vv = np.zeros_like(uu)
    ww = np.zeros_like(uu)
    vis = np.zeros([num_times, num_channels,
                       num_baselines, num_pols], dtype='c8')
    for t in range(num_times):
        time_stamp = mjd_20210921 * 86400.0 + t
        uu[:] = uv_new[t,:,0]
        vv[:] = uv_new[t,:,1]
        ww[:] = uv_new[t,:,2]
        for c in range(num_channels):
            for b in range(num_baselines):
                vis[t, c, b, :] = ms_new[t,b]
        start_row = t * num_baselines
        ms.write_coords(start_row, num_baselines, uu, vv, ww,exposure_sec, interval_sec, time_stamp)
        ms.write_vis(start_row, 0, num_channels, num_baselines, vis[t, ...])

ii="%04d" % i
freqs = 1.06e8+np.arange(901) * 1.e5 # Hz
imagepath='/scratch/snx3000/rsharma/gleam_point_source_ms/image_'+ii
path_out = '/scratch/snx3000/rsharma/gleam_point_source_ms/'
root_name='GLEAM_point_sources'
path_telescope = '/scratch/snx3000/rsharma/telescope.tm'
#path_telescope = '/users/mibianco/codes/sdc3/data/files/telescope.tm'
run_name =path_out+root_name
filename_point=run_name+'_'+ii+'.MS'
iono_fits='/scratch/snx3000/rsharma/atmo/screen_4h_i0_ch'+"%03d" % i+'.fits'
r0, sampling = 7e3, 100.0
sky = SkyModel()
inner_sky,outter_sky=get_gleam(root_name)
sky.add_point_sources(inner_sky)
#sky.add_point_sources(inner_sky);sky.add_point_sources(outter_sky)
telescope = Telescope.read_from_file(path_telescope)
t_start = datetime(2021, 9, 21, 14, 12, 40, 0) # HA between -2h to +2h, obs start at '2021-09-21 14:12:40.1'
t_day = timedelta(hours=4, minutes=0, seconds=0, milliseconds=0)
t_int = timedelta(seconds=10)
nr_tsteps = int(t_day.total_seconds() / t_int.total_seconds())
nr_days_obs = int(t_day.total_seconds() / t_day.total_seconds())
print(' Simulating %d days observation\n time steps: %d' %(nr_days_obs, nr_tsteps))
# name of the outputfile to change manually to the time selected
observation_settings = Observation(phase_centre_ra_deg=0,mode="Tracking",
                                        phase_centre_dec_deg=-30,
                                        start_date_and_time=t_start,
                                        start_frequency_hz=freqs[i],
                                        number_of_channels=1,
                                        number_of_time_steps=nr_tsteps,
                                        length=t_day)
simulation = InterferometerSimulation(
            ms_file_path=filename_point,
            vis_path=run_name+'_'+ii+'.vis',
            channel_bandwidth_hz=1e6,
            time_average_sec=1,)
simulation = InterferometerSimulation(ms_file_path=filename_point,
                                       vis_path=run_name+'_'+ii+'.vis',
                                        use_gpus=True, use_dask=False, channel_bandwidth_hz=1e5,noise_enable=False,ionosphere_fits_path=iono_fits,ionosphere_screen_type="External",ionosphere_screen_height_km=r0, ionosphere_screen_pixel_size_m=sampling,ionosphere_isoplanatic_screen=True)
#if(os.path.isfile(run_name+'_'+kk+'.vis')==False):
visibilities = simulation.run_simulation(telescope, sky, observation_settings)

