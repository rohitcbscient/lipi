import numpy as np
import os
import pickle
from karabo.imaging.imager_oskar import OskarDirtyImager
from karabo.imaging import imager_wsclean
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from karabo.simulation.visibility import Visibility
from karabo.imaging.imager_oskar import OskarDirtyImager, OskarDirtyImagerConfig
from karabo.simulator_backend import SimulatorBackend
from astropy.io import fits
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib as mpl
#mpl.use('Tkagg')
from IPython.core.magic import register_cell_magic
@register_cell_magic
def skip(line, cell):
    return


#telescope_path = '/data/rohit/skao_repo/meerkat.tm'
#telescope_path = '/data/rohit/skao_repo/mwa.phase1.tm'
#telescope_path = '/data/rohit/skao_repo/ska1mid.tm'
#telescope_path = '/data/rohit/skao_repo/ska1low.tm'
#telescope_path = '/data/rohit/skao_repo/ska-mid-AAstar.tm'
#telescope_path = '/data/rohit/skao_repo/ska-mid-AA2.tm'
#telescope_path = '/data/rohit/skao_repo/ska-mid-AA1.tm'
#telescope_path = '/data/rohit/skao_repo/ska-low-AAstar.tm'
#telescope_path = '/data/rohit/skao_repo/ska-low-AA2.tm'
#telescope_path = '/data/rohit/skao_repo/ska-low-AA1.tm'
#telescope_path = '/media/rohit/sdata/ska-solar-files/mwa.phase1.tm'
#telescope_path = '/media/rohit/sdata/ska-solar-files/ska1low.tm'
#telescope_name='skamid_AA2_full_snap0_100MHz-2-2000GHz'
telescope_name='skalow_full_point_one_snap0_100MHz-2-2000GHz'
#telescope_name='meerkat_point_rand_snap0_100MHz-2-2000GHz'
#telescope_name='mwa_full_snap0_100MHz-2-2000GHz'
#dtime=datetime(2000, 1, 1, 10, 0, 00, 0) # MeerKAT/ SKA-mid
dtime=datetime(2000, 1, 1, 4, 0, 00, 0) # MWA / SKA-low
hour_=0
minutes_=1
noise_enable_ = False
enable_array_beam=False
point_source_bool = False # tag 1 / tag 0 is ideal case
random_points = False # tag 2
nchan = 1
ntchan = 1
if(random_points==True):
    size_ = 100
else:
    size_ = 1
point_flux = 100.e4

ran = range(18)


layout=np.loadtxt(telescope_path+'/layout.txt')
nant=len(layout)
nb=int(nant*(nant-1)*0.5)
print('Number of Baselines:',nb)
base_length=[0]*nb
k=0
for i in range(nant):
    for j in range(i,nant):
        if(i!=j):
            base_length[k] = np.sqrt((layout[i][0]-layout[j][0])**2 + (layout[i][1]-layout[j][1])**2)
            k=k+1
base_length=np.array(base_length)
print('Maximum Baseline',base_length.max())
#----------------------------
freq_list = np.arange(1,20)*100 # in MHz
beamsize_arr = 3.e8/base_length.max()/(freq_list*1.e6)
beamsize_arr_arcsec = beamsize_arr*180/np.pi*3600
f,ax=plt.subplots(1,1)
ax.plot(freq_list,beamsize_arr_arcsec,'o-')
ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel('Max. Resolution (arcsec)')
ax.set_title('SKA-mid')
plt.close()


path = '/data/rohit/ska-low-sim/'
for i in ran:
    start_frequency_hz_ = freq_list[i]*1.e6
    beam_size_arcsec = 3.e8/start_frequency_hz_/base_length.max()*180/np.pi*3600
    print('Maximum Baseline',base_length.max(), 'Beam (arcsec)',beam_size_arcsec)
    print("Frequency (MHz): ",start_frequency_hz_)
    aa=fits.open('/data/rohit/20151203_240MHz_psimas.fits')
    ms_filename = path+"solar_"+"freq_"+str(int(freq_list[i]))+'_'+telescope_name+".ms"
    vis_filename = path+"solar_"+"freq_"+str(int(freq_list[i]))+'_'+telescope_name+".vis"
    solar_map=aa[0].data;solar_map_jy=solar_map/np.nanmax(solar_map)*20*1.e4*(start_frequency_hz_/2.4e8)**1
    ra_sun_center=249.141666667;dec_sun_center=21.986 #16 34 34.52 -21 59 09.7
    ra_grid,dec_grid=np.meshgrid((np.arange(256)-128)*22.5/3600.,(np.arange(256)-128)*22.5/3600.)
    ra_grid=ra_grid+ra_sun_center;dec_grid=dec_grid+dec_sun_center
    idx=np.where(solar_map>0.1*np.nanmax(solar_map))
    sky_model_ra=ra_grid[idx];sky_model_dec=dec_grid[idx];flux=solar_map_jy[idx]
    # Simulation starts
    sky = SkyModel()
    sky_data = np.array([sky_model_ra, sky_model_dec, flux,np.zeros(len(flux)), \
        np.zeros(len(flux)),np.zeros(len(flux)), np.zeros(len(flux)),np.zeros(len(flux)), \
    np.zeros(len(flux)),np.zeros(len(flux)), np.zeros(len(flux)),np.zeros(len(flux))]).T
    sky_data=sky_data[0:16000,:]
    # Add Point Soure
    add_source = np.zeros((size_,12))
    if(point_source_bool):
        add_source[0][0] = ra_sun_center # RA
        add_source[0][1] = dec_sun_center # DEC
        add_source[0][2] = point_flux # Flux
        sky_data1 = np.vstack((sky_data,add_source))
    elif(random_points):
        ra_point_array = np.random.uniform(low=ra_sun_center-0.3, high=ra_sun_center+0.3, size=size_)
        dec_point_array = np.random.uniform(low=dec_sun_center-0.3, high=dec_sun_center+0.3, size=size_)
        flux_array = np.random.uniform(low=1.e1, high=1.e2, size=size_)
        add_source[:,0] = ra_point_array # RA
        add_source[:,1] = dec_point_array # DEC
        add_source[:,2] = flux_array # Flux
        sky_data1 = np.vstack((sky_data,add_source))
    else:
        sky_data1 = sky_data
    print('Maximum flux density (SFU):',sky_data1[:,2].max()/1.e4)
    print('Number of Sources:',len(sky_data1))
    sky.add_point_sources(sky_data1)
    np.savetxt(path+"source_"+telescope_name+".txt", add_source.T)
    #telescope_name="SKA1MID"
    #telescope=Telescope.get_SKA1_LOW_Telescope()
    backend=SimulatorBackend.OSKAR
    #telescope = Telescope.constructor(telescope_name, backend=backend)
    telescope=Telescope.read_OSKAR_tm_file(telescope_path)
    telescope.read_OSKAR_tm_file(telescope_path)
    simulation = InterferometerSimulation(vis_path=vis_filename, ms_file_path=ms_filename,
        channel_bandwidth_hz=1, time_average_sec=10, noise_enable=noise_enable_, use_gpus=True,
        noise_seed="time", noise_freq="Range", noise_rms="Range",
        noise_start_freq=1.e9, noise_inc_freq=1.e8, noise_number_freq=24,
        noise_rms_start=0, noise_rms_end=0, enable_numerical_beam=enable_array_beam,
        enable_array_beam=enable_array_beam)
    observation = Observation(mode='Tracking',phase_centre_ra_deg=ra_sun_center, start_date_and_time=dtime,
        length=timedelta(hours=hour_, minutes=minutes_, seconds=0, milliseconds=0),
        phase_centre_dec_deg=dec_sun_center, number_of_time_steps=ntchan,
        start_frequency_hz=start_frequency_hz_, frequency_increment_hz=1,
        number_of_channels=nchan, )
    visibility = simulation.run_simulation(telescope, sky, observation, backend=backend)
    npix_per_beam = 3
    imgsize=4096
    cellsize_arcsec=1 #beam_size_arcsec/npix_per_beam
    cellsize_rad=cellsize_arcsec/3600*np.pi/180 # in rad
    print('Cellsize:',cellsize_arcsec,'Beam (arcsec):',beam_size_arcsec,'Max length:',base_length.max())
    print('Field of View (deg):',imgsize*cellsize_arcsec/3600.)
    dirty_imager = OskarDirtyImager(
        OskarDirtyImagerConfig(
            imaging_npixel=imgsize,
            imaging_cellsize=cellsize_rad,
        ))
    dirty_oskar_img = path+"solar_"+"freq_"+str(int(freq_list[i]))+'_'+telescope_name+"_oskar.fits"
    dirty_image = dirty_imager.create_dirty_image(visibility,output_fits_path=dirty_oskar_img)
    dirty_wsclean_img = path+"solar_"+"freq_"+str(int(freq_list[i]))+'_'+telescope_name+"_wsclean.fits"
    img_cmd = 'wsclean \
            -size '+str(imgsize)+' '+str(imgsize)+' \
            -name '+dirty_wsclean_img+' \
            -scale '+str(cellsize_rad)+'rad -niter 25000 -mgain 0.8 \
            -weight uniform\
            -maxuv-l 9000 -minuv-l 10\
            -channels-out '+str(nchan)+' '+ms_filename
    print(img_cmd)
    try:
        restored = imager_wsclean.create_image_custom_command(command=img_cmd)
    except:
        pass # doing nothing on exception



