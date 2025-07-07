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
import sys
from karabo_simutils import *
@register_cell_magic
def skip(line, cell):
    return



#----- Define your Telescopes for the Simulation
skao_repo_tel_path = '/data/rohit/skao_repo/'
tel_all = get_telescope(skao_repo_tel_path) # 0 is SKA, and 1 is SKA-precursor
telescope_list = ['skalow','skamid','mwa','meerkat']
ska_aa_list = ['AA0.5','AA1','AA2','AAs','Full']
mid_bands=['1','3']
nt=1; tel = telescope_list[nt] # Put in your telescope in `tel`

#----- Define your Frequency Array for the Simulation
freq_samples = 10
ska_freq_param = ska_frequency_params(bandwidth_samples=freq_samples)
if(tel[0:3]=='ska' ):
    if(tel=='skamid' or tel=='meerkat'):
        bands=mid_bands
        midbands_list=['1','2','3','4','5a','5b']
        freq_array=[0]*len(mid_bands);i=0
        for m in mid_bands:
            bandsidx_ = midbands_list.index(m)
            freq_array[i] = ska_freq_param['mid'][0][bandsidx_]
            i+=1
        freq_array = np.array(freq_array).flatten()
    else:
        freq_array = ska_freq_param['low'][0][0]
    
#----- Define the Observations
if(tel=='skamid' or tel=='meerkat'):
    dtime=datetime(2000, 1, 1, 10, 0, 00, 0) # MeerKAT/ SKA-mid
if(tel=='skalow' or tel=='mwa'):
    dtime=datetime(2000, 1, 1, 4, 0, 00, 0) # MWA / SKA-low
hour_=0
minutes_=1
noise_enable_ = False
enable_array_beam=False

#-------- Define SkyModel Parameters
skymodel_path='/data/rohit/20151203_240MHz_psimas.fits'
ra_sun_center=249.141666667;dec_sun_center=21.986 #16 34 34.52 -21 59 09.7
skymodel_cellsize=22.5 # in arcsec
sm_threshold = 0.1
sm_fscale = 1
add_sm = 'random' # 'point' / 'random'
point_source_bool = False # tag 1 / tag 0 is ideal case
point_flux = 100.e4
random_points = False # tag 2
if(add_sm=='random'):
    size_ = 100
else:
    size_ = 1
randsize = 100
flux_rand_low = 10 # Jy
flux_rand_high = 100 # Jy
angle_rand = 0.3 # in degrees
source_cut=16000
save_sources = True

#-------- Define Image Parameters
#npix_per_beam = 3
imgsize=4096
cellsize_arcsec=1 #beam_size_arcsec/npix_per_beam
cellsize_rad=cellsize_arcsec/3600*np.pi/180 # in rad
imager_str = 'wsclean'
niter = 18000
maxuv = 9000 # in lambda units
minuv = 10 # in lambda units
nchan = 1
ntchan = 1
weight = 'uniform'
    
for naa in range(len(ska_aa_list)): # loop over array assemblies 
#for naa in range(1): # loop over array assemblies 
    telescope_path = tel_all[0][telescope_list[nt]][naa] # 0 --> length of 4 where 0 is SKA-full and 4 is AA0.5
    bl_array,maxbl,medbl=get_telescope_config(telescope=telescope_path,info=True) # bl -> Baseline
    
    for j in range(len(freq_array)): # Loop over the frequencies
    #for j in range(1): # Loop over the frequencies
        start_frequency_hz_ = freq_array[j]*1.e6
        #----------- Define strings and names
        path_ = '/data/rohit/ska-solar-sim-repo/'+tel+'/'
        sm_type = 'rand'
        tel_str = telescope_path.split('repo/')[1].split('.t')[0]
        sm_str = 'forward-20151203-'+sm_type
        freq_str = 'freq-'+str(np.round(freq_array[j],2))
        obs_str = 'snap-1-min-1-channel'
        img_str = 'img-cellsize-'+str(cellsize_arcsec)+'-imsize-'+str(imgsize)
        vis_str = path_+'solar_'+tel_str+'_'+sm_str+'_'+freq_str+'_'+obs_str+'_sim'
        sm_save_str = path_+'solar_'+tel_str+'_'+sm_str+'_'+freq_str+'_'+obs_str+'_skymodel'
        img_str = vis_str+'_'+img_str
        
        #------------- Get Skymodel
        sky_data1,solar_map,solar_map_jy=get_solar_skydata(skymodel_path,sm_save_str, start_frequency_hz_,\
            ra_sun_center,dec_sun_center, skymodel_cellsize,sm_fscale,\
            sm_threshold, save_sources,source_cut,\
            add_sm, randsize, point_flux, flux_rand_low, flux_rand_high,angle_rand)
        
        #------------ Simulation starts
        sky = SkyModel()
        sky.add_point_sources(sky_data1)
        backend=SimulatorBackend.OSKAR
        telescope=Telescope.read_OSKAR_tm_file(telescope_path)
        telescope.read_OSKAR_tm_file(telescope_path)
        simulation = InterferometerSimulation(vis_path=vis_str+'.vis', ms_file_path=vis_str+'.ms',
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
                
        #--------- Get the images
        image=make_images(visibility,imager_str,imgsize,img_str,cellsize_rad,\
                        niter,weight,maxuv,minuv,nchan,vis_str)


sys.exit()



#----------------------------
#freq_list = freq_array
#beamsize_arr = 3.e8/base_length.max()/(freq_list*1.e6)
#beamsize_arr_arcsec = beamsize_arr*180/np.pi*3600
#f,ax=plt.subplots(1,1)
#ax.plot(freq_list,beamsize_arr_arcsec,'o-')
#ax.set_xlabel('Frequency (MHz)')
#ax.set_ylabel('Max. Resolution (arcsec)')
#ax.set_title('SKA-mid')
#plt.close()


for i in ran:
    start_frequency_hz_ = freq_list[i]*1.e6
    beam_size_arcsec = 3.e8/start_frequency_hz_/base_length.max()*180/np.pi*3600
    print('Maximum Baseline',base_length.max(), 'Beam (arcsec)',beam_size_arcsec)
    print("Frequency (MHz): ",start_frequency_hz_)
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






