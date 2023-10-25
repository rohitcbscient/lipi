import casacore.tables as ct
from oskar.measurement_set import MeasurementSet
import oskar
import os
import numpy as np
import tempfile
import matplotlib.pyplot as plt
import numpy as np
import pytest
from karabo.imaging.imager import Imager
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from karabo.sourcedetection.evaluation import SourceDetectionEvaluation
from karabo.sourcedetection.result import (
    PyBDSFSourceDetectionResult,
    SourceDetectionResult,
)
from karabo.test.conftest import NNImageDiffCallable, TFiles
from astropy.io import fits
import bdsf
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.coordinates import Angle
import astropy.units as u
from datetime import datetime, timedelta
import tools21cm as t2c
import sys
from astropy.wcs import WCS

def do_bdsf(imagepath,output_model,freqs,idx_f):
    img=fits.open(imagepath+'_I.fits')
    head=img[0].header
    data=np.abs(img[0].data[0])
    head['CDELT3']=1
    uvw_max=73508.6 # max  
    b_deg=(3.e8/freqs[idx_f])/uvw_max*180/np.pi
    head['BMIN']= b_deg # in Deg 2.8 arcsec for 120 km baseline
    head['BMAJ']=b_deg
    out_imagepath=imagepath+'_mod_I.fits'
    fits.writeto(out_imagepath,data,head,overwrite=True)
    out_imagedata=fits.open(out_imagepath)[0].data
    img_bdsf = bdsf.process_image(out_imagepath, beam=(10/3600.,10/3600.,0),thresh_pix=5,thresh_isl=3,thresh='hard')
    bdsf_source_list='/scratch/snx3000/rsharma/bdsf_source_list.fits'
    img_bdsf.write_catalog(outfile=bdsf_source_list,catalog_type='gaul',format='fits', clobber=True)
    sl=fits.open(bdsf_source_list)
    xpix,ypix= sl[1].data['Xposn'],sl[1].data['Yposn'];xpix=xpix.astype(int);ypix=ypix.astype(int)
    RA=sl[1].data['RA'];DEC=sl[1].data['DEC']
    sl_peak=data[ypix,xpix]#sl[1].data['Peak_flux']
    RA[np.where(RA>200)]=RA[np.where(RA>200)]-360
    sl_bmaj=sl[1].data['Maj']
    sl_bmin=sl[1].data['Min']
    sl_bpa=sl[1].data['PA']
    gauss_n_size=RA.shape[0]
    sky_data=np.hstack((RA.reshape(gauss_n_size,1),DEC.reshape(gauss_n_size,1),sl_peak.reshape(gauss_n_size,1),np.zeros((gauss_n_size,9))))
    sky_data[:,9]=sl_bmaj
    sky_data[:,10]=sl_bmin
    sky_data[:,11]=sl_bpa
    np.savetxt(X=sky_data,fname=output_model)
    w = WCS(head) 
    #xpix,ypix,_=w.wcs_world2pix(RA,DEC,0,0)
    return sky_data,data,xpix,ypix

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
    num_times=1 # 1440 number of time channels and 10 sec integration
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
        uu[:] = uv_new[t,:,0];vv[:] = uv_new[t,:,1];ww[:] = uv_new[t,:,2]
        #uu[:] = uv_new[:,0];vv[:] = uv_new[:,1];ww[:] = uv_new[:,2]
        for c in range(num_channels):
            for b in range(num_baselines):
                vis[t, c, b, :] = ms_new[t,b]
                #vis[t, c, b, :] = ms_new[b]
        start_row = t * num_baselines
        ms.write_coords(start_row, num_baselines, uu, vv, ww,exposure_sec, interval_sec, time_stamp)
        ms.write_vis(start_row, 0, num_channels, num_baselines, vis[t, ...])

def plot_bdsf(data_img,xpix,ypix):
	f,ax=plt.subplots(1,1)
	im=ax.imshow(data_img,aspect='auto',origin='lower',vmin=-1.e-3,vmax=1.e-3)
	ax.plot(xpix,ypix,'o',color='red')
	f.colorbar(im)
	plt.show()

#i = int(os.environ['SLURM_ARRAY_TASK_ID'])+585
#i=0
for i in range(868,887):
	ii="%04d" % i
	freqs = 1.06e8+np.arange(901) * 1.e5 # Hz
	#os.system('rm -rf /scratch/snx3000/rsharma/pybdsf_tests/*')
	filename_ska='/scratch/snx3000/rsharma/ska_data_reduced/ZW3_IFRQ_for_karabo_'+ii+'.MS'
	filename_point0='/scratch/snx3000/rsharma/gleam_point_source_ms/GLEAM_point_sources_'+ii+'.MS'
	skadc=ct.table(filename_ska)
	point=ct.table(filename_point0)

	skadc_data=skadc.getcol('DATA')
	skadc_uvw = skadc.getcol("UVW")
	ms_new_ska=skadc_data[:,0,0]
	#-----------------------------------------
	ska_imagepath='/scratch/snx3000/rsharma/residual1_pybdsf/ska_images/ska_image_ch'+ii
	ska_read=MeasurementSet.open(filename_ska,readonly=True)
	num_baselines_=int(ska_read.num_stations*(ska_read.num_stations-1)/2.)
	uu,vv,ww=ska_read.read_coords(start_row=0,num_baselines=num_baselines_)
	ska_vis=ska_read.read_vis(start_row=0,start_channel=0,num_channels=1,num_baselines=num_baselines_)
	imager = oskar.Imager()
	imager.fov_deg = 4             # 0.1 degrees across.
	imager.image_size = 2048          # 256 pixels across.
	imager.set_vis_frequency(180e6)  # 100 MHz, single channel data.
	imager.set(output_root=ska_imagepath)
	imager.update(uu, vv, ww, ska_vis[0,:,0])
	ska_image = imager.finalise(return_images=1)['images'][0]
	print('SKA Image Maximum:',np.max(ska_image))
	freqs = 1.06e8+np.arange(901) * 1.e5 # Hz

	#--- Make Point Image-----
	imagepath='/scratch/snx3000/rsharma/residual1_pybdsf/res_ch'+ii
	do_point_source_det=1
	nk=5;ms_new_ska_array=[0]*(nk+1)
	ms_new_ska_array[0]=ms_new_ska
	uvlim=1.e12 # in mts
	for k in range(nk):
	        kk="%02d" % k
	        print('Iteration: '+kk)
	        if((k>0) & (do_point_source_det==1)):
	            print('Doing pybdf.....')
	            output_model='/scratch/snx3000/rsharma/residual1_pybdsf/pybdsf_skymodel/pybdsf_skymodel_ch'+ii+'_'+kk+'.txt'
	            bdsf_data,data_img,xpix,ypix= do_bdsf(imagepath,output_model,freqs,i)
	            #print(bdsf_data[:,2])
	            #f,ax=plt.subplots(1,1)
	            #im=ax.imshow(data_img,aspect='auto',origin='lower',vmin=-1.e-3,vmax=1.e-3)
	            #ax.plot(xpix,ypix,'o',color='red')
	            #f.colorbar(im)
	            #plt.show()
		#---------------- Do Point Source Simulation-------------------------
	        imagepath='/scratch/snx3000/rsharma/residual1_pybdsf/res_images/res_image_ch'+ii+'_'+"%02d" % k
	        path_out = '/scratch/snx3000/rsharma/residual1_pybdsf/res_images/'
	        root_name='point_source_iteration_ch'+ii
	        path_telescope = '/store/ska/sk014/dataset_sdc3/inputs/telescope.tm'
	        run_name =path_out+root_name
	        filename_point=run_name+'_'+kk+'.MS'
	        iono_fits='/scratch/snx3000/rsharma/atmo/screen_4h_i0_ch'+"%03d" % i+'.fits'
	        r0, sampling = 7e3, 100.0
	        inner_sky,outter_sky=get_gleam(root_name)
	        sky = SkyModel()
	        #sky_test=np.zeros((1,13));sky_test[0][0]=-1;sky_test[0][1]=-31;sky_test[0][2]=1
	        #sky.add_point_sources(sky_test)
	        if(k==0):
	               sky.add_point_sources(inner_sky);sky.add_point_sources(outter_sky)
	        else:
	               sky.add_point_sources(bdsf_data)
	        telescope = Telescope.read_from_file(path_telescope)
	        t_start = datetime(2021, 9, 21, 14, 12, 40, 0) # HA between -2h to +2h, obs start at '2021-09-21 14:12:40.1'
	        t_obs = timedelta(hours=4, minutes=0, seconds=0, milliseconds=0)
	        #t_obs = timedelta(hours=0, minutes=1, seconds=0, milliseconds=0)
	        observation_settings = Observation(phase_centre_ra_deg=0,mode="Tracking",
	                                        phase_centre_dec_deg=-30,
	                                        start_date_and_time=t_start,
	                                        start_frequency_hz=freqs[i],
	                                        number_of_channels=1,
	                                        number_of_time_steps=1440,
	                                        length=t_obs)
	        simulation = InterferometerSimulation(ms_file_path=filename_point,
	                                        vis_path=run_name+'_'+kk+'.vis',
	                                        use_gpus=True, use_dask=False,
	                                        channel_bandwidth_hz=1e5, enable_numerical_beam=False,enable_array_beam=False, 	noise_enable=False,ionosphere_fits_path=iono_fits,ionosphere_screen_type="External",ionosphere_screen_height_km=r0,ionosphere_screen_pixel_size_m=sampling,ionosphere_isoplanatic_screen=True)
	        #if(os.path.isfile(run_name+'_'+kk+'.vis')==False):
	        print('Simulating Point Visibilities..'+filename_point)
	        visibilities = simulation.run_simulation(telescope, sky, observation_settings)
		#--------------- Subtractig the Point Source Vis
	        point1=ct.table(filename_point)
	        point_data=point1.getcol('DATA');point_uvw = point1.getcol("UVW")
	        ms_new_point=point_data[:,0,0]
	        print('Point Source Max:',ms_new_point.max())
	        ms_new_res=ms_new_ska_array[k]-ms_new_point
	        filename_res='/scratch/snx3000/rsharma/residual1_pybdsf/ms/residual_sdc3point_ch'+ii+'_'+kk+'.MS'
	        write_ms(filename_res,ms_new_res,skadc_uvw,130816)
	        ms_new_ska_array[k+1]=ms_new_res
	        #------- Image Residuals
	        imager = oskar.Imager()
	        imager.fov_deg = 4             # 0.1 degrees across.
	        imager.image_size = 2048          # 256 pixels across.
	        imager.set_output_root(imagepath)
	        imager.set_vis_frequency(freqs[i])  # 100 MHz, single channel data.
	        imager.update(uu, vv, ww, ms_new_res)
	        res_image = imager.finalise(return_images=1)['images'][0]
	
from astropy.io import fits
import numpy as np
plot_test=0
ch='000'
if(plot_test):
	ska_img=fits.open('../ska_images/ska_image_ch0'+ch+'_I.fits');ska_image=ska_img[0].data[0]
	img0=fits.open('res_image_ch0'+ch+'_00_I.fits');d0=img0[0].data[0]
	img1=fits.open('res_image_ch0'+ch+'_02_I.fits');d1=img1[0].data[0]
	img2=fits.open('res_image_ch0'+ch+'_04_I.fits');d2=img2[0].data[0]
	f,((ax0,ax1),(ax2,ax3))=plt.subplots(2,2,sharex=True,sharey=True)
	ax0.imshow(d0,origin='lower',cmap='coolwarm',extent=[-2,2,-32,-28],aspect='auto',vmin=-0.1,vmax=0.1)
	ax0.set_xlabel('RA (deg)');ax1.set_xlabel('RA (deg)');ax0.set_ylabel('Dec (deg)');ax1.set_ylabel('Dec (deg)')
	ax1.imshow(d1,origin='lower',cmap='coolwarm',extent=[-2,2,-32,-28],aspect='auto',vmin=-0.1,vmax=0.1)
	ax2.imshow(d2,origin='lower',cmap='coolwarm',extent=[-2,2,-32,-28],aspect='auto',vmin=-0.1,vmax=0.1)
	ax3.imshow(ska_image,origin='lower',cmap='coolwarm',extent=[-2,2,-32,-28],aspect='auto',vmin=-0.1,vmax=0.1)
	ax0.set_title('Residual0')
	ax1.set_title('Residual1')
	ax2.set_title('Residual2')
	ax3.set_title('SKA Image')
	plt.show()
