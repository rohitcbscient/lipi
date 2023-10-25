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
from rascil.apps import rascil_imager

def write_ms(filename,ms_new,skadc_uvw,n,ref_freq_hz):
    ms_new=ms_new.reshape(1440,n)
    uv_new=skadc_uvw.reshape(1440,n,3)
    num_stations=512;num_channels=1;num_pols=1;freq_inc_hz=1.e5
    num_baselines=int(num_stations*(num_stations-1)*0.5)
    num_times=1440 #number of time channels and 10 sec integration
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

def plot_ms(imagepath,uu,vv,ww,res_vis):
	imager = oskar.Imager()
	imager.fov_deg = 4             # 0.1 degrees across.
	imager.image_size = 2048          # 256 pixels across.
	imager.set_output_root(imagepath)
	imager.set_vis_frequency(freqs[i])  # 100 MHz, single channel data.
	imager.update(uu, vv, ww, res_vis[0,:,0])
	res_image = imager.finalise(return_images=1)['images'][0]
	return res_image

n=5
path='/scratch/snx3000/rsharma/galactic_fg_test/'
filename0=path+'residual_sdc3point_ch0000_04.MS'
filename_ave=path+'residual_sdc3point_ch_ave'+str(n)+'_04'
dc=ct.table(filename0)
uvw =dc.getcol("UVW")
ms_ave_data=dc.getcol("DATA")[:,0,0]
for k in range(1,n):
	kk="%04d" % k
	filename_low=path+'residual_sdc3point_ch'+kk+'_04.MS'
	#filename_low='/scratch/snx3000/rsharma/ska_data_reduced/ZW3_IFRQ_for_karabo_'+kk+'.MS'
	dc=ct.table(filename_low)
	dc_data=dc.getcol('DATA')
	ms_ave_data=np.nanmean(np.vstack((ms_ave_data,dc_data[:,0,0])),axis=0)
filename_high=path+'residual_sdc3point_ch0700_04.MS'
ms_diff=ms_data

freqs = 1.06e8+np.arange(901) * 1.e5 # Hz
freq_ave=freqs[0:n].mean()
write_ms(filename_ave,ms_data,uvw,130816,freq_ave)

image_ave_path=path+'image_ave_sdc3point_ch_ave'+str(n)+'_04'
res_read=MeasurementSet.open(filename_ave,readonly=True)
num_baselines_=int(res_read.num_stations*(res_read.num_stations-1)/2.)
uu,vv,ww=res_read.read_coords(start_row=0,num_baselines=num_baselines_)
vis=res_read.read_vis(start_row=0,start_channel=0,num_channels=1,num_baselines=num_baselines_)
ave_image=plot_ms(image_ave_path,uu,vv,ww,vis)

#-----------------------------
class Args:
	ingest_msname='residual_sdc3point_ch_ave5_04.MS'
	mode='ical';logfile='rascil_logfile';use_dask="False"
	ingest_dd= [0, 1, 2, 3];ingest_vis_nchan=1;imaging_npixel=2048
	imaging_cellsize=1.e-5;ingest_chan_per_vis =1
	clean_nmajor=5
	ingest_average_vis=True
	performance_file='rascil_performance_file'
	imaging_pol='XX';imaging_weighting='natural';imaging_gaussian_taper=False
	imaging_context='ng';imaging_flat_sky=False
	clean_algorithm='hogbom';clean_beam=[5.e-5,5.e-5,0];clean_scales=[0,6,10];clean_niter=0;clean_fractional_threshold=0;clean_threshold=0
	clean_component_threshold=0;clean_gain=0.1;clean_component_method='hogbom'
	clean_nmoment=5;clean_psf_support=3;clean_restored_output='rascil_output';clean_facets=1;clean_overlap=32;clean_taper='tukey';clean_restore_facets=1
	clean_restore_overlap=32;clean_restore_taper='tukey';imaging_dft_kernel='gpu_raw';perform_flagging=False;flagging_strategy_name="generic-default.lua"
	imaging_ng_threads=4
	imaging_w_stacking=True

	

args=Args()
rascil_imager.imager(args)
#-----------------------------


from astropy.io import fits
import numpy as np
plot_test=1
if(plot_test):
	ska_img=fits.open('../residual1_pybdsf/ska_images/ska_image_ch0002_I.fits');ska_image=ska_img[0].data[0]
	img0=fits.open(image_ave_path+'_I.fits');d0=img0[0].data[0]
	img1=fits.open('/scratch/snx3000/rsharma/residual1_pybdsf/res_images/res_image_ch0000_04_I.fits');d1=img1[0].data[0]
	img2=fits.open('/scratch/snx3000/rsharma/residual1_pybdsf/res_images/res_image_ch0700_04_I.fits');d2=img2[0].data[0]
	f,((ax0,ax1),(ax2,ax3))=plt.subplots(2,2,sharex=True,sharey=True)
	ax0.imshow(d0,origin='lower',cmap='coolwarm',extent=[-2,2,-32,-28],aspect='auto',vmin=-0.1,vmax=0.1)
	ax0.set_xlabel('RA (deg)');ax1.set_xlabel('RA (deg)');ax0.set_ylabel('Dec (deg)');ax1.set_ylabel('Dec (deg)')
	ax1.imshow(d1,origin='lower',cmap='coolwarm',extent=[-2,2,-32,-28],aspect='auto',vmin=-0.1,vmax=0.1)
	ax2.imshow(d1-d0,origin='lower',cmap='coolwarm',extent=[-2,2,-32,-28],aspect='auto',vmin=-0.1,vmax=0.1)
	ax3.imshow(ska_image,origin='lower',cmap='coolwarm',extent=[-2,2,-32,-28],aspect='auto',vmin=-0.1,vmax=0.1)
	ax0.set_title('Average')
	ax1.set_title('Channel 000')
	ax2.set_title('Channel diff')
	ax3.set_title('SKA Image')
	plt.show()



