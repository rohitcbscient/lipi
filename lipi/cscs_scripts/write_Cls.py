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

def plot_bdsf(data_img,xpix,ypix):
	f,ax=plt.subplots(1,1)
	im=ax.imshow(data_img,aspect='auto',origin='lower',vmin=-1.e-3,vmax=1.e-3)
	ax.plot(xpix,ypix,'o',color='red')
	f.colorbar(im)
	plt.show()

ilist=[753]
#i = int(os.environ['SLURM_ARRAY_TASK_ID'])
i=ilist[0]
#i=0
ii="%04d" % i
freqs = 1.06e8+np.arange(901) * 1.e5 # Hz
#os.system('rm -rf /scratch/snx3000/rsharma/pybdsf_tests/*')
#--------------------------------------------------------------------
filename_ska='/scratch/snx3000/rsharma/ska_data_reduced/ZW3_IFRQ_for_karabo_'+ii+'.MS'
#filename_ska='/scratch/snx3000/mibianco/output_sdc3/dataLC_130923/ms/lc_256_train_130923_i0_dTnoisegainiongfpoint_ch'+ii[1:]+'_4h1d_256.MS'
#--------------------------------------------------------------------
skadc=ct.table(filename_ska)
skadc_data_=skadc.getcol('DATA')
skadc_uvw = skadc.getcol("UVW")
skadc_data=skadc_data_[:,0,0]
uv_radial=np.sqrt(skadc_uvw[:,0]**2 + skadc_uvw[:,1]**2)
#uv_phase=np.arccos(skadc_uvw[:,0]/np.sqrt(skadc_uvw[:,0]**2+skadc_uvw[:,1]**2))
n_uvbins=50
uvradial_bins=np.logspace(np.log10(uv_radial.min()),np.log10(uv_radial.max()),n_uvbins)
uvradial_num= np.histogram(uv_radial,bins=uvradial_bins)[0]
del_uvradial_bins=uvradial_bins[1:]-uvradial_bins[:-1]
l=2*np.pi*uvradial_bins
idx=[0]*n_uvbins;Vabs=[0]*n_uvbins;Vcross=[0]*n_uvbins
idx[0]=np.where((uvradial_bins[0]<uv_radial)&(uv_radial<uvradial_bins[1]))[0]
for i in range(1,n_uvbins):
	print(i)
	idx[i]=np.where((uvradial_bins[i]<uv_radial)&(uv_radial<uvradial_bins[i+1]))[0]
	Vabs[i]=np.mean(np.abs(skadc_data[idx[i]]))
	Vcross[i]=np.mean(skadc_data[idx[i]])*np.mean(np.conj(skadc_data[idx[i-1]]))






































