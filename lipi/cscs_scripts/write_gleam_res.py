import casacore.tables as ct
import numpy as np
import sys
import os
from oskar.measurement_set import MeasurementSet
import oskar

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

def plot_ms(imagepath,uu,vv,ww,res_vis):
	imager = oskar.Imager()
	imager.fov_deg = 4             # 0.1 degrees across.
	imager.image_size = 2048          # 256 pixels across.
	imager.set_output_root(imagepath)
	imager.set_vis_frequency(freqs[i])  # 100 MHz, single channel data.
	imager.update(uu, vv, ww, res_vis[0,:,0])
	res_image = imager.finalise(return_images=1)['images'][0]
	return res_image

i = int(os.environ['SLURM_ARRAY_TASK_ID'])+0
#i=0
freqs = 1.06e8+np.arange(901) * 1.e5 # Hz
ii="%04d" % i
skams='/scratch/snx3000/rsharma/ska_data_reduced/ZW3_IFRQ_for_karabo_'+ii+'.MS'
pointms='/scratch/snx3000/rsharma/gleam_point_source_ms/GLEAM_point_sources_'+ii+'.MS'
skadc=ct.table(skams)
point=ct.table(pointms)

skadc_data=skadc.getcol('DATA')
skadc_uvw = skadc.getcol("UVW")
ms_new_ska=skadc_data[:,0,0]
ska_read=MeasurementSet.open(skams,readonly=True)
num_baselines_=int(ska_read.num_stations*(ska_read.num_stations-1)/2.)
ska_uu,ska_vv,ska_ww=ska_read.read_coords(start_row=0,num_baselines=num_baselines_)
ska_vis=ska_read.read_vis(start_row=0,start_channel=0,num_channels=1,num_baselines=num_baselines_)
ska_imagepath='/scratch/snx3000/rsharma/residual0_gleam/images/ska_'+ii
ska_image=plot_ms(ska_imagepath,ska_uu,ska_vv,ska_ww,ska_vis)

point_data=point.getcol('DATA')
point_uvw = point.getcol("UVW")
ms_new_point=point_data[:,0,0]
point_read=MeasurementSet.open(pointms,readonly=True)
num_baselines_=int(point_read.num_stations*(point_read.num_stations-1)/2.)
point_uu,point_vv,point_ww=ska_read.read_coords(start_row=0,num_baselines=num_baselines_)
point_vis=point_read.read_vis(start_row=0,start_channel=0,num_channels=1,num_baselines=num_baselines_)
point_imagepath='/scratch/snx3000/rsharma/residual0_gleam/images/point_'+ii
point_image=plot_ms(point_imagepath,point_uu,point_vv,point_ww,point_vis)


filename_res='/scratch/snx3000/rsharma/residual0_gleam/residual0_gleam_'+ii+'.MS'
ms_new_res=ms_new_ska-ms_new_point
write_ms(filename_res,ms_new_res,skadc_uvw,130816)
#---
uvlim=1000
filename_res_uvlim='/scratch/snx3000/rsharma/residual0_gleam/residual0_gleam_uvlim'+str(uvlim)+'_'+ii+'.MS'
uvdist=np.sqrt(ska_uu**2 + ska_vv**2)
idx_low=np.where(uvdist<uvlim)[0]
ms_new_res1=ms_new_res.reshape(1440,130816)[:,idx_low].flatten()
res1_uvw=skadc_uvw.reshape(1440,130816,3)[:,idx_low,:].flatten()
write_ms(filename_res_uvlim,ms_new_res1,res1_uvw,idx_low.shape[0])

#----------------- Make Residual Image-----------------------
imagepath='/scratch/snx3000/rsharma/residual0_gleam/images/residual0_gleam_'+ii
res_read=MeasurementSet.open(filename_res,readonly=True)
num_baselines_=int(res_read.num_stations*(res_read.num_stations-1)/2.)
uu,vv,ww=res_read.read_coords(start_row=0,num_baselines=num_baselines_)
res_vis=res_read.read_vis(start_row=0,start_channel=0,num_channels=1,num_baselines=num_baselines_)
res_image=plot_ms(imagepath,uu,vv,ww,res_vis)

#----------------- Residual small baselines-------------------
imagepath_uvlim='/scratch/snx3000/rsharma/residual0_gleam/images/residual0_gleam_uv'+str(uvlim)+'_'+ii
res_read_uvlim=MeasurementSet.open(filename_res_uvlim,readonly=True)
num_baselines_uvlim=int(res_read_uvlim.num_stations*(res_read_uvlim.num_stations-1)/2.)
uu_uvlim,vv_uvlim,ww_uvlim=res_read_uvlim.read_coords(start_row=0,num_baselines=num_baselines_uvlim)
res_vis_uvlim=res_read_uvlim.read_vis(start_row=0,start_channel=0,num_channels=1,num_baselines=num_baselines_uvlim)
res_image_uvlim=plot_ms(imagepath_uvlim,uu_uvlim,vv_uvlim,ww_uvlim,res_vis_uvlim)


