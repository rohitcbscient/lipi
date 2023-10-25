import oskar
from oskar.measurement_set import MeasurementSet
res='residual_sdc3point_ch750_4h1d.MS'
ska='sdc3point_ch750_4h1d.MS'
#----------------------------------------
res_read=MeasurementSet.open(res,readonly=True)
num_baselines_=int(res_read.num_stations*(res_read.num_stations-1)/2.)
uu,vv,ww=res_read.read_coords(start_row=0,num_baselines=num_baselines_)
res_vis=res_read.read_vis(start_row=0,start_channel=0,num_channels=1,num_baselines=num_baselines_)
imager = oskar.Imager()
imager.fov_deg = 4             # 0.1 degrees across.
imager.image_size = 2048          # 256 pixels across.
imager.set_vis_frequency(150e6)  # 100 MHz, single channel data.
imager.update(uu, vv, ww, res_vis)
res_image = imager.finalise(return_images=1)



#-----------------------------------------
ska_read=MeasurementSet.open(ska,readonly=True)
num_baselines_=int(ska_read.num_stations*(ska_read.num_stations-1)/2.)
uu,vv,ww=ska_read.read_coords(start_row=0,num_baselines=num_baselines_)
ska_vis=ska_read.read_vis(start_row=0,start_channel=0,num_channels=1,num_baselines=num_baselines_)
imager = oskar.Imager()
imager.fov_deg = 4             # 0.1 degrees across.
imager.image_size = 2048          # 256 pixels across.
imager.set_vis_frequency(150e6)  # 100 MHz, single channel data.
imager.update(uu, vv, ww, ska_vis)
ska_image = imager.finalise(return_images=1)

