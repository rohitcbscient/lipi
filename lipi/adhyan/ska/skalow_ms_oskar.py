#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os, oskar

from astropy.io import fits
from astropy.time import Time, TimeDelta
import astropy.wcs as awcs
import astropy.coordinates as coord
import astropy.units as u
import astropy.constants as cst

path_in = '/project/c31/SKA_low_images/eos_fits/'
path_out = '/scratch/snx3000/mibianco/test_skalow/'

# --- Sky model ---
filename = 'EOS202MHz_21cm'
with fits.open('%s.fits' %(path_in+filename), mode="readonly", memmap=True) as hdulist:
    hdr = hdulist[0].header
    data = hdulist[0].data
    WCS = awcs.WCS(hdr)

# calculate lat, lon base on the fits information
N_1, N_2 = hdr['NAXIS1'], hdr['NAXIS2']
idx_1 = np.arange(1, N_1 + 1).reshape(N_1, 1)
idx_2 = np.arange(1, N_2 + 1).reshape(1, N_2)

d_1, d_2 = WCS.celestial.all_pix2world(idx_1, idx_2, 1)  # [deg]
if WCS.wcs.lng == 0:
    lon, lat = d_1, d_2
else:
    lat, lon = d_1, d_2

# create sky model
sky_coord = np.vstack((lon[np.newaxis, ...], lat[np.newaxis, ...])).reshape(2,N_1*N_2) # lon, lat [deg]
sky_data = np.hstack((sky_coord.T, data.flatten()[..., np.newaxis]))
sky = oskar.Sky.from_array(sky_data, 'single')

# Field of view
FoV = hdr['CDELT1']*N_1          # [deg]
RA, DEC = hdr['CRVAL1']*u.deg, hdr['CRVAL2']*u.deg
sky_detected_sources = sky.create_copy()
sky_detected_sources.filter_by_radius(0.0, FoV, RA.value, DEC.value)
print("Number of sources in sky model: %d"  %sky.num_sources)
print("Percent of sources detected in sky model: %.2f %%" %(100*sky_detected_sources.num_sources/sky.num_sources))

# Observation starts such that the phase centre crosses the meridian at obs_length/2
obs_length = (6.*u.h).to(u.s)
t = Time('J2000', scale='utc', location=(116.7644482, -26.82472208))
obs_start = t - obs_length/2
obs_start.format = 'iso'
obs_start = '-'.join(obs_start.value.split()[0].split('-')[::-1])+' '+obs_start.value.split()[1]

# --- Params ---
params = {  
    "simulator": {
        'double_precision': True,
        "use_gpus": True,
        'max_sources_per_chunk': 23000},
    "observation" : {
        "num_channels": 1,
        "start_frequency_hz": hdr['CRVAL3'],
        "phase_centre_ra_deg": RA.value,
        "phase_centre_dec_deg": DEC.value,
        "start_time_utc": obs_start, 
        "length": obs_length.value,
        "num_time_steps": 2160},
    "telescope": {

        "input_directory": "/project/c31/SKA_low_images/ska1low.tm"},
    "interferometer": {
        "oskar_vis_filename": "%s.vis" %(path_out+filename),
        "ms_filename": path_out+filename,
        "channel_bandwidth_hz": 1e5,
        "time_average_sec": 10.0}
}
np.save(path_out+'params.dict.npy', params)

settings = oskar.SettingsTree("oskar_sim_interferometer")
settings.from_dict(params)
settings["simulator/double_precision"] = False

# Run the simulation
sim = oskar.Interferometer(settings=settings)
sim.set_output_measurement_set(path_out+filename)
sim.set_sky_model(sky)
sim.run()

# Make Image
imager = oskar.Imager('single')
imager.set(fov_deg=FoV, image_size=N_1)
imager.set(input_file="%s.vis" %(path_out+filename), output_root=path_out+filename)
output = imager.run(return_images=1)
image = output["images"][0]

# Plot
im = plt.imshow(image, cmap="jet")
plt.gca().invert_yaxis()
plt.colorbar(im)
plt.savefig("%s.png" %imager.output_root)
plt.close("all")
