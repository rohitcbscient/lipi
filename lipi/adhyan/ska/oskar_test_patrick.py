import matplotlib.pyplot as plt
import numpy as np
import oskar



#ra dec flux specindex shape majoraxis minoraxis PA

sky_data = np.array([
        [150.24433333334338, 2.4918888888888913, 0.00014356435643564358, None, None, None, 100000000.0, 0.0, None, 0.32960396039603956, 0.20267326732673266, 52.67326732673268],
        [150.24822222223213, 2.4932777777777804, 0.0003466336633663367, None, None, None, 100000000.0, 0.0, None, 0.6381188118811881, 0.34405940594059403, 348.33483348334835]]
)

params = {
    "simulator": {
        "use_gpus": False
    },
    "observation" : {
        "num_channels": 5,
        "start_frequency_hz": 100e6,
        "frequency_inc_hz": 20e6,
        "phase_centre_ra_deg": 150,
        "phase_centre_dec_deg": 2,
        "num_time_steps": 24,
        "start_time_utc": "01-01-2000 12:00:00.000",
        "length": "1:00:00.000"
    },
    "telescope": {
        "input_directory": "/home/rohit/simulations/meerKat/telescope.tm"
    },
    "interferometer": {
        "oskar_vis_filename": "/home/rohit/simulations/meerKat/example_pat.vis",
        "ms_filename": "",
        "channel_bandwidth_hz": 1e6,
        "time_average_sec": 10
    }
}
settings = oskar.SettingsTree("oskar_sim_interferometer")
settings.from_dict(params)

# Set the numerical precision to use.
precision = "single"
if precision == "single":
    settings["simulator/double_precision"] = False

# Create a sky model containing three sources from a numpy array.
#sky_data = np.array([
#        [20.0, -30.0, 1, 0, 0, 0, 100.0e6, -0.7, 0.0, 0,   0,   0],
#        [20.0, -30.5, 3, 2, 2, 0, 100.0e6, -0.7, 0.0, 600, 50,  45],
#        [20.5, -30.5, 3, 0, 0, 2, 100.0e6, -0.7, 0.0, 700, 10, -10]])

print(sky_data.shape)


sky = oskar.Sky.from_array(sky_data, precision)  # Pass precision here.

# Set the sky model and run the simulation.
sim = oskar.Interferometer(settings=settings)
sim.set_sky_model(sky)
sim.run()

# Make an image 4 degrees across and return it to Python.
# (It will also be saved with the filename "example_I.fits".)
imager = oskar.Imager(precision)
imager.set(fov_deg=4, image_size=2048)
imager.set(input_file="/home/rohit/simulations/meerKat/example_pat.vis", output_root="/home/rohit/simulations/meerKat/example")
#imager.set(image_type='PSF') #to get the psf
output = imager.run(return_images=1)
image = output["images"][0] #this is the dirty image


# Render the image using matplotlib and save it as a PNG file.
im = plt.imshow(image, cmap="jet")
plt.gca().invert_yaxis()
plt.colorbar(im)
plt.show()
