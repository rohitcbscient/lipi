########################### Program for point source subtraction from the SKA data ##################################

# Here we removed the point sources interatively using GLEAM catalog and pybdsf (See the for loop below)
# In first iteration, we used GLEAM catalog to remove the point sources
# Since not all point sources are subtracted just by GLEAM, we employ pybdsf to search for point sources in the image
# The point sources were converted into skymodel and subtrated from the residuals in each step (See the for loop below)
# Author: Rohit Sharma @ 2023 (rohitcbscient@gmail.com)

import casacore.tables as ct
from oskar.measurement_set import MeasurementSet
import oskar
import os
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from astropy.io import fits
import bdsf
from astropy.coordinates import Angle
import astropy.units as u
from datetime import datetime, timedelta
from astropy.wcs import WCS


def do_bdsf(imagepath, output_model, freqs, idx_f):
    """
    Function to do Pybdsf on the inout image and return the sky model in the form of a catalog
    ----------------------
    Inputs:
    1. imagepath: input fitsfile image
    2. output_model: output skymodel filename
    3. freqs: freqeuncy array
    4. idx_f: Index of the channel number
    -----------------------
    Outputs:
    1. sky_data: Karabo/oskar sky model array
    2. data: Data for the point sources
    3. xpix: X pixel position of the sources
    4. ypix: Y pixel position of the sources
    """
    img = fits.open(imagepath + "_I.fits")
    head = img[0].header
    data = np.abs(img[0].data[0])
    head["CDELT3"] = 1
    uvw_max = 73508.6  # max
    b_deg = (3.0e8 / freqs[idx_f]) / uvw_max * 180 / np.pi
    head["BMIN"] = b_deg  # in Deg 2.8 arcsec for 120 km baseline
    head["BMAJ"] = b_deg
    out_imagepath = imagepath + "_mod_I.fits"
    fits.writeto(out_imagepath, data, head, overwrite=True)
    # ------ PYBDSF Begins here
    img_bdsf = bdsf.process_image(
        out_imagepath,
        beam=(10 / 3600.0, 10 / 3600.0, 0),
        thresh_pix=5,
        thresh_isl=6,
        thresh="hard",
    )
    bdsf_source_list = (
        "/scratch/snx3000/rsharma/bdsf_source_list.fits"  # Fits catalog from the pybdsf
    )
    # Writing the fits source catalog
    img_bdsf.write_catalog(
        outfile=bdsf_source_list, catalog_type="gaul", format="fits", clobber=True
    )
    # Reading the catalog and putting into karabo/oskar skymodel array format
    sl = fits.open(bdsf_source_list)
    xpix, ypix = sl[1].data["Xposn"], sl[1].data["Yposn"]
    xpix = xpix.astype(int)
    ypix = ypix.astype(int)
    RA = sl[1].data["RA"]
    DEC = sl[1].data["DEC"]
    sl_peak = data[ypix, xpix]
    RA[np.where(RA > 200)] = RA[np.where(RA > 200)] - 360
    sl_bmaj = sl[1].data["Maj"]
    sl_bmin = sl[1].data["Min"]
    sl_bpa = sl[1].data["PA"]
    gauss_n_size = RA.shape[0]
    sky_data = np.hstack(
        (
            RA.reshape(gauss_n_size, 1),
            DEC.reshape(gauss_n_size, 1),
            sl_peak.reshape(gauss_n_size, 1),
            np.zeros((gauss_n_size, 9)),
        )
    )
    sky_data[:, 9] = sl_bmaj
    sky_data[:, 10] = sl_bmin
    sky_data[:, 11] = sl_bpa
    # Saving the skymodel below
    np.savetxt(X=sky_data, fname=output_model)
    return sky_data, data, xpix, ypix


def get_gleam(root_name):
    """
    Returns the GLEAM skymodel inner and outer sky catalog array from the ascii text file
    """
    # add GLEAM point sources foreground
    root_name += "gleam"
    path_point = (
        "/store/ska/sk014/dataset_sdc3/inputs/dataLC_256_train_090523_test/frg/exgf/"
    )
    gleam_data = np.loadtxt(path_point + "rohit_sdc3cat_skymodel_4deg.txt")
    gleam_data[:, 2] = gleam_data[:, 2]
    # create inner & outter sky
    ra_wrap, dec_wrap = (
        Angle([gleam_data[:, 0], gleam_data[:, 1]], unit="degree")
        .wrap_at(180 * u.deg)
        .deg
    )
    # select inner sky sources
    inner_mask = np.sqrt(ra_wrap**2 + (dec_wrap + 30) ** 2) <= 2
    inner_sky = gleam_data[inner_mask]
    # select outter sky sources
    outter_mask = np.sqrt(ra_wrap**2 + (dec_wrap + 30) ** 2) > 2
    outter_sky = gleam_data[outter_mask]
    outter_sky[:, 2] *= 1e-3
    return inner_sky, outter_sky


def write_ms(filename, ms_new, skadc_uvw, n, ref_freq_hz):
    """
    Write the measurement sets for given UVW data and corresponding visibilities
    Input:
    1. filename: MS output filename
    2. ms_new: visibilities to be planted
    3. skadc_uvw: UVW array
    4. ref_freq_hz: Reference frequency in Hz
    """
    # Below are the specifics operation related to ska data
    ms_new = ms_new.reshape(1440, n)
    uv_new = skadc_uvw.reshape(1440, n, 3)
    num_stations = 512
    num_channels = 1
    num_pols = 1
    freq_inc_hz = 1.0e5
    num_baselines = int(num_stations * (num_stations - 1) * 0.5)
    num_times = 1440  # number of time channels and 10 sec integration
    num_pols = 1
    num_channels = 1
    # Create an empty MS---------------------------------------------------------------------
    ms = oskar.MeasurementSet.create(
        filename, num_stations, num_channels, num_pols, ref_freq_hz, freq_inc_hz
    )
    ra_rad = 0.0
    dec_rad = -30.0
    exposure_sec = 1.0
    interval_sec = 1.0
    ms.set_phase_centre(ra_rad, dec_rad)
    mjd_20210921 = 59478.592071770836
    uu = np.zeros([num_baselines])
    vv = np.zeros_like(uu)
    ww = np.zeros_like(uu)
    vis = np.zeros([num_times, num_channels, num_baselines, num_pols], dtype="c8")
    # Write the visibilities row by row in times -----------------------------------------------
    for t in range(num_times):
        time_stamp = mjd_20210921 * 86400.0 + t
        uu[:] = uv_new[t, :, 0]
        vv[:] = uv_new[t, :, 1]
        ww[:] = uv_new[t, :, 2]
        # uu[:] = uv_new[:,0];vv[:] = uv_new[:,1];ww[:] = uv_new[:,2]
        for c in range(num_channels):
            for b in range(num_baselines):
                vis[t, c, b, :] = ms_new[t, b]
                # vis[t, c, b, :] = ms_new[b]
        start_row = t * num_baselines
        ms.write_coords(
            start_row, num_baselines, uu, vv, ww, exposure_sec, interval_sec, time_stamp
        )
        ms.write_vis(start_row, 0, num_channels, num_baselines, vis[t, ...])


def plot_bdsf(data_img, xpix, ypix):
    """
    Plotting the pybdsf results
    Inputs:
    1. data_img: Image array
    2. xpix, ypix: X and Y position of the point sources
    """
    f, ax = plt.subplots(1, 1)
    im = ax.imshow(data_img, aspect="auto", origin="lower", vmin=-1.0e-3, vmax=1.0e-3)
    ax.plot(xpix, ypix, "o", color="red")
    f.colorbar(im)
    plt.show()


# ------------------------------------------- Main Program --------
i = 753  # Choose the channel number for the simulation
ilist = [i]
# i = int(os.environ["SLURM_ARRAY_TASK_ID"]) # Argument for the SLURM task on CSCS
ii = "%04d" % i
freqs = 1.06e8 + np.arange(901) * 1.0e5  # Hz
# os.system('rm -rf /scratch/snx3000/rsharma/pybdsf_tests/*') #  For entire clean up
# --------------------------------------------------------------------
# Define SKA filename and path
filename_ska = (
    "/scratch/snx3000/rsharma/ska_data_reduced/ZW3_IFRQ_for_karabo_" + ii + ".MS"
)
# Define Point Source filename and path
filename_point0 = (
    "/scratch/snx3000/rsharma/gleam_point_source_ms/GLEAM_point_sources_" + ii + ".MS"
)
# ---- Reading the visibilities and UVW coordinates from SKA and point sources
skadc = ct.table(filename_ska)
point = ct.table(filename_point0)
skadc_data = skadc.getcol("DATA")
skadc_uvw = skadc.getcol("UVW")
ms_new_ska = skadc_data[:, 0, 0]

# Make SKA  image for sanity checks
make_ska_image = 0
ska_imagepath = "/scratch/snx3000/rsharma/residual1_pybdsf/ska_images/ska_image_ch" + ii
if make_ska_image:
    ska_read = MeasurementSet.open(filename_ska, readonly=True)
    num_baselines_ = int(ska_read.num_stations * (ska_read.num_stations - 1) / 2.0)
    uu, vv, ww = ska_read.read_coords(start_row=0, num_baselines=num_baselines_)
    ska_vis = ska_read.read_vis(
        start_row=0, start_channel=0, num_channels=1, num_baselines=num_baselines_
    )
    imager = oskar.Imager()
    imager.fov_deg = 4  # 0.1 degrees across.
    imager.image_size = 2048  # 256 pixels across.
    imager.set_vis_frequency(180e6)  # 100 MHz, single channel data.
    imager.set(output_root=ska_imagepath)
    imager.update(uu, vv, ww, ska_vis[0, :, 0])
    ska_image = imager.finalise(return_images=1)["images"][0]
    print("SKA Image Maximum:", np.max(ska_image))

# ----------------------------------- Program starts -----------------------------------------------
# Make the full frequency array
freqs = 1.06e8 + np.arange(901) * 1.0e5  # Hz
imagepath_point = "/scratch/snx3000/rsharma/residual1_pybdsf/point_ch" + ii
# ``do_point_source_det" allows you simulate the GLEAM point in case you want to testing and sanity purposes
do_point_source_det = 1
# nk is the total number of iterations / iteration 0 is from the GLEAM catalog, while iterations >0 are pybdsf point model substractions
nk = 5
# Array initialisation
ms_new_ska_array = [0] * (nk + 1)
ms_new_ska_array[
    0
] = ms_new_ska  # Since we want to subtract everything from SKA dattttta, the index 0 is the SKA data
uvlim = 1.0e12  # in mts # Upper limit on UVdist in case you want to limit the UV range, we give it very large value to make use of all UV points

# ------- Loop for the subtraction starts -------------------------------------------------------
for k in range(nk):
    kk = "%02d" % k
    print("Iteration: " + kk)
    # Below in if statement we simulate the pybdsf sky model
    if (k > 0) & (do_point_source_det == 1):
        print("Doing pybdf.....")
        output_model = (
            "/scratch/snx3000/rsharma/residual1_pybdsf/pybdsf_skymodel/pybdsf_skymodel_ch"
            + ii
            + "_"
            + kk
            + ".txt"
        )  # Output_model is the path of the pybdsf skymodel
        bdsf_data, data_img, xpix, ypix = do_bdsf(
            imagepath, output_model, freqs, i
        )  # Note that imagepath will be defined for k>0 iterations
        check_image = 0
        if check_image:
            print(bdsf_data[:, 2])
            f, ax = plt.subplots(1, 1)
            im = ax.imshow(
                data_img, aspect="auto", origin="lower", vmin=-1.0e-3, vmax=1.0e-3
            )
            ax.plot(xpix, ypix, "o", color="red")
            f.colorbar(im)
            plt.show()
    # ---------------- Do Point Source Simulation-------------------------
    imagepath = (
        "/scratch/snx3000/rsharma/residual1_pybdsf/res_images/res_images_ch"
        + ii
        + "_"
        + "%02d" % k
    )  # Path for the output images iteration k and channel i
    path_out = (
        "/scratch/snx3000/rsharma/residual1_pybdsf/res_images/"  # Output Image path
    )
    root_name = "point_source_iteration_ch" + ii  # Point source filenames
    path_telescope = "/store/ska/sk014/dataset_sdc3/inputs/telescope.tm"  # Path to the telecope directory
    run_name = path_out + root_name  # Path for output point sources
    filename_point = run_name + "_" + kk + ".MS"  # Same as above
    iono_fits = (
        "/scratch/snx3000/rsharma/atmo/screen_4h_i0_ch" + "%03d" % i + ".fits"
    )  # Path for the ionosphere file for channel i
    r0, sampling = 7e3, 100.0  # Ionosphere parameters
    # ----------------Sky Model Defination ------------------------------------
    inner_sky, outter_sky = get_gleam(
        root_name
    )  # Path for the skymodel from GLEAM catalog
    sky = SkyModel()
    # Below the choice for the 0th iteration is made for the simulation of the GLEAM catalog otherwise it will pybdsf sky model
    if k == 0:
        sky.add_point_sources(inner_sky)
        sky.add_point_sources(outter_sky)
    else:
        sky.add_point_sources(bdsf_data)
    # ------ Telescope Defination -----------------------------------------------
    telescope = Telescope.read_from_file(path_telescope)
    t_start = datetime(
        2021, 9, 21, 14, 12, 40, 0
    )  # HA between -2h to +2h, obs start at '2021-09-21 14:12:40.1'
    t_obs = timedelta(hours=4, minutes=0, seconds=0, milliseconds=0)
    # t_obs = timedelta(hours=0, minutes=1, seconds=0, milliseconds=0)
    # ---------Observation Defination ---------------------------------------------
    observation_settings = Observation(
        phase_centre_ra_deg=0,
        mode="Tracking",
        phase_centre_dec_deg=-30,
        start_date_and_time=t_start,
        start_frequency_hz=freqs[i],
        number_of_channels=1,
        number_of_time_steps=1440,
        length=t_obs,
    )
    # ----------- Interferometer Defination ---------------------------------------
    simulation = InterferometerSimulation(
        ms_file_path=filename_point,
        use_gpus=True,
        use_dask=False,
        channel_bandwidth_hz=1e5,
        enable_numerical_beam=False,
        enable_array_beam=False,
        noise_enable=False,
        ionosphere_fits_path=iono_fits,
        ionosphere_screen_type="External",
        ionosphere_screen_height_km=r0,
        ionosphere_screen_pixel_size_m=sampling,
        ionosphere_isoplanatic_screen=True,
    )
    print("Simulating Point Visibilities.." + filename_point)
    visibilities = simulation.run_simulation(telescope, sky, observation_settings)
    # --------------- Subtractig the Point Source Vis----------------------------
    point1 = ct.table(filename_point)
    point_data = point1.getcol("DATA")
    point_uvw = point1.getcol("UVW")
    ms_new_point = point_data[:, 0, 0]
    ms_new_res = ms_new_ska_array[k] - ms_new_point
    print("Point Source Max:", ms_new_point.max(), " Residual Shape:", ms_new_res.shape)
    filename_res = (
        "/scratch/snx3000/rsharma/residual1_pybdsf/ms/residual_sdc3point_ch"
        + ii
        + "_"
        + kk
        + ".MS"
    )
    write_ms(filename_res, ms_new_res, skadc_uvw, 130816, freqs[i])
    ms_new_ska_array[k + 1] = ms_new_res
    # ------- Image Residuals ----------------------------------------------------
    imager = oskar.Imager()
    imager.fov_deg = 4  # 0.1 degrees across.
    imager.image_size = 2048  # 256 pixels across.
    imager.set_output_root(imagepath)
    imager.set_vis_frequency(freqs[i])  # 100 MHz, single channel data.
    imager.update(uu, vv, ww, ms_new_res)
    res_image = imager.finalise(return_images=1)["images"][0]
    # ------- Point Source Image -------------------------------------------------
    imager = oskar.Imager()
    imager.fov_deg = 4  # 0.1 degrees across.
    imager.image_size = 2048  # 256 pixels across.
    imager.set_output_root(imagepath_point)
    imager.set_vis_frequency(freqs[i])  # 100 MHz, single channel data.
    imager.update(uu, vv, ww, ms_new_point)
    res_image = imager.finalise(return_images=1)["images"][0]
    os.system(
        "rm -rf " + filename_point
    )  # Delete the point source MS to save the disk space
os.system("rm -rf /users/rsharma/core.*")  # Clear the home memory area

# ----------- Visualisation test below -------------------------------------------
plot_test = 0
ch = "801"
if plot_test:
    ska_img = fits.open("../ska_images/ska_image_ch0" + ch + "_I.fits")
    ska_image = ska_img[0].data[0]
    point_img = fits.open("../ska_images/ska_image_ch0" + ch + "_I.fits")
    point_image = point_img[0].data[0]
    img0 = fits.open("res_image_ch0" + ch + "_00_I.fits")
    d0 = img0[0].data[0]
    img1 = fits.open("res_image_ch0" + ch + "_02_I.fits")
    d1 = img1[0].data[0]
    img2 = fits.open("res_image_ch0" + ch + "_04_I.fits")
    d2 = img2[0].data[0]
    ska_std = ska_image[np.where(ska_image < 5 * np.std(ska_image))].std()
    d0 = d0 / np.std(d0) * ska_image[np.where(ska_image < 5 * np.std(ska_image))].std()
    d1 = d1 / np.std(d1) * ska_image[np.where(ska_image < 5 * np.std(ska_image))].std()
    d2 = d2 / np.std(d2) * ska_image[np.where(ska_image < 5 * np.std(ska_image))].std()
    f, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, sharex=True, sharey=True)
    ax0.imshow(
        d0,
        origin="lower",
        cmap="coolwarm",
        extent=[-2, 2, -32, -28],
        aspect="auto",
        vmin=-0.1,
        vmax=0.1,
    )
    ax0.set_xlabel("RA (deg)")
    ax1.set_xlabel("RA (deg)")
    ax0.set_ylabel("Dec (deg)")
    ax1.set_ylabel("Dec (deg)")
    ax1.imshow(
        d1,
        origin="lower",
        cmap="coolwarm",
        extent=[-2, 2, -32, -28],
        aspect="auto",
        vmin=-0.1,
        vmax=0.1,
    )
    ax2.imshow(
        d2,
        origin="lower",
        cmap="coolwarm",
        extent=[-2, 2, -32, -28],
        aspect="auto",
        vmin=-0.1,
        vmax=0.1,
    )
    ax3.imshow(
        ska_image,
        origin="lower",
        cmap="coolwarm",
        extent=[-2, 2, -32, -28],
        aspect="auto",
        vmin=-0.1,
        vmax=0.1,
    )
    ax0.set_title("Residual0")
    ax1.set_title("Residual2")
    ax2.set_title("Residual4")
    ax3.set_title("SKA Image")
    plt.show()
