import casacore.tables as ct
from oskar.measurement_set import MeasurementSet
import oskar
import matplotlib.pyplot as plt
import numpy as np
from karabo.imaging.imager import Imager
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from astropy.io import fits
import bdsf
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.coordinates import Angle
import astropy.units as u
from datetime import datetime, timedelta
import tools21cm as t2c
import sys

matplotlib.use("TkAgg")
from astropy.wcs import WCS


def do_bdsf(imagepath, output_model, idx_f):
    img = fits.open(imagepath + "_I.fits")
    head = img[0].header
    data = np.abs(img[0].data[0])
    head["CDELT3"] = 1
    head["BMIN"] = 10 / 3600.0  # in Deg 2.8 arcsec for 120 km baseline
    head["BMAJ"] = 10 / 3600.0
    out_imagepath = imagepath + "_mod_I.fits"
    fits.writeto(out_imagepath, data, head, overwrite=True)
    img_bdsf = bdsf.process_image(
        out_imagepath,
        beam=(10 / 3600.0, 10 / 3600.0, 0),
        thresh_pix=3,
        thresh_isl=3,
        thresh="hard",
    )
    bdsf_source_list = "/scratch/snx3000/rsharma/bdsf_source_list.fits"
    img_bdsf.write_catalog(
        outfile=bdsf_source_list, catalog_type="gaul", format="fits", clobber=True
    )
    sl = fits.open(bdsf_source_list)
    xpix, ypix = sl[1].data["Xposn"], sl[1].data["Yposn"]
    xpix = xpix.astype(int)
    ypix = ypix.astype(int)
    RA = sl[1].data["RA"]
    DEC = sl[1].data["DEC"]
    sl_peak = data[ypix, xpix]  # sl[1].data['Peak_flux']
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
    np.savetxt(X=sky_data, fname=output_model)
    w = WCS(head)
    # xpix,ypix,_=w.wcs_world2pix(RA,DEC,0,0)
    return sky_data, data, xpix, ypix


def get_gleam(root_name):
    # add GLEAM point sources foreground
    root_name += "gleam"
    path_point = (
        "/store/ska/sk014/dataset_sdc3/inputs/dataLC_256_train_090523_test/frg/exgf/"
    )
    gleam_data = np.loadtxt(path_point + "rohit_sdc3cat_skymodel_4deg.txt")
    gleam_data[:, 2] = gleam_data[:, 2]
    # gleam_data = np.loadtxt(path_point+'gleamcat_skymodel_4deg.txt')
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


def write_ms(filename, ms_new, skadc_uvw, n):
    ms_new = ms_new.reshape(1440, n)
    uv_new = skadc_uvw.reshape(1440, n, 3)
    num_stations = 512
    num_channels = 1
    num_pols = 1
    ref_freq_hz = 1.8e8
    freq_inc_hz = 1.0e5
    num_baselines = n  # int(num_stations*(num_stations-1)*0.5)
    num_times = 1440  # 1440 number of time channels and 10 sec integration
    num_pols = 1
    num_channels = 1
    ms = oskar.MeasurementSet.create(
        filename, num_stations, num_channels, num_pols, ref_freq_hz, freq_inc_hz
    )
    ra_rad = 0.0
    dec_rad = -30.0
    exposure_sec = 1.0
    interval_sec = 1.0
    ms.set_phase_centre(ra_rad, dec_rad)
    mjd_20210921 = 59478.592071770836
    # UTC 2021-09-21T14:12:35.1
    # Modified Julian Day (MJD)   59478.592071770836
    # Julian Day (JD) 2459479.0920717707
    uu = np.zeros([num_baselines])
    vv = np.zeros_like(uu)
    ww = np.zeros_like(uu)
    vis = np.zeros([num_times, num_channels, num_baselines, num_pols], dtype="c8")
    for t in range(num_times):
        time_stamp = mjd_20210921 * 86400.0 + t
        uu[:] = uv_new[t, :, 0]
        vv[:] = uv_new[t, :, 1]
        ww[:] = uv_new[t, :, 2]
        for c in range(num_channels):
            for b in range(num_baselines):
                vis[t, c, b, :] = ms_new[t, b]
        start_row = t * num_baselines
        ms.write_coords(
            start_row, num_baselines, uu, vv, ww, exposure_sec, interval_sec, time_stamp
        )
        ms.write_vis(start_row, 0, num_channels, num_baselines, vis[t, ...])


def plot_bdsf(data_img, xpix, ypix):
    f, ax = plt.subplots(1, 1)
    im = ax.imshow(data_img, aspect="auto", origin="lower", vmin=-1.0e-3, vmax=1.0e-3)
    ax.plot(xpix, ypix, "o", color="red")
    f.colorbar(im)
    plt.show()


# -------------------- Simulate MS files--------------------------------
iter_ska = 0
if iter_ska:
    for i in range(295, 901):
        ii = "%04d" % i
        skadc = ct.table("/store/ska/sk01/sdc3-new/MS/ZW3_IFRQ_" + ii + ".ms")
        skadc_data = skadc.getcol("DATA")
        skadc_uvw = skadc.getcol("UVW")
        filename_ska = (
            "/scratch/snx3000/rsharma/ska_data_reduced/ZW3_IFRQ_for_karabo_"
            + ii
            + ".MS"
        )
        ms_new_ska = skadc_data[:, 0, 0]  # /np.nanmax(np.abs(skadc_data))
        write_ms(filename_ska, ms_new_ska, skadc_uvw, 130816)

write_ska_uvw = 0
if write_ska_uvw:
    for i in range(277, 901):
        ii = "%04d" % i
        skadc = ct.table("/store/ska/sk01/sdc3-new/MS/ZW3_IFRQ_" + ii + ".ms")
        skadc_u = skadc.getcol("UVW")[:, 0].reshape(1440, 130816)
        skadc_v = skadc.getcol("UVW")[:, 1].reshape(1440, 130816)
        skadc_w = skadc.getcol("UVW")[:, 2].reshape(1440, 130816)
        np.save(
            "/scratch/snx3000/rsharma/ska_data_reduced/SKADC_U_" + ii + ".npy", skadc_u
        )
        np.save(
            "/scratch/snx3000/rsharma/ska_data_reduced/SKADC_V_" + ii + ".npy", skadc_v
        )
        np.save(
            "/scratch/snx3000/rsharma/ska_data_reduced/SKADC_W_" + ii + ".npy", skadc_w
        )

skadc = ct.table("/scratch/snx3000/rsharma/pybdsf_tests/ZW3_IFRQ_for_karabo_0750.MS")
point = ct.table("/scratch/snx3000/rsharma/pybdsf_tests/GLEAM_point_sources_0750.MS")

skadc_data = skadc.getcol("DATA")
# point_data=point.getcol('DATA')
skadc_uvw = skadc.getcol("UVW")
# point_uvw = point.getcol("UVW")
# ms1=MeasurementSet.open('/nas08-data02/rohit/SDC3/ZW3_IFRQ_0750.ms',readonly=True)
# ms1_uvw=ms1.read_column(column='UVW',start_row=0,num_rows=188375040)
# ms1_point_data=ms1.read_column(column='DATA',start_row=0,num_rows=188375040)[:,0,0]
write_ska = 0
write_point = 0
write_res = 0
filename_ska = "/scratch/snx3000/rsharma/pybdsf_tests/ZW3_IFRQ_for_karabo_0750.MS"
ms_new_ska = skadc_data[:, 0, 0]  # /np.nanmax(np.abs(skadc_data))
print("SKA Image max:", ms_new_ska.max())
if write_ska:
    write_ms(filename_ska, ms_new_ska, skadc_uvw, 130816)
# -----------------------------------------
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
imager.update(uu, vv, ww, ska_vis[0, :, 0])
ska_image = imager.finalise(return_images=1)["images"][0]
print("SKA Image Maximum:", np.max(ska_image))

# filename_point='/scratch/snx3000/rsharma/subvis_tests/point_sdc3point_ch750_4h1d.MS'
# ms_new_point=point_data[:,0,0]/np.nanmax(np.abs(point_data))
# if(write_ska):
#    write_ms(filename_point,ms_new_point,skadc_uvw,130816)

# filename_res='/scratch/snx3000/rsharma/subvis_tests/residual_sdc3point_ch750_4h1d.MS'
# ms_new_res=ms_new_ska-ms_new_point
# if(write_res):
#    write_ms(filename_res,ms_new_res,skadc_uvw,130816)
# ------- Image read -----
idx_f = 750
freqs = 1.06e8 + np.arange(901) * 1.0e5  # Hz
z = t2c.nu_to_z(freqs[idx_f] * 1e-6)

# --- Make Point Image-----
imagepath = "/scratch/snx3000/rsharma/pybdsf_tests/image_ch750"
make_image = 0
if make_image:
    point_read = MeasurementSet.open(filename_point, readonly=True)
    num_baselines_ = int(point_read.num_stations * (point_read.num_stations - 1) / 2.0)
    uu, vv, ww = point_read.read_coords(start_row=0, num_baselines=num_baselines_)
    point_vis = point_read.read_vis(
        start_row=0, start_channel=0, num_channels=1, num_baselines=num_baselines_
    )
    imager = oskar.Imager()
    imager.set(output_root=imagepath)
    imager.fov_deg = 4  # 0.1 degrees across.
    imager.image_size = 2048  # 256 pixels across.
    imager.set_vis_frequency(ref_hz=180e6)  # 100 MHz, single channel data.
    imager.update(uu, vv, ww, point_vis[0, :, 0])
    point_image = imager.finalise(return_images=1)["images"][0]
do_point_source_det = 1
nk = 3
ms_new_ska_array = [0] * (nk + 1)
ms_new_ska_array[0] = ms_new_ska
uvlim = 6000000  # in mts
for k in range(nk):
    kk = "%02d" % k
    print("Iteration: " + kk)
    if (k > 5) & (do_point_source_det == 1):
        print("Doing pybdf.....")
        output_model = (
            "/scratch/snx3000/rsharma/pybdsf_tests/pybdsf_skymodel_" + kk + ".txt"
        )
        bdsf_data, data_img, xpix, ypix = do_bdsf(imagepath, output_model, idx_f)
        print(bdsf_data[:, 2])
        f, ax = plt.subplots(1, 1)
        im = ax.imshow(
            data_img, aspect="auto", origin="lower", vmin=-1.0e-3, vmax=1.0e-3
        )
        ax.plot(xpix, ypix, "o", color="red")
        f.colorbar(im)
        plt.show()
    # ---------------- Do Point Source Simulation-------------------------
    imagepath = "/scratch/snx3000/rsharma/pybdsf_tests/res_image_ch750_" + "%02d" % k
    path_out = "/scratch/snx3000/rsharma/pybdsf_tests/"
    root_name = "point_source_iteration"
    path_telescope = "/store/ska/sk014/dataset_sdc3/inputs/telescope.tm"
    run_name = path_out + root_name
    filename_point = run_name + "_" + kk + ".MS"
    iono_fits = (
        "/scratch/snx3000/mibianco/output_sdc3/dataLC_256_train_090523/atmo/screen_4h_i0_ch"
        + str(idx_f)
        + ".fits"
    )
    r0, sampling = 7e3, 100.0
    # inner_sky,outter_sky=get_gleam(root_name)
    sky = SkyModel()
    sky_test = np.zeros((1, 13))
    sky_test[0][0] = -1
    sky_test[0][1] = -31
    sky_test[0][2] = 1
    sky.add_point_sources(sky_test)
    # if(k==0):
    #       sky.add_point_sources(inner_sky);sky.add_point_sources(outter_sky)
    # else:
    #       sky.add_point_sources(bdsf_data)
    telescope = Telescope.read_from_file(path_telescope)
    t_start = datetime(
        2021, 9, 21, 14, 12, 40, 0
    )  # HA between -2h to +2h, obs start at '2021-09-21 14:12:40.1'
    t_obs = timedelta(hours=4, minutes=0, seconds=0, milliseconds=0)
    # t_obs = timedelta(hours=0, minutes=1, seconds=0, milliseconds=0)
    t_day = t_obs
    t_int = timedelta(seconds=10)
    nr_tsteps = int(t_day.total_seconds() / t_int.total_seconds())
    nr_days_obs = int(t_obs.total_seconds() / t_day.total_seconds())
    print(" Simulating %d days observation\n time steps: %d" % (nr_days_obs, nr_tsteps))
    observation_settings = Observation(
        phase_centre_ra_deg=0,
        mode="Tracking",
        phase_centre_dec_deg=-30,
        start_date_and_time=t_start,
        start_frequency_hz=freqs[idx_f],
        number_of_channels=1,
        number_of_time_steps=nr_tsteps,
        length=t_day,
    )
    simulation = InterferometerSimulation(
        ms_file_path=filename_point,
        vis_path=run_name + "_" + kk + ".vis",
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
    # if(os.path.isfile(run_name+'_'+kk+'.vis')==False):
    print("Simulating Point Visibilities.." + filename_point)
    visibilities = simulation.run_simulation(telescope, sky, observation_settings)
    # --------------- Subtractig the Point Source Vis
    point1 = ct.table(filename_point)
    point_data = point1.getcol("DATA")
    point_uvw = point1.getcol("UVW")
    ms_new_point = point_data[:, 0, 0]
    print("Point Source Max:", ms_new_point.max())
    ms_new_res = ms_new_ska_array[k] - ms_new_point
    filename_res = (
        "/scratch/snx3000/rsharma/pybdsf_tests/residual_sdc3point_ch750_" + kk + ".MS"
    )
    write_ms(filename_res, ms_new_res, skadc_uvw, 130816)
    ms_new_ska_array[k + 1] = ms_new_res
    imager = oskar.Imager()
    imager.fov_deg = 4  # 0.1 degrees across.
    imager.image_size = 2048  # 256 pixels across.
    imager.set_output_root(imagepath)
    imager.set_vis_frequency(180e6)  # 100 MHz, single channel data.
    imager.update(uu, vv, ww, ms_new_res)
    res_image = imager.finalise(return_images=1)["images"][0]
    get_uvlim = 0
    if get_uvlim:
        # ----------------- Residual small baselines-------------------
        uvdist = np.sqrt(uu**2 + vv**2)
        idx = np.where(uvdist > uvlim)[0]
        idx_low = np.where(uvdist < uvlim)[0]
        ms_new_res1 = ms_new_res.reshape(1440, 130816)[:, idx_low].flatten()
        res1_uvw = skadc_uvw.reshape(1440, 130816, 3)[:, idx_low, :].flatten()
        filename_res1 = (
            "/scratch/snx3000/rsharma/pybdsf_tests/residual_sdc3point_ch750_uvlim"
            + str(uvlim)
            + "_"
            + kk
            + ".MS"
        )
        # write_ms(filename_res1,ms_new_res1,res1_uvw,idx_low.shape[0])
        res_read1 = MeasurementSet.open(filename_res1, readonly=True)
        num_baselines_ = int(
            res_read1.num_stations * (res_read1.num_stations - 1) / 2.0
        )
        uu, vv, ww = res_read1.read_coords(start_row=0, num_baselines=num_baselines_)
        res_vis1 = res_read1.read_vis(
            start_row=0, start_channel=0, num_channels=1, num_baselines=num_baselines_
        )
        imager = oskar.Imager()
        imager.fov_deg = 4  # 0.1 degrees across.
        imager.image_size = 2048  # 256 pixels across.
        imager.set_output_root(imagepath)
        imager.set_vis_frequency(180e6)  # 100 MHz, single channel data.
        imager.update(uu, vv, ww, res_vis1[0, :, 0])
        res_image1 = imager.finalise(return_images=1)["images"][0]
        print("Residual Image Maximum (UV cut):", np.max(res_image1))

img0 = fits.open("res_image_ch750_00_I.fits")
d0 = img0[0].data[0]
img1 = fits.open("image_ch750_01_I.fits")
d1 = img1[0].data[0]
img2 = fits.open("image_ch750_02_I.fits")
d2 = img2[0].data[0]
f, (ax0, ax1, ax2) = plt.subplots(3, 1, sharex=True, sharey=True)
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
ax0.set_title("Residual0")
ax1.set_title("Residual1")
ax2.set_title("Residual2")
plt.show()


sys.exit()

# -------------------------------------------
print("Reading MS..")
# ----------------------------------------
point_read = MeasurementSet.open(filename_point, readonly=True)
num_baselines_ = int(point_read.num_stations * (point_read.num_stations - 1) / 2.0)
uu, vv, ww = point_read.read_coords(start_row=0, num_baselines=num_baselines_)
point_vis = point_read.read_vis(
    start_row=0, start_channel=0, num_channels=1, num_baselines=num_baselines_
)
imager = oskar.Imager()
imager.fov_deg = 4  # 0.1 degrees across.
imager.image_size = 2048  # 256 pixels across.
imager.set_vis_frequency(ref_hz=180e6)  # 100 MHz, single channel data.
imager.update(uu, vv, ww, point_vis[0, :, 0])
point_image = imager.finalise(return_images=1)["images"][0]
print("Point Source Image Maximum:", np.max(point_image))

matplotlib.use("TkAgg")
uvdist = np.sqrt(uu**2 + vv**2)
idx = np.where(uvdist > 1000)[0]
idx_low = np.where(uvdist < 2000)[0]
fact = point_vis[0, :, 0][idx].imag.std() / ska_vis[0, :, 0][idx].imag.std()

# ----------------------------------------
res_read = MeasurementSet.open(filename_res, readonly=True)
num_baselines_ = int(res_read.num_stations * (res_read.num_stations - 1) / 2.0)
uu, vv, ww = res_read.read_coords(start_row=0, num_baselines=num_baselines_)
res_vis = res_read.read_vis(
    start_row=0, start_channel=0, num_channels=1, num_baselines=num_baselines_
)
imager = oskar.Imager()
imager.fov_deg = 4  # 0.1 degrees across.
imager.image_size = 2048  # 256 pixels across.
imager.set_vis_frequency(180e6)  # 100 MHz, single channel data.
imager.update(uu, vv, ww, res_vis[0, :, 0])
res_image = imager.finalise(return_images=1)["images"][0]
print("Residual Image Maximum:", np.max(res_image))

# -----------------------------------------
res_read1 = MeasurementSet.open(filename_res1, readonly=True)
num_baselines_ = int(res_read1.num_stations * (res_read1.num_stations - 1) / 2.0)
uu, vv, ww = res_read1.read_coords(start_row=0, num_baselines=num_baselines_)
res_vis1 = res_read1.read_vis(
    start_row=0, start_channel=0, num_channels=1, num_baselines=num_baselines_
)
imager = oskar.Imager()
imager.fov_deg = 4  # 0.1 degrees across.
imager.image_size = 2048  # 256 pixels across.
imager.set_vis_frequency(180e6)  # 100 MHz, single channel data.
imager.update(uu, vv, ww, res_vis1[0, :, 0])
res_image1 = imager.finalise(return_images=1)["images"][0]
print("Residual Image Maximum (UV cut):", np.max(res_image1))


# ---------------------------------------------
sigma = 5
mask = res_image1 * 0
mask[np.where(res_image1 > sigma * np.nanstd(res_image1))] = np.nan
mask[np.where(res_image1 < -1 * sigma * np.nanstd(res_image1))] = np.nan
res_image_masked = mask + res_image1

f, (ax0, ax1, ax2) = plt.subplots(1, 3, sharex=True, sharey=True)
ax0.plot(uvdist, np.abs(ska_vis[0, :, 0]), ".")
ax1.plot(uvdist, np.abs(res_vis[0, :, 0]), ".")
ax2.plot(uvdist, np.abs(point_vis[0, :, 0]), ".")
ax0.set_xscale("log")
ax1.set_xscale("log")
ax2.set_xscale("log")
ax0.set_xlabel("UV distance (m)")
ax1.set_xlabel("UV distance (m)")
ax2.set_xlabel("UV distance (m)")
ax0.set_ylabel("Amplitude")
ax0.set_title("SKA")
ax1.set_title("Residual")
ax2.set_title("Point")
plt.show()

f, ax = plt.subplots(1, 1)
ax.plot(uu.real, vv.real, ".", color="red")
ax.plot(-uu, -vv, ".", color="red")
ax.set_xlabel("U (m)")
ax.set_ylabel("V (m)")
ax.set_title("SKA-low")
plt.show()

f, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, sharex=True, sharey=True)
ax0.imshow(
    ska_image,
    origin="lower",
    cmap="coolwarm",
    extent=[-2, 2, -32, -28],
    aspect="auto",
    vmin=-0.1,
    vmax=0.1,
)
ax0.set_title("SKA Data")
im1 = ax1.imshow(
    res_image1,
    origin="lower",
    cmap="coolwarm",
    extent=[-2, 2, -32, -28],
    aspect="auto",
    vmin=-0.1,
    vmax=0.1,
)
ax1.set_title("Residual")
ax2.imshow(
    point_image,
    origin="lower",
    cmap="coolwarm",
    extent=[-2, 2, -32, -28],
    aspect="auto",
    vmin=-0.1,
    vmax=0.1,
)
ax2.set_title("Point Source")
ax3.imshow(
    res_image_masked,
    origin="lower",
    extent=[-2, 2, -32, -28],
    cmap="coolwarm",
    aspect="auto",
    vmin=-0.1,
    vmax=0.1,
)
ax3.set_title("Residual Masked 5-sigma")
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
f.colorbar(im1, cax=cax, orientation="vertical")
ax3.set_xlabel("RA (deg)")
ax2.set_xlabel("RA (deg)")
ax2.set_ylabel("Dec (deg)")
ax0.set_ylabel("Dec (deg)")
plt.show()
