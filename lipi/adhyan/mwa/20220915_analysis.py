import pickle
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from sunpy.map import Map
import itertools
import astropy.units as u
from scipy.ndimage.interpolation import rotate
from surya.utils import model as mdl
from surya.utils import main as ut
from scipy.io import readsav
import pfsspy
from pfsspy import tracing
from astropy.coordinates import SkyCoord
import astropy.constants as const
import astropy.units as u
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime


baseline = ["000-002", "000-005", "000-007", "002-005", "002-007", "005-007"]
chlist = [
    61,
    65,
    69,
    73,
    81,
    86,
    91,
    96,
    101,
    113,
    120,
    127,
    134,
    142,
    150,
    158,
    167,
    187,
    226,
]
flux = [0] * len(chlist)
freq = [0] * len(chlist)
dflux = [0] * len(chlist)
i = 0
for ch in chlist:
    f2 = pickle.load(
        open("/data/20220915_MWA/pickle/flux_V1_20220915_ch" + str(ch) + ".p", "rb"),
        encoding="latin1",
    )
    flux[i] = np.mean(f2[17][3][0], axis=0)
    freq[i] = ch * 1.28
    dflux[i] = flux[i] - flux[i][0]
    i = i + 1
freq = np.array(freq)
flux = np.array(flux)
dflux = np.array(dflux)

# ------- Newkirk Model
height_Mm = [0] * len(freq)
height_Mm_harmonic = [0] * len(freq)
for i in range(len(freq)):
    height_Mm_harmonic[i] = (float(mdl.nk_freq2r(freq[i] / 2.0, 1)[0]) - 1) * 6.99e2
    if freq[i] < 240:
        height_Mm[i] = (float(mdl.nk_freq2r(freq[i], 1)[0]) - 1) * 6.99e2
    else:
        height_Mm[i] = 15  # 15 Mm threshold


f, ax = plt.subplots(1, 1)
im = ax.imshow(
    flux,
    aspect="auto",
    origin="lower",
    interpolation=None,
    vmin=0,
    vmax=40,
    cmap="jet",
)
# f.colorbar(label='NCCF')
ax.set_xticks([0, 50, 100, 150, 200, 250, 300, 350])
ax.set_xticklabels(
    [
        "03:29:58",
        "03:38:18",
        "03:46:38",
        "03:54:58",
        "04:03:18",
        "04:11:38",
        "04:19:58",
        "04:28:18",
    ]
)
ax.set_yticks(np.arange(len(chlist)))
ax.set_yticklabels(np.round(freq, 2))
ax.set_xlabel("Time (2022-09-15 HH:MM:SS UT)")
ax.set_ylabel("Frequency (MHz)")
plt.colorbar(im, ax=ax, label="Flux Density(SFU)")
plt.show()

f, ax = plt.subplots(1, 1)
ax.plot(flux[0:5].mean(axis=0), label="78-110 MHz")
ax.plot(flux[5:10].mean(axis=0), label="110-153 MHz")
ax.plot(flux[10:15].mean(axis=0), label="153-202 MHz")
ax.plot(flux[15:19].mean(axis=0), label="202-240 MHz")
ax.set_xticks([0, 50, 100, 150, 200, 250, 300, 350])
ax.set_xticklabels(
    [
        "03:29:58",
        "03:38:18",
        "03:46:38",
        "03:54:58",
        "04:03:18",
        "04:11:38",
        "04:19:58",
        "04:28:18",
    ]
)
ax.set_xlabel("Time (2022-09-15 HH:MM:SS UT)")
ax.set_ylabel("Flux Density (SFU)")
ax.legend()
plt.show()

f, ax = plt.subplots(1, 1)
ax.plot(freq[:-2], flux.mean(axis=1)[:-2], "o-", label="Total mean")
ax.plot(freq[:-2], flux[:, 268:275].mean(axis=1)[:-2], "o-", label="Radio burst")
ax.plot(
    freq[:-2],
    dflux[:, 268:275].mean(axis=1)[:-2],
    "o-",
    label="Radio burst (subtracted)",
)
ax.legend()
ax.set_xlabel("Frequency (MHz)")
ax.set_ylabel("Flux Density (SFU)")
plt.show()

write_data = 0
if write_data:
    pa = 24
    img = [0] * len(chlist)
    i = 0
    freq_mwa = [0] * len(chlist)
    Tb = [0] * len(chlist)
    S = [0] * len(chlist)
    bmin = [0] * len(chlist)
    bmaj = [0] * len(chlist)
    for c in chlist:
        print("Channel: ", c)
        imglist = sorted(
            glob.glob(
                "/sdata/20220915_MWA/images/20220915_ch"
                + str(c)
                + ".pol.I.time.*.image_image.FITS"
            )
        )[0:380]
        j = 0
        img[i] = [0] * len(imglist)
        Tb[i] = [0] * len(imglist)
        sumdata = [0] * len(imglist)
        S[i] = [0] * len(imglist)
        for im in imglist:
            aa = fits.open(im)
            data = aa[0].data[0][0][924:1124, 924:1124]
            std = np.nanstd(aa[0].data[0][0])
            data[np.where(data < 5 * std)] = np.nan
            img[i][j] = data
            sumdata[j] = np.nansum(data)
            img_ = img[i][j]
            img_[np.isnan(img_)] = 0
            padX = [img_.shape[0] - 99, 99]
            padY = [img_.shape[1] - 99, 99]
            imgP = np.pad(img_, [padY, padX], "constant")
            img[i][j] = rotate(imgP, angle=pa, reshape=False, mode="constant")
            img[i][j] = img[i][j][padY[0] : -padY[1], padX[0] : -padX[1]]
            j = j + 1
        bmin[i] = aa[0].header["BMIN"] * 3600
        bmaj[i] = aa[0].header["BMAJ"] * 3600
        sumdata = np.array(sumdata)
        minid = np.where(sumdata == np.nanmin(sumdata))[0][0]
        fact = flux[i] / np.nansum(img[i][0]) * (freq[i] / 108.0) ** 2  # Flux Setting
        for k in range(len(imglist)):
            S[i][k] = img[i][k] * fact[k]
            Tb[i][k] = (
                1224 * S[i][k] * 1.0e7 / (freq[i] ** 2 * 1.0e-6 * bmin[i] * bmaj[i])
            )
            k = k + 1
        freq_mwa[i] = np.round(aa[0].header["CRVAL3"] / 1.0e6, 1)
        i = i + 1
    Tb = np.array(Tb)
    S = np.array(S)
    img = np.array(img)
    Tbmax = np.nanmax(Tb, axis=(2, 3))
    freq_mwa = np.array(freq_mwa)
    bmin = np.array(bmin)
    bmaj = np.array(bmaj)
    pickle.dump(
        [freq_mwa, Tb, Tbmax, S, bmin, bmaj, img],
        open("/sdata/20220915_MWA/Tb_20220915.p", "wb"),
    )


chan = np.arange(14)
freq_mwa, Tb, Tbmax, S, bmin, bmaj, img = pickle.load(
    open("/data/20220915_MWA/Tb_20220915.p", "rb")
)
Tb[:, :, 0:70, :] = np.nan
Tb[:, :, 140:, :] = np.nan
Tb[:, :, :, 0:70] = np.nan
Tb[:, :, :, 140:] = np.nan
# Tb[np.where(Tb<1.e3)]=np.nan#;Tb[np.where(Tb>1.e7)]=np.nan
# region1=img[:,:,105:130,105:130];region1_max=np.nanmax(region1,axis=(2,3))
# region1=img[:,:,105:120,110:135];region1_max=np.nanmax(region1,axis=(2,3))
region1_10s = img[:, :, 115:124, 115:124]
region1_max_10s = np.nanmax(region1_10s, axis=(2, 3))
Tb_region1_10s = Tb[:, :, 115:124, 115:124]
Tb_region1_max_10s = np.nanmax(Tb_region1_10s, axis=(2, 3))
Tb_region1_mean_10s = np.nanmean(Tb_region1_10s, axis=(2, 3))
region2_10s = img[:, :, 95:112, 115:128]
region2_max_10s = np.nanmax(region2_10s, axis=(2, 3))
Tb_region2_10s = Tb[:, :, 95:112, 115:128]
Tb_region2_max_10s = np.nanmax(Tb_region2_10s, axis=(2, 3))
Tb_region2_mean_10s = np.nanmean(Tb_region2_10s, axis=(2, 3))
Tb_region3_10s = Tb[:, :, 86:90, 101:103]
Tb_region3_max_10s = np.nanmax(Tb_region3_10s, axis=(2, 3))
Tb_region3_mean_10s = np.nanmean(Tb_region3_10s, axis=(2, 3))
tmwa_10s = 12608 + np.arange(Tb_region1_10s.shape[1]) * 10

Tb_region1_10s_full = np.zeros(Tb.shape)
Tb_region1_10s_full[:, :, 115:124, 115:124] = Tb[:, :, 115:124, 115:124]

Tb_region2_10s_full = np.zeros(Tb.shape)
Tb_region2_10s_full[:, :, 95:112, 115:128] = Tb[:, :, 95:112, 115:128]


# region1=img;region1_max=np.nanmean(region1,axis=(2,3))
# region1_mean=np.nanmean(region1,axis=(2,3));Tb_region1_mean=np.nanmean(Tb_region1,axis=(2,3))

freq_mwa = np.array(freq_mwa)
# for i in range(len(chlist)):
#    region1_mean[i]=region1_mean[i]-region1_mean[i].min()

Tb_pixel0 = Tb[:, :, 98, 120]
Tb_pixel0[Tb_pixel0 < 5.0e3] = np.nan
f, ax = plt.subplots(1, 1)
j = [34, 36, 37, 38, 40, 41, 42]
cl = ["g", "r", "b", "k", "orange", "brown", "magenta"]
al = [0.5, 1.0, 1.0, 0.5, 0.5, 1.0, 0.5]
ll = [
    "03:35:40",
    "03:36:00",
    "03:36:10",
    "03:36:20",
    "03:36:40",
    "03:36:50",
    "03:37:00",
]
for i in range(7):
    ax.plot(
        freq, Tb_pixel0[:, j[i]] / 1.0e6, "o-", color=cl[i], alpha=al[i], label=ll[i]
    )
ax.set_ylabel("$T_B$ (MK)")
ax.set_xlabel("Frequency (MHz)")
ax.legend(loc=2)
ax.set_ylim(0, 10)
plt.show()

# ---------------------------- new

# chlist_img=[61,65,69,77,81,86,91,96,101,107,113,120,127,134,142,150,167,177,187,210,226]
dump_mwa_Tb = 0
if dump_mwa_Tb:
    chlist_img = [
        61,
        65,
        69,
        81,
        86,
        91,
        96,
        101,
        107,
        120,
        127,
        134,
        142,
        150,
        167,
        177,
        187,
    ]
    Tbmax_r1 = [0] * len(chlist_img)
    Tbmax_r2 = [0] * len(chlist_img)
    Smax_r1 = [0] * len(chlist_img)
    Smax_r2 = [0] * len(chlist_img)
    bmin = [0] * len(chlist_img)
    bmaj = [0] * len(chlist_img)
    freq_mwa = [0] * len(chlist_img)
    Tb_mean_region1 = [0] * len(chlist_img)
    Tb_mean_region2 = [0] * len(chlist_img)
    Tb_r2_1 = [0] * len(chlist_img)
    Tb_r1_1 = [0] * len(chlist_img)
    xc_max_reg1 = [0] * len(chlist_img)
    xc_max_reg2 = [0] * len(chlist_img)
    yc_max_reg1 = [0] * len(chlist_img)
    yc_max_reg2 = [0] * len(chlist_img)
    idx_r1 = [0] * len(chlist_img)
    idx_r2 = [0] * len(chlist_img)
    reg = "reg1"
    for i in range(len(chlist_img)):
        c = chlist_img[i]
        print("Channel: " + str(c))
        freq_mwa[i], Tb, Tbmax_r1[i], S, bmin[i], bmaj[i], img = pickle.load(
            open(
                "/data/20220915_MWA/Tb_20220915_" + str(reg) + "_ch" + str(c) + ".p",
                "rb",
            ),
            encoding="latin1",
        )
        print(Tb.shape)
        Smax_r1[i] = np.nanmax(S, axis=(1, 2))
        Tb_r1 = Tb[:, 115:128, 95:122]
        # Tb_r1 = Tb[:,119,116]
        Tb_r1_1[i] = Tb_r1 * 0
        xc_max_reg1[i] = [0] * len(Tb_r1)
        yc_max_reg1[i] = [0] * len(Tb_r1)
        idx_r1[i] = [0] * len(Tb_r1)
        for k in range(len(Tb_r1)):
            if (
                any(
                    map(
                        lambda x: x == 0,
                        np.array(np.where(Tb_r1[k] == np.nanmax(Tb_r1[k]))).flatten(),
                    )
                )
                == False
            ):
                Tb_r1_1[i][k] = Tb_r1[k]
                xc_max_reg1[i][k], yc_max_reg1[i][k] = np.where(
                    Tb[k][115:128, 95:122] == np.nanmax(Tb[k][115:128, 95:122])
                )
            else:
                Tbmax_r1[i][k] = np.nan
                Tb_r1_1[i][k] = np.nan
        Tbmax_r1[i][Tbmax_r1[i] > 1.0e10] = np.nan
        Tbmax_r1[i][Tbmax_r1[i] < 1.0e5] = np.nan
        Tb_r1_1[i][Tb_r1_1[i] > 1.0e10] = np.nan
        Tb_r1_1[i][Tb_r1_1[i] < 1.0e5] = np.nan
        idx_r1[i] = np.arange(len(Tbmax_r1[i]))[np.isfinite(Tbmax_r1[i])]
        # Tb_mean_region1[i]=np.nanmean(Tb_r1_1[i][idx_r1[i]],axis=(1,2))
        Tb_mean_region1[i] = Tb[:, 119, 116]

    reg = "reg2"
    for i in range(len(chlist_img)):
        c = chlist_img[i]
        print("Channel: " + str(c))
        freq_mwa[i], Tb, Tbmax_r2[i], S, bmin[i], bmaj[i], img = pickle.load(
            open(
                "/data/20220915_MWA/Tb_20220915_" + str(reg) + "_ch" + str(c) + ".p",
                "rb",
            ),
            encoding="latin1",
        )
        Smax_r2[i] = np.nanmax(S, axis=(1, 2))
        Tb_r2 = Tb[:, 99:106, 116:123]
        # Tb_r2=Tb[:,102,120]
        Tb_r2_1[i] = Tb_r2 * 0
        idx_r2[i] = [0] * len(Tb_r2)
        xc_max_reg2[i] = [0] * len(Tb_r2)
        yc_max_reg2[i] = [0] * len(Tb_r2)
        for k in range(len(Tb_r2)):
            if (
                any(
                    map(
                        lambda x: x == 0,
                        np.array(np.where(Tb_r2[k] == np.nanmax(Tb_r2[k]))).flatten(),
                    )
                )
                == False
            ):
                Tb_r2_1[i][k] = Tb_r2[k]
                xc_max_reg2[i][k], yc_max_reg2[i][k] = np.where(
                    Tb[k][115:125, 112:122] == np.nanmax(Tb[k][115:125, 112:122])
                )
                idx_r2[i][k] = k
        Tbmax_r2[i][Tbmax_r2[i] > 1.0e10] = np.nan
        Tbmax_r2[i][Tbmax_r2[i] < 1.0e5] = np.nan
        Tb_r2_1[i][Tb_r2_1[i] > 1.0e10] = np.nan
        Tb_r2_1[i][Tb_r2_1[i] < 1.0e5] = np.nan
        idx_r2[i] = np.arange(len(Tbmax_r2[i]))[np.isfinite(Tbmax_r2[i])]
        # Tb_mean_region2[i]=np.nanmean(Tb_r2_1[i][idx_r2[i]],axis=(1,2))
        Tb_mean_region2[i] = Tb[:, 102, 120]

    freq_mwa = np.array(freq_mwa)
    Tbmax_r1 = np.array(Tbmax_r1)
    Tbmax_r2 = np.array(Tbmax_r2)
    Smax_r1 = np.array(Smax_r1)
    Smax_r2 = np.array(Smax_r2)
    bmin = np.array(bmin)
    bmaj = np.array(bmaj)
    Tb_mean_region2 = np.array(Tb_mean_region2)
    Tb_mean_region1 = np.array(Tb_mean_region1)
    Tb_r1_1 = np.array(Tb_r1_1)
    Tb_r2_1 = np.array(Tb_r2_1)
    pickle.dump(
        [
            freq_mwa,
            Tbmax_r1,
            Tbmax_r2,
            Tb_mean_region1,
            Tb_mean_region2,
            Smax_r1,
            Smax_r2,
            bmin,
            bmaj,
            Tb_r1_1,
            Tb_r2_1,
            idx_r1,
            idx_r2,
        ],
        open("/data/20220915_MWA/20220915_Tbmax.p", "wb"),
    )

# freq_mwa,Tbmax_reg1,Smax_reg1,bmin,bmaj = pickle.load(open('/data/20220915_MWA/20220915_reg1_Tbmax.p','rb'),encoding='latin1')
# freq_mwa,Tbmax_reg2,Smax_reg2,bmin,bmaj = pickle.load(open('/data/20220915_MWA/20220915_reg2_Tbmax.p','rb'),encoding='latin1')
(
    freq_mwa,
    Tbmax_r1,
    Tbmax_r2,
    Tb_mean_region1,
    Tb_mean_region2,
    Smax_r1,
    Smax_r2,
    bmin,
    bmaj,
    Tb_r1_1,
    Tb_r2_1,
    idx_r1,
    idx_r2,
) = pickle.load(open("/data/20220915_MWA/20220915_Tbmax.p", "rb"), encoding="latin1")
tmwa = 12608 + np.arange(Tb_mean_region1.shape[1]) * 0.5

# ----------------- HMI & some AIA
aiamap = Map(
    "/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220915_EUV/aia.lev1.171A_2022-09-15T03_30_09.35Z.image_lev1.fits"
)
hmifile = sorted(
    glob.glob("/sdata/20220915_hmi/hmi.m_45s.2022.09.12_03_01_30_TAI.magnetogram.fits")
)
hmimap = Map(hmifile)
norm = colors.Normalize(vmin=-0.1, vmax=0.1)
hmimap.plot_settings["norm"] = norm
hmimap.plot()
plt.show()

# -------------------------- FORWARD analysis
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel

axy = np.array(([850, 845], [220, 356]))
rxy = np.array(([900, 1300], [1000, 300]))
axy_skycoord = SkyCoord(
    axy[0] * u.arcsec, axy[1] * u.arcsec, frame=aiamap.coordinate_frame
)
rxy_skycoord = SkyCoord(
    rxy[0] * u.arcsec, rxy[1] * u.arcsec, frame=aiamap.coordinate_frame
)
aiahead = fits.open(
    "/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220915_EUV/aia.lev1.171A_2022-09-15T03_30_09.35Z.image_lev1.fits"
)[0].header
w = WCS(aiahead)
axy_xp, axy_yp = skycoord_to_pixel(axy_skycoord, w, origin=0, mode="all")

axy_p = 2048 + axy / 0.5
rxy_p = 2048 + rxy / 0.5
axy_p_fwd = 128 + axy * (1 / 18.75)
rxy_p_fwd = 128 + rxy * (1 / 18.75)


from scipy.io import readsav

rt = readsav("/data/Dropbox/STIX-MWA/20220915/forward/RT_params.sav")
r = rt["r3dall"]
tempobs = rt["tempall"]
densobs = rt["densall"]
brobs = rt["brall"]
bthobs = rt["bthall"]
bphobs = rt["bphall"]
bobs = np.sqrt(brobs * brobs + bthobs * bthobs + bphobs * bphobs)
taur = rt["taur"]
taul = rt["taul"]
tau_fwd = (taul + taur) * 0.5
dtaur = rt["dtaur"]
dtaul = rt["dtaul"]

t100 = readsav("/data/Dropbox/STIX-MWA/20220915/forward/widget_temp_100MHz.sav")
t100_v = readsav(
    "/data/Dropbox/STIX-MWA/20220915/forward/widget_temp_100MHz_stokesV.sav"
)
t200 = readsav("/data/Dropbox/STIX-MWA/20220915/forward/widget_temp_200MHz.sav")
fwdI_100 = t100["stokesstruct"]["I"][0]
fwdI_200 = t200["stokesstruct"]["I"][0]
fwdV_100 = t100_v["stokesstruct"]["V"][0]

print(
    "R1/ne",
    np.mean(densobs, axis=2)[int(rxy_p_fwd[0][0]), int(rxy_p_fwd[0][1])] / 1.0e8,
    "Plasma frequency:",
    np.sqrt(np.mean(densobs, axis=2)[int(rxy_p_fwd[0][0]), int(rxy_p_fwd[0][1])])
    * 9000
    / 1.0e6,
)
print(
    "R2/ne",
    np.mean(densobs, axis=2)[int(rxy_p_fwd[1][0]), int(rxy_p_fwd[1][1])] / 1.0e8,
    "Plasma frequency:",
    np.sqrt(np.mean(densobs, axis=2)[int(rxy_p_fwd[1][0]), int(rxy_p_fwd[1][1])])
    * 9000
    / 1.0e6,
)
print("R1/br", np.mean(brobs, axis=2)[int(rxy_p_fwd[0][0]), int(rxy_p_fwd[0][1])])
print("R2/br", np.mean(brobs, axis=2)[int(rxy_p_fwd[1][0]), int(rxy_p_fwd[1][1])])


f, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, sharex=True, sharey=True)
im0 = ax0.imshow(
    brobs.mean(axis=2),
    origin="lower",
    vmin=-3,
    vmax=3,
    extent=[-2400, 2400, -2400, 2400],
    cmap="coolwarm",
)
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="5%", pad=0.05)
f.colorbar(im0, cax=cax, orientation="vertical", label="B$_r$")
im1 = ax1.imshow(
    tempobs.mean(axis=2),
    origin="lower",
    vmin=1.0e5,
    vmax=2.0e6,
    extent=[-2400, 2400, -2400, 2400],
)
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
f.colorbar(im1, cax=cax, orientation="vertical", label="$T_e$ (K)")
im2 = ax2.imshow(
    np.log10(densobs.mean(axis=2)),
    origin="lower",
    vmin=6,
    vmax=9,
    extent=[-2400, 2400, -2400, 2400],
)
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
f.colorbar(im2, cax=cax, orientation="vertical", label="log($n_e$ cm$^{-3}$)")
im3 = ax3.imshow(
    bthobs.mean(axis=2),
    origin="lower",
    vmin=-3,
    vmax=3,
    cmap="coolwarm",
    extent=[-2400, 2400, -2400, 2400],
)
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="5%", pad=0.05)
f.colorbar(im3, cax=cax, orientation="vertical", label="B$_{\\theta}$")
ax3.set_xlabel("Solar-X (arcsec)")
ax2.set_xlabel("Solar-X (arcsec)")
ax2.set_ylabel("Solar-Y (arcsec)")
ax0.set_ylabel("Solar-Y (arcsec)")
plt.show()

f, ax = plt.subplots(1, 1)
im = ax.imshow(
    np.log10(tau_fwd.mean(axis=2)),
    origin="lower",
    cmap="coolwarm",
    vmin=-2,
    vmax=2,
    extent=[-2400, 2400, -2400, 2400],
)
f.colorbar(im, label="log($\\tau$)")
ax.set_xlabel("Solar-X (arcsec)")
ax.set_ylabel("Solar-Y (arcsec)")
plt.show()

f, (ax, ax0, ax1) = plt.subplots(1, 3, sharex=True, sharey=True)
im = ax.imshow(
    fwdI_100,
    origin="lower",
    cmap="coolwarm",
    vmin=1.0e4,
    vmax=2.0e6,
    extent=[-2400, 2400, -2400, 2400],
)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
f.colorbar(im, cax=cax, orientation="vertical", label="$T_B$ (K)")
ax.set_xlabel("Solar-X (arcsec)")
ax.set_ylabel("Solar-Y (arcsec)")
im0 = ax0.imshow(
    fwdV_100,
    origin="lower",
    cmap="coolwarm",
    vmin=1.0e3,
    vmax=5.0e4,
    extent=[-2400, 2400, -2400, 2400],
)
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="5%", pad=0.05)
f.colorbar(im0, cax=cax, orientation="vertical", label="$T_B$ (K)")
ax0.set_xlabel("Solar-X (arcsec)")
im1 = ax1.imshow(
    fwdV_100 / fwdI_100 * 100,
    origin="lower",
    cmap="coolwarm",
    vmin=0.01,
    vmax=10,
    extent=[-2400, 2400, -2400, 2400],
)
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
f.colorbar(im1, cax=cax, orientation="vertical", label="V/I (%)")
ax1.set_xlabel("Solar-X (arcsec)")
plt.show()


# ------------------------- Carrington HMI Synoptic maps & PFSS
# Radial HMI file read

do_pfss = 0
if do_pfss:  # Migrated to server1062 due to memory issues
    print(datetime.datetime.now())
    hmimap_radial_ = Map("/data/20220915_MWA/pfss/hmi.Synoptic_Mr.2262.fits")
    hmimap_radial_data = hmimap_radial_.data
    hmimap_radial_data[np.isnan(hmimap_radial_data)] = 0
    hmimap_radial = Map(hmimap_radial_data, hmimap_radial_.meta)
    nrho = 50
    rss = 1.5
    pfss_in = pfsspy.Input(hmimap_radial, nrho, rss)
    pfss_out = pfsspy.pfss(pfss_in)
    tracer = tracing.FortranTracer()
    r = 1.2 * const.R_sun
    hp_lon = np.linspace(800, 900, 10) * u.arcsec
    hp_lat = np.linspace(300, 400, 10) * u.arcsec
    # Make a 2D grid from these 1D points
    lon, lat = np.meshgrid(hp_lon, hp_lat)
    seeds = SkyCoord(lon.ravel(), lat.ravel(), frame=hmimap.coordinate_frame)
    field_lines = tracer.trace(seeds, pfss_out)
    field_line = field_lines[0]
    B = field_line.b_along_fline
    r = field_line.coords.radius
    print(datetime.datetime.now())
    fc = [0] * 4
    i = 0
    for fline in field_lines:
        fc[i] = fline.coords
        i = i + 1
    pickle.dump(fc, open("/data/20220915_MWA/pfss/hmi_2262_fc.p", "rb"))

field_lines = pickle.load(open("/data/20220915_MWA/pfss/hmi_2262_fc.p", "rb"))

m = aiamap
fig = plt.figure()
ax = plt.subplot(projection=m)
m.plot()
plt.colorbar()
for fline in field_lines:
    ax.plot_coord(fline, color="red", linewidth=1)
    # Set the axes limits. These limits have to be in pixel values
    # ax.set_xlim(0, 180)
    # ax.set_ylim(45, 135)
    ax.set_title("Photospheric field and traced field lines")
plt.show()

plt.plot((r[144, 164] - 1) * 7.0e2, bobs[144, 164], "o-")
plt.xlabel("Radial Distance (Mm)")
plt.ylabel("|B|(Gauss)")
plt.show()

plt.plot((r[144, 164] - 1) * 7.0e2, densobs[144, 164], "o-")
plt.xlabel("Radial Distance (Mm)")
plt.ylabel("|B|(Gauss)")
plt.show()

# ----------------------- PFSS from GONG
# gong_map = Map("/data/20220915_MWA/pfss/mrzqs220915t0004c2262_289.fits")
gong_map = Map("/data/20220915_MWA/pfss/mrzqs220911t1104c2262_336.fits")
hmimap_exp = Map(
    "/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220911_hmi/hmi.m_45s.2022.09.11_03_31_30_TAI.magnetogram.fits"
)
# aiamap_exp=Map("/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220911_hmi/")
nrho = 100
rss = 1.5
pfss_in = pfsspy.Input(gong_map, nrho, rss)
pfss_out = pfsspy.pfss(pfss_in)
tracer = tracing.FortranTracer()
r = 1.2 * const.R_sun
hp_lon = np.linspace(100, 800, 20) * u.arcsec
hp_lat = np.linspace(-400, 500, 20) * u.arcsec
# Make a 2D grid from these 1D points
lon, lat = np.meshgrid(hp_lon, hp_lat)
gg_lon = np.linspace(40, 160, 20) * u.deg
gg_lat = np.linspace(0, 50, 20) * u.deg
# Make a 2D grid from these 1D points
# lon, lat = np.meshgrid(gg_lon, gg_lat)
# seeds = SkyCoord(lon.ravel(), lat.ravel(), frame=gong_map.coordinate_frame)
seeds = SkyCoord(lon.ravel(), lat.ravel(), frame=hmimap_exp.coordinate_frame)
field_lines = tracer.trace(seeds, pfss_out)
field_line = field_lines[0]
B = field_line.b_along_fline
r = field_line.coords.radius

from sunpy.coordinates import frames

fhp = frames.Helioprojective(obstime="20220915T03:30:00")


from sunpy.coordinates import RotatedSunFrame

durations = np.array([5]) * u.day
diffrot_point = SkyCoord(
    RotatedSunFrame(base=fline.coords.helioprojective, duration=durations)
)


f = plt.figure()
ax = plt.subplot(projection=pfss_in.map)
gong_map.plot()
for fline in field_lines:
    ax.plot_coord(fline.coords, color="black", linewidth=1)
plt.colorbar()
plt.show()

plt.imshow(
    hmimap_exp.data[::-1, ::-1],
    origin="lower",
    extent=[-1024, 1024, -1024, 1024],
    vmin=-20,
    vmax=20,
)
plt.show()

m = hmimap_exp
fhp = frames.HeliographicCarrington(obstime="2022-09-15T23:30:00")
fig = plt.figure()
ax = plt.subplot(projection=m)
m.plot()
plt.colorbar()
for fline in field_lines[0:10]:
    f1 = fline.coords  # .transform_to(frame=fhp)
    ax.plot_coord(f1, color="red", linewidth=1)
    ax.set_title("Photospheric field and traced field lines")
plt.show()

# ---------------------- STIX X-ray

stix = readsav(
    "/data/Dropbox/STIX-MWA/20220915/STIX_data/stix_lightcurves_spec_8-s_Earth-UT_20220915.sav"
)
stix_data = stix["stix_lcstr"]
stix_data0 = stix_data["data"][0]
tstix_low = (
    stix_data["ut"][0][:, 0] - stix_data["ut"][0][:, 0][0] + 11022
)  # Start time 03:03:42
tstix_high = (
    stix_data["ut"][0][:, 1] - stix_data["ut"][0][:, 1][0] + 11041
)  # Start time 03:04:01

# plt.imshow(region1_max,aspect='auto',origin='lower',interpolation=None)

# plt.imshow(img[0].mean(axis=0),aspect='auto',origin='lower',interpolation=None)

stix_img2 = fits.open(
    "/data/Dropbox/STIX-MWA/20220915/STIX_data/stix-thermal-image_flare2_peak2_2022-09-15.fits"
)
stix_img2_head = stix_img2[0].header
stix_img2_map = Map(
    "/data/Dropbox/STIX-MWA/20220915/STIX_data/stix-thermal-image_flare2_peak2_2022-09-15.fits"
)
stix_wcs = WCS(stix_img2_head)
aia_wcs = WCS(aiahead)
stix_loc = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=stix_img2_map.coordinate_frame)
aia_loc = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=aiamap.coordinate_frame)


Tb_ds = np.mean(Tb[:, :, 120:150, 140:170], axis=(2, 3))
for i in range(5):
    Tb_ds[i] = Tb_ds[i]  # -Tb_ds[i][0]
plt.imshow(Tb_ds, origin="lower", aspect="auto", cmap="jet", vmin=1.0e3, vmax=1.0e5)
plt.show()


f, ax = plt.subplots(1, 1)
ax1 = ax.twinx()
c = ["k", "g", "magenta", "cyan", "brown"]
# fact=[7,1,0.5,0.5,0.4,1,1,1,1,1,1,1,1]
i = 3
ax1.plot(
    tmwa,
    Tb_ds[i] * 5.0e-4,
    "o-",
    linewidth=1,
    markersize=1,
    label="MWA (" + str(freq_mwa[i]) + " MHz)",
)
# ax1.plot(tmwa,region1_mean[i]*fact[i],'o-',linewidth=1,markersize=1,label='MWA ('+str(freq_mwa[i])+' MHz)',color=c[i])
# ax1.plot(tmwa,(region1_mean[2]-17)*0.2,'o-',linewidth=1,markersize=3,label='MWA ('+str(freq_mwa[2])+' MHz)',color='green')
# ax1.plot(tmwa,region1_mean[0:].mean(axis=0),'o-',linewidth=1,markersize=3,label='MWA ($\\nu$-Average)',color='k')
ax.plot(tstix_low, stix_data0[:, 0], "-", label="STIX (6-7 keV)", color="r")
ax.plot(tstix_high, stix_data0[:, 1], "-", label="STIX (16-22 keV)", color="b")
ax.set_xticks(tstix_low[::50])
ax.set_xticklabels(
    [
        "03:03:42",
        "03:14:12",
        "03:21:14",
        "03:28:16",
        "03:36:20",
        "03:44:58",
        "03:55:11",
        "04:03:29",
        "04:14:13",
    ]
)
ax.set_ylabel("Amplitude")
ax.set_xlabel("Time (HH:MM:SS)")
ax.legend()
ax1.legend(loc=4)
ax.set_yscale("log")
ax.set_ylabel("STIX Count Flux (counts/s/cm$^2$/keV)")
ax.set_xlabel("Time (HH:MM:SS UT) (Start Time: 03:03:42 UT)")
ax1.set_ylabel("Amplitude")
# ax1.set_ylim(0.5,4)
plt.show()

# Tb_region1_max[:,272:273] = np.nan

mway = np.nanmax(Tb_region2_mean[0:2], axis=0) * 1.0e-5
f, ax = plt.subplots(1, 1)
ax1 = ax.twinx()
c = ["k", "g", "magenta", "cyan", "brown"]
fact = [7, 1, 0.5, 0.5, 0.4]
i = 5
ax1.plot(
    tmwa,
    mway,
    "o-",
    linewidth=1,
    markersize=1,
    label="MWA (<"
    + str(freq_mwa[i])
    + " MHz) | Harmonic Height"
    + str(np.round(height_Mm_harmonic[i], 0))
    + " Mm",
    color=c[0],
)
# ax1.plot(tmwa,(region1_mean[2]-17)*0.2,'o-',linewidth=1,markersize=3,label='MWA ('+str(freq_mwa[2])+' MHz)',color='green')
# ax1.plot(tmwa,region1_mean[0:].mean(axis=0),'o-',linewidth=1,markersize=3,label='MWA ($\\nu$-Average)',color='k')
ax.plot(tstix_low, stix_data0[:, 0], "-", label="STIX (6-7 keV)", color="r")
ax.plot(tstix_high, stix_data0[:, 1], "-", label="STIX (16-22 keV)", color="b")
ax.set_xticks(tstix_low[::50])
ax.set_xticklabels(
    [
        "03:03:42",
        "03:14:12",
        "03:21:14",
        "03:28:16",
        "03:36:20",
        "03:44:58",
        "03:55:11",
        "04:03:29",
        "04:14:13",
    ]
)
ax1.set_ylabel("T$_B$ ($\\times 10^5$ K)")
ax.set_xlabel("Time (HH:MM:SS)")
ax.legend()
ax1.legend(loc=2)
ax.set_yscale("log")
ax.set_ylabel("STIX Count Flux (counts/s/cm$^2$/keV)")
ax.set_xlabel("Time (HH:MM:SS UT) (Start Time: 03:03:42 UT)")
ax1.set_ylim(0.1, 0.6)  # Chan 0
plt.show()

mway = np.nanmax(Tb_region2_mean[15:], axis=0) * 1.0e-5
f, ax = plt.subplots(1, 1)
ax1 = ax.twinx()
c = ["k", "g", "magenta", "cyan", "brown"]
fact = [7, 1, 0.5, 0.5, 0.4]
i = 15
ax1.plot(
    tmwa,
    mway,
    "o-",
    linewidth=1,
    markersize=1,
    label="MWA (>"
    + str(freq_mwa[i])
    + " MHz) | Harmonic Height"
    + str(np.round(height_Mm_harmonic[i], 0))
    + " Mm",
    color=c[0],
)
# ax1.plot(tmwa,(region1_mean[2]-17)*0.2,'o-',linewidth=1,markersize=3,label='MWA ('+str(freq_mwa[2])+' MHz)',color='green')
# ax1.plot(tmwa,region1_mean[0:].mean(axis=0),'o-',linewidth=1,markersize=3,label='MWA ($\\nu$-Average)',color='k')
ax.plot(tstix_low, stix_data0[:, 0], "-", label="STIX (6-7 keV)", color="r")
ax.plot(tstix_high, stix_data0[:, 1], "-", label="STIX (16-22 keV)", color="b")
ax.set_xticks(tstix_low[::50])
ax.set_xticklabels(
    [
        "03:03:42",
        "03:14:12",
        "03:21:14",
        "03:28:16",
        "03:36:20",
        "03:44:58",
        "03:55:11",
        "04:03:29",
        "04:14:13",
    ]
)
ax1.set_ylabel("T$_B$ ($\\times 10^5$ K)")
ax.set_xlabel("Time (HH:MM:SS)")
ax.legend()
ax1.legend(loc=2)
ax.set_yscale("log")
ax.set_ylabel("STIX Count Flux (counts/s/cm$^2$/keV)")
ax.set_xlabel("Time (HH:MM:SS UT) (Start Time: 03:03:42 UT)")
ax1.set_ylim(0.25, 10.5)  # Chan 0
plt.show()


# r2_y=Tb_region2_max[:,14:16].mean(axis=1);r1_y=Tb_region1_max[:,258:261].mean(axis=1)
r2_y = Tb_region2_max[:, 15]
r1_y = Tb_region1_max[:, 259]
r1_y0 = Tb_region1_max[:, 20]
r2_y0 = Tb_region2_max[:, 0]
r2_y[10] = np.nan
r2_y[11] = np.nan
r2_y[-1] = np.nan
r1_y[-1] = np.nan
diff_r1 = r1_y - r1_y0
diff_r2 = r2_y - r2_y0

f, ax = plt.subplots(2, 1)
ax0 = ax[0]
ax1 = ax[1]
ax0.plot(freq[0:5], diff_r1[0:5] / 1.0e3, "o-", label="Region 1 (04:16:30)")
ax1.plot(freq[13:], diff_r2[13:] / 1.0e6, "o-", label="Region 2 (03:32:30)")
ax0.set_xlabel("Frequency (MHz)")
ax1.set_xlabel("Frequency (MHz)")
ax0.set_ylabel("$T_B$ (kK)")
ax1.set_ylabel("$T_B$ (MK)")
ax0.legend()
ax1.legend()
plt.show()

f, ax = plt.subplots(1, 1)
ax.plot(freq, diff_r2 / 1.0e6, "o-", label="Region 2")
ax.plot(freq, diff_r1 / 1.0e6, "o-", label="Region 1")
ax.legend()
ax.set_xlabel("Frequency (MHz)")
ax.set_ylabel("$T_B$ (MK)")
ax.set_ylim(0, 0.6)
plt.show()

f, ax = plt.subplots(1, 1)
ax.plot(freq, r2_y / 1.0e6, "o-", label="Region 2 (03:32:30)")
ax.plot(freq, r1_y / 1.0e6, "o-", label="Region 1 (04:16:30)")
# ax.plot(freq,r2_y0/1.e6,'o-',label='Region 2 (03:30:10)')
# ax.plot(freq,r1_y0/1.e6,'o-',label='Region 1 (03:33:20)')
ax.legend()
ax.set_xlabel("Frequency (MHz)")
ax.set_ylabel("$T_B$ (MK)")
ax.set_ylim(1.0e-3, 10)
ax.set_yscale("log")
plt.show()

plt.plot(region1_max[0], "o-", label="Channel " + str(chan[i]))
plt.legend()
plt.show()

region1_max1 = region1_max * 1.0
for j in range(5):
    region1_max1[j] = region1_max[j] - region1_max[j][0]

f, ax = plt.subplots(1, 1)
ax.imshow(
    region1_max1,
    aspect="auto",
    origin="lower",
    cmap="jet",
    interpolation=None,
    vmin=0,
    vmax=3,
)
ax.set_xticks([0, 50, 100, 150, 200, 250, 300, 350])
ax.set_xticklabels(
    [
        "03:30:00",
        "03:38:20",
        "03:46:40",
        "03:55:00",
        "04:03:20",
        "04:11:40",
        "04:20:00",
        "04:28:20",
    ]
)
ax.set_yticks(np.arange(5))
ax.set_yticklabels(freq[chan])
plt.show()

# euv094_list=sorted(glob.glob('/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220830/20220830_EUV/*.94A*.fits'))
# euv171_list=sorted(glob.glob('/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220830/20220830_EUV/*.171A*.fits'))

# -------------------- GOES
goes = readsav("/data/Dropbox/STIX-MWA/20220915/20220915_goes.sav")
gf0540 = goes["YCLEAN"][1]
gf1080 = goes["YCLEAN"][0]
gtime = goes["tarray"] + 10800

# ------------------ AIA Analysis
euv094_list = sorted(
    glob.glob("/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220915_EUV/*.94A*.fits")
)[0:400]
euv171_list = sorted(
    glob.glob("/media/rohit/Seagate_Expansion_Drive/MWA-STIX/20220915_EUV/*.171A*.fits")
)[0:400]


aiamap171_diff = [0] * (len(euv171_list) - 10)
aiamap171 = [0] * len(euv171_list)
t171 = [0] * len(euv094_list)
t171_str = [0] * len(euv094_list)
count171_r1 = [0] * len(euv171_list)
count171_r2 = [0] * len(euv171_list)
count171_r1_max = [0] * len(euv171_list)
count171_r2_max = [0] * len(euv171_list)
for i in range(len(aiamap171_diff)):
    print("171:", i)
    aia1 = Map(euv171_list[i + 10])
    aia_ref = Map(euv171_list[i])
    aiamap171_diff[i] = Map(aia1.data - aia_ref.data, aia1.meta)
    # aiamap171[i] = aia1
    t171_str[i] = aia1.meta["date-obs"].split("T")[1]
    tt = t171_str[i].split(":")
    t171[i] = int(tt[0]) * 3600 + int(tt[1]) * 60 + float(tt[2])
    count171_r1[i] = aiamap171_diff[i].data[2560:2660, 3420:3440].mean()
    count171_r2[i] = aiamap171_diff[i].data[2300:2440, 3400:3550].mean()
    count171_r1_max[i] = aiamap171_diff[i].data[2560:2660, 3420:3440].max()
    count171_r2_max[i] = aiamap171_diff[i].data[2300:2440, 3400:3550].max()
    plot_aia = 0
    if plot_aia:
        ii = "%04d" % i
        fig = plt.figure(figsize=(20, 10))
        ax0 = fig.add_subplot(121, projection=aia1)
        p0 = aiamap171[i].plot(axes=ax0, aspect="auto", cmap="coolwarm")
        ax0.set_ylim(2200, 3000)
        ax0.set_xlim(3200, 4000)
        ax1 = fig.add_subplot(122, projection=aia1)
        p1 = aiamap171_diff[i].plot(
            axes=ax1, aspect="auto", cmap="coolwarm", vmin=-400, vmax=800
        )
        ax1.set_ylim(2200, 3000)
        ax1.set_xlim(3200, 4000)
        plt.savefig(
            "/data/20220915_MWA_old/20220915_pngs/aia171_pngs/aia171_"
            + str(ii)
            + ".png",
            dpi=60,
        )
        plt.close()

pickle.dump(
    [
        t171,
        t171_str,
        count171_r1,
        count171_r2,
        count171_r1_max,
        count171_r2_max,
    ],
    open("/data/20220915_MWA_old/20220915_171AA_timeseries.p", "wb"),
)

aiamap094_diff = [0] * (len(euv094_list) - 10)
aiamap094 = [0] * len(euv094_list)
t094 = [0] * len(euv094_list)
t094_str = [0] * len(euv094_list)
count094_r1 = [0] * len(euv094_list)
count094_r2 = [0] * len(euv094_list)
count094_r1_max = [0] * len(euv094_list)
count094_r2_max = [0] * len(euv094_list)
for i in range(len(aiamap094_diff)):
    print("94:", i)
    aia1 = Map(euv094_list[i + 10])
    aia_ref = Map(euv094_list[i])
    aiamap094_diff[i] = Map(aia1.data - aia_ref.data, aia1.meta)
    # aiamap094[i] = aia1
    t094_str[i] = aia1.meta["date-obs"].split("T")[1]
    tt = t094_str[i].split(":")
    t094[i] = int(tt[0]) * 3600 + int(tt[1]) * 60 + float(tt[2])
    count094_r1[i] = aiamap094_diff[i].data[2560:2660, 3420:3440].mean()
    count094_r2[i] = aiamap094_diff[i].data[2300:2440, 3400:3550].mean()
    count094_r1_max[i] = aiamap094_diff[i].data[2560:2660, 3420:3440].max()
    count094_r2_max[i] = aiamap094_diff[i].data[2300:2440, 3400:3550].max()
    plot_aia = 0
    if plot_aia:
        ii = "%04d" % i
        fig = plt.figure(figsize=(20, 10))
        ax0 = fig.add_subplot(121, projection=aia1)
        p0 = aiamap094[i].plot(axes=ax0, aspect="auto", cmap="coolwarm")
        ax0.set_ylim(2200, 3000)
        ax0.set_xlim(3200, 4000)
        ax1 = fig.add_subplot(122, projection=aia1)
        p1 = aiamap094_diff[i].plot(
            axes=ax1, aspect="auto", cmap="coolwarm", vmin=-400, vmax=800
        )
        ax1.set_ylim(2200, 3000)
        ax1.set_xlim(3200, 4000)
        plt.savefig(
            "/data/20220915_MWA_old/20220915_pngs/aia094_pngs/aia094_"
            + str(ii)
            + ".png",
            dpi=60,
        )
        plt.close()

pickle.dump(
    [t094, t094_str, count094_r1, count094_r2, count094_r1_max, count094_r2_max],
    open("/data/20220915_MWA_old/20220915_094AA_timeseries.p", "wb"),
)

(
    t094,
    t094_str,
    count094_r1,
    count094_r2,
    count094_r1_max,
    count094_r2_max,
) = pickle.load(
    open("/data/20220915_MWA_old/20220915_094AA_timeseries.p", "rb"),
)
(
    t171,
    t171_str,
    count171_r1,
    count171_r2,
    count171_r1_max,
    count171_r2_max,
) = pickle.load(
    open("/data/20220915_MWA_old/20220915_171AA_timeseries.p", "rb"),
)

f, ax = plt.subplots(1, 1)
ax.plot(
    tmwa_10s,
    np.nanmax(Tb_region2_max_10s[0:2, :], axis=0),
    "o-",
    label=str(freq_mwa[0]),
)
ax.plot(
    tmwa_10s,
    np.nanmax(Tb_region2_max_10s[2:4, :], axis=0),
    "o-",
    label=str(freq_mwa[2]),
)
ax.plot(
    tmwa_10s,
    np.nanmax(Tb_region2_max_10s[4:6, :], axis=0),
    "o-",
    label=str(freq_mwa[4]),
)
ax.plot(
    tmwa_10s,
    np.nanmax(Tb_region2_max_10s[6:8, :], axis=0),
    "-",
    label=str(freq_mwa[6]),
)
ax1 = ax.twinx()
ax1.plot(
    tstix_low,
    stix_data0[:, 1],
    ".-",
    color="k",
    label="STIX (16-22 keV)",
)
ax.set_ylim(1.0e3, 1.20e5)
ax1.set_ylim(0, 2.5)
ax.legend()
ax1.legend(loc=2)
ax1.set_xlim(3600 * 3, 3600 * 4.5)
ax1.set_xticks([12000, 12500, 13000, 13500, 14000])
ax1.set_xticklabels(
    [
        "03:20:00",
        "03:28:20",
        "03:36:40",
        "03:45:00",
        "03:53:20",
    ]
)
ax.set_xlabel("Time (HH:MM:SS UT)")
ax.set_ylabel("T$_B$ (K)")
ax.set_xlim(12608, 14000)
plt.show()


f, (ax, ax1) = plt.subplots(2, 1, sharex=True)
ax.plot(
    t094[0:370],
    count094_r1_max[0:370] / np.max(count094_r1_max),
    "-",
    label="A1 (94 $\AA$)",
)
ax.plot(
    t094[0:370],
    count094_r2_max[0:370] / np.max(count094_r2_max),
    "-",
    label="A2 (94 $\AA$)",
)
# ax.plot(
#    t094[:-10],
#    count094_r2_max[:-10] / np.max(count094_r2_max),
#    "-",
#    label="Region 2 (94 $\AA$)",
# )
# ax.plot(
#    t171[:-10],
#    count171_r1_max[:-10] / np.max(count171_r1_max),
#    "-",
#    label="Region 1 (171 $\AA$)",
# )
ax.plot(
    t171[:-10],
    count171_r2_max[:-10] / np.max(count171_r2_max),
    "-",
    label="A2 (171 $\AA$)",
)
ax.plot(gtime, gf0540 / np.max(gf0540), label="GOES $0.5-4.0 \AA$")
ax.plot(gtime, gf1080 / np.max(gf1080), label="GOES $1.0-8.0 \AA$")
ax2 = ax1.twinx()
yp = np.nanmax(Tb_region2_max_10s[0:5, :], axis=0)
yp[yp < 3.0e4] = np.nan
yp[yp > 7.0e4] = np.nan
ax2.plot(
    tmwa_10s,
    yp,
    "-",
    label="(R2) MWA ("
    + str(np.round(freq_mwa[0], 2))
    + "-"
    + str(np.round(freq_mwa[5], 2))
    + ") MHz",
    color="k",
)
yp1 = np.nanmax(Tb_region1_max_10s[13:18, :], axis=0) / 200.0 + 2.8e4
yp1[np.isinf(yp1)] = np.nan
ax2.plot(
    tmwa_10s,
    yp1,
    "-",
    label="(R1) MWA ("
    + str(np.round(freq_mwa[13], 2))
    + "-"
    + str(np.round(freq_mwa[16], 2))
    + ") MHz",
    color="brown",
)
ax1.plot(
    tstix_low,
    stix_data0[:, 0] / np.nanmax(stix_data0[:, 0]),
    "-",
    label="STIX (6-7 keV)",
)
ax1.plot(
    tstix_low,
    stix_data0[:, 1] / np.nanmax(stix_data0[:, 1]),
    "-",
    label="STIX (16-22 keV)",
)
ax.legend()
ax1.legend(loc=3)
ax2.legend(loc=2)
ax2.set_ylabel("$T_B$ (K)")
ax2.set_ylim(3.0e4, 7.0e4)
ax.set_ylabel("Amplitude")
ax1.set_ylabel("Amplitude")
ax1.set_xlabel("Time (HH:MM:SS UT)")
ax1.set_xlim(3600 * 3, 3600 * 4.5)
ax1.set_xticks(
    [11000, 11500, 12000, 12500, 13000, 13500, 14000, 14500, 15000, 15500, 16000]
)
ax1.set_xticklabels(
    [
        "03:03:20",
        "03:11:40",
        "03:20:00",
        "03:28:20",
        "03:36:40",
        "03:45:00",
        "03:53:20",
        "04:01:40",
        "04:10:00",
        "04:18:20",
        "04:26:40",
    ]
)
plt.show()

f, ax = plt.subplots(1, 1)
ax.plot(np.nanmean(Tb, axis=(0, 2, 3)), "o-", label="Average")
ax.plot(np.nanmax(Tb_region1_max_10s[13:18], axis=0), "o-", label="Region1")
ax.plot(np.nanmax(Tb_region2_max_10s[0:5, :], axis=0), "o-", label="Region2")
ax.legend()
ax.set_ylim(1.0e4, 1.0e7)
plt.show()

f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
im1 = ax1.imshow(
    np.nanmean(Tb[13:18, 200:210], axis=(0, 1)),
    origin="lower",
    vmin=1.0e4,
    vmax=1.0e6,
    cmap="jet",
    extent=[-5000, 5000, -5000, 5000],
)
im2 = ax2.imshow(
    np.nanmax(Tb[0:5, 0:60], axis=(0, 1)),
    origin="lower",
    vmin=1.0e4,
    vmax=5.0e4,
    cmap="jet",
    extent=[-5000, 5000, -5000, 5000],
)
ax1.add_patch(
    plt.Circle((0, 0), 16 * 60, color="k", alpha=0.8, fill=False, linewidth=5)
)
ax2.add_patch(
    plt.Circle((0, 0), 16 * 60, color="k", alpha=0.8, fill=False, linewidth=5)
)
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
f.colorbar(im1, cax=cax, orientation="vertical")
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
f.colorbar(im2, cax=cax, orientation="vertical", label="$T_B$(K)")
ax1.set_xlim(-1500, 1500)
ax2.set_xlim(-1500, 1500)
ax1.set_ylim(-1500, 1500)
ax2.set_ylim(-1500, 1500)
ax1.set_ylabel("Solar-Y (arcsec)")
ax1.set_xlabel("Solar-X (arcsec)")
ax2.set_xlabel("Solar-X (arcsec)")
plt.show()


r0 = region2_10s[0][54]
r1 = region2_10s[2][51]
r2 = region2_10s[3][47]
r0_max_ = np.where(r0 == np.max(r0))
r1_max_ = np.where(r1 == np.max(r1))
r2_max_ = np.where(r2 == np.max(r2))
r0_max = (np.array([r0_max_[0][0] + 115, r0_max_[1][0] + 95]) - 100) * 50
r1_max = (np.array([r1_max_[0][0] + 115, r1_max_[1][0] + 95]) - 100) * 50
r2_max = (np.array([r2_max_[0][0] + 115, r2_max_[1][0] + 95]) - 100) * 50
print(r0_max, r1_max, r2_max)


time_mwa = np.arange(380) * 10 + 3 * 3600 + 30 * 60
i = 0
time_aia171 = [0] * len(euv171_list)
for e in euv171_list:
    l = euv171_list[i].split("T")[2].split("Z")[0].split("_")
    time_aia171[i] = int(l[0]) * 3600 + int(l[1]) * 60 + float(l[2])
    i = i + 1
time_aia171 = np.array(time_aia171)


import matplotlib.patches as patches

a = 00 / 0.6
b = 00 / 0.6
stix_loc1 = np.array([1590 + 100, 455 - 100]) * 0.6
stix_loc1_pix = stix_loc1 / 0.6 + 2096
euv094_file = euv094_list[150]
map94 = Map(euv094_file)
map94_fits = fits.open(euv094_file)
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection=map94)
map94.plot(axes=ax)
ax.contour(
    Tb_region2_10s_full[0][54],
    [2.2e4, 2.4e4],
    extent=[-6285 + a, 10381 + b, -6285 + a, 10381 + b],
    colors="cyan",
)
ax.contour(
    Tb_region2_10s_full[2][51],
    [4.0e4, 4.2e4],
    extent=[-6285 + a, 10381 + b, -6285 + a, 10381 + b],
    colors="green",
)
ax.contour(
    Tb_region2_10s_full[3][47],
    [2.7e4, 2.8e4],
    extent=[-6285 + a, 10381 + b, -6285 + a, 10381 + b],
    colors="red",
)
ax.text(
    3100,
    2000,
    "Freq: 78 MHz | 03:37:54",
    style="italic",
    bbox={"facecolor": "cyan", "alpha": 0.5, "pad": 10},
)
ax.text(
    3100,
    2100,
    "Freq: 88 MHz | 03:38:39",
    style="italic",
    bbox={"facecolor": "green", "alpha": 0.5, "pad": 10},
)
ax.text(
    3100,
    2200,
    "Freq: 104 MHz | 03:39:09",
    style="italic",
    bbox={"facecolor": "red", "alpha": 0.5, "pad": 10},
)
patch = patches.Ellipse(
    stix_loc1_pix, 40 / 0.6, 20 / 0.6, 0, fc="none", ls="solid", ec="orange", lw="3."
)
ax.add_patch(patch)
plt.show()

for j in range(30, 50):
    idx = ut.find_predecessor(time_aia171, time_mwa[j])
    euv171_file = euv171_list[idx[0]]
    map171 = Map(euv171_file)
    aiahead1 = []
    aiahead1 = map171.meta.copy()
    aiahead1["naxis1"] = 200
    aiahead1["naxis2"] = 200
    aiahead1["CRPIX1"] = 99
    aiahead1["CRPIX2"] = 99
    aiahead1["CRVAL1"] = 0
    aiahead1["CRVAL2"] = 0
    aiahead1["CDELT1"] = 50
    aiahead1["CDELT2"] = 50
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection=map171)
    # lev_ct = [[1.0e5, 2.0e5, 3.0e5], [3.0e4,4.0e4], [2.0e4, 3.0e4], [3.e5, 5.e5], [2.e4, 3.e5]]
    lev_ct = [[1.0e6], [1.0e5], [3.0e4, 5.0e4], [3.0e5], [3.0e5]]
    cc = ["r", "g", "b", "k", "brown"]
    for k in range(5):
        mwadata = Tb[k][j]
        mwadata1 = mwadata * 0
        mwadata1[90:160, 90:160] = mwadata[90:160, 90:160]
        mwadata1[mwadata1 < 0.5] = 0
        mwamap = Map(mwadata, aiahead1)
        mwamap1 = Map(mwadata1, aiahead1)
        ii = "%04d" % j
        p0 = map171.plot(axes=ax, aspect="auto", vmin=0, vmax=800)
        # frac_r1 = [0.05,0.1,0.12,0.2,0.3,0.4,0.5]
        for f in lev_ct[k]:
            # lev_r1 = np.nanmax(mwamap.data) * f
            c1 = mwamap1.contour(level=f * u.ct)
            if len(c1) != 0:
                ax.plot_coord(c1[0], color=cc[k])
        ax.text(
            -1400,
            -700 + k * 500,
            "Freq: " + str(np.round(freq[chan[k]], 1)) + " MHz",
            style="italic",
            bbox={"facecolor": cc[k], "alpha": 0.5, "pad": 10},
        )
    ax.set_ylim(-1600, 5000)
    ax.set_xlim(-1600, 5000)
    plt.savefig("/home/rohit/20220915/20220915_pngs/img_" + str(ii) + ".png")
    plt.close()

for j in range(len(time_mwa)):
    idx = ut.find_predecessor(time_aia171, time_mwa[j])
    euv171_file = euv171_list[idx[0]]
    map171 = Map(euv171_file)
    aiahead1 = []
    aiahead1 = map171.meta.copy()
    aiahead1["naxis1"] = 200
    aiahead1["naxis2"] = 200
    aiahead1["CRPIX1"] = 99
    aiahead1["CRPIX2"] = 99
    aiahead1["CRVAL1"] = 0
    aiahead1["CRVAL2"] = 0
    aiahead1["CDELT1"] = 50
    aiahead1["CDELT2"] = 50
    fig = plt.figure(figsize=(8, 20))
    ax0 = fig.add_subplot(211, projection=map171)
    ax1 = fig.add_subplot(212)
    lev_ct = [[1.0e5], [1.0e5], [3.0e4, 5.0e4], [3.0e5], [3.0e5]]
    cc = ["r", "g", "b", "k", "brown"]
    for k in range(5):
        mwadata = Tb[k][j]
        mwadata1 = mwadata * 0
        mwadata1[90:160, 110:170] = mwadata[90:160, 110:170]
        mwadata1[mwadata1 < 0.5] = 0
        mwamap = Map(mwadata, aiahead1)
        mwamap1 = Map(mwadata1, aiahead1)
        ii = "%04d" % j
        p0 = map171.plot(axes=ax0, aspect="auto", vmin=0, vmax=800)
        # frac_r1 = [0.05,0.1,0.12,0.2,0.3,0.4,0.5]
        for f in lev_ct[k]:
            # lev_r1 = np.nanmax(mwamap.data) * f
            c1 = mwamap1.contour(level=f * u.ct)
            if len(c1) != 0:
                ax0.plot_coord(c1[0], color=cc[k])
        ax0.text(
            5400,
            -700 + k * 500,
            "Freq: " + str(np.round(freq[chan[k]], 1)) + " MHz",
            style="italic",
            bbox={"facecolor": cc[k], "alpha": 0.5, "pad": 10},
        )
    ax0.set_ylim(-1600, 7000)
    ax0.set_xlim(-1600, 7000)
    idx1 = ut.find_predecessor(tstix_high, time_mwa[j])[0]
    ax1.plot(
        tstix_high[idx1 - 100 : idx1 + 100] - tstix_high[0],
        stix_data0[:, 1][idx1 - 100 : idx1 + 100],
        "o-",
        label="STIX (16-22 keV)",
        color="b",
    )
    # ax1.set_xticklabels(['03:03:42','03:14:12','03:21:14','03:28:16','03:36:20','03:44:58','03:55:11','04:03:29','04:14:13'])
    ax1.axvline(x=time_mwa[j] - tstix_high[0], color="k")
    ax1.set_ylabel("STIX Count Flux (counts/s/cm$^2$/keV)")
    ax1.set_xlabel("Time (HH:MM:SS) / Start time: 03:03:42")
    ax1.legend()
    plt.savefig("/home/rohit/20220915/20220915_pngs/img_ts_" + str(ii) + ".png")
    plt.close()


idx = ut.find_predecessor(time_aia171, time_mwa[j])
euv171_file = euv171_list[idx[0]]
map171 = Map(euv171_file)
aiahead1 = []
aiahead1 = map171.meta.copy()
aiahead1["naxis1"] = 200
aiahead1["naxis2"] = 200
aiahead1["CRPIX1"] = 99
aiahead1["CRPIX2"] = 99
aiahead1["CRVAL1"] = 0
aiahead1["CRVAL2"] = 0
aiahead1["CDELT1"] = 50
aiahead1["CDELT2"] = 50
xlmwa, xrmwa, ylmwa, yrmwa = -5000, 5000, -5000, 5000
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection=map171)
lev_ct = [
    [7.0e4, 1.0e5],
    [1.0e5, 2.0e5, 3.0e5],
    [1.0e5, 3.0e5],
    [8.0e4, 3.0e5],
    [2.0e3, 8.0e4, 1.0e5, 1.5e5, 2.0e5, 3.0e6],
]
cc = ["r", "k", "b", "g", "brown"]
for k in range(5):
    k1 = k + 13
    mwadata = np.nanmean(Tb[k1], axis=0)
    mwadata[np.isnan(mwadata)] = 0
    # mwadata1=mwadata*0;mwadata1[90:160,90:160] = mwadata[90:160,90:160]
    # mwadata1[mwadata1<0.5] = 0
    mwamap = Map(mwadata, aiahead1)  # ;mwamap1 = Map(mwadata1, aiahead1)
    ii = "%04d" % j
    p0 = map171.plot(axes=ax, aspect="auto", vmin=0, vmax=800)
    lev1 = np.array([20, 40]) * u.percent
    mwamap.draw_contours(
        axes=ax,
        levels=lev1,
        colors=cc[k],
        linewidths=3,
        extent=[xlmwa, xrmwa, ylmwa, yrmwa],
    )
    # frac_r1 = [0.05,0.1,0.12,0.2,0.3,0.4,0.5]
    # for f in lev_ct[k]:
    #    #lev_r1 = np.nanmax(mwamap.data) * f
    #    c1 = mwamap.contour(level=f * u.ct)
    #    if(len(c1)!=0):
    #        ax.plot_coord(c1[0], color=cc[k])
    ax.text(
        -1400,
        -700 + k * 500,
        "Freq: " + str(np.round(freq[k1], 1)) + " MHz",
        style="italic",
        bbox={"facecolor": cc[k], "alpha": 0.5, "pad": 10},
    )
ax.set_ylim(-1600, 7000)
ax.set_xlim(-1600, 7000)
for fline in field_lines:
    ax.plot_coord(fline.coords, color="red", linewidth=1)
plt.show()


j = 36
idx = ut.find_predecessor(time_aia171, time_mwa[j])
euv171_file = euv171_list[idx[0]]
map171 = Map(euv171_file)
aiahead1 = []
aiahead1 = map171.meta.copy()
aiahead1["naxis1"] = 200
aiahead1["naxis2"] = 200
aiahead1["CRPIX1"] = 99
aiahead1["CRPIX2"] = 99
aiahead1["CRVAL1"] = 0
aiahead1["CRVAL2"] = 0
aiahead1["CDELT1"] = 50
aiahead1["CDELT2"] = 50
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection=map171)
lev_ct = [[3.0e4], [3.0e4], [3.0e4, 5.0e4], [3.0e5, 8.0e5], [3.0e5, 8.0e5]]
# lev_ct = [[1.0e5, 0.5e5, 1.0e5], [2.0e4,3.0e4], [2.0e4, 3.0e4], [2.e5, 3.e5], [2.e3, 3.e5, 5.e5]]
cc = ["r", "g", "b", "k", "brown"]
for k in range(5):
    mwadata = np.nanmean(Tb[k][j - 2 : j + 2], axis=0)
    mwadata1 = mwadata * 0
    mwadata1[90:160, 90:160] = mwadata[90:160, 90:160]
    mwadata1[mwadata1 < 0.5] = 0
    mwamap = Map(mwadata, aiahead1)
    mwamap1 = Map(mwadata1, aiahead1)
    ii = "%04d" % j
    p0 = map171.plot(axes=ax, aspect="auto", vmin=0, vmax=800)
    # frac_r1 = [0.05,0.1,0.12,0.2,0.3,0.4,0.5]
    for f in lev_ct[k]:
        # lev_r1 = np.nanmax(mwamap.data) * f
        c1 = mwamap1.contour(level=f * u.ct)
        if len(c1) != 0:
            ax.plot_coord(c1[0], color=cc[k])
    ax.text(
        -1400,
        -700 + k * 500,
        "Freq: " + str(np.round(freq[chan[k]], 1)) + " MHz",
        style="italic",
        bbox={"facecolor": cc[k], "alpha": 0.5, "pad": 10},
    )
ax.set_ylim(-1600, 7000)
ax.set_xlim(-1600, 7000)
plt.show()

level = [0.6]
cc = ["r", "g", "b", "k", "brown"]
f, ax = plt.subplots(1, 1)
for i in range(5):
    ax.contour(
        Tb[i].mean(axis=0) / np.nanmax(Tb[i].mean(axis=0)),
        levels=level,
        colors=cc[i],
        origin="lower",
    )
pc = plt.Circle((130, 130), 19.2, fill=False)
ax.add_artist(pc)
plt.show()
