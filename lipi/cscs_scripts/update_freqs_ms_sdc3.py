"""
Program for updating reference frequency in MS (Correcting from last runs!) 
"""

import casacore.tables as ct
import numpy as np
import os
import oskar
import astropy.units as u


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


# i = int(os.environ['SLURM_ARRAY_TASK_ID'])+0
# Channel number is given by i, which will be parallelised in SLURM batch file
i = 0
ii = "%04d" % i
# Frequency array
freqs = 1.06e8 + np.arange(901) * 1.0e5  # Hz
# Input file
ms_file = (
    "/scratch/snx3000/rsharma/residual1_pybdsf/ms/residual_sdc3point_ch" + ii + "_02.MS"
)
# Output file
ms_file_updated = (
    "/scratch/snx3000/rsharma/residual1_pybdsf/ms/residual_sdc3point_freq_updated_ch"
    + ii
    + "_02.MS"
)
# Read the MS file
ms = ct.table(ms_file, readonly=True)
ms_vis = ms.getcol("DATA")[:, 0, 0]
ms_uvw = ms.getcol("UVW")
freq = ms.SPECTRAL_WINDOW.getcol("CHAN_FREQ").squeeze() * u.Hz
n = 130816  # Total number of baselines
# Write the MS file
write_ms(
    ms_file_updated, ms_vis, ms_uvw, n, freqs[i]
)  # The frequency (freqs[i]) is  updated in MS file
os.system("rm -rf " + ms_file)  # Remove old MS to save disk space
