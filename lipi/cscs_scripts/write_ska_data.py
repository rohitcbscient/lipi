import casacore.tables as ct
import numpy as np
import os
import oskar


i = int(os.environ["SLURM_ARRAY_TASK_ID"]) + 122


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


ii = "%04d" % i
skadc = ct.table("/store/ska/sk01/sdc3-new/MS/ZW3_IFRQ_" + ii + ".ms")
skadc_data = skadc.getcol("DATA")
skadc_uvw = skadc.getcol("UVW")
filename_ska = (
    "/scratch/snx3000/rsharma/ska_data_reduced/ZW3_IFRQ_for_karabo_" + ii + ".MS"
)
ms_new_ska = skadc_data[:, 0, 0]  # /np.nanmax(np.abs(skadc_data))
write_ms(filename_ska, ms_new_ska, skadc_uvw, 130816)
