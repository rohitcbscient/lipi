import os   
import numpy as np, gc, sys
from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime, timedelta

i = int(os.environ['SLURM_ARRAY_TASK_ID'])+750


def simulate_TEC(screens, screen_width_metres, bmax, rate, sampling, num_times, frequency, fits_filename):
    """
    Inputs:
    - screen_width_metres (float): self explanatory
    - r0 (float): Scale size in km
    - bmax (float): sub-aperture size in km
    - sampling (float): meter per pixel
    - speed (float): speed of the layer in m/s (default for 150 km/s)
    - rate (float): the inverse frame rate in 1/s
    - alpha (float): evolve the screen slovely
    - num_times (float): length of the observation in sec
    - frequency (float): in Hz
    - fits_filename (string): where to store the simulated screen
    """ 
    m = int(bmax / sampling)  # Pixels per sub-aperture
    n = int(screen_width_metres / bmax)  # Sub-apertures across the screen
    num_pix = n * m
    pscale = screen_width_metres / (n * m)  # Pixel scale (in m/pixel).
    print("Number of pixels %d, pixel size %.3f m" % (num_pix, pscale))
    print("Field of view %.1f (m)" % (num_pix * pscale))
    # Convert to TEC
    phase2tec = -frequency / 8.44797245e9 * 1e-2 # TODO: Here the factor 1e-2 is introduced by the SDC3, and it represent the outcome of a successful DD calibration.
    w = WCS(naxis=4)
    w.naxis = 4
    w.wcs.cdelt = [pscale, pscale, 1.0/rate, 1.0]
    w.wcs.crpix = [num_pix // 2 + 1, num_pix // 2 + 1, num_times // 2 + 1, 1.0]
    w.wcs.ctype = ['XX', 'YY', 'TIME', 'FREQ']
    w.wcs.crval = [0.0, 0.0, 0.0, frequency]
    data = phase2tec*screens.sum(axis=0)[None, ...]
    fits.writeto(filename=fits_filename, data=data, header=w.to_header(), overwrite=True)
    del data, screens
    gc.collect()
    return 0

t_start = datetime(2021, 9, 21, 14, 12, 40, 0)
t_obs = timedelta(hours=4, minutes=0, seconds=0, milliseconds=0)
t_day = timedelta(hours=4, minutes=0)
t_int = timedelta(seconds=10)
nr_tsteps = int(t_day.total_seconds() / t_int.total_seconds())
nr_days_obs = int(t_obs.total_seconds() / t_day.total_seconds())

ii="%03d" % i
freqs = 1.06e8+np.arange(901) * 1.e5 # Hz
freq=freqs[i]
path_in = '/scratch/snx3000/rsharma/'
screen_file = path_in+'atmo/screen_4h_i0.npy'
my_screen = np.load(file=screen_file, mmap_mode='r')
iono_fits = screen_file.replace('.npy', '_ch%03d.fits' %i)
r0, sampling = 7e3, 100.0
if not (os.path.exists(iono_fits)):
	simulate_TEC(screens=my_screen, screen_width_metres=200e3, rate=0.1, bmax=20e3, sampling=sampling, num_times=nr_tsteps, frequency=freq, fits_filename=iono_fits)

