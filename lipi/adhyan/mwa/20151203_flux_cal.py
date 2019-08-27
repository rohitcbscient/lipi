
from surya.mwa import mwa_fluxcal as mfl
from surya.mpis import do_parallel as dp
import glob


mwa_phase=1
DS_DIR = '/media/rohit/MWA/20151203_HerA/'                   # Path to directory containing DS files
BEAM_DIR = '/media/rohit/MWA/MWA_BEAMS/'             # Path to directory containing beam files
WORKING_DIR = '/media/rohit/MWA/20151203_HerA/'              # Path to working directory
HASLAM_DIR = '/media/rohit/MWA/20151203_HerA/haslam/'

DS_LIST=glob.glob('1*.p')
#DS='1133329472_177MHz_T000-008.DS.dat.p'

# MWA Coordinates

array_lon = 116.6708            # MWA core coordinates (in degrees)
array_lat = -26.7033
array_elevation = 377.83

# Coordinates of source
star_ra=252.784
star_dec=4.992

rec_path='/home/i4ds1807205/scripts/flux_calibration/Trec.p'
grd_path='/home/i4ds1807205/scripts/flux_calibration/Tpickup.p'
spidx_path='/home/i4ds1807205/scripts/flux_calibration/haslam_spec_gal_guzman.p'
ifsun=0
#### FOR ONE FILE
#mfl.tsun_computation(DS,DS_DIR,BEAM_DIR,WORKING_DIR,HASLAM_DIR,mwa_phase,array_lon,array_lat,array_elevation,rec_path,grd_path,spidx_path,ifsun,star_ra,star_dec)
#### FOR MPI
dp.MPI_MAIN_FLUX(mfl,DS_LIST,DS_DIR,BEAM_DIR,WORKING_DIR,HASLAM_DIR,mwa_phase,array_lon,array_lat,array_elevation,rec_path,grd_path,spidx_path,ifsun,star_ra,star_dec)
