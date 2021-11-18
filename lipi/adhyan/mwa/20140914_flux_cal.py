
from surya.mwa import mwa_fluxcal as mfl
from surya.mpis import do_parallel as dp
import glob


mwa_phase=1
DS_DIR = '/media/rohit/MWA/20140914/pickle/'                   # Path to directory containing DS files
BEAM_DIR = '/media/rohit/MWA/MWA_BEAMS/'             # Path to directory containing beam files
WORKING_DIR = '/media/rohit/MWA/20140914/pickle/'              # Path to working directory
HASLAM_DIR = '/media/rohit/MWA/haslam/'

DS_LIST=sorted(glob.glob('1*169-*.p'))
#DS='1133329472_177MHz_T000-008.DS.dat.p'
#for i in range(len(DS_LIST)):
#    os.system('mv '+DS_LIST[i]+' '+DS_LIST[i].split('-%b')[0]+'_'+DS_LIST[i].split('-%b')[1])

# MWA Coordinates

array_lon = 116.6708            # MWA core coordinates (in degrees)
array_lat = -26.7033
array_elevation = 377.83

# Coordinates of source
star_ra=171.72
star_dec=3.569
rec_path='/home/i4ds1807205/scripts/flux_calibration/Trec.p'
grd_path='/home/i4ds1807205/scripts/flux_calibration/Tpickup.p'
spidx_path='/home/i4ds1807205/scripts/flux_calibration/haslam_spec_gal_guzman.p'
ifsun=1
#### FOR ONE FILE
for DS in DS_LIST:
    mfl.tsun_computation(DS,DS_DIR,BEAM_DIR,WORKING_DIR,HASLAM_DIR,mwa_phase,array_lon,array_lat,array_elevation,rec_path,grd_path,spidx_path,ifsun,star_ra,star_dec)
#### FOR MPI
#dp.MPI_MAIN_FLUX(mfl,DS_LIST,DS_DIR,BEAM_DIR,WORKING_DIR,HASLAM_DIR,mwa_phase,array_lon,array_lat,array_elevation,rec_path,grd_path,spidx_path,ifsun,star_ra,star_dec)
