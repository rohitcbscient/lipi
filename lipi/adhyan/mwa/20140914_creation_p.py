from surya.mwa import get_pickle as gp
from astropy.io import fits
import glob

#############################################################################################################################################
INDATA_DIR = '/nas08-data02/rohit/20140914/187-188/raw_data/'                # Path of directory where MS files are sitting
OUTDATA_DIR = '/nas08-data02/rohit/20140914/pickle'       # Path of directory to which the output DS files are directed
DB_DIR  = '/nas08-data02/rohit/20140914/187-188/raw_data/'           # Path of directory where database of MWA events is located
WORKING_DIR = '/nas08-data02/rohit/20140914/187-188/raw_data/'               # Path of working directory

POL_LIST = ['XX','YY']          # List of polarisation for computing autocorrelations
CPOL_LIST = ['XX','XY','YX','YY']       # List of polarisations for computing crosscorrelations
tilenames_list = ['Tile011MWA', 'Tile021MWA', 'Tile022MWA','Tile023MWA']       # List of tile names
#Long-baslines
#tilenames_list = ['Tile118MWA', 'Tile121MWA', 'Tile142MWA','Tile164MWA']       # List of tile names
MWA_PHASE=1 # Put 1 or 2
#tilenames_list = ['Tile074MWA', 'Tile075MWA', 'Tile076MWA','Tile077MWA']        # List of tile names (PHASE-II)

if(MWA_PHASE==2):
      tile1=[62,60,26,25,0,16]
      tile2=[63,61,28,27,7,17]
      baseline_list=['62-63', '60-61', '26-28', '25-27', '0-7', '16-17']
########
metalist=sorted(glob.glob('/nas08-data02/rohit/20140914/meta/*.metafits'));m=len(metalist);azi_pointing=[0]*m;ele_pointing=[0]*m
for i in range(m):
    header=fits.open(metalist[i])
    azi_pointing[i]=header[0].header['AZIMUTH']
    ele_pointing[i]=header[0].header['ALTITUDE']
#######
MS_LIST=sorted(glob.glob('*.ms'))
for i in range(m):
    gp.main(INDATA_DIR,OUTDATA_DIR,DB_DIR,WORKING_DIR,[MS_LIST[i]],POL_LIST,CPOL_LIST,MWA_PHASE,azi_pointing[i],ele_pointing[i],tilenames_list,ms,qa,tb,casalog)

