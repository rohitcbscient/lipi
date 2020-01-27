from surya.mwa import get_pickle as gp
from astropy.io import fits

#############################################################################################################################################
INDATA_DIR = '/media/rohit/MWA/20151203_HerA/'                # Path of directory where MS files are sitting
OUTDATA_DIR = '/media/rohit/MWA/20151203_HerA/'       # Path of directory to which the output DS files are directed
DB_DIR  = '/media/rohit/MWA/20151203_HerA/'           # Path of directory where database of MWA events is located
WORKING_DIR = '/media/rohit/MWA/20151203_HerA/'               # Path of working directory

POL_LIST = ['XX','YY']          # List of polarisation for computing autocorrelations
CPOL_LIST = ['XX','XY','YX','YY']       # List of polarisations for computing crosscorrelations
tilenames_list = ['Tile011MWA', 'Tile021MWA', 'Tile022MWA','Tile023MWA']       # List of tile names
MWA_PHASE=1 # Put 1 or 2
#tilenames_list = ['Tile074MWA', 'Tile075MWA', 'Tile076MWA','Tile077MWA']        # List of tile names (PHASE-II)

if(MWA_PHASE==2):
      tile1=[62,60,26,25,0,16]
      tile2=[63,61,28,27,7,17]
      baseline_list=['62-63', '60-61', '26-28', '25-27', '0-7', '16-17']
########
metafits='/media/rohit/MWA/20151203_HerA/1133329472.metafits'
header=fits.open(metafits)
azi_pointing=header[0].header['AZIMUTH']
ele_pointing=header[0].header['ALTITUDE']
#######

MS_LIST=['1133329472_177MHz.MS']
gp.main(INDATA_DIR,OUTDATA_DIR,DB_DIR,WORKING_DIR,MS_LIST,POL_LIST,CPOL_LIST,MWA_PHASE,azi_pointing,ele_pointing,tilenames_list,ms,qa,tb,casalog)

