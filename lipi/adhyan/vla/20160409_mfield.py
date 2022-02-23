import numpy as np
from astropy.io import fits
from surya.utils import Bextrap
from astropy import units as u


ff='/media/rohit/VLA/20160409_EUV/hmi.M_720s.20160409_183417.E18N10CR.CEA.NAS.sav'
bsmap,bcarrmap=Bextrap.get_gxs_sav2hpp(ff,'2016-04-09T18:34:17.00')
file_name_radio='/media/rohit/VLA/paraview/EUV_loop_radio.csv'
x,y,z,bx,by,bz=Bextrap.get_fieldlines('/media/rohit/VLA/paraview/EUV_loop_radio.csv')
x,y,z,bx,by,bz,b_hp,b_proj=Bextrap.transform_fieldlines(x,y,z,bx,by,bz,'2016/04/09T18:45:00',bsmap[0].wcs,bcarrmap[0])

#####
