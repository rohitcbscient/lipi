from astropy.io import fits
import numpy as np
import glob

make_model=0
if(make_model):
    #aa=fits.open('/nas08-data02/rohit/20140914/calibrators/cal/hyd.fits')
    aa=fits.open('/nas08-data02/rohit/20140914/calibrators/cal_cen/cenA.fits')
    h=aa[0].header;d=aa[0].data
    #h['CRVAL3']=2.4e8;d[np.where(d<4.e-3)]=0
    h['CRVAL3']=1.606e8;d[np.where(d<0.8)]=0
    fits.writeto('/nas08-data02/rohit/20140914/calibrators/cal_cen/cenA.fits_161MHz.fits', d, h)

do_applycal=0
if(do_applycal):
    listms=sorted(glob.glob('*.ms'))
    for ll in listms:
        print ll
        applycal(vis=ll,field="",spw="",intent="",selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",msselect="",docallib=False,callib="",gaintable=['1094715672_187-188.cal.p'],gainfield=[],interp=[],spwmap=[],calwt=[True],parang=False,applymode="calonly",flagbackup=True)

merge_ms=0
if(merge_ms):
    allms=sorted(glob.glob('*.ms'))


