from astropy.io import fits
import glob
import pickle
from sunpy.map import Map
import matplotlib.pyplot as plt

listTb=sorted(glob.glob('/media/rohit/MWA/20140914/Tb*.p'))
aia=fits.open('/media/rohit/MWA/20140914/EUV/aia.lev1.171A_2014-09-14T03_06_11.34Z.image_lev1.fits')
freq=[108,161,240]
for j in range(3):
    aa=pickle.load(open(listTb[j],'rb'));Tb=aa[0];time=aa[9]
    for i in range(len(aa[0])):
        td='2014-09-14T0'+time[i]
        h=aia[0].header;d=aia[0].data 
        h['TELESCOP']='sdo';h['INSTRUME']='';h['WAVELNTH']=94;h['WAVEUNIT']='angstrom';h['CDELT2']=50.0;h['DATE-OBS']=td;h['T_REC']=td;h['DATE']=td;h['T_OBS']=td;h['ISPPKTIM']=td
        h['CDELT1']=50.0;h['CRPIX1']=49;h['CRPIX2']=49
        if((i==0) & (j==0)):
            h.remove('OSCNMEAN');h.remove('OSCNRMS')
        fits.writeto('/media/rohit/MWA/20140914/mwa_maps/mwa_'+str(freq[j])+'_'+str(time[i])+'.fits',Tb[i],h)
#mm=Map('test.fits');am=Map('/media/rohit/MWA/20140914/EUV/aia.lev1.171A_2014-09-14T03_06_11.34Z.image_lev1.fits')
#mm.plot();plt.show()


