import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
from scipy.io import readsav
from sunpy.map import Map
import os

sav94d=sorted(glob.glob('/sdata/fits/*.94A*ddiff.sav'))
sav171d=sorted(glob.glob('/sdata/fits/*.171A*ddiff.sav'))
sav193d=sorted(glob.glob('/sdata/fits/*.193A*ddiff.sav'))
sav211d=sorted(glob.glob('/sdata/fits/*.211A*ddiff.sav'))

sav94b=sorted(glob.glob('/sdata/fits/*.94A*bdiff.sav'))
sav171b=sorted(glob.glob('/sdata/fits/*.171A*bdiff.sav'))
sav193b=sorted(glob.glob('/sdata/fits/*.193A*bdiff.sav'))
sav211b=sorted(glob.glob('/sdata/fits/*.211A*bdiff.sav'))

fit94=sorted(glob.glob('/sdata/fits/*.94A*.fits'))
fit171=sorted(glob.glob('/sdata/fits/*.171A*.fits'))
fit131=sorted(glob.glob('/sdata/fits/*.131A*.fits'))
fit193=sorted(glob.glob('/sdata/fits/*.193A*.fits'))
fit211=sorted(glob.glob('/sdata/fits/*.211A*.fits'))

def convert_fits(sav94b,sav94d):
    for i in range(len(sav94d)):
        print(i,sav94d[i])
        #os.system('rm -rf '+sav94d[i].split('_dd')[0].split('.fit')[0]+'.ddiff.fits')
        savd=readsav(sav94d[i])['dd2']
        hdu = fits.PrimaryHDU(savd)
        mapd=fits.open(sav94d[i].split('_dd')[0]);mapd[0].verify('fix')
        hdu.header=mapd[0].header;hdulist = fits.HDUList([hdu])
        hdulist.writeto(sav94d[i].split('_dd')[0].split('.fit')[0]+'.ddiff.fits',overwrite=True)
        #os.system('rm -rf '+sav94d[i])
        #---------------
        #os.system('rm -rf '+sav94b[i].split('_bd')[0].split('.fit')[0]+'.bdiff.fits')
        savb=readsav(sav94b[i])['dd1']
        hdu = fits.PrimaryHDU(savb)
        mapb=fits.open(sav94b[i].split('_bd')[0]);mapb[0].verify('fix')
        hdu.header=mapb[0].header;hdulist = fits.HDUList([hdu])
        hdulist.writeto(sav94b[i].split('_bd')[0].split('.fit')[0]+'.bdiff.fits',overwrite=True)
        #os.system('rm -rf '+sav94b[i])

def convert_sunpy(sav94b,sav94d):
    for i in range(len(sav94d)):
        print(sav94d[i])
        mapd=Map(sav94d[i].split('_dd')[0])
        savd=readsav(sav94d[i])['dd2']
        mapd1=Map(savd,mapd.meta)
        #---------------
        mapb=Map(sav94b[i].split('_bd')[0])
        savb=readsav(sav94b[i])['dd1']
        mapb1=Map(savb,mapb.meta)
        np.save(sav94b[i]+'.sunpy',[mapb1,mapd1])

convert_fits(sav171b,sav171d)
convert_fits(sav193b,sav193d)
convert_fits(sav211b,sav211d)

convert_sunpy(sav94b,sav94d)
convert_sunpy(sav193b,sav193d)

