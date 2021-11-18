from astropy.io import fits
import os

# Create a black image
path='/media/rohit/MWA/20190405_PSP_MWA/sun_1238472032/192MHz/'
imgpath=path+'point_sun.fits'
#dir_='J2000 01h06m15.387231s +09d46m47.60506s';f=107.5
#dir_='J2000 00h54m37.5s +05d50m36.2s';f=107.5
dir_='J2000 00h54m37.5s +05d50m36.2s';f=191.8
aa=fits.open(imgpath)
#aa[0].data=aa[0].data*0;aa[0].header['CRVAL3']=240.6e6
aa[0].data=aa[0].data*0;aa[0].header['CRVAL3']=f*1.e6
os.system('rm '+path+'/blank.fits')
aa.writeto(path+'blank.fits')
cl.addcomponent(flux=1.0, fluxunit='Jy', polarization='Stokes',freq=str(f)+'MHz', dir=dir_, shape='gaussian', majoraxis='20arcsec', \
                minoraxis='20arcsec', positionangle='0deg')
os.system('rm -rf gaussian_model.cl')
cl.rename('gaussian_model.cl')

os.system('rm -rf '+path+'blank.im')
importfits(fitsimage=path+'blank.fits',imagename=path+'blank.im')
ia.open(path+'blank.im')
ia.modify(cl.torecord(),subtract=False)

cl.close()
cl.done()
ia.close()
