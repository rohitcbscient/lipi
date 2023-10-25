import os
import tempfile
import matplotlib.pyplot as plt
import numpy as np
import pytest

from karabo.imaging.imager import Imager
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from karabo.sourcedetection.evaluation import SourceDetectionEvaluation
from karabo.sourcedetection.result import (
    PyBDSFSourceDetectionResult,
    SourceDetectionResult,
)
from karabo.test.conftest import NNImageDiffCallable, TFiles
from astropy.io import fits
import bdsf
import matplotlib
from astropy.wcs import WCS

#matplotlib.use('TkAgg')
def do_bdsf_old(imagepath,output_model,idx_f):
        img=fits.open(imagepath+'_I.fits')
        head=img[0].header
        data=img[0].data[0]
        head['CDELT3']=1
        head['BMIN']=10/3600. # in Deg 2.8 arcsec for 120 km baseline
        head['BMAJ']=10/3600.
        out_imagepath=imagepath+'_mod.fits'
        fits.writeto(out_imagepath,data,head,overwrite=True)
        out_imagedata=fits.open(out_imagepath)[0].data
        img_bdsf = bdsf.process_image(out_imagepath, beam=(10/3600.,10/3600.,0),thresh_pix=5,thresh='hard')
        nisland=len(img_bdsf.islands);gauss_list=[0]*nisland
        gauss_n_size=nisland#len(img_bdsf.gaussians)
        gauss_position=[0]*gauss_n_size
        gauss_position_pix=[0]*gauss_n_size
        gauss_peak_flux=[0]*gauss_n_size
        gauss_size=[0]*gauss_n_size
        for i in range(gauss_n_size):
                img_gauss=img_bdsf.islands[i]
                gauss_position_pix[i]=np.mean(img_gauss.border,axis=1)
                gauss_position[i]=(-gauss_position_pix[i]+1024)*abs(head['CDELT1'])
                gauss_position[i][1]=-30-gauss_position[i][1]
                gauss_peak_flux[i]=img_gauss.max_value
                gauss_size[i]=(2*np.sqrt(img_gauss.border[0].shape[0]/np.pi)*0.667)*abs(head['CDELT1'])
        gauss_position=np.array(gauss_position)
        gauss_position_pix=np.array(gauss_position_pix)
        gauss_peak_flux=np.array(gauss_peak_flux)
        gauss_size=np.array(gauss_size)
        sky_data=np.hstack((gauss_position,gauss_peak_flux.reshape(gauss_n_size,1),np.zeros((gauss_n_size,8))))
        sky_data[:,9]=gauss_size;sky_data[:,10]=gauss_size
        np.savetxt(X=sky_data,fname=output_model)
        return sky_data,gauss_position,gauss_position_pix,gauss_peak_flux,gauss_size

def do_bdsf(imagepath,output_model):
    img=fits.open(imagepath)
    head=img[0].header
    data=img[0].data[0]
    head['CDELT3']=1
    head['BMIN']=10/3600. # in Deg 2.8 arcsec for 120 km baseline
    head['BMAJ']=10/3600.
    out_imagepath='/scratch/snx3000/rsharma/sdc3point_ch750_4h1d_I_mod.fits'
    fits.writeto(out_imagepath,data,head,overwrite=True)
    out_imagedata=fits.open(out_imagepath)[0].data
    img_bdsf = bdsf.process_image(out_imagepath, beam=(10/3600.,10/3600.,0),thresh_pix=1,thresh_isl=2,thresh='hard')
    bdsf_source_list='/scratch/snx3000/rsharma/bdsf_source_list.fits'
    img_bdsf.write_catalog(outfile=bdsf_source_list,catalog_type='gaul',format='fits', clobber=True)
    sl=fits.open(bdsf_source_list)
    RA=sl[1].data['RA'];DEC=sl[1].data['DEC'];sl_peak=sl[1].data['Peak_flux']
    RA[np.where(RA>200)]=RA[np.where(RA>200)]-360
    sl_bmaj=sl[1].data['Maj']
    sl_bmin=sl[1].data['Min']
    sl_bpa=sl[1].data['PA']
    gauss_n_size=RA.shape[0]
    sky_data=np.hstack((RA.reshape(gauss_n_size,1),DEC.reshape(gauss_n_size,1),sl_peak.reshape(gauss_n_size,1),np.zeros((gauss_n_size,9))))
    sky_data[:,9]=sl_bmaj
    sky_data[:,10]=sl_bmin
    sky_data[:,11]=sl_bpa
    np.savetxt(X=sky_data,fname=output_model)
    w = WCS(head) 
    #xpix,ypix,_=w.wcs_world2pix(RA,DEC,0,0)
    xpix,ypix= sl[1].data['Xposn'],sl[1].data['Yposn']
    return sky_data,data,xpix,ypix


beam_path='/store/ska/sk014/sdc3/station_beam.fits'
beam=fits.open(beam_path);beam_head=beam[0].header;beam_data=beam[0].data
#imagepath='/scratch/snx3000/rsharma/sdc3point_ch750_4h1d_I.fits'
imagepath='/scratch/snx3000/rsharma/subvis_tests/image_sdc3point_ch750_4h1d_I.fits'
output_model='/scratch/snx3000/rsharma/sdc3point_ch750_4h1d_I_mod_bdsf_skymodel.txt'


sky_data,data,xpix,ypix = do_bdsf(imagepath,output_model)
f,ax=plt.subplots(1,1)
im=ax.imshow(data,aspect='auto',origin='lower',vmin=-1.e-3,vmax=1.e-3)
ax.plot(xpix,ypix,'o',color='red')
f.colorbar(im)
plt.show()

