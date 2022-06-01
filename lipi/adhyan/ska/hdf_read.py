import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import healpy as hp
import healpy.visufunc
from astropy.io import fits
from astropy.wcs import WCS


def h5py_dataset_iterator(g, prefix=''):
    for key, item in g.items():
        path = '{}/{}'.format(prefix, key)
        if isinstance(item, h5.Dataset): # test for dataset
            yield (path, item)
        elif isinstance(item, h5.Group): # test for group (go down)
            yield from h5py_dataset_iterator(item, path)


filename='/home/rohit/simulations/ethz/galaxy_fullPol_512.h5'
with h5.File(filename, 'r') as f:
    for (path, dset) in h5py_dataset_iterator(f):
        print(path, dset)


with h5.File('/home/rohit/simulations/ethz/RadioSims/map_0.h5', 'r') as f:
    print(f.keys())
    DS1=f['map'][0,3,:]
    print(DS1.shape)
    out=hp.mollview(DS1,cmap='jet',return_projected_map=True)
    plt.imshow(out);plt.colorbar()
    #plt.savefig('/home/rohit/simulations/ethz/RadioSims/HP_f10_p3.png',dpi=300)
    plt.show()


with h5.File('/home/rohit/simulations/ethz/galaxy_fullPol_512.h5', 'r') as f:
    for (path, dset) in h5py_dataset_iterator(f):
        print(path, dset)
    print(f.keys())
    mapp=f['map'][:]
    DS1=f['map'][0,3,:]
    print(DS1.shape)
    out=hp.mollview(DS1,cmap='jet',return_projected_map=True)
    plt.imshow(out);plt.colorbar()
    #plt.savefig('/home/rohit/simulations/ethz/RadioSims/HP_f10_p3.png',dpi=300)
    plt.show()

group='index_map'
data = h5.File('/home/rohit/simulations/ethz/galaxy_fullPol_512.h5', 'r')
for dset in data[group].keys():
    print(dset)
    freq1=data[group]['freq'][:]
    pixel=data[group]['pixel'][:]

freqlist=freq.dtype.names
freq=freq1['centre']
width=freq1['width']

mapp_I=mapp[:,0,:]
mapp_Q=mapp[:,1,:]
mapp_U=mapp[:,2,:]
mapp_V=mapp[:,3,:]

ra,dec=hp.pix2ang(ipix=pixel,nside=512,lonlat=True)

out=hp.mollview(mapp_I[-1],cmap='jet',return_projected_map=True)
plt.show()


aa=fits.open('/home/rohit/simulations/MIGHTEE/MIGHTEE_Continuum_Early_Science_COSMOS_Level1.fits')
mightee_head=aa[0].header;mightee_data=aa[0].data
mightee_cat_head=aa[1].header;mightee_cat_ra=aa[1].data['RA'];mightee_cat_dec=aa[1].data['DEC']

#aa=fits.open('/home/rohit/simulations/MIGHTEE/MIGHTEE_Continuum_Early_Science_XMMLSS_r0p0_circ.fits')
#mightee_img_head=aa[0].header;mightee_img_data=aa[0].data

wcs=WCS(mightee_head)
wcsimg=WCS(mightee_img_head)
img_skycoord=wcsimg.pixel_to_world(0,0,0,0)

plt.subplot(projection=wcsimg,slices=('x','y',0,0))
plt.imshow(mightee_img_data,origin='lower',aspect='auto',vmin=1.e-7,vmax=1.e-4)
plt.show()

wcscat=WCS(mightee_cat_head)
plt.subplot(projection=wcscat,slices=('x','y'))
plt.plot(mightee_cat_ra,mightee_cat_dec,'o',markersize=1,color='k')
plt.show()

