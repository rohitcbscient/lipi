import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import healpy as hp
import healpy.visufunc



with h5.File('/home/rohit/simulations/ethz/RadioSims/map_0.h5', 'r') as f:
    print(f.keys())
    DS1=f['map'][0,3,:]
    print(DS1.shape)
    out=hp.mollview(DS1,cmap='jet',return_projected_map=True)
    plt.imshow(out);plt.colorbar()
    plt.savefig('/home/rohit/simulations/ethz/RadioSims/HP_f10_p3.png',dpi=300)
