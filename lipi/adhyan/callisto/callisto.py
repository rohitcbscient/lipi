from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import glob

ooty_lst=sorted(glob.glob('OOTY*58.fit'))[:-3]
alaska_lst=sorted(glob.glob('ALASKA*58.fit'))[2:-2]
def get_ecallisto_data(lst):
    d=[0]*len(lst)
    for i in range(len(lst)):
        f=ooty_lst[i]
        a=fits.open(f)
        d[i]=a[0].data
    d=np.array(d).swapaxes(0,1)
    nd=np.array(d).shape
    d=np.array(d).reshape((nd[0],nd[1]*nd[2]))
    d1=d-np.array(list(d[:,0])*d.shape[1]).reshape((d.shape[1],d.shape[0])).swapaxes(0,1)
    return d1,d

data_ooty=get_ecallisto_data(ooty_lst)[0]
data_alas=get_ecallisto_data(alaska_lst)[0]

