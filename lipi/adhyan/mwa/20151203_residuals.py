import numpy as np
from astropy.io import fits
import glob

freq_list=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
maxi=[0]*len(freq_list)
std=[0]*len(freq_list)
i=0
for f in freq_list:
    print f
    l=sorted(glob.glob('1*'+f+'*.FITS'))[0:20]
    maxi[i]=[0]*len(l)
    std[i]=[0]*len(l)
    j=0
    for ll in l:
        ff=fits.open(ll)
        d=ff[0].data[0][0]
        maxi[i][j]=np.max(d)
        std[i][j]=np.std(d[np.where(d==maxi[i][j])[0][0]-10:np.where(d==maxi[i][j])[0][0]+10,np.where(d==maxi[i][j])[1][0]-10:np.where(d==maxi[i][j])[1][0]+10])
        j=j+1
    i=i+1

maxi=np.array(maxi)
std=np.array(std)

