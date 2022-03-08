import numpy as np
from astropy.io import fits
import glob

flist=['084-085','093-094','103-104','113-114','125-126','139-140','153-154','169-170','187-188']
freq=[109,121,134,147,162,180,198,218,241]

std_res_ons=[0]*len(flist)
std_res_ofs=[0]*len(flist)
mean=[0]*len(flist)
std=[0]*len(flist)
i=0
for ff in flist:
    ll=glob.glob('1133139776-%b'+ff+'.*.image.FITS')
    ll_res=glob.glob('1133139776-%b'+ff+'.*.image.FITS')
    mean[i]=[0]*len(ll)
    std_res_ons[i]=[0]*len(ll)
    std_res_ofs[i]=[0]*len(ll)
    data=[0]*len(ll)
    data_res=[0]*len(ll)
    j=0
    for l in ll:
        f=fits.open(l)
        data[j]=f[0].data[0][0]
        fres=fits.open(ll_res[j])
        data_res[j]=fres[0].data[0][0]
        std_res_ons[i][j]=np.std(data_res[j][310:330,623:640])
        std_res_ofs[i][j]=np.std(data_res[j])
        j=j+1
    mean[i]=np.array(data).mean(axis=0)
    std[i]=np.array(data).std(axis=0)
    print freq[i], mean[i].max(), np.median(std_res_ons[i]),  np.median(std_res_ofs[-1])
    i=i+1


# 623, 330
# 640, 310







