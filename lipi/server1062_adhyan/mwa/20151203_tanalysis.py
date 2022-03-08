import numpy as np
import matplotlib.pyplot as plt
import pickle

def scale_161MHz(d):
    fact=d[1650]/d[1750]
    d[1712:2297]=fact*d[1712:2297]
    fact1=d[2250]/d[2350]
    d[2297:]=fact1*d[2297:]
    return d

mwa_ar1xl,mwa_ar1xr,mwa_ar1yl,mwa_ar1yr=57,62,40,48
mwa_qsxl,mwa_qsxr,mwa_qsyl,mwa_qsyr=49,53,42,47
mwa_chxl,mwa_chxr,mwa_chyl,mwa_chyr=43,53,52,58

Tb_file='Tb_20151203_169-170.p'
res_file='res_20151203_169-170.p'

img=pickle.load(open(Tb_file,'r'))
res=pickle.load(open(res_file,'r'))

Tb=np.array(img[0])
flag=np.where(np.array(img[2])>1)[0]

Tbmax=[0]*Tb.shape[0]
for i in range(Tb.shape[0]):
    Tbmax[i]=np.max(Tb[i])

for i in range(100):
    for j in range(100):
        Tb[:,i,j]=scale_161MHz(Tb[:,i,j])
        #Tb[:,i,j][] 

Tbmax=np.array(Tbmax)
Tbmax[flag]=np.nan
Tbmax[1712:2297]=1.0472*Tbmax[1712:2297]

