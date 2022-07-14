import pickle
import glob
import numpy as np
import matplotlib.pyplot as plt
import itertools

plist=sorted(glob.glob("*T011-012*.p"))
plist.remove(plist[14]);plist.remove(plist[7]);plist.remove(plist[20]);ncorr=[0]*len(plist);tmid=[0]*len(plist)
for i in range(len(plist)):
    aa=pickle.load(open(plist[i],'rb'),encoding='latin1')
    ncorr[i]=aa[5][0].reshape(48,24)
    tmid[i]=aa[12].split(' ')[1]
ncorr=np.array(ncorr).swapaxes(0,1).reshape(48,24*len(plist))
#ncorr=list(itertools.chain(*ncorr))
plt.imshow(ncorr,origin='lower',aspect="auto",vmin=0.1,vmax=0.5,cmap='jet')
for i in range(len(plist)):
    plt.axvline(i*24,color='white')
plt.xticks(np.arange(len(plist))*24,tmid);plt.xlabel('Time (25th June 2022 HH:MM:SS UT)')
plt.show()


