import numpy as np
import matplotlib.pyplot as plt

specfile='/media/rohit/VLA/20160409/sun_L_20160409T1844-1846UT.50ms.cal.ms.dspec.median.npz'
data=np.load(specfile) 
tim_=data['tim'];freq=data['freq'];ds_=data['spec']
tim=[0]*8
for i in range(8):
    tim[i]=np.hstack((tim_[2089*i:2089*(i+1)],tim_[16712+i*311:16712+(i+1)*311]))
tim=np.array(tim)
ds_=ds_[:,0,:,:]
ds=[0]*2
for i in range(2):
    ds[i]=[0]*64
    for j in range(64):
        ds[i][j]=[0]*8
        for k in range(8):
            ds[i][j][k]=np.hstack((ds_[i,j,2089*k:2089*(k+1)],ds_[i,j,16712+k*311:16712+(k+1)*311]))

ds=np.array(ds)
ds_LL=ds[0].swapaxes(0,1).reshape(512,2400)
ds_RR=ds[1].swapaxes(0,1).reshape(512,2400)
freq=freq.swapaxes(0,1).flatten()/1.e9
plt.imshow(ds_LL,aspect='auto',interpolation='None',origin=0,vmin=20,vmax=100)
#plt.imshow(ds_RR,aspect='auto',interpolation='None',origin=0,vmin=5,vmax=60)
plt.xticks(np.arange(2400)[::240],['18:44:00','18:44:12','18:44:24','18:44:36','18:44:48','18:45:00','18:45:12','18:45:24','18:45:36','18:45:48','18:46:00'])
plt.yticks(np.arange(512)[::50],freq[::50])
plt.xlabel('Time (HH:MM:SS UT)');plt.ylabel('Frequency (GHz)')
plt.show()
