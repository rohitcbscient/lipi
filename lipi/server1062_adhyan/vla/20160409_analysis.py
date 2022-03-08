import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
from scipy.signal import correlate

def autocorrelation (x) :
    """
    Compute the autocorrelation of the signal, based on the properties of the
    power spectral density of the signal.
    """
    xp = x-np.mean(x)
    f = np.fft.fft(xp)
    p = np.array([np.real(v)**2+np.imag(v)**2 for v in f])
    pi = np.fft.ifft(p)
    return np.real(pi)[:x.size/2]/np.sum(xp**2),p


spwlist=[0,1,2,3,4,5,6,7]
freqcntr=[1057,1185,1313,1441,1557,1685,1813,1941]

sumall=[0]*len(spwlist)
maxall=[0]*len(spwlist)
i=0
for s in spwlist:
    print 'SPW: '+str(s)
    filelist=sorted(glob.glob('*.spw.'+str(s)+'_16-31*.FITS'))
    sumall[i]=[0]*len(filelist)
    maxall[i]=[0]*len(filelist)
    for j in range(len(filelist)):
        f=fits.open(filelist[j])
        sumall[i][j]=np.sum(f[0].data)
        maxall[i][j]=np.max(f[0].data)
    i=i+1

sumall=np.array(sumall)
maxall=np.array(maxall)*100/848

yl=[[20,250],[50,150],[20,60],[20,40],[15,50],[20,40],[20,40],[20,40]]
f,ax=plt.subplots(nrows=8,sharex=True)
for i in range(8):
    ax[i].plot(maxall[i],'-',color='k',label='SPW='+str(i))
    ax[i].set_ylim(yl[i])
    ax[i].legend()
ax[4].set_ylabel('Flux (SFU)')
ax[-1].set_xlabel('Time (s)')
plt.show()

maxcor=autocorrelation(maxall[0])
maxauto=correlate(maxall[0],maxall[0],mode='full')

for i in range(8):
    plt.hist(maxall[i]-np.mean(maxall[i]),bins=100,histtype='step',label='SPW='+str(i))
plt.legend()
plt.xlabel('Flux (SFU)')
plt.show()

