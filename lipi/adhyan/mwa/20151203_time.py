import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
import sys


filelist=sorted(glob.glob('*.fts'))[6:]
sum0=[0]*len(filelist)
sum1=[0]*len(filelist)
sum2=[0]*len(filelist)
header=[0]*len(filelist)
c2=[0]*len(filelist)
d0=fits.open(filelist[0])
time_c2=[0]*len(filelist)
for i in range(len(filelist)):
    print filelist[i]
    d=fits.open(filelist[i])
    diff=d[0].data-d0[0].data
    c2[i]=diff
    sum0[i]=np.mean(diff)
    sum1[i]=np.mean(diff[400:600,600:850])
    sum2[i]=np.mean(diff[420:520,200:300])
    header[i]=d[0].header
    time_c2[i]=header[i]['TIME-OBS']

plot_lasco=1
if(plot_lasco):
    for i in range(len(c2)):
        plt.imshow(c2[i],origin=0,aspect='equal',vmin=-200,vmax=200)
        plt.colorbar()
        plt.show()

plot_mean_c2=1
if(plot_mean_c2):
    plt.plot(np.arange(len(sum1))*12,sum1,'o-',label='NP')
    plt.plot(np.arange(len(sum1))*12,sum2,'o-',label='SP')
    plt.legend()
    plt.xlabel('Time (mins)');plt.ylabel('Mean')
    plt.show()

sys.exit()
#################
filelist=sorted(glob.glob('aia*.fits'))
aia193data=[0]*len(filelist)
sum0=[0]*len(filelist)
sum1=[0]*len(filelist)
sum2=[0]*len(filelist)
header=[0]*len(filelist)
aia193=[0]*len(filelist)
d0=fits.open(filelist[0])
#for i in range(len(filelist)):
for i in range(50):
    i+120
    print filelist[i]
    d=fits.open(filelist[i])
    aia193[i]=d[0].data
    diff=d[0].data-d0[0].data
    #aia193data[i]=np.sum(d[0].data[3000:,1000:3000])
    sum0[i]=np.mean(diff)
    sum1[i]=np.mean(diff[2800:3700,1000:3000])
    sum2[i]=np.mean(diff[250:600,1700:2400])
    header[i]=d[0].header



