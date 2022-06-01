import numpy as np
import matplotlib.pyplot as plt
import glob
from sunpy.map import Map
import pickle

############333
tmwa,tsubmwa,Tbsubmax,xcsub90,ycsub90,maxsubX,maxsubY,pa_angle=pickle.load(open('/media/rohit/MWA/20140914/Tb_centroid.p','rb'))

############### AIA

ddiff_list=sorted(glob.glob('/sdata/fits/running_diff/*171*ddiff.fits'))
bdiff_list=sorted(glob.glob('/sdata/fits/running_diff/*171*bdiff.fits'))
linex=np.arange(3100,3500);liney=np.arange(1270,1670)[::-1]
line_check=1
if line_check:
    mapp=Map(ddiff_list[0]);d=mapp.data
    mapp.plot()
    plt.plot(linex,liney,color='k')
    plt.show()

intd171=[0]*len(ddiff_list)
time171=[0]*len(ddiff_list)
for i in range(len(ddiff_list)):
    mapp=Map(ddiff_list[i]);d=mapp.data;intd171[i]=[0]*len(linex)
    time171[i]=ddiff_list[i].split('T')[1].split('Z')[0]
    for j in range(len(linex)):
        intd171[i][j]=d[liney[j]-5:liney[j]+5,linex[j]-5:linex[j]+5].mean()
intd171=np.array(intd171).swapaxes(0,1)[:,::-1]

intb171=[0]*len(bdiff_list)
for i in range(len(bdiff_list)):
    mapp=Map(bdiff_list[i]);d=mapp.data
    intb171[i]=d[liney,linex]
intb171=np.array(intb171).swapaxes(0,1)[:,::-1]

f,ax=plt.subplots(4,1,sharex=True);ax0=ax[0];ax1=ax[1];ax2=ax[2];ax3=ax[3]
ax0.imshow(intd171[:,::-1],aspect='auto',origin='lower',vmin=-60,vmax=60,cmap='coolwarm')
ax0.set_xticks(np.arange(1311)[::200]);ax0.set_xticklabels(time171[::200]);ax0.set_xlabel('Time (HH:MM:SS)')
ax0.set_yticks(np.arange(400)[::50]);ax0.set_yticklabels(np.arange(400)[::50]*726/1.e3);ax0.set_ylabel('Radial Distance (Mm)')
ax1.plot(225+np.arange(len(Tbsubmax[0])),Tbsubmax[0]/1.e6,'o-');ax1.set_yscale('log')
ax2.plot(225+np.arange(len(ycsub90[0])),ycsub90[0],'o-');ax2.set_ylim([-600,0])
ax3.plot(225+np.arange(len(ycsub90[0])),pa_angle,'o')
ax3.set_xlabel('Time (HH:MM:SS) UT');ax3.set_ylabel('P.A. (degrees)');ax2.set_ylabel('Y-Coordinate (arcsec)');ax1.set_ylabel('$T_B$ (MK)')
plt.show()


