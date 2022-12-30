import numpy as np
import matplotlib.pyplot as plt
import glob
from sunpy.map import Map
import pickle
from surya.utils import main as ut
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
from matplotlib.patches import Ellipse


############333
tmwa,tsubmwa,Tbsubmax,xcsub90,ycsub90,maxsubX,maxsubY,pa_angle=pickle.load(open('/media/rohit/MWA/20140914/Tb_centroid.p','rb'))

############### AIA

ddiff_list=sorted(glob.glob('/sdata/fits/running_diff/*171*ddiff.fits'))
bdiff_list=sorted(glob.glob('/sdata/fits/running_diff/*171*bdiff.fits'))
linex=np.arange(3100,3500);liney=np.arange(1270,1670)[::-1]
mapp=Map(ddiff_list[0]);d=mapp.data
p2w=mapp.pixel_to_world(linex*u.pix,liney*u.pix);linex_arcsec=p2w.Tx.value;liney_arcsec=p2w.Ty.value
liner_arcsec=np.sqrt(linex_arcsec**2 + liney_arcsec**2)
rcsub90=np.sqrt(xcsub90**2 + ycsub90**2)
line_check=1
if line_check:
    mapp=Map(ddiff_list[0]);d=mapp.data
    mapp.plot(vmin=-800,vmax=1000)
    plt.plot(linex,liney,color='red',linewidth=3)
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


y1=ut.find_nearest(liner_arcsec,rcsub90[-1][80])[0]
y2=ut.find_nearest(liner_arcsec,rcsub90[-1][-1])[0]
liner_plot=[0]*rcsub90.shape[1]
for i in range(rcsub90.shape[1]):
    liner_plot[i]=ut.find_nearest(liner_arcsec,rcsub90[-1][i])[0]

c = SkyCoord(xcsub90*u.arcsec, ycsub90*u.arcsec, frame=mapp.coordinate_frame)
w2p=mapp.world_to_pixel(c)
linex_pix=w2p.x.value;liney_pix=w2p.y.value
liner_pix=np.sqrt(linex_pix**2 + liney_pix**2)
speedx=np.linspace(185,630,10)
speedy=np.linspace(130,350,10)[::-1]

speedupx=np.linspace(131,1300,10)
speedupy=np.linspace(131,247,10)

f,ax0=plt.subplots(1,1)
ax0.imshow(intd171[:,::-1],aspect='auto',origin='lower',vmin=-60,vmax=60,cmap='coolwarm')
ax0.set_xticks(np.arange(1311)[::200]);ax0.set_xticklabels(time171[::200])
ax0.set_yticks(np.arange(400)[::50]);ax0.set_yticklabels(np.round(liner_arcsec[::50],1))
#ax0.plot(225+np.arange(len(Tbsubmax[0]))[80:],liner_pix[-1][80:],'o-',color='yellow',linewidth=0.5,markersize=2)
ax1=ax0.twinx()
ax1.plot(158+np.arange(len(Tbsubmax[0]))[80:],rcsub90[-1][80:],'o-',color='k',linewidth=1,markersize=2,label='240 MHz Type-IV radial location')
ax1.set_ylim(liner_arcsec[0],liner_arcsec[-1])
ax0.axvline(x=562,linestyle='--',color='k',label='FH-Transition')
ax0.axvline(x=217,linestyle='-.',color='r',label='Type-IV Start')
ax0.axvline(x=1056,linestyle='-.',color='g',label='Type-IV Ends')
ax0.plot(speedx,speedy,'-',color='cyan',linewidth=3,label='Slope: -0.05"/s $\\approx$ -36 km/s')
ax0.plot(speedupx,speedupy,'-',color='lime',linewidth=3,label='Slope: 0.01"/s $\\approx$ 7.2 km/s')
#ax0.axhline(y=y1,color='r',label='240 MHz Location at Start')
#ax0.axhline(y=y2,color='g',label='240 MHz Location at End')
ax0.legend(loc=1);ax1.legend(loc=2);ax0.set_ylabel('Radial Coordinate (arcsec)');ax0.set_xlabel('Time (HH:MM:SS UT)')
e1 = Ellipse(xy=(240, 315), width=100, height=70, edgecolor='magenta', fc='None', lw=2);ax0.add_patch(e1)
e1 = Ellipse(xy=(205, 250), width=80, height=50, edgecolor='magenta', fc='None', lw=2);ax0.add_patch(e1)
plt.show()

f,ax=plt.subplots(4,1,sharex=True);ax0=ax[0];ax1=ax[1];ax2=ax[2];ax3=ax[3]
ax0.imshow(intd171[:,::-1],aspect='auto',origin='lower',vmin=-60,vmax=60,cmap='coolwarm')
ax0.set_xticks(np.arange(1311)[::200]);ax0.set_xticklabels(time171[::200]);ax0.set_xlabel('Time (HH:MM:SS)')
ax0.set_yticks(np.arange(400)[::50]);ax0.set_yticklabels(np.arange(400)[::50]*726/1.e3);ax0.set_ylabel('Radial Distance (Mm)')
ax1.plot(225+np.arange(len(Tbsubmax[0])),Tbsubmax[0]/1.e6,'o-');ax1.set_yscale('log')
ax2.plot(225+np.arange(len(ycsub90[0])),ycsub90[0],'o-');ax2.set_ylim([-600,0])
ax3.plot(225+np.arange(len(ycsub90[0])),pa_angle,'o')
ax3.set_xlabel('Time (HH:MM:SS) UT');ax3.set_ylabel('P.A. (degrees)');ax2.set_ylabel('Y-Coordinate (arcsec)');ax1.set_ylabel('$T_B$ (MK)')
plt.show()


f,ax=plt.subplots(1,1)
ax.plot(rhessi_spec[1].data['RATE'][0],'o-',label='RATE')
ax.plot(rhessi_spec[1].data['PHMODEL'][0],'o-',label='PHMODEL')
ax.plot(rhessi_spec[1].data['CONVFAC'][0],'o-',label='CONVFAC')
#ax.plot(rhessi_spec[1].data['RESIDUAL'][0],label='RESIDUAL')
ax.plot(rhessi_spec[1].data['CONVFAC'][0]+rhessi_spec[1].data['PHMODEL'][0],'o-',label='CONVFAC+PHMODEL')
ax.set_xscale('log');ax.set_yscale('log');ax.legend()
plt.show()
